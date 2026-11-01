/* Production stage driver: accepted flags, MPI progress and selective PF. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>

static int mode, updates, pf_calls, maps;
int g_rmhd_failsafe_low_order;
static const int cells[5][3]={{3,3,3},{3,3,4},{3,4,3},{4,3,3},{4,3,4}};

void printLog(const char *fmt, ...) {(void)fmt;}
void print(const char *fmt, ...) {(void)fmt;}
void LogFileFlush(void) {}
void PlutoError(int bad, char *msg) {if(bad){puts(msg);abort();}}
#ifdef PARALLEL
int AL_Exchange(void *buffer,int size) {(void)buffer;(void)size;return 0;}
#endif

static void reject(Data *d,int cell)
{
  const int *q=cells[cell];
  if(prank==0) d->Uc[q[0]][q[1]][q[2]][MX1]=100.0;
}

void UpdateStage(Data *d,Data_Arr u,Data_Arr b,double **f,double dt,
                 timeStep *ts,Grid *g)
{
  int i,j,k;
  (void)u;(void)b;(void)f;(void)dt;(void)g;
  assert(ts->cmax[0]==7.0 && ts->invDt_hyp==9.0);
  DOM_LOOP(k,j,i) {
    assert(d->Uc[k][j][i][MX1]==0.0);
    assert(!(d->flag[k][j][i]&FLAG_CONS2PRIM_FAIL));
  }
  if(updates==0) {
    reject(d,0);
    if(mode==1 || mode==3 || mode==4) reject(d,1);
    if(mode==3 || mode==4) reject(d,2);
  } else if(mode==1) {
    reject(d,1); /* no new seeds: stop after one retry */
  } else if(mode==2) {
    reject(d,1); /* same count but a new seed: one more correction needed */
  } else if(mode==3) {
    reject(d,1);reject(d,2);reject(d,3);
  } else if(mode==4 && updates==1) {
    reject(d,1);reject(d,3);
  } else if(mode==4 && updates==2) {
    reject(d,4);
  } else if(mode==5) {
    reject(d,updates); /* expanding seeds reach the configured retry cap */
  }
  d->Uc[5][5][5][MX1]=updates==0?1.0:2.0;
  if(mode==0 && updates>0 && prank==0)
    d->flag[3][3][3]|=FLAG_CONS2PRIM_FAIL|FLAG_PRESSURE_FIX_FAIL|
                      FLAG_NEGATIVE_DENSITY|FLAG_NEGATIVE_ENERGY;
  updates++;
  ts->cmax[0]=-1;ts->invDt_hyp=-2;
}

void UserDefCFloors(Data *d,Grid *g)
{
  (void)g;
  assert(d->Uc[5][5][5][MX1]==1.0); /* cache restore precedes floor hook */
}

int RMHD_PressureFix(double *u,double *v)
{
  assert(!RMHD_PCP_Check(u)); /* never PF a criterion-accepted seed */
  pf_calls++;
  u[ENG]=200.0;
  v[RHO]=1.0;v[PRS]=0.1;v[VX1]=v[VX2]=v[VX3]=0.0;
  return 0;
}

int ConsToPrim(double **u,double **v,int beg,int end,uint16_t *flags)
{
  int i,n;
  for(i=beg;i<=end;i++) {
    assert(RMHD_PCP_Check(u[i]));
    assert(!(flags[i]&FLAG_CONS2PRIM_FAIL));
    NVAR_LOOP(n) v[i][n]=0.0;
    v[i][RHO]=1.0;v[i][PRS]=0.1;
  }
  return 0;
}

int RMHD_PCP_ConsToPrim3D_CheckModified(Data_Arr u,Data_Arr v,
                                      uint16_t ***flags,RBox *box,Grid *g)
{
  int i,j,k,n;
  (void)box;(void)g;
  maps++;
  DOM_LOOP(k,j,i) {
    assert(RMHD_PCP_Check(u[k][j][i]));
    assert(!(flags[k][j][i]&FLAG_CONS2PRIM_FAIL));
    NVAR_LOOP(n) v[n][k][j][i]=0.0;
    v[RHO][k][j][i]=1.0;v[PRS][k][j][i]=0.1;
  }
  return 0;
}

int main(int argc,char **argv)
{
  Data d={0};Grid g={0};timeStep ts={0};RBox box;
  int i,j,k,n;
#ifdef PARALLEL
  int size;
  MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&prank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);assert(size==2);
#endif
  mode=atoi(argv[1]);
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2;IEND=JEND=KEND=5;
  g_gamma=5.0/3.0;g_intStage=1;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);
  d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.flag=ARRAY_3D(8,8,8,uint16_t);
  g.dV=ARRAY_3D(8,8,8,double);
  ts.cmax=ARRAY_1D(8,double);ts.cmax[0]=7;ts.invDt_hyp=9;
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) {d.Uc[k][j][i][n]=0;d.Vc[n][k][j][i]=0;}
    d.Uc[k][j][i][RHO]=1;d.Uc[k][j][i][ENG]=3;
    d.Vc[RHO][k][j][i]=1;d.Vc[PRS][k][j][i]=0.1;
    d.flag[k][j][i]=FLAG_MINMOD;g.dV[k][j][i]=1;
  }
  RBoxDefine(IBEG,IEND,JBEG,JEND,KBEG,KEND,CENTER,&box);
  assert(!RMHD_PCP_RunStage(&d,&ts,0.25,NULL,NULL,0,0.5,0.5,&box,&g));
  assert(PCP_STAGE_LOCAL_PASSES==2 || PCP_STAGE_LOCAL_PASSES==4);
  assert(updates==(mode==5?PCP_STAGE_LOCAL_PASSES+1:
                   mode==4?MIN(4,PCP_STAGE_LOCAL_PASSES+1):
                   mode==2 || mode==3?3:2));
  assert(pf_calls==(prank!=0?0:mode==0 || (mode==4 && PCP_STAGE_LOCAL_PASSES>=3)?0:mode==3?3:1));
  assert(maps==(mode==0 || (mode==4 && PCP_STAGE_LOCAL_PASSES>=3)));
  assert(!g_rmhd_failsafe_low_order);
  assert(d.Uc[5][5][5][MX1]==1.0);
  DOM_LOOP(k,j,i) {
    assert(d.flag[k][j][i]&FLAG_MINMOD);
    assert(!(d.flag[k][j][i]&(FLAG_LF|FLAG_CONS2PRIM_FAIL|FLAG_PRESSURE_FIX_FAIL|
                             FLAG_NEGATIVE_DENSITY|FLAG_NEGATIVE_ENERGY)));
  }
  printf("PASS stage progress mode=%d rank=%d updates=%d pf=%d\n",
         mode,prank,updates,pf_calls);
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
}

