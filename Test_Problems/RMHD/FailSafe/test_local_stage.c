/* Real RK driver and stage rollback, with controlled flux/CT/C2P failures.
 * The separate CT test exercises actual shared-edge locality. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>
static int calls[4], mode, target, maps, recoveries;
void UserDefCFloors(Data *d, Grid *grid)
{
  (void)grid;
  if (mode==1 && g_rmhd_failsafe_low_order)
    d->flag[5][5][8]&=~FLAG_CONS2PRIM_FAIL;
}
void printLog(const char *fmt, ...) { va_list a; va_start(a,fmt); vprintf(fmt,a); va_end(a); }
void Where(int i, Grid *grid) {(void)i;(void)grid;}
void print(const char *fmt, ...) {(void)fmt;}
void PlutoError(int b,char *s) {if(b){puts(s);abort();}}
void LogFileFlush(void){}
void Boundary(const Data *d,int side,Grid *g){(void)d;(void)side;(void)g;}
void FlagShock(const Data *d,Grid *g){(void)d;(void)g;}
void PrimToCons3D(Data_Arr v,Data_Arr u,RBox *b,Grid *g){(void)v;(void)u;(void)b;(void)g;}
#ifdef PARALLEL
/* Two x-decomposed ranks, three ghost cells, six active cells. Transport
 * real MPI data, emulating the flag-only ArrayLib exchange used by FailSafe.
 */
int AL_Exchange(void *buffer, int size)
{
  uint16_t remote[12*12*12];
  uint16_t *local = buffer;
  int i,j,k;
  (void)size;
  MPI_Sendrecv(local,12*12*12,MPI_UNSIGNED_SHORT,1-prank,0,
               remote,12*12*12,MPI_UNSIGNED_SHORT,1-prank,0,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  for (k=0; k<12; k++) for (j=0; j<12; j++) for (i=0; i<3; i++) {
    if (prank == 0) local[(k*12+j)*12+9+i] = remote[(k*12+j)*12+3+i];
    else           local[(k*12+j)*12+i] = remote[(k*12+j)*12+6+i];
  }
  return 0;
}
#endif

void UpdateStage(Data *d,Data_Arr u,Data_Arr b,double **f,double dt,timeStep *Dts,Grid *g)
{
 int i,j,k,n,pass=calls[g_intStage]++;
 (void)f;(void)g;
 assert(dt==0.25);
 assert(Dts->cmax[0]==(g_intStage==1?10:20));
 if(g_intStage==target && pass) {
   assert(g_rmhd_failsafe_low_order);
   if(prank==0) {
     assert(d->flag[5][5][8]&FLAG_LF);
   }
#ifdef PARALLEL
   else { assert(d->flag[5][5][2]&FLAG_LF); }
#endif
   /* One seed, no reconstruction buffer, including across the MPI face. */
   TOT_LOOP(k,j,i) {
     assert(!(d->flag[k][j][i]&(FLAG_FLAT|FLAG_MINMOD|FLAG_HLL)));
     assert(!!(d->flag[k][j][i]&FLAG_LF) ==
            (k==5 && j==5 && i==(prank==0?8:2)));
   }
 }
 /* Restore must recover entry primitives, not failed CT or mapper output. */
 assert(d->Vc[RHO][5][5][5]==(g_intStage==1?1:2));
 DOM_LOOP(k,j,i) {u[k][j][i][RHO]+=0.01; DIM_LOOP(n)b[n][k][j][i]+=0.125;}
 Dts->cmax[0]=20; Dts->invDt_hyp=30;
 if(g_intStage==target && prank==0 && (mode>=2 || pass==0)) {
   u[5][5][8][MX1]=100; /* D/E positive: full candidate check must reject. */
   if(mode==6)b[IDIR][5][5][8]=NAN;
   if(mode==7)u[5][5][8][RHO]=-100;
 }
}
int RMHD_PressureFix(double *u,double *v) {(void)u;(void)v;return 1;}
int ConsToPrim(double **u,double **v,int beg,int end,uint16_t *flags)
{(void)u;(void)v;(void)beg;(void)end;(void)flags;return 1;}
void CT_AverageStaggeredFields(const Data *d,int correct,RBox *b,Grid *g)
{
 (void)correct;(void)b;(void)g;
 /* Simulate CT internal mapper touching V even on a rejected candidate. */
 d->Vc[RHO][5][5][5]=-99;
}
int RMHD_PCP_ConsToPrim3D_CheckModified(Data_Arr u,Data_Arr v,uint16_t ***flags,RBox *b,Grid *g)
{
 int i,j,k,n,bad=0;
 (void)u;(void)b;(void)g;
 maps++;
 DOM_LOOP(k,j,i) NVAR_LOOP(n) v[n][k][j][i]=n==RHO?2:0;
 if(g_intStage==target && g_rmhd_failsafe_low_order && mode==1 && prank==0) {
   flags[5][5][8]|=FLAG_CONS2PRIM_FAIL;bad=1;
 }
#ifdef PARALLEL
 int global;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);bad=global;
#endif
 if(g_intStage==target && mode==1 && calls[g_intStage]>1 && !bad)
   flags[5][5][8]&=~FLAG_CONS2PRIM_FAIL;
 return bad;
}
int ConsToPrim3D(Data_Arr u,Data_Arr v,uint16_t ***flags,RBox *b,Grid *g)
{
 int i,j,k,n,bad=0;
 (void)b;(void)g;
 assert(mode>=2 && mode!=6 && mode!=7);
 assert(calls[target]==PCP_STAGE_LOCAL_PASSES+1);
 assert(++recoveries==1);
 DOM_LOOP(k,j,i) {
   assert(!(flags[k][j][i]&(FLAG_CONS2PRIM_FAIL|FLAG_PRESSURE_FIX_FAIL)));
   NVAR_LOOP(n)v[n][k][j][i]=(n==RHO?2:n==PRS?0.1:0);
 }
 if(prank==0) {
   assert(u[5][5][8][MX1]>1); /* Recover the final candidate, not entry U. */
   u[5][5][8][ENG]=200;
   if(mode==2)bad=1; /* Error code without flags. */
   if(mode==4)v[VX1][5][5][8]=NAN; /* Success return with unusable V. */
   if(mode==5)flags[5][5][8]|=FLAG_PRESSURE_FIX_FAIL;
 }
#ifdef PARALLEL
 int global;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);bad=global;
#endif
 return bad;
}
int main(int argc,char **argv)
{
 Data d={0};Grid grid={0};timeStep ts={0};int i,j,k,n;
#ifdef PARALLEL
 int size;MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&prank);MPI_Comm_size(MPI_COMM_WORLD,&size);assert(size==2);
#endif
 mode=argc>1?atoi(argv[1]):1;target=argc>2?atoi(argv[2]):2;
 NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=12;
 IBEG=JBEG=KBEG=3;IEND=JEND=KEND=8;g_dt=.25;g_gamma=5.0/3;
 d.Uc=ARRAY_4D(12,12,12,NVAR,double);d.Vc=ARRAY_4D(NVAR,12,12,12,double);
 d.Vs=ARRAY_4D(3,12,12,12,double);d.flag=ARRAY_3D(12,12,12,uint16_t);
 grid.dV=ARRAY_3D(12,12,12,double);
 ts.cmax=ARRAY_1D(12,double);for(i=0;i<12;i++)ts.cmax[i]=10;
 TOT_LOOP(k,j,i){
   grid.dV[k][j][i]=1.0;
   NVAR_LOOP(n){d.Uc[k][j][i][n]=n==RHO?1:n==ENG?3:0;d.Vc[n][k][j][i]=n==RHO?1:0;}
   DIM_LOOP(n)d.Vs[n][k][j][i]=0;d.flag[k][j][i]=0;
 }
 int err=AdvanceStep(&d,&ts,&grid);
 if(mode==2 || mode>=4){
   assert(FAILSAFE==YES); /* FAILSAFE=NO must terminate before returning. */
   assert(err);assert(calls[target]==PCP_STAGE_LOCAL_PASSES+1);
   assert(recoveries==(mode==6 || mode==7 ? 0:1));
 }
 else {
   assert(!err);assert(calls[target]==(mode==0?2:3));
   assert(recoveries==(mode==3));
   assert(calls[target==1?2:1]==1);
   assert(fabs(d.Uc[5][5][5][RHO]-1.01)<1e-14);
   DOM_LOOP(k,j,i)assert(!(d.flag[k][j][i]&(FLAG_LF|FLAG_FLAT|FLAG_HLL|FLAG_MINMOD|FLAG_CONS2PRIM_FAIL)));
 }
 assert(g_rmhd_failsafe_low_order==0);
 printf("PASS local stage integrator=%d mode=%d target=%d calls=%d,%d,%d rank=%d\n",TIME_STEPPING,mode,target,calls[1],calls[2],calls[3],prank);
#ifdef PARALLEL
 MPI_Finalize();
#endif
 return 0;
}
