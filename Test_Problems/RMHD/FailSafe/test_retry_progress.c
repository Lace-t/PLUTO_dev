/* Fixed-dt stagnation, growing masks, rollback and rank-local failures. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>

static int mode, calls;
void printLog(const char *fmt, ...)
{ va_list a; va_start(a,fmt); vprintf(fmt,a); va_end(a); }
void print(const char *fmt, ...) { (void)fmt; }
void Where(int i, Grid *g) { (void)i; (void)g; }
void LogFileFlush(void) { fflush(stdout); }
void PlutoError(int bad,char *msg) { if (bad) { puts(msg); abort(); } }
#ifdef PARALLEL
int AL_Exchange(void *buffer,int size)
{
  uint16_t other[8*8*8], *local=buffer;
  int i,j,k;
  (void)size;
  MPI_Sendrecv(local,8*8*8,MPI_UINT16_T,1-prank,0,
               other,8*8*8,MPI_UINT16_T,1-prank,0,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  for(k=0;k<8;k++) for(j=0;j<8;j++) for(i=0;i<2;i++) {
    if(prank==0) local[(k*8+j)*8+6+i]=other[(k*8+j)*8+2+i];
    else local[(k*8+j)*8+i]=other[(k*8+j)*8+4+i];
  }
  return 0;
}
#endif

int Integrate(Data *d,timeStep *ts,Grid *g)
{
  int i,j,k,n, owner=0;
  int extra=(mode==2 && RMHD_FAILSAFE_LOCAL_PASSES>1);
  double expected=ldexp(1.0,-MAX(0,calls-1-extra));
  (void)g;
  assert(g_dt==expected && g_time==0.0 && g_stepNumber==42);
  assert(ts->cmax[0]==7.0 && ts->invDt_hyp==9.0);
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) {
      assert(d->Uc[k][j][i][n]==1+n+0.01*i);
      assert(d->Vc[n][k][j][i]==2+n+0.01*j);
      d->Uc[k][j][i][n]=-123;
      d->Vc[n][k][j][i]=-456;
    }
    DIM_LOOP(n) {
      assert(d->Vs[n][k][j][i]==3+n+0.01*k);
      d->Vs[n][k][j][i]=-789;
    }
    d->flag[k][j][i]=0;
  }
  ts->cmax[0]=-1; ts->invDt_hyp=-2;
  calls++;
  if(mode==0 || (mode!=4 && g_dt==0.5)) return 0;
#ifdef PARALLEL
  if(mode==2 && calls>1) owner=1;
#endif
  if(prank==owner && mode!=3)
    d->flag[3][3][mode==2 && calls>1 ? 4:3] |= FLAG_CONS2PRIM_FAIL;
  return prank==owner;
}

static void CheckExit(void)
{
  if(mode==4) assert(calls==RMHD_FAILSAFE_MAX_RETRIES+2);
}

int main(int argc,char **argv)
{
  Data d={0}; Grid g={0}; timeStep ts={0};
  int i,j,k,n;
#ifdef PARALLEL
  int ranks;
  MPI_Init(&argc,&argv); MPI_Comm_rank(MPI_COMM_WORLD,&prank);
  MPI_Comm_size(MPI_COMM_WORLD,&ranks); assert(ranks==2);
#endif
  assert(argc==2); mode=atoi(argv[1]);
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2; IEND=JEND=KEND=5;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);
  d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.Vs=ARRAY_4D(3,8,8,8,double);
  d.flag=ARRAY_3D(8,8,8,uint16_t);
  ts.cmax=ARRAY_1D(8,double);
  for(n=0;n<3;n++) g.x[n]=ARRAY_1D(8,double);
  ts.cmax[0]=7; ts.invDt_hyp=9;
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) { d.Uc[k][j][i][n]=0; d.Vc[n][k][j][i]=0; }
    DIM_LOOP(n)d.Vs[n][k][j][i]=0;
    d.flag[k][j][i]=0;
  }

  /* First and newly growing masks retry; repeated seeds retain the rejected
   * candidate and pass accumulated seeds to the outer fallback. */
  RMHD_PCP_BeginStage(&d,&ts);
  if(prank==0)d.flag[3][3][3]=FLAG_CONS2PRIM_FAIL;
  assert(RMHD_PCP_RetryStage(&d,&ts));
#ifdef PARALLEL
  if(prank==1)d.flag[4][4][4]|=FLAG_CONS2PRIM_FAIL;
#else
  d.flag[4][4][4]|=FLAG_CONS2PRIM_FAIL;
#endif
  assert(RMHD_PCP_RetryStage(&d,&ts));
  if(prank==0)d.flag[3][3][3]|=FLAG_CONS2PRIM_FAIL;
  d.Uc[3][3][3][ENG]=12345;
  assert(!RMHD_PCP_RetryStage(&d,&ts));
  assert(d.Uc[3][3][3][ENG]==12345);
  RMHD_PCP_EndStage(&d);
  assert(!g_rmhd_failsafe_low_order);

  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) { d.Uc[k][j][i][n]=1+n+0.01*i; d.Vc[n][k][j][i]=2+n+0.01*j; }
    DIM_LOOP(n)d.Vs[n][k][j][i]=3+n+0.01*k;
    d.flag[k][j][i]=0;
  }
  g_dt=1; g_time=0; g_stepNumber=42;
  atexit(CheckExit);
  FailSafe(&d,&ts,&g);
  assert(mode!=4);
  assert(calls==(mode==0 ? 1 : 3+(mode==2 && RMHD_FAILSAFE_LOCAL_PASSES>1)));
  assert(!g_rmhd_failsafe_low_order);
  assert(d.Uc[3][3][3][ENG]==-123 && ts.cmax[0]==-1);
  puts("PASS stage growth/stagnation, full-step early halving, rollback, commit");
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
}
