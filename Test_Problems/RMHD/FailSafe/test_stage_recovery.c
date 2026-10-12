/* Real ordinary C2P + PressureFix after exhausting production FOFC retries. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>
void printLog(const char *fmt, ...)
{ va_list a;va_start(a,fmt);vprintf(fmt,a);va_end(a); }
void Where(int i,Grid *g) {}
void PlutoError(int bad,char *msg) {if(bad)abort();}
void LogFileFlush(void) {fflush(stdout);}
static double *direct_fix_cell;
int __real_ConsToPrim(double **,double **,int,int,uint16_t *);
int __wrap_ConsToPrim(double **u,double **v,int beg,int end,uint16_t *flags)
{
  for(int i=beg;i<=end;i++) assert(u[i]!=direct_fix_cell);
  return __real_ConsToPrim(u,v,beg,end,flags);
}
int main(void)
{
  Data d={0};Grid grid={0};timeStep ts={0};RBox box;
  int i,j,k,n,pass;
  double v[NVAR]={0},u[NVAR],saved[NVAR];
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2;IEND=JEND=KEND=5;
  g_gamma=5.0/3.0;g_smallPressure=1.e-8;g_smallDensity=1.e-10;
  g_intStage=2;g_dt=0.25;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.Vs=ARRAY_4D(3,8,8,8,double);d.flag=ARRAY_3D(8,8,8,uint16_t);
  ts.cmax=ARRAY_1D(8,double);
  grid.dV=ARRAY_3D(8,8,8,double);
  TOT_LOOP(k,j,i) grid.dV[k][j][i]=1.0;
  RBoxDefine(IBEG,IEND,JBEG,JEND,KBEG,KEND,CENTER,&box);
  v[RHO]=1;v[PRS]=0.1;v[VX1]=0.2;v[BX1]=0.3;
  PrimToConsLoc(v,u);
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n){d.Uc[k][j][i][n]=u[n];d.Vc[n][k][j][i]=v[n];}
    DIM_LOOP(n)d.Vs[n][k][j][i]=v[BX1+n];
    d.flag[k][j][i]=FLAG_MINMOD; /* Pre-existing flags must survive. */
  }
#if PCP_FOFC == YES
  RMHD_PCP_BeginStage(&d,&ts);
  /* Repair the first rejected candidate without any FOFC retry. */
  d.Uc[3][3][3][ENG]=-1;
  direct_fix_cell=d.Uc[3][3][3];
  assert(RMHD_PCP_CheckStage(&d,&grid));
  assert(!RMHD_PCP_TryPressureFix(&d,1,&box,&grid));
  assert(!g_rmhd_failsafe_low_order);
  assert(!RMHD_PCP_FinishStage(&d,0,&box,&grid));
  RMHD_PCP_BeginStage(&d,&ts);
  for(pass=0;pass<=PCP_STAGE_LOCAL_PASSES;pass++) {
    d.Uc[3][3][3][ENG]=-1; /* Force the native pressure-fix branch. */
    assert(RMHD_PCP_CheckStage(&d,&grid));
    int retry=RMHD_PCP_RetryStage(&d,&ts);
    if (pass==0) assert(retry);
    if (!retry) break; /* No new seeds: identical FOFC pass is exhausted. */
  }
#else
  d.Uc[3][3][3][ENG]=-1;
  assert(RMHD_PCP_CheckStage(&d,&grid));
#endif
  memcpy(saved,d.Uc[3][3][3],sizeof saved);
  direct_fix_cell=d.Uc[3][3][3];
  assert(!RMHD_PCP_TryPressureFix(&d,1,&box,&grid));
  assert(!RMHD_PCP_FinishStage(&d,0,&box,&grid));
  assert(d.Vc[PRS][3][3][3]==g_smallPressure);
  assert(d.Uc[3][3][3][ENG]>0);
  NVAR_LOOP(n)if(n!=ENG)assert(d.Uc[3][3][3][n]==saved[n]);
  DOM_LOOP(k,j,i) {
    assert(d.flag[k][j][i]&FLAG_MINMOD);
    assert(!(d.flag[k][j][i]&(FLAG_LF|FLAG_CONS2PRIM_FAIL|FLAG_PRESSURE_FIX_FAIL)));
  }
  assert(!g_rmhd_failsafe_low_order);
  /* The next stage must not inherit the exhausted retry mask. */
#if PCP_FOFC == YES
  RMHD_PCP_BeginStage(&d,&ts);
#endif
  assert(!RMHD_PCP_CheckStage(&d,&grid));
#if PCP_FOFC == YES
  RMHD_PCP_EndStage(&d);
#endif
#if FAILSAFE == YES
  /* Unrepairable density must reach outer failsafe unchanged. */
#if PCP_FOFC == YES
  RMHD_PCP_BeginStage(&d,&ts);
#endif
  d.Uc[3][3][3][RHO]=-1;
  memcpy(saved,d.Uc[3][3][3],sizeof saved);
  assert(RMHD_PCP_CheckStage(&d,&grid));
  assert(RMHD_PCP_TryPressureFix(&d,1,&box,&grid));
  assert(!memcmp(saved,d.Uc[3][3][3],sizeof saved));
  assert(d.flag[3][3][3]&FLAG_CONS2PRIM_FAIL);
  /* Failed inverse must restore the candidate, including healthy cells. */
#if PCP_FOFC == YES
  RMHD_PCP_EndStage(&d);
  RMHD_PCP_BeginStage(&d,&ts);
#endif
  d.Uc[3][3][3][RHO]=1;
  d.Uc[3][3][3][MX1]=1.e300;
  memcpy(saved,d.Uc[3][3][3],sizeof saved);
  assert(RMHD_PCP_CheckStage(&d,&grid));
  assert(RMHD_PCP_TryPressureFix(&d,1,&box,&grid));
  assert(!memcmp(saved,d.Uc[3][3][3],sizeof saved));
  assert(d.flag[3][3][3]&FLAG_PRESSURE_FIX_FAIL);
#endif
  printf("PASS real stage PressureFix failsafe=%d reduced_energy=%d\n",FAILSAFE,RMHD_REDUCED_ENERGY);
  return 0;
}
