/* Real RK + consistent CT + native recovery. Count full-domain C2P calls
 * and verify accepted primitives against an independent native inversion. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>
static int full_maps, calls[4], fail_stage;
void printLog(const char *fmt, ...)
{ va_list a; va_start(a,fmt); vprintf(fmt,a); va_end(a); }
void print(const char *fmt, ...) { (void)fmt; }
void Where(int i,Grid *g) { (void)i; (void)g; }
void LogFileFlush(void) { fflush(stdout); }
void PlutoError(int bad,char *msg) { if(bad) { puts(msg); abort(); } }
void Boundary(const Data *d,int side,Grid *g) { (void)d; (void)side; (void)g; }
void FlagShock(const Data *d,Grid *g) { (void)d; (void)g; }
int __real_ConsToPrim3D(Data_Arr,Data_Arr,uint16_t ***,RBox *,Grid *);
int __wrap_ConsToPrim3D(Data_Arr u,Data_Arr v,uint16_t ***f,RBox *b,Grid *g)
{ full_maps++; return __real_ConsToPrim3D(u,v,f,b,g); }
int __real_RMHD_PCP_ConsToPrim3D_CheckModified(Data_Arr,Data_Arr,uint16_t ***,RBox *,Grid *);
int __wrap_RMHD_PCP_ConsToPrim3D_CheckModified(Data_Arr u,Data_Arr v,uint16_t ***f,RBox *b,Grid *g)
{ full_maps++; return __real_RMHD_PCP_ConsToPrim3D_CheckModified(u,v,f,b,g); }

void UpdateStage(Data *d,Data_Arr u,Data_Arr b,double **flux,double dt,timeStep *ts,Grid *g)
{
  int i,j,k;
  (void)d; (void)flux; (void)ts; (void)g;
  calls[g_intStage]++;
  DOM_LOOP(k,j,i) { u[k][j][i][ENG]+=0.01*dt; u[k][j][i][MX1]+=0.001*dt; }
  TOT_LOOP(k,j,i)b[IDIR][k][j][i]+=0.003*dt;
  /* Stay negative after every RK convex combination, not just stage 1. */
  if(g_intStage==fail_stage && calls[g_intStage]==1) u[3][3][3][RHO]=-100;
}

int main(int argc,char **argv)
{
  Data d={0}; Grid g={0}; timeStep ts={0};
  int i,j,k,n,s,stages=TIME_STEPPING==RK3 ? 3:2;
  double initial[NVAR]={0};
  assert(argc==2); fail_stage=atoi(argv[1]);
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2; IEND=JEND=KEND=5;
  g_gamma=5.0/3; g_smallDensity=g_smallPressure=1.e-12; g_dt=0.1;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double); d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.Vs=ARRAY_4D(3,8,8,8,double); d.flag=ARRAY_3D(8,8,8,uint16_t);
  ts.cmax=ARRAY_1D(8,double);
  initial[RHO]=1; initial[PRS]=0.2; initial[VX1]=0.3; initial[VX2]=-0.1;
  initial[BX1]=0.7; initial[BX2]=0.2; initial[BX3]=-0.4;
#if NTRACER > 0
  initial[TRC]=0.37;
#endif
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n)d.Vc[n][k][j][i]=initial[n];
    PrimToConsLoc(initial,d.Uc[k][j][i]);
    DIM_LOOP(n)d.Vs[n][k][j][i]=initial[BX1+n];
    d.flag[k][j][i]=0;
  }
  assert(AdvanceStep(&d,&ts,&g)==0);
  assert(full_maps==0);
  for(s=1;s<=stages;s++)assert(calls[s]==1+(s==fail_stage));
  DOM_LOOP(k,j,i) {
    double u[NVAR],v[NVAR]; uint16_t flag=0;
    NVAR_LOOP(n) { u[n]=d.Uc[k][j][i][n]; v[n]=d.Vc[n][k][j][i]; }
    assert(!ConsToPrimLoc(u,v,&flag));
    assert(!(flag&FLAG_CONS2PRIM_FAIL));
    NVAR_LOOP(n) assert(fabs(v[n]-d.Vc[n][k][j][i])<1.e-9*(1+fabs(v[n])));
    assert(!(d.flag[k][j][i]&(FLAG_CONS2PRIM_FAIL|FLAG_LF)));
  }
  assert(!g_rmhd_failsafe_low_order);
  puts("PASS real RK/CT primitive reuse, independent inversion, retry and mask cleanup");
  return 0;
}
