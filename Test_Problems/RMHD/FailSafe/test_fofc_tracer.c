/* Replay complete scalar fluxes; recompute only corrected LLF faces. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
int g_rmhd_failsafe_low_order;
void printLog(const char *fmt,...) {(void)fmt;}
void PlutoError(int bad,char *msg) {if(bad){puts(msg);abort();}}
void PrimToConsLoc(double *v,double *u)
{int n;NVAR_LOOP(n)u[n]=v[n];}
void LF_Solver(const Sweep *s,int beg,int end,double *cmax,Grid *g)
{
  int i,n;(void)g;
  for(i=beg;i<=end;i++) {
    NFLX_LOOP(n)s->flux[i][n]=0;
    s->flux[i][RHO]=i%2?-3.0:2.0;cmax[i]=1;
  }
}
int main(void)
{
  Data d={0};timeStep ts={0};Sweep s={0};Grid g={0};
  int i,j,k,n,dir;double cached[8];
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2;IEND=JEND=KEND=5;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.flag=ARRAY_3D(8,8,8,uint16_t);ts.cmax=ARRAY_1D(8,double);
  s.stateC.v=ARRAY_2D(8,NVAR,double);
  s.stateL.v=ARRAY_2D(8,NVAR,double);s.stateR.v=ARRAY_2D(8,NVAR,double);
  s.stateL.u=ARRAY_2D(8,NVAR,double);s.stateR.u=ARRAY_2D(8,NVAR,double);
  s.flux=ARRAY_2D(8,NVAR,double);s.press=ARRAY_1D(8,double);
  s.flag=ARRAY_1D(8,uint16_t);
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n){d.Uc[k][j][i][n]=0;d.Vc[n][k][j][i]=0;}
    d.flag[k][j][i]=0;
  }
  for(i=0;i<8;i++)ts.cmax[i]=1;
  RMHD_PCP_BeginStage(&d,&ts);
  DIM_LOOP(dir) {
    for(i=0;i<8;i++) {
      NVAR_LOOP(n){s.flux[i][n]=0;s.stateL.v[i][n]=s.stateR.v[i][n]=0;}
      s.flux[i][RHO]=i%2?-1.0:1.0;
      s.stateL.v[i][TRC]=0.1+0.01*i;
      s.stateR.v[i][TRC]=0.6+0.01*i;
      s.press[i]=0;
    }
    AdvectFlux(&s,1,5,&g);
    for(i=1;i<=5;i++) {
      cached[i]=s.flux[i][TRC];
      RMHD_PCP_CacheSweepFlux(dir,3,3,i,&s,ts.cmax,i);
    }
  }
  d.flag[3][3][3]=FLAG_CONS2PRIM_FAIL;
  assert(RMHD_PCP_RetryStage(&d,&ts));
  DIM_LOOP(dir) {
    for(i=0;i<8;i++) {
      s.flag[i]=i==3?FLAG_LF:0;
      NVAR_LOOP(n){s.stateC.v[i][n]=0;s.stateL.v[i][n]=s.stateR.v[i][n]=NAN;s.flux[i][n]=NAN;}
      s.stateC.v[i][TRC]=0.3+0.02*i;
    }
    for(i=1;i<=5;i++)RMHD_PCP_LoadSweepFlux(dir,3,3,i,&s,ts.cmax,i);
    RMHD_FailSafe_CorrectFluxes(&s,1,5,ts.cmax,&g);
    for(i=1;i<=5;i++) {
      double expected=i==2?2.0*s.stateC.v[i][TRC]:
                      i==3?-3.0*s.stateC.v[i+1][TRC]:cached[i];
      assert(s.flux[i][TRC]==expected);
      if(i!=2 && i!=3)assert(isnan(s.stateL.v[i][TRC]));
    }
    double delta=0;
    for(i=2;i<=5;i++)delta+=s.flux[i-1][TRC]-s.flux[i][TRC];
    assert(fabs(delta-(s.flux[1][TRC]-s.flux[5][TRC]))<1.e-14);
  }
  RMHD_PCP_EndStage(&d);
  puts("PASS tracer replay: cached untouched faces, fresh signed LLF scalars, conservation");
  return 0;
}

