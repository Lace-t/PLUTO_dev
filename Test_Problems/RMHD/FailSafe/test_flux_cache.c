/* A retry restores poisoned entry speeds; cached faces must restore CFL too. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>

int g_rmhd_failsafe_low_order;
void printLog(const char *fmt, ...) { (void)fmt; }
void PlutoError(int bad, char *msg) { if (bad) { puts(msg); abort(); } }

int main(void)
{
  Data d = {0};
  timeStep ts = {0};
  Sweep s = {0};
  int i,j,k,n,dir,q;
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2; IEND=JEND=KEND=5;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);
  d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.flag=ARRAY_3D(8,8,8,uint16_t);
  ts.cmax=ARRAY_1D(8,double);
  s.flux=ARRAY_2D(8,NVAR,double);
  s.press=ARRAY_1D(8,double);
  TOT_LOOP(k,j,i) {
    d.flag[k][j][i]=0;
    NVAR_LOOP(n) {d.Uc[k][j][i][n]=0; d.Vc[n][k][j][i]=0;}
  }
  for(q=0;q<8;q++) ts.cmax[q]=NAN;
  RMHD_PCP_BeginStage(&d,&ts);
  DIM_LOOP(dir) for(q=1;q<=5;q++) {
    NVAR_LOOP(n) s.flux[q][n]=100*dir+10*q+n;
    s.press[q]=dir+0.1*q;
    ts.cmax[q]=0.2+0.1*dir+0.01*q;
    RMHD_PCP_CacheSweepFlux(dir,3,3,q,&s,ts.cmax,q);
  }
  d.flag[3][3][3]=FLAG_CONS2PRIM_FAIL;
  assert(RMHD_PCP_RetryStage(&d,&ts));
  assert(RMHD_PCP_UseFluxCache());
  for(q=0;q<8;q++) assert(isnan(ts.cmax[q]));
  DIM_LOOP(dir) for(q=1;q<=5;q++) {
    ts.cmax[q]=NAN;
    s.press[q]=NAN;
    NVAR_LOOP(n) s.flux[q][n]=NAN;
    RMHD_PCP_LoadSweepFlux(dir,3,3,q,&s,ts.cmax,q);
    assert(ts.cmax[q]==0.2+0.1*dir+0.01*q);
    assert(s.press[q]==dir+0.1*q);
    NVAR_LOOP(n) assert(s.flux[q][n]==100*dir+10*q+n);
  }
  RMHD_PCP_EndStage(&d);
  puts("PASS flux cache: retry restores each face's flux, pressure and CFL speed");
  return 0;
}
