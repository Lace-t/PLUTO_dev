#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>
#include <string.h>

void printLog(const char *fmt, ...) { (void)fmt; }
void Where(int i, Grid *grid) { (void)i; (void)grid; }
void PlutoError(int bad, char *msg) { if (bad) { puts(msg); abort(); } }
void LogFileFlush(void) {}

static void close_to(double a, double b) {
  assert(fabs(a-b) <= 1.e-8*fmax(1.e-20,fmax(fabs(a),fabs(b))));
}

int main(void)
{
  double u[NVAR]={0},v[NVAR]={0},saved[NVAR],pv[NVAR],mapped[NVAR];
  double b[3]={0.7,-0.4,0.2},delta[4],vel[3]={0.3,0.4,0};
  uint16_t flag=0;
  int n;
  g_gamma=5.0/3.0; g_smallPressure=1.e-12; g_smallDensity=1.e-11;
  /* Real remote rank-432 failure: old v^2 bisection exhausted 50 iterations. */
  {
    double remote[NVAR]={3.2600346692726572e-8,-5.5745491947143587e-5,
      -8.5558907718439635e-5,-5.1096013189279363e-8,0.0016656814710510795,
      0.002830913330261386,-0.00012672759120328919,0.00010754176560916415};
    double recovered[NVAR]={0},check[NVAR]={0};
    assert(!RMHD_PressureFix(remote,recovered));
    PrimToConsLoc(recovered,check);
    assert(recovered[PRS]==g_smallPressure);
    for(n=MX1;n<=MX3;n++)assert(fabs(remote[n]-check[n])<1.e-7*1.e-4);
  }
  /* Fixed-pressure inversions across Gamma and magnetic orientations. */
  for(int orientation=0;orientation<3;orientation++)
  for(int level=0;level<5;level++) {
    double gammas[5]={1,10,50,1000,10000};
    double initial[NVAR]={0}, cons[NVAR]={0},recovered[NVAR]={0};
    initial[RHO]=1.e-11; initial[PRS]=g_smallPressure;
    initial[VX1]=sqrt(1-1/(gammas[level]*gammas[level]));
    initial[BX1+orientation]=1.e-3;
    PrimToConsLoc(initial,cons);
    assert(!RMHD_PressureFix(cons,recovered));
    close_to(recovered[PRS],g_smallPressure);
    assert(fabs(recovered[VX1]-initial[VX1])<1.e-12);
    assert(fabs(recovered[RHO]-initial[RHO])<1.e-7*initial[RHO]);
  }
  /* Raw FV input captured on remote rank409, before any mapper mutates U.
   * Native recovery succeeds at p=3.2e-15, below the existing 1e-12 floor. */
  {
    double raw[NVAR]={3.9556726842597231e-9,-4.3020288269765366e-8,
      -2.0609819537442577e-7,4.7456939986752012e-7,-0.00030282605811526669,
      0.00046287248977126014,0.00033663252356898874,5.4750141969760301e-7};
    double target[3]={-0.00030314534597928775,0.00046312861045027316,
                      0.00033691882962382468}, recovered[NVAR]={0}, change[4];
    uint16_t f=0; double D=raw[RHO];
#if RMHD_REDUCED_ENERGY == NO
    raw[ENG]+=D;
#endif
    assert(!CT_RMHD_CorrectCell(raw,recovered,target,&f,change));
    assert(raw[RHO]==D && recovered[PRS]>=g_smallPressure);
    assert(RMHD_PCP_Check(raw));
  }
  v[RHO]=1.2; v[PRS]=0.3;
  for(n=0;n<3;n++){v[VX1+n]=vel[n];v[BX1+n]=0.1*(n+1);}
  v[TRC]=0.37;
  PrimToConsLoc(v,u); memcpy(saved,u,sizeof u);
  assert(!CT_RMHD_CorrectCell(u,v,b,&flag,delta));
  assert(u[RHO]==saved[RHO] && u[TRC]==saved[TRC]);
  close_to(v[RHO],1.2); close_to(v[PRS],0.3);
  for(n=0;n<3;n++){close_to(v[VX1+n],vel[n]);assert(u[BX1+n]==b[n]);}
  assert(RMHD_PCP_Check(u)); PrimToConsLoc(v,mapped);
  for(n=0;n<NVAR;n++)close_to(u[n],mapped[n]);
  for(n=0;n<3;n++)close_to(delta[n],u[MX1+n]-saved[MX1+n]);
  close_to(delta[3],u[ENG]-saved[ENG]);

  /* Independent electromagnetic stress/energy identity at fixed fluid state. */
  double oldb2=0.14,newb2=0.69,oldvb=0.11,newvb=0.05;
  for(n=0;n<3;n++)close_to(delta[n],(newb2-oldb2)*vel[n]
                              -newvb*b[n]+oldvb*0.1*(n+1));
  close_to(delta[3],0.5*(1.0+0.25)*(newb2-oldb2)
                              -0.5*(newvb*newvb-oldvb*oldvb));

  /* Recovery may repair energy: final state must still be consistent. */
  u[ENG]=0.0; flag=0;
  assert(!CT_RMHD_CorrectCell(u,v,b,&flag,delta));
  assert(u[RHO]==saved[RHO] && u[TRC]==saved[TRC]);
  assert(RMHD_PCP_Check(u)); PrimToConsLoc(v,mapped);
  for(n=0;n<NVAR;n++)close_to(u[n],mapped[n]);

  /* Invalid data and preexisting recovery failure must be atomic. */
  double positive_D=u[RHO];
  for(int test=0;test<4;test++){
    flag=0;
    if(test==0)u[RHO]=-1;
    if(test==1){u[RHO]=positive_D;b[1]=NAN;}
    if(test==2){b[1]=1.e300;}
    if(test==3){b[1]=-0.4;flag=FLAG_CONS2PRIM_FAIL;}
    memcpy(saved,u,sizeof u);memcpy(pv,v,sizeof v);
    assert(CT_RMHD_CorrectCell(u,v,b,&flag,delta));
    assert(flag&FLAG_CONS2PRIM_FAIL);
    assert(!memcmp(saved,u,sizeof u) && !memcmp(pv,v,sizeof v));
    for(n=0;n<4;n++)assert(delta[n]==0);
  }
  puts("CT correction: mapping, electromagnetic identity, repair, mass/scalar invariance and atomic failure passed");
  return 0;
}
