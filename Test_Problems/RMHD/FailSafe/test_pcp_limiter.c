/* Exercise the real fixed-field limiter, including independent rest-state
 * energy bounds and randomized relativistic fluid states. */
#include "pluto.h"
#include <assert.h>
#include <float.h>
#include <stdint.h>

static void SetEnergy(double *u, double E)
{
#if RMHD_REDUCED_ENERGY == YES
  u[ENG] = E-u[RHO];
#else
  u[ENG] = E;
#endif
}

static double FullEnergy(const double *u)
{
#if RMHD_REDUCED_ENERGY == YES
  return u[ENG]+u[RHO];
#else
  return u[ENG];
#endif
}

static void Rest(double *u, double D, double E, double B)
{
  memset(u,0,NVAR*sizeof(*u));
  u[RHO]=D; u[BX1]=B; SetEnergy(u,E);
}

static uint32_t random_state=123456789;
static double Uniform(void)
{
  random_state=1664525u*random_state+1013904223u;
  return (double)random_state/4294967296.0;
}

static void PhysicalState(double *u)
{
  double rho=0.1+10.0*Uniform(), p=0.1+10.0*Uniform();
  double v[3], B[3], v2=0.0, B2=0.0, vB=0.0, W2, xi;
  int n;
  memset(u,0,NVAR*sizeof(*u));
  for (n=0;n<3;n++) {
    v[n]=0.8*(Uniform()-0.5);
    B[n]=20.0*(Uniform()-0.5);
    v2+=v[n]*v[n]; B2+=B[n]*B[n]; vB+=v[n]*B[n];
  }
  W2=1.0/(1.0-v2);
  xi=(rho+4.0*p)*W2; /* Gamma=4/3, rho h=rho+4p. */
  u[RHO]=rho*sqrt(W2);
  for (n=0;n<3;n++) {
    u[MX1+n]=(xi+B2)*v[n]-vB*B[n];
    u[BX1+n]=B[n];
  }
  SetEnergy(u,xi-p+0.5*(B2*(1.0+v2)-vB*vB));
}

int main(void)
{
  double anchor[NVAR], u[NVAR], saved[NVAR], trial[NVAR];
  int i,n,status,limited=0;
  const int fluid[]={RHO,MX1,MX2,MX3,ENG};

  Rest(anchor,1.0,7.0,2.0);
  memcpy(u,anchor,sizeof(u));
  assert(RMHD_PCP_Limit(anchor,u)==0);
  assert(memcmp(u,anchor,sizeof(u))==0);

  /* Psi-only and simultaneous q/Psi violations: independently, a rest
   * state with D=1 and B=2 requires full E>3. */
  for (i=0;i<2;i++) {
    Rest(u,1.0,i ? -1.0 : 1.5,2.0);
    assert(!RMHD_PCP_Check(u));
    assert(RMHD_PCP_Limit(anchor,u)==1);
    assert(RMHD_PCP_Check(u));
    assert(fabs(FullEnergy(u)-3.0)<1.e-10);
    assert(u[BX1]==2.0 && u[BX2]==0.0 && u[BX3]==0.0);
  }

  /* Density preprocessing must preserve full E, including reduced E. */
  Rest(anchor,1.0,20.0,1.0);
  Rest(u,-1.0,20.0,1.0);
  assert(RMHD_PCP_Limit(anchor,u)==1);
  assert(u[RHO]>8.0*DBL_EPSILON);
  assert(fabs(FullEnergy(u)-20.0)<=2.0*DBL_EPSILON*20.0);
  assert(RMHD_PCP_Check(u));

  /* A physical center need not remain physical with the face field;
   * rejection must leave the candidate untouched. */
  Rest(anchor,1.0,7.0,2.0);
  Rest(u,1.0,1.5,8.0);
  memcpy(saved,u,sizeof(u));
  assert(RMHD_PCP_Limit(anchor,u)==-1);
  assert(memcmp(saved,u,sizeof(u))==0);
  u[BX2]=NAN;
  memcpy(saved,u,sizeof(u));
  assert(RMHD_PCP_Limit(anchor,u)==-1);
  assert(memcmp(saved,u,sizeof(u))==0);
  assert(RMHD_PCP_Limit(anchor,NULL)==-1);

  for (i=0;i<6000;i++) {
    double E,theta,max_delta=0.0;
    int pivot=RHO;
    PhysicalState(anchor);
    assert(RMHD_PCP_Check(anchor));
    memcpy(u,anchor,sizeof(u));
    E=FullEnergy(u);
    u[RHO]*=-0.5+2.0*Uniform();
    for (n=0;n<3;n++) u[MX1+n]*=0.5+5.0*Uniform();
    SetEnergy(u,E*(-0.5+2.0*Uniform()));
    memcpy(saved,u,sizeof(u));
    memcpy(trial,u,sizeof(u));
    if (trial[RHO]<nextafter(8.0*DBL_EPSILON,INFINITY)) {
      trial[RHO]=nextafter(8.0*DBL_EPSILON,INFINITY);
      SetEnergy(trial,FullEnergy(saved));
    }
    status=RMHD_PCP_Limit(anchor,u);
    assert(status>=0);
    assert(RMHD_PCP_Check(u));
    for (n=0;n<NVAR;n++) {
      if (n==RHO || n==MX1 || n==MX2 || n==MX3 || n==ENG) continue;
      assert(memcmp(&saved[n],&u[n],sizeof(double))==0);
    }
    if (status==0) {
      assert(memcmp(saved,u,sizeof(u))==0);
      continue;
    }
    limited++;
    for (n=0;n<5;n++) {
      double delta=fabs(trial[fluid[n]]-anchor[fluid[n]]);
      if (delta>max_delta) { max_delta=delta; pivot=fluid[n]; }
    }
    assert(max_delta>0.0);
    theta=(u[pivot]-anchor[pivot])/(trial[pivot]-anchor[pivot]);
    assert(theta>=-1.e-12 && theta<=1.0+1.e-12);
    for (n=0;n<5;n++) {
      int f=fluid[n];
      double expected=anchor[f]+theta*(trial[f]-anchor[f]);
      assert(fabs(u[f]-expected)<1.e-10*fmax(1.0,fabs(expected)));
    }
  }
  assert(limited>1000);
  printf("PASS PCP limiter: rest bounds, density/full energy, atomic failure, "
         "fixed B, common contraction; %d/6000 random states limited\n",limited);
  return 0;
}
