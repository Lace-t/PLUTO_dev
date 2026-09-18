/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RMHD using a pressure fix.

  Fix p to a small value. For IDEAL EOS, bracket the enthalpy variable W
  with extended-precision residuals and verify the recovered momentum.
  For other EOS retain the original velocity-square bisection of Eq (A3).
  Once the root has been found, we recompute total energy E.
  Return 0 if succesful, 1 otherwise.

  \authors A. Mignone \n
           C. Zanni

  \date   Jul 03, 2021
*/ 
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <float.h>

#if EOS == IDEAL
/* At fixed p, W = D*Gamma + gamma/(gamma-1)*p*Gamma^2.
 * Both |v(W)|^2 and 1/Gamma(W)^2 decrease with W. Bracket their
 * sum minus one in W, not v^2: near light speed adjacent doubles in v^2
 * can straddle the root without satisfying the old absolute 1e-9 test.
 * Long-double intermediates resolve the small 1/Gamma^2 contribution. */
static long double PressureResidualW(long double W, long double D,
                                    long double a, long double B2,
                                    long double S, const double *u)
{
  long double invg=(D+sqrtl(D*D+4*a*W))/(2*W);
  long double v2=0;
  for(int n=0;n<3;n++) {
    long double v=((long double)u[MX1+n]+S*u[BX1+n]/W)/(W+B2);
    v2+=v*v;
  }
  return (v2-1)+invg*invg;
}

static int PressureFixIdealW(double *u, double *v)
{
  long double D=u[RHO],a=(long double)g_gamma/(g_gamma-1)*g_smallPressure;
  long double B2=0,S=0,lo,hi,mid,f;
  double out[NVAR]={0}, mapped[NVAR]={0},v2=0;
  int n,it;
  if (!isfinite(u[RHO]) || D<=0 || !isfinite((double)a) || a<=0) return 1;
  for(n=0;n<3;n++) {
    if(!isfinite(u[MX1+n]) || !isfinite(u[BX1+n])) return 1;
    B2+=(long double)u[BX1+n]*u[BX1+n];
    S+=(long double)u[MX1+n]*u[BX1+n];
  }
  lo=D+a;
  hi=fmaxl(2*lo, sqrtl((long double)u[MX1]*u[MX1]
                     +(long double)u[MX2]*u[MX2]
                     +(long double)u[MX3]*u[MX3])+B2+lo);
  for(it=0;it<128;it++) {
    f=PressureResidualW(hi,D,a,B2,S,u);
    if(!isfinite(f))return 1;
    if(f<=0)break;
    hi*=2;
  }
  if(it==128)return 2;
  for(it=0;it<256;it++) {
    mid=lo+(hi-lo)/2;
    if(mid==lo || mid==hi || hi-lo<=16*LDBL_EPSILON*hi)break;
    f=PressureResidualW(mid,D,a,B2,S,u);
    if(!isfinite(f))return 1;
    if(f>0)lo=mid;else hi=mid;
  }
  if(it==256)return 2;
  mid=lo+(hi-lo)/2;
  for(n=0;n<3;n++) {
    out[VX1+n]=(double)(((long double)u[MX1+n]+S*u[BX1+n]/mid)/(mid+B2));
    out[BX1+n]=u[BX1+n];
    v2+=out[VX1+n]*out[VX1+n];
  }
  if(!isfinite(v2) || v2>=1)return 1;
  out[RHO]=u[RHO]*sqrt(1-v2); out[PRS]=g_smallPressure;
  PrimToConsLoc(out,mapped);
  if(!isfinite(mapped[ENG]))return 1;
  /* Reject a poorly resolved ultra-relativistic inverse rather than silently
   * changing momentum. Normalization uses the momentum/enthalpy scale. */
  for(n=0;n<3;n++)
    if(fabsl((long double)mapped[MX1+n]-u[MX1+n])>1.e-7L*(mid+B2))return 2;
  for(n=0;n<8;n++)v[n]=out[n];
  u[ENG]=mapped[ENG];
#if ENTROPY_SWITCH
  u[ENTR]=g_smallPressure/sqrt(1-v2)/pow(out[RHO],g_gamma-1);
#endif
  return 0;
}
#endif

#define MAX_ITER 50
/* ********************************************************************* */
int RMHD_PressureFix(double *u, double *v)
/*!
 *
 * \return Error codes are:
 * - 0 = success
 * - 1 = v^2 > 1
 * - 2 = too many iterations
 *********************************************************************** */
{
#if EOS == IDEAL
  return PressureFixIdealW(u,v);
#else
  int    k, done=0;
  double v2,  f, W, W2, dW, S2_W2;
  double fmin, fmax, v2min, v2max;
  double lor, plor, scrh, w_1;

  double D  = u[RHO];
  double S  = u[MX1]*u[BX1] + u[MX2]*u[BX2] + u[MX3]*u[BX3];
  double m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
  double B2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];
  double S2 = S*S;

  double p = g_smallPressure; 

/* ------------------------------------------------
    Use secant algorithm
   ------------------------------------------------ */

/*
  v2max = 1.0-1.e-8;
  v2c = 0.95;
  fc  = VelocitySquareFunc(v2c, par);
  v2  = 0.96;
  for (k = 1; k < MAX_ITER; k++){
    f   = VelocitySquareFunc(v2, par);
print ("%d,  v2 = %12.6e  f = %12.6e\n",k,v2,f);
    if (done == 1) break;
    dW  = (v2 - v2c)/(f - fc)*f;
    v2c = v2; fc = f;
    v2 -= dW;
    v2 = MIN(v2max,v2);
    v2 = MAX(v2, 0.0);
    if (fabs(f) < 1.e-9) done = 1;
  }

  if (v2 >= 1.0 || k >= MAX_ITER) {
    FILE *fp;
    fp = fopen("pressure_fix.dat","w");
    for (v2 = 0.99; v2 < 1.0; v2 += 1.e-5){
      f = VelocitySquareFunc (v2, par);
      fprintf (fp, "%12.6e  %12.6e\n", v2, f);
    }
    fclose(fp);
    printLog ("! PressureFix: too many iter while fixing p , v^2 = %f\n", v2);
     
    return (1);
  }
*/

/* --------------------------------------------------------
   2. Use bisection algorithm to solve f(v^2) = 0
   -------------------------------------------------------- */

/* ----------------------------------------------
   2a. Set lower bound for v^2
   ---------------------------------------------- */

  v2min = v2 = 0.0;
  lor  = 1.0/sqrt(1.0 - v2);
  plor = p*lor;
  #if EOS == IDEAL
  W = (D + plor*g_gamma/(g_gamma - 1.0))*lor;
  #elif EOS == TAUB
  W = (2.5*plor + sqrt(2.25*plor*plor + D*D))*lor;
  #endif  
  W2  = W*W;

  f  = (S2*(2.0*W + B2) + m2*W2)/((W + B2)*(W + B2)*W2);
  f -= v2;
  fmin = f;

/* ----------------------------------------------
   2b. Set upper bound for v^2
   ---------------------------------------------- */

  v2max = v2 = 1.0-1.e-9;
  lor  = 1.0/sqrt(1.0 - v2);
  plor = p*lor;
  #if EOS == IDEAL
  W = (D + plor*g_gamma/(g_gamma - 1.0))*lor;
  #elif EOS == TAUB
  W = (2.5*plor + sqrt(2.25*plor*plor + D*D))*lor;
  #endif  
  W2  = W*W;

  f  = (S2*(2.0*W + B2) + m2*W2)/((W + B2)*(W + B2)*W2);
  f -= v2;
  fmax = f;

/* ----------------------------------------------
   2c. Start iterating
   ---------------------------------------------- */

  for (k = 1; k < MAX_ITER; k++){

    v2 = 0.5*(v2min + v2max);

    lor  = 1.0/sqrt(1.0 - v2);
    plor = p*lor;
    #if EOS == IDEAL
    W = (D + plor*g_gamma/(g_gamma - 1.0))*lor;
    #elif EOS == TAUB
    W = (2.5*plor + sqrt(2.25*plor*plor + D*D))*lor;
    #endif  
    W2  = W*W;

    f  = (S2*(2.0*W + B2) + m2*W2)/((W + B2)*(W + B2)*W2);
    f -= v2;

    if (f*fmin > 0.0){
      v2min = v2; fmin = f;
    }else{
      v2max = v2; fmax = f;
    }
    if (fabs(f) < 1.e-9) break;
  }

  if (v2 >= 1.0)      return 1;
  if (k  == MAX_ITER) return 2;
  
/* ----------------------------------------------
   3. Complete conversion
   ---------------------------------------------- */

  v[RHO] = D/lor;
  v[PRS] = p;

  w_1  = 1.0/(W + B2);
  scrh = S/W;
  v[VX1] = w_1*(u[MX1] + scrh*u[BX1]);
  v[VX2] = w_1*(u[MX2] + scrh*u[BX2]);
  v[VX3] = w_1*(u[MX3] + scrh*u[BX3]);

/* ----------------------------------------------
   4. For consistency, redefine energy, proper
      density and entropy
   ---------------------------------------------- */

  S2_W2   = S2/W2;
#if RMHD_REDUCED_ENERGY == YES
  u[ENG] = W - D - p + 0.5*(1.0 + v2)*B2 - 0.5*S2_W2;
#else
  u[ENG] = W - p + 0.5*(1.0 + v2)*B2 - 0.5*S2_W2;
#endif

/* -- Recompute entropy consistently -- */

#if ENTROPY_SWITCH
{
  double rho = v[RHO];
  double th  = v[PRS]/rho;  
  #if EOS == IDEAL
  u[ENTR] = p*lor/pow(rho,g_gamma-1);
  #elif EOS == TAUB
  u[ENTR] = p*lor/pow(rho,2.0/3.0)*(1.5*th + sqrt(2.25*th*th + 1.0));
  #endif
}
#endif

  return 0;  /* -- success -- */
#endif
}

#undef MAX_ITER
