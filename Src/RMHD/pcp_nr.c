/* ///////////////////////////////////////////////////////////////////// */
/*
   Admissibility checks, conservative-state limiting, and the ideal-EOS
   physical-constraint-preserving primitive recovery for RMHD.

   The Wu admissibility and reconstruction-limiter routines are independent
   of the normal Mignone-McKinney recovery path and can also support stage
   limiting.  Only the Cai primitive-recovery solver is selected by PCP_NR.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include <float.h>

#define PCP_MAX_ITER   60
#define PCP_LIMIT_ITER 48
#define PCP_TOL        (64.0*DBL_EPSILON)
#define PCP_ROOT_TOL   1.e-10
#define PCP_ABS_FLOOR  1.e-14

static double PCPFullEnergy (const double *u)
{
#if RMHD_REDUCED_ENERGY == YES
  return u[ENG] + u[RHO];
#else
  return u[ENG];
#endif
}

static int PCPFinite (double x)
{
  return isfinite(x);
}

/* Wu et al. G_0 admissibility test.  The small relative margins prevent a
   state on a floating-point boundary from being passed to the NR solver. */
int RMHD_PCP_Check (const double *u)
{
  double D, E, m2, B2, tau, phi, rad, lhs, rhs, psi_d, psi_scale;
  long double psi;
  double scale, eps_q, eps_E, a, Etest, q_d;

  if (u == NULL) return 0;
  D   = u[RHO];
  E   = PCPFullEnergy(u);
  m2  = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
  B2  = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];
  tau = u[MX1]*u[BX1] + u[MX2]*u[BX2] + u[MX3]*u[BX3];

  if (!PCPFinite(D) || !PCPFinite(E) || !PCPFinite(m2) ||
      !PCPFinite(B2) || !PCPFinite(tau) || D <= 0.0 || E <= 0.0 ||
      m2 < 0.0 || B2 < 0.0) return 0;

  a = D*D + m2;
  if (!PCPFinite(a) || a < 0.0) return 0;
  scale = fmax(1.0, fabs(E));
  eps_q = fmax(PCP_ABS_FLOOR, PCP_TOL*scale);
  /* Most states are well inside q > 0.  Avoid an unconditional long-double
     square root in this hot checker; use extended precision only when the
     double result is close enough to the boundary for cancellation to
     matter. */
  q_d = E - sqrt(a);
  if (PCPFinite(q_d)) {
    if (q_d > eps_q) {
      /* Fast path: q has a clear positive margin. */
    } else if (q_d < -eps_q) {
      return 0;
    } else {
      long double q = (long double)E - sqrtl((long double)a);
      if (!isfinite((double)q) || q <= (long double)eps_q) return 0;
    }
  } else {
    return 0;
  }

  a = E*E - D*D - m2;
  if (!PCPFinite(a) || a < 0.0) return 0;
  /* Evaluate Psi on E - epsilon_E, as in the definition of G_epsilon.
     This avoids comparing Psi, which has a different physical scale, to an
     energy-sized tolerance. */
  /* The Psi margin is applied to the fluid-energy scale.  Scaling it with
     the total energy would reject cold, strongly magnetized states merely
     because their magnetic energy is large. */
  eps_E = fmax(PCP_ABS_FLOOR, PCP_TOL*fmax(1.0, D));
  Etest = E - eps_E;
  if (!PCPFinite(Etest) || Etest <= 0.0) return 0;
  a = B2 - Etest;
  phi = sqrt(a*a + 3.0*(Etest*Etest - D*D - m2));
  rad = phi + a;
  if (!PCPFinite(phi) || !PCPFinite(rad) || rad < 0.0) return 0;
  /* Most states are well inside the admissible set.  Resolve those with a
     double-precision fast path; reserve long double for the narrow region
     where the Psi subtraction is vulnerable to cancellation. */
  lhs = (phi - 2.0*a)*sqrt(rad);
  rhs = sqrt(fmax(0.0, 13.5*(D*D*B2 + tau*tau)));
  if (PCPFinite(lhs) && PCPFinite(rhs)) {
    psi_d = lhs - rhs;
    psi_scale = fmax(1.0, rhs);
    if (psi_d > PCP_TOL*psi_scale) return 1;
    if (psi_d < -PCP_TOL*psi_scale) return 0;
  }
  /* Psi is a difference of two nearly equal large terms in magnetically
     dominated states.  Evaluate only this checker expression in extended
     precision; the hot NR path remains double precision. */
  {
    long double ldD = (long double)D, ldE = (long double)Etest;
    long double ldB2 = (long double)B2, ldm2 = (long double)m2;
    long double ldtau = (long double)tau;
    long double lda = ldB2 - ldE;
    long double ldrad = lda*lda + 3.0L*(ldE*ldE - ldD*ldD - ldm2);
    long double ldphi, ldradp;
    if (!(isfinite((double)ldrad) && ldrad >= 0.0L)) return 0;
    ldphi = sqrtl(ldrad);
    ldradp = ldphi + lda;
    if (!(isfinite((double)ldphi) && isfinite((double)ldradp) && ldradp >= 0.0L)) return 0;
    psi = (ldphi - 2.0L*lda)*sqrtl(ldradp)
        - sqrtl(fmaxl(0.0L, 13.5L*(ldD*ldD*ldB2 + ldtau*ldtau)));
  }
  return isfinite((double)psi) && psi >= 0.0L;
}

/* Scale only (D,m,E), keeping the face magnetic field untouched.  This is
   important for staggered CT: both states on a face retain the same B_n. */
int RMHD_PCP_Limit (const double *u0, double *u)
{
  int k, n;
  double theta, lo, hi, mid;
  double cand[NVAR], anchor[NVAR];

  /* Return values: -1 = no admissible anchor, 0 = unchanged, 1 = limited. */
  if (u == NULL) return -1;
  if (RMHD_PCP_Check(u)) return 0;
  if (u0 == NULL || !RMHD_PCP_Check(u0)) return -1;

  /* Keep the face magnetic field fixed for CT.  The actual bisection anchor
     is therefore the center fluid state paired with that face field, which
     must itself be admissible before convex limiting is meaningful. */
  for (n = 0; n < NVAR; n++) anchor[n] = u0[n];
  anchor[BX1] = u[BX1];
  anchor[BX2] = u[BX2];
  anchor[BX3] = u[BX3];
  if (!RMHD_PCP_Check(anchor)) return -1;

  //左右界面的重构状态相对于单元中心参考状态的偏差
  for (n = 0; n < NVAR; n++) cand[n] = u[n];
  lo = 0.0;
  hi = 1.0;
  for (k = 0; k < PCP_LIMIT_ITER; k++) {
    mid = 0.5*(lo + hi);
    for (n = 0; n < NVAR; n++) cand[n] = u[n];
    cand[RHO] = anchor[RHO] + mid*(u[RHO] - anchor[RHO]);
    cand[MX1] = anchor[MX1] + mid*(u[MX1] - anchor[MX1]);
    cand[MX2] = anchor[MX2] + mid*(u[MX2] - anchor[MX2]);
    cand[MX3] = anchor[MX3] + mid*(u[MX3] - anchor[MX3]);
    cand[ENG] = anchor[ENG] + mid*(u[ENG] - anchor[ENG]);
    if (RMHD_PCP_Check(cand)) lo = mid;
    else                      hi = mid;
  }
  theta = lo;
  u[RHO] = anchor[RHO] + theta*(u[RHO] - anchor[RHO]);
  u[MX1] = anchor[MX1] + theta*(u[MX1] - anchor[MX1]);
  u[MX2] = anchor[MX2] + theta*(u[MX2] - anchor[MX2]);
  u[MX3] = anchor[MX3] + theta*(u[MX3] - anchor[MX3]);
  u[ENG] = anchor[ENG] + theta*(u[ENG] - anchor[ENG]);
  if (!RMHD_PCP_Check(u)) return -1;
  return 1;
}

#if PCP_NR == YES
/* Produce a finite primitive state after an unsuccessful recovery.  The
   density uses the supplied lab-frame density when possible, so mass is not
   needlessly discarded; zero velocity makes it a valid rest-frame state. */
void RMHD_PCP_SetFallback (const double *u, double *v)
{
  int n;
  double rho;

  for (n = 0; n < NVAR; n++) v[n] = 0.0;
  rho = (u != NULL && isfinite(u[RHO])) ? u[RHO] : g_smallDensity;
  v[RHO] = fmax(g_smallDensity, rho);
  v[PRS] = g_smallPressure;
  if (u != NULL) {
    v[BX1] = isfinite(u[BX1]) ? u[BX1] : 0.0;
    v[BX2] = isfinite(u[BX2]) ? u[BX2] : 0.0;
    v[BX3] = isfinite(u[BX3]) ? u[BX3] : 0.0;
  }
}

#endif

#if (PCP_NR == YES) || (PCP_STAGE_LIMITING == YES)
/* Repair a trial reconstructed state and keep its primitive/conservative
   forms synchronized.  The caller discards failed trials, leaving PLUTO's
   original reconstruction for the normal update/C2P/failsafe paths. */
static int PCPRepairState (const double *center_u, double *face_u,
                           double *face_v)
{
  uint16_t local_flag = 0;
  int status, n;

  status = RMHD_PCP_Limit(center_u, face_u);
  if (status < 0) return -1;
  if (status == 0) return 0;

  /* The limiter changes conservative variables.  Recompute primitives before
     the Riemann solver sees the state; otherwise v and u would describe two
     different states and create a non-physical flux. */
  if (ConsToPrimLoc(face_u, face_v, &local_flag) != 0
      || (local_flag & (FLAG_CONS2PRIM_FAIL | FLAG_NEGATIVE_PRESSURE))
      || !RMHD_PCP_Check(face_u)) return -1;
  /* PressureFix is allowed in the ordinary cell C2P path, but is not a
     successful convex repair: it replaces the energy at fixed floor pressure.
     Reject it here even when ConsToPrimLoc reports successful recovery. */
  for (n = 0; n < NVAR; n++) {
    if (!isfinite(face_u[n]) || !isfinite(face_v[n])) return -1;
  }
  if (face_v[RHO] <= 0.0 || face_v[PRS] <= 0.0) return -1;
  return 1;
}

/* Limit both states at each interface around their respective cell-center
   states while retaining the magnetic components supplied by CT.  The
   operation is atomic per interface: if either side cannot be repaired, both
   trials are discarded and the original PLUTO reconstruction is left intact
   for the ordinary update/C2P/failsafe machinery. */
void RMHD_PCP_LimitStates (const Sweep *sweep, int beg, int end)
{
  int i, statusL, statusR, checkL, checkR;
  double uL[NVAR], uR[NVAR], vL[NVAR], vR[NVAR];
  double centerL[NVAR], centerR[NVAR];
  const State *stateC = &(sweep->stateC);
  const State *stateL = &(sweep->stateL);
  const State *stateR = &(sweep->stateR);

  /* beg/end are reconstructed CELL indices.  The paired interfaces are
     beg...end-1: stateR[i] comes from cell i+1 (vm = stateR.v - 1).
     Neither stateR[end] nor stateL[beg-1] is reconstructed by this call. */
  for (i = beg; i < end; i++) {
    checkL = RMHD_PCP_Check(stateL->u[i]);
    checkR = RMHD_PCP_Check(stateR->u[i]);
    if (checkL && checkR) continue;

    memcpy(uL, stateL->u[i], sizeof(uL));
    memcpy(uR, stateR->u[i], sizeof(uR));
    memcpy(vL, stateL->v[i], sizeof(vL));
    memcpy(vR, stateR->v[i], sizeof(vR));

    /* RK UpdateStage fills stateC.v, not stateC.u.  Build anchors on demand
       from the current stage's center primitives, never stale sweep storage.
       All States callers have restored 3-velocities before this hook. */
    statusL = statusR = 0;
    if (!checkL) {
      PrimToConsLoc(stateC->v[i], centerL);
      statusL = PCPRepairState(centerL, uL, vL);
    }
    if (!checkR) {
      PrimToConsLoc(stateC->v[i+1], centerR);
      statusR = PCPRepairState(centerR, uR, vR);
    }

    if (statusL < 0 || statusR < 0) {
      /* Trial arrays are discarded; original u/v were never modified. */
      if (sweep->flag != NULL) {
        sweep->flag[i] |= FLAG_PCP_LIMIT_FAIL;
        sweep->flag[i + 1] |= FLAG_PCP_LIMIT_FAIL;
      }
    } else {
      if (statusL > 0) {
        memcpy(stateL->u[i], uL, sizeof(uL));
        memcpy(stateL->v[i], vL, sizeof(vL));
      }
      if (statusR > 0) {
        memcpy(stateR->u[i], uR, sizeof(uR));
        memcpy(stateR->v[i], vR, sizeof(vR));
      }
    }
  }
}
#endif


// PCP_NR

typedef struct {
  double f, df, lor, rho, p, v2;
} PCPEval;

static int PCPEvaluate (double xi, double D, double E, double m2,
                        double B2, double tau, double beta1,
                        double sm, double g0, PCPEval *out)
{
  double alpha1, alpha2, beta2, eta, eta2;
  double invw2, invlor, lor, phi_a, vterm;

  if (xi <= 0.0 || !PCPFinite(xi)) return 0;
  eta = xi + B2;
  if (eta <= 0.0 || !PCPFinite(eta)) return 0;
  eta2 = eta*eta;
  if (!PCPFinite(eta2) || eta2 <= 0.0) return 0;
  alpha1 = B2 - E;
  alpha2 = B2 - sm;
  /* beta1 = (m dot Bhat)^2 is computed once with a scaled magnetic norm.
     Do not reconstruct it as (m dot B)^2/B^2 here: a subnormal B^2 loses
     relative accuracy and spuriously violates beta1 <= m^2. */
  beta2 = m2 - beta1;
  /* Factor the small difference of inverse squares.  Direct subtraction
     loses relative accuracy when xi >> B^2, which is common in magnetized
     high-Lorentz states. */
  invw2 = ((xi + alpha2)*(eta + sm))/eta2
        - beta1*B2*(2.0*xi + B2)/(eta2*xi*xi);
  if (!PCPFinite(invw2) || invw2 <= 0.0) return 0;
  invlor = sqrt(invw2);
  lor = 1.0/invlor;
  if (!PCPFinite(lor) || lor < 1.0) return 0;
  out->lor = lor;
  out->v2 = 1.0 - invw2;
  if (!PCPFinite(out->v2) || out->v2 < 0.0 || out->v2 >= 1.0) return 0;

  vterm = xi*invw2 - D*invlor;
  out->rho = D/lor;
  out->p = g0*vterm;
  phi_a = -(beta1/(xi*xi*xi) + beta2/(eta*eta*eta));
  out->f = xi - g0*vterm - 0.5*(B2*invw2 + tau*tau/(xi*xi)) + alpha1;
  out->df = 1.0 + B2*phi_a + tau*tau/(xi*xi*xi)
          - g0*(invw2 - 2.0*xi*phi_a + D*phi_a/invlor);
  return PCPFinite(out->f) && PCPFinite(out->df) && PCPFinite(out->rho)
      && PCPFinite(out->p);
}

static double PCPCubicRoot (double E, double D, double B2, double tau)
{
  double a = B2 - E;
  double c = -0.5*(B2*D*D + tau*tau);
  double p = -a*a/3.0;
  double q = 2.0*a*a*a/27.0 + c;
  double disc = 0.25*q*q + p*p*p/27.0;
  double shift = -a/3.0;

  if (disc >= 0.0) {
    double s = sqrt(disc);
    return cbrt(-0.5*q + s) + cbrt(-0.5*q - s) + shift;
  } else {
    double r = 2.0*sqrt(-p/3.0);
    double arg = (3.0*q/(2.0*p))*sqrt(-3.0/p);
    double ang;
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;
    ang = acos(arg)/3.0;
    return r*cos(ang) + shift;
  }
}

/* Cai et al. one-dimensional PCP Newton recovery for an ideal EOS. */
int RMHD_PCP_NR (const double *u, double *v)
{
#if EOS != IDEAL
  (void)u; (void)v;
  return 1;
#else
  int k;
  int n, j;
  double D, E, m2, B2, tau, base, bad, good, bscale, bpar, beta1;
  double sm, g0;
  double lo, hi, flo, fhi, xi, xnew, xcard;
  double v2, lor_check, rel, ucheck[NVAR];
  PCPEval ev, ev_good;

  if (u == NULL || v == NULL || !RMHD_PCP_Check(u)) return 1;
  D = u[RHO];
  E = PCPFullEnergy(u);
  m2 = u[MX1]*u[MX1] + u[MX2]*u[MX2] + u[MX3]*u[MX3];
  sm = sqrt(fmax(0.0, m2));
  g0 = (g_gamma - 1.0)/g_gamma;
  B2 = u[BX1]*u[BX1] + u[BX2]*u[BX2] + u[BX3]*u[BX3];
  tau = u[MX1]*u[BX1] + u[MX2]*u[BX2] + u[MX3]*u[BX3];
  /* Normalize by the largest component BEFORE forming the norm.  Even
     hypot(B) has poor relative accuracy when B itself is subnormal.
     The scaled direction remains O(1) down to the smallest representable
     field, retaining the hydrodynamic limit without clipping/changing B. */
  bscale = fmax(fabs(u[BX1]), fmax(fabs(u[BX2]), fabs(u[BX3])));
  beta1 = 0.0;
  if (bscale > 0.0) {
    double b1 = u[BX1]/bscale, b2 = u[BX2]/bscale, b3 = u[BX3]/bscale;
    double bnorm = sqrt(b1*b1 + b2*b2 + b3*b3);
    bpar = (u[MX1]*b1 + u[MX2]*b2 + u[MX3]*b3)/bnorm;
    beta1 = bpar*bpar;
    if (!PCPFinite(beta1) || beta1 > m2*(1.0 + 128.0*DBL_EPSILON)) return 1;
    beta1 = fmin(m2, beta1);
  }
  /* xi = rho h W^2 is always >= D.  Starting at |m| is not a valid
     lower bound: magnetic stresses can make |m| much larger than xi. */
  base = fmax(D, PCP_ABS_FLOOR);
  bad = base;
  good = base*(1.0 + 128.0*DBL_EPSILON) + PCP_ABS_FLOOR;
  for (k = 0; k < 64; k++) {
    if (PCPEvaluate(good, D, E, m2, B2, tau, beta1, sm, g0, &ev) && ev.p > 0.0) {
      ev_good = ev;
      break;
    }
    bad = good;
    good *= 2.0;
    if (!PCPFinite(good)) return 1;
  }
  if (k == 64) return 1;
  /* Refine the lower edge of the physical interval instead of jumping from
     max(D,|m|) directly to a factor-of-two larger xi. */
  for (j = 0; j < 64; j++) {
    double mid = 0.5*(bad + good);
    if (PCPEvaluate(mid, D, E, m2, B2, tau, beta1, sm, g0, &ev) && ev.p > 0.0) {
      good = mid;
      ev_good = ev;
    } else {
      bad = mid;
    }
    /* Once the endpoints are adjacent in floating point, another
       evaluation cannot refine the physical interval. */
    if (good - bad <= 8.0*DBL_EPSILON*fmax(1.0, fabs(good))) break;
  }
  lo = good;
  ev = ev_good;
  if (!PCPFinite(ev.f) || !PCPFinite(ev.df) || ev.p <= 0.0) return 1;
  flo = ev.f;
  /* The Cardano estimate is used both to enlarge the bracket and as the
     Newton starting point.  Compute it once per recovery. */
  xcard = PCPCubicRoot(E, D, B2, tau);
  hi = fmax(lo*2.0, xcard);
  if (!PCPFinite(hi) || hi <= lo) hi = lo*2.0;
  for (k = 0; k < 32; k++) {
    if (!PCPEvaluate(hi, D, E, m2, B2, tau, beta1, sm, g0, &ev)) return 1;
    fhi = ev.f;
    if (fhi > 0.0) break;
    hi *= 2.0;
    if (!PCPFinite(hi)) return 1;
  }
  if (k == 32 || flo > PCP_ROOT_TOL*fmax(1.0, fabs(lo))) return 1;

  xi = xcard;
  if (!PCPFinite(xi) || xi <= lo || xi >= hi) xi = 0.5*(lo + hi);
  for (k = 0; k < PCP_MAX_ITER; k++) {
    if (!PCPEvaluate(xi, D, E, m2, B2, tau, beta1, sm, g0, &ev)) return 1;
    if (fabs(ev.f) <= PCP_ROOT_TOL*fmax(1.0, fabs(xi))) break;
    if (ev.f < 0.0) lo = xi; else hi = xi;
    xnew = (fabs(ev.df) > PCP_TOL) ? xi - ev.f/ev.df : NAN;
    if (!PCPFinite(xnew) || xnew <= lo || xnew >= hi) xnew = 0.5*(lo + hi);
    xi = xnew;
  }
  if (k == PCP_MAX_ITER) return 1;
  if (ev.p <= 0.0 || ev.rho <= 0.0 || ev.v2 >= 1.0) return 1;

  for (n = 0; n < NVAR; n++) v[n] = 0.0;
  v[RHO] = ev.rho;
  v[PRS] = ev.p;
  {
    double den = xi + B2;
    double fac = tau/xi;
    v[VX1] = (u[MX1] + fac*u[BX1])/den;
    v[VX2] = (u[MX2] + fac*u[BX2])/den;
    v[VX3] = (u[MX3] + fac*u[BX3])/den;
  }
  v[BX1] = u[BX1];
  v[BX2] = u[BX2];
  v[BX3] = u[BX3];
  v2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];
  if (!PCPFinite(v[VX1]) || !PCPFinite(v[VX2]) || !PCPFinite(v[VX3])
      || !PCPFinite(v2) || v2 >= 1.0) return 1;
  lor_check = 1.0/sqrt(1.0 - v2);
  rel = fabs(lor_check - ev.lor)/fmax(1.0, lor_check);
  /* The explicit velocity formula and the stable Wcal expression can differ
     by several ulps after cancellation in magnetically dominated states.
     Keep this tighter than the conservative round-trip check, but do not
     require bitwise agreement. */
  if (!PCPFinite(lor_check) || rel > 1.e-10) return 1;

  /* Verify that the recovered primitive state maps back to the supplied
     conservative state.  This catches cancellation in extreme magnetized
     or ultra-relativistic states before they enter the flux calculation. */
  PrimToConsLoc(v, ucheck);
  rel = fmax(1.0, fabs(E));
  if (fabs(PCPFullEnergy(ucheck) - E) > 1.e-8*rel) return 1;
  rel = fmax(1.0, fabs(D));
  if (fabs(ucheck[RHO] - D) > 1.e-8*rel) return 1;
  for (n = MX1; n <= MX3; n++) {
    rel = fmax(1.0, fabs(u[n]));
    if (fabs(ucheck[n] - u[n]) > 1.e-8*rel) return 1;
  }
  return 0;
#endif
}

#undef PCP_MAX_ITER
#undef PCP_LIMIT_ITER
#undef PCP_TOL
#undef PCP_ROOT_TOL
#undef PCP_ABS_FLOOR
