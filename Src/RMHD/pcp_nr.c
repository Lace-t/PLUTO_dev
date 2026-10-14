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

#define PCP_MAX_ITER   500
#define PCP_LIMIT_ITER 48
#define PCP_TOL        (8.0*DBL_EPSILON)
#define PCP_ROOT_TOL   1.e-14
#define PCP_ABS_FLOOR  (8.0*DBL_EPSILON)

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
static int PCPCheck (const double *u)
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
      !PCPFinite(B2) || !PCPFinite(tau) || D <= PCP_ABS_FLOOR || E <= 0.0 ||
      m2 < 0.0 || B2 < 0.0) return 0;

  a = D*D + m2;
  if (!PCPFinite(a) || a < 0.0) return 0;
  scale = fmax(1.0, fabs(E));
  /* A small rounding-scale margin protects the q subtraction. */
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
  /* Evaluate the safety-set energy shift at rounding-error scale. */
  eps_E = fmax(PCP_ABS_FLOOR, PCP_TOL*fmax(1.0, fabs(E)));
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


int RMHD_PCP_Check (const double *u)
{
  return PCPCheck(u);
}

/* Scale only (D,m,E), keeping the face magnetic field untouched.  This is
   important for staggered CT: both states on a face retain the same B_n. */
int RMHD_PCP_Limit (const double *u0, double *u)
{
  int k, n;
  double theta, lo, hi, mid, epsD;
  double cand[NVAR], anchor[NVAR], stage[NVAR];

  /* Return values: -1 = failed (input unchanged), 0 = unchanged, 1 = limited. */
  if (u == NULL) return -1;
  if (RMHD_PCP_Check(u)) return 0;
  if (u0 == NULL || !RMHD_PCP_Check(u0)) return -1;
  for (n = 0; n < NVAR; n++) if (!PCPFinite(u[n])) return -1;

  /* Keep the face magnetic field fixed for CT.  The actual bisection anchor
     is therefore the center fluid state paired with that face field, which
     must itself be admissible before convex limiting is meaningful. */
  for (n = 0; n < NVAR; n++) anchor[n] = u0[n];
  anchor[BX1] = u[BX1];
  anchor[BX2] = u[BX2];
  anchor[BX3] = u[BX3];
  if (!RMHD_PCP_Check(anchor)) return -1;

  /* Adapted to CT: first enforce the density floor at fixed full energy,
     then enforce q and Psi together along one fixed-B fluid-state segment. */
  for (n = 0; n < NVAR; n++) cand[n] = u[n];
  epsD = nextafter(PCP_ABS_FLOOR, INFINITY);
  if (u[RHO] < epsD) {
    theta = (anchor[RHO]-epsD)/(anchor[RHO]-u[RHO]);
    if (!PCPFinite(theta)) return -1;
    theta = fmax(0.0, fmin(1.0, theta));
    cand[RHO] = epsD;
  }
#if RMHD_REDUCED_ENERGY == YES
  /* ENG stores E-D: changing only D must leave the full energy unchanged. */
  cand[ENG] = u[ENG] + (u[RHO] - cand[RHO]);
#endif
  memcpy(stage, cand, sizeof(stage));

  /* The former q and Psi stages scaled the same fields toward the same
     anchor.  Combine them into one search using the full checker, which
     rejects failed q states before evaluating Psi.  Keep the admissible
     endpoint; neither B nor the other unscaled fields change. */
  if (!RMHD_PCP_Check(cand)) {
    lo=0.0; hi=1.0;
    for (k=0; k<PCP_LIMIT_ITER; k++) {
      mid=0.5*(lo+hi);
      cand[RHO]=anchor[RHO]+mid*(stage[RHO]-anchor[RHO]);
      cand[MX1]=anchor[MX1]+mid*(stage[MX1]-anchor[MX1]);
      cand[MX2]=anchor[MX2]+mid*(stage[MX2]-anchor[MX2]);
      cand[MX3]=anchor[MX3]+mid*(stage[MX3]-anchor[MX3]);
      cand[ENG]=anchor[ENG]+mid*(stage[ENG]-anchor[ENG]);
      if (RMHD_PCP_Check(cand)) lo=mid; else hi=mid;
    }
    theta=lo;
    cand[RHO]=anchor[RHO]+theta*(stage[RHO]-anchor[RHO]);
    cand[MX1]=anchor[MX1]+theta*(stage[MX1]-anchor[MX1]);
    cand[MX2]=anchor[MX2]+theta*(stage[MX2]-anchor[MX2]);
    cand[MX3]=anchor[MX3]+theta*(stage[MX3]-anchor[MX3]);
    cand[ENG]=anchor[ENG]+theta*(stage[ENG]-anchor[ENG]);
  }
  if (!RMHD_PCP_Check(cand)) return -1;
  memcpy(u, cand, sizeof(cand));
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

#if (PCP_NR == YES) || (PCP_FOFC == YES)
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
  double invw2, invlor, lor, phi_a, vterm, v2;

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
  /* Resolve small velocities as a sum of nonnegative terms.  Computing
     v2 = 1 - invw2 instead can make a resting state spuriously negative.
     At high velocity, retain the factored inverse Lorentz expression to
     avoid cancellation in 1 - v2. */
  v2 = (beta1/xi)/xi + (beta2/eta)/eta;
  if (PCPFinite(v2) && v2 >= 0.0 && v2 <= 0.5) {
    invw2 = 1.0 - v2;
  } else {
    invw2 = ((xi + alpha2)*(eta + sm))/eta2
          - beta1*B2*(2.0*xi + B2)/(eta2*xi*xi);
    v2 = 1.0 - invw2;
  }
  if (!PCPFinite(invw2) || invw2 <= 0.0) return 0;
  invlor = sqrt(invw2);
  lor = 1.0/invlor;
  if (!PCPFinite(lor) || lor < 1.0) return 0;
  out->lor = lor;
  out->v2 = v2;
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
  /* Cai et al. (2.11)-(2.12).  cbrt is the real cube root, including
     negative arguments; the trigonometric branch avoids complex roots. */
  double a = B2 - E;
  double a0 = -0.5*(B2*D*D + tau*tau);
  double delta = 27.0*a0 + 4.0*a*a*a;

  if (delta > 0.0) {
    double arg = 1.0 + 13.5*a0/(a*a*a);
    double theta = acos(fmax(-1.0, fmin(1.0, arg)));
    return -a/3.0*(1.0 - 2.0*cos(theta/3.0 - acos(-1.0)/3.0));
  } else {
    double X1 = a*a*a + 13.5*a0;
    double X2 = 1.5*sqrt(fmax(0.0, 3.0*a0*delta));
    return -(a + cbrt(X1 + X2) + cbrt(X1 - X2))/3.0;
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
  int n, nosc = 0;
  double D, E, m2, B2, tau, bscale, bpar, beta1;
  double sm, g0;
  double xi, xnew, alpha1, phi, previous_f = 0.0, step;
  double v2, lor_check, rel, ucheck[NVAR];
  PCPEval ev;

  /* Stage/limiter admission is separate from recovery. Allow NR to try
     states inside the checker's safety layer; validate the recovered state. */
  if (u == NULL || v == NULL) return 1;
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
  /* Paper (2.9)-(2.10): use xi_d when F(xi_d) <= 0, otherwise xi_c.
     The proven safe initial value replaces the former bracket search. */
  alpha1 = B2 - E;
  phi = sqrt(alpha1*alpha1 + 3.0*(E*E - D*D - m2));
  xi = (phi - 2.0*alpha1)/3.0;
  if (!PCPEvaluate(xi, D, E, m2, B2, tau, beta1, sm, g0, &ev)) return 1;
  if (ev.f > 0.0) {
    xi = PCPCubicRoot(E, D, B2, tau);
    if (!PCPEvaluate(xi, D, E, m2, B2, tau, beta1, sm, g0, &ev)) return 1;
  }
  for (k = 0; k < PCP_MAX_ITER; k++) {
    if (ev.p <= 0.0 || ev.rho <= 0.0 || ev.df <= 0.0) return 1;
    /* Algorithm A.1: unmodified Newton update, absolute step tolerance,
       and stop after more than three residual sign changes (round-off
       oscillation).  Compare signs without an underflow-prone product. */
    xnew = xi - ev.f/ev.df;
    if (!PCPFinite(xnew) || xnew <= 0.0) return 1;
    if ((previous_f < 0.0 && ev.f > 0.0) ||
        (previous_f > 0.0 && ev.f < 0.0)) nosc++;
    previous_f = ev.f;
    step = fabs(xnew - xi);
    xi = xnew;
    /* Re-evaluate at the OUTPUT iterate: do not reuse rho,p at xi_n
       together with velocities calculated at xi_{n+1}. */
    if (!PCPEvaluate(xi, D, E, m2, B2, tau, beta1, sm, g0, &ev)) return 1;
    if (step <= PCP_ROOT_TOL || nosc > 3) break;
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
  if (!PCPFinite(lor_check)) return 1;

  /* Paper (2.13): calculate the final thermodynamics from the recovered
     velocity and the final xi, so D = rho W and p = P(rho,h) use the
     same Lorentz factor as the returned primitive vector.  Do not require
     agreement with the internal Wcal to a fixed relative tolerance:
     reconstructing W from explicit velocities amplifies round-off by O(W^2).
     Validate the returned state with the physical and conservative round-trip
     checks below instead. */
  v[RHO] = D/lor_check;
  v[PRS] = g0*(xi*(1.0 - v2) - v[RHO]);
  if (!PCPFinite(v[RHO]) || !PCPFinite(v[PRS]) ||
      v[RHO] <= 0.0 || v[PRS] <= 0.0) return 1;

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
