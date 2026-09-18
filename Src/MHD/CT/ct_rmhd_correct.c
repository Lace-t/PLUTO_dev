#include "pluto.h"

#if PHYSICS == RMHD
/* Transactional CA2'-style correction. All recovery, including the native
 * PressureFix policy, takes place on private copies with the original FV B.
 * Rebuilding with CT B preserves the recovered fluid primitive state.
 * D and passive conserved scalars are kept bit-for-bit. Unlike a conservative
 * flux update this intentionally changes momentum/energy; delta reports the
 * TOTAL change, including recovery and round-trip projection, for accounting.
 * Failure changes only flags. The caller still installs CT B, and the normal
 * fail-safe path rejects the flagged stage and rolls back the entire step. */
int CT_RMHD_CorrectCell(double *u, double *v, const double *b,
                       uint16_t *flag, double *delta)
{
  double work[NVAR], prim[NVAR], mapped[NVAR], v2;
  uint16_t trial = *flag;
  int n;
  for (n=0; n<4; n++) delta[n] = 0.0;
  NVAR_LOOP(n) {
    work[n] = u[n]; prim[n] = v[n]; mapped[n] = u[n];
    if (!isfinite(u[n])) goto fail;
  }
  for (n=0; n<3; n++) if (!isfinite(b[n])) goto fail;
  if (u[RHO] <= 0.0 || (trial & FLAG_CONS2PRIM_FAIL)) goto fail;
  if (ConsToPrimLoc(work, prim, &trial)) goto fail;
  /* The energy solver may succeed with p below the case's existing floor.
   * Apply the SAME pressure-floor recovery at fixed D and momentum before
   * switching B, rather than rejecting a positive but sub-floor state or
   * raising p at fixed velocity (which would add uncontrolled momentum). */
  if (prim[PRS] < g_smallPressure && RMHD_PressureFix(work, prim)) {
    trial |= FLAG_PRESSURE_FIX_FAIL;
    goto fail;
  }
  NVAR_LOOP(n) if (!isfinite(prim[n])) goto fail;
  v2 = DOT_PRODUCT(prim+VX1, prim+VX1);
  if (!(v2 < 1.0) || prim[RHO] <= 0.0 || prim[PRS] <= 0.0) goto fail;

  /* Reconcile density with the explicit recovered velocity, even when the
   * native recovery's intermediate Lorentz factor differs by roundoff. */
  prim[RHO] = u[RHO]*sqrt(1.0-v2);
  for (n=0; n<3; n++) prim[BX1+n] = b[n];
  PrimToConsLoc(prim, mapped);
  mapped[RHO] = u[RHO];
  NVAR_LOOP(n) if (!isfinite(mapped[n])) goto fail;
  if (!RMHD_PCP_Check(mapped)) goto fail;

  for (n=0; n<3; n++) {
    delta[n] = mapped[MX1+n]-u[MX1+n];
    u[MX1+n] = mapped[MX1+n]; u[BX1+n] = b[n];
  }
  delta[3] = mapped[ENG]-u[ENG];
  u[ENG] = mapped[ENG];
  NVAR_LOOP(n) v[n] = prim[n];
  *flag = trial;
  return 0;
fail:
  *flag = trial | FLAG_CONS2PRIM_FAIL;
  return 1;
}
#endif
