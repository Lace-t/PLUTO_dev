/* Conservative, CT-compatible a-posteriori PCP stage limiting.
 *
 * Both candidates start from the SAME FE input. A single MPI-global theta
 * blends every conservative component and every staggered face. Thus shared
 * flux/edge corrections cancel, and CT divergence is preserved to roundoff.
 * The LF+CT low-order state is CHECKED, never presumed PCP from Wu's theorem.
 * If it fails, FailSafe restores the entire step and first retries with
 * local order reduction at the original dt, then with smaller time steps.
 * See robust/PCP_NR_PLUTO_Integration.md, section 4.9.
 * Xufan Hu 2026.9
 */
#include "pluto.h"

int g_pcp_low_order = 0;

/* Construct the fluid states used by a failsafe retry after States() has
 * allocated its ordinary sweep storage.  Normal B is copied once from the
 * staggered face shared by both sides of each interface.  This is the
 * compatible first-order FV input for the LF flux and CT edge EMF; it never
 * replaces a face field independently. */
void RMHD_FailSafe_StageStates(const Sweep *s, int beg, int end)
{
  int i, n;
  for (i = beg; i <= end; i++) {
    NVAR_LOOP(n) {
      s->stateL.v[i][n] = s->stateC.v[i][n];
      s->stateR.v[i][n] = s->stateC.v[i+1][n];
    }
#ifdef STAGGERED_MHD
    s->stateL.v[i][BXn] = s->stateR.v[i][BXn] = s->Bn[i];
#endif
    PrimToConsLoc(s->stateL.v[i], s->stateL.u[i]);
    PrimToConsLoc(s->stateR.v[i], s->stateR.u[i]);
  }
}

/* Check the actual RK average after its RK combination and CT averaging,
 * before the mapper can change U. No primitive recovery or state repair is
 * performed here. All ranks must make the same decision before returning
 * from AdvanceStep, since the following C2P contains MPI collectives. */
int RMHD_PCP_CheckStage(Data *d, Grid *grid)
{
  int i,j,k,n, local_cells = 0, deferred_cells = 0;
  int counts[2];
  (void)grid;
  DOM_LOOP(k,j,i) {
    double *u = d->Uc[k][j][i];
    int full_bad = !RMHD_PCP_Check(u);
    int cell_bad = full_bad;
#if RMHD_FAILSAFE_PCP_STRICT == NO
    if (full_bad) {
      double E = u[ENG];
#if RMHD_REDUCED_ENERGY == YES
      E += u[RHO];
#endif
      /* Necessary conditions with no arbitrary absolute energy floor.
         A positive-energy state failing only the remaining full check is
         left to the existing mapper/PressureFix policy, explicitly counted
         below. This compatibility mode is NOT a full PCP acceptance gate. */
      cell_bad = !isfinite(E) || u[RHO] <= 0.0 || u[ENG] <= 0.0
              || !(E > hypot(u[RHO],hypot(u[MX1],hypot(u[MX2],u[MX3]))));
    }
#endif
    NVAR_LOOP(n) if (!isfinite(u[n])) cell_bad = 1;
#ifdef STAGGERED_MHD
    /* Both faces bounding the active cell must be finite. */
    DIM_EXPAND(
      if (!isfinite(d->Vs[IDIR][k][j][i]) ||
          !isfinite(d->Vs[IDIR][k][j][i-1])) cell_bad = 1;,
      if (!isfinite(d->Vs[JDIR][k][j][i]) ||
          !isfinite(d->Vs[JDIR][k][j-1][i])) cell_bad = 1;,
      if (!isfinite(d->Vs[KDIR][k][j][i]) ||
          !isfinite(d->Vs[KDIR][k-1][j][i])) cell_bad = 1;
    )
#endif
    if (cell_bad) {
      d->flag[k][j][i] |= FLAG_CONS2PRIM_FAIL;
      if (u[RHO] <= 0.0) d->flag[k][j][i] |= FLAG_NEGATIVE_DENSITY;
      if (u[ENG] <= 0.0) d->flag[k][j][i] |= FLAG_NEGATIVE_ENERGY;
      local_cells++;
    } else if (full_bad) deferred_cells++;
  }
  if (local_cells) printLog("RMHD_STAGE_REJECT: step=%ld stage=%d rank=%d "
                   "cells=%d dt=%.4g pre_c2p=YES\n",
                   g_stepNumber,g_intStage,prank,local_cells,g_dt);
  counts[0] = local_cells;
  counts[1] = deferred_cells;
#ifdef PARALLEL
  {
    int local_counts[2] = {local_cells,deferred_cells};
    MPI_Allreduce(local_counts,counts,2,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
#endif
  if (prank == 0 && (counts[0] || counts[1])){
    WARNING(printLog("RMHD_STAGE_CHECK: step=%ld stage=%d rejected=%d "
             "deferred_to_c2p=%d strict=%d\n",
             g_stepNumber,g_intStage,counts[0],counts[1],RMHD_FAILSAFE_PCP_STRICT);
            )
    }

  return counts[0] != 0;
}

#if PCP_STAGE_LIMITING == YES
#if EOS != IDEAL
#error "PCP_STAGE_LIMITING requires IDEAL EOS"
#endif
#if ((TIME_STEPPING != RK2) && (TIME_STEPPING != RK3))
#error "PCP_STAGE_LIMITING requires FAILSAFE=YES and SSP RK2/RK3"
#endif
#if ENTROPY_SWITCH != NO || PARABOLIC_FLUX != NO \
 || FORCED_TURB != NO || INCLUDE_LES != NO || DUST_FLUID != NO \
 || RING_AVERAGE > 1
#error "PCP_STAGE_LIMITING: this source/auxiliary update is not supported"
#endif
#if (UPDATE_VECTOR_POTENTIAL == YES) && !defined(STAGGERED_MHD)
#error "PCP_STAGE_LIMITING needs staggered MHD when vector potential is updated"
#endif
#if defined(CHOMBO) || defined(FARGO) || defined(SHEARINGBOX) || defined(HIGH_ORDER)
#error "PCP_STAGE_LIMITING: AMR/orbital/high-order quadrature paths need separate validation"
#endif

static int GlobalBad(int bad)
{
#ifdef PARALLEL
  int all;
  MPI_Allreduce(&bad, &all, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return all;
#else
  return bad;
#endif
}

static int PhysicalPrimitive(const double *v)
{
  int n;
  double v2 = 0.0;
  NVAR_LOOP(n) if (!isfinite(v[n])) return 0;
  for (n = VX1; n <= VX3; n++) v2 += v[n]*v[n];
  return v[RHO] > 0.0 && v[PRS] > 0.0 && v2 < 1.0;
}

/* CandidateBad() is called before the global reduction.  Dump the local
 * evidence first; an MPI abort can otherwise erase the only rank that saw
 * the non-finite or non-admissible cell.  Files are deliberately rank-local
 * and bounded so this remains usable on the production-sized mesh. */
static void PCPDumpBadCells (Data *d, Grid *grid, const char *phase)
{
  char fname[256];
  FILE *fp;
  int i, j, k, n, nbad = 0, truncated = 0;

  snprintf(fname, sizeof(fname),
           "pcp_bad_cells_%s_step%08ld_stage%02d_rank%04d.dat",
           phase, g_stepNumber, g_intStage, prank);
  fp = fopen(fname, "w");
  if (fp == NULL) return;
  fprintf(fp, "# rank=%d step=%ld stage=%d time=%.4g dt=%.4g phase=%s\n",
          prank, g_stepNumber, g_intStage, g_time, g_dt, phase);
  fprintf(fp, "# columns: i j k x1 x2 x3 flags v2 "
              "U[rho mx1 mx2 mx3 bx1 bx2 bx3 eng] "
              "V[rho vx1 vx2 vx3 bx1 by1 bz1 prs]"
#ifdef STAGGERED_MHD
              " Bface[x1 x2 x3]"
#endif
              "\n");

  DOM_LOOP(k,j,i) {
    double v2 = 0.0;
    int bad = (d->flag[k][j][i] &
               (FLAG_CONS2PRIM_FAIL | FLAG_NEGATIVE_PRESSURE |
                FLAG_NEGATIVE_ENERGY | FLAG_NEGATIVE_DENSITY |
                FLAG_PRESSURE_FIX_FAIL | FLAG_PCP_LIMIT_FAIL |
                FLAG_PCP_NR_FAIL)) != 0 ||
              !isfinite(d->Uc[k][j][i][RHO]) ||
              !isfinite(d->Uc[k][j][i][ENG]) ||
              d->Uc[k][j][i][RHO] <= 0.0 ||
              d->Uc[k][j][i][ENG] <= 0.0;
    NVAR_LOOP(n) {
      if (!isfinite(d->Uc[k][j][i][n]) ||
          !isfinite(d->Vc[n][k][j][i])) bad = 1;
    }
    for (n = VX1; n <= VX3; n++) v2 += d->Vc[n][k][j][i]*d->Vc[n][k][j][i];
    if (!(d->Vc[RHO][k][j][i] > 0.0) ||
        !(d->Vc[PRS][k][j][i] > 0.0) || !(v2 < 1.0)) bad = 1;
#ifdef STAGGERED_MHD
    DIM_LOOP(n) if (!isfinite(d->Vs[n][k][j][i])) bad = 1;
#endif
    if (!bad) continue;
    if (nbad >= 1024) { truncated = 1; continue; }
    fprintf(fp, "%d %d %d %.4g %.4g %.4g 0x%04x %.4g",
            i, j, k, grid->x[IDIR][i], grid->x[JDIR][j], grid->x[KDIR][k],
            (unsigned)d->flag[k][j][i], v2);
    for (n = RHO; n <= ENG; n++) fprintf(fp, " %.4g", d->Uc[k][j][i][n]);
    for (n = RHO; n <= PRS; n++) fprintf(fp, " %.4g", d->Vc[n][k][j][i]);
#ifdef STAGGERED_MHD
    DIM_LOOP(n) fprintf(fp, " %.4g", d->Vs[n][k][j][i]);
#endif
    fputc('\n', fp);
    nbad++;
  }
  fprintf(fp, "# bad_cells=%d%s\n", nbad,
          truncated ? " (truncated at 1024)" : "");
  fflush(fp);
  fclose(fp);
}

/* This is ordinary first-order reconstruction, NOT a fabricated safe state.
 * A bad high-order face is downgraded before the Riemann solver can see it.
 * This retains the original interface limiter and extends its failure path
 * only in the opt-in stage-limited mode. */
void RMHD_PCP_StageStates(const Sweep *s, int beg, int end)
{
  int i, n;

  /* UpdateStage calls this with [nbeg-1,nend], the complete interface
   * range consumed by the Riemann solver.  PLUTO stores stateR[i] as the
   * left state of cell i+1, hence stateC[i+1] is required at the final
   * interface.  stateL/stateR also receive the normal staggered field for
   * the current sweep direction (BX1, BX2, or BX3 via BXn). */
  if (end < beg) return;
  for (i = beg; i <= end; i++) {
    if (g_pcp_low_order || (g_rmhd_failsafe_low_order &&
        ((s->flag[i] | s->flag[i+1]) & FLAG_FLAT))
        || !RMHD_PCP_Check(s->stateL.u[i])
        || !RMHD_PCP_Check(s->stateR.u[i])
        || !PhysicalPrimitive(s->stateL.v[i])
        || !PhysicalPrimitive(s->stateR.v[i])) {
      RMHD_FailSafe_StageStates(s, i, i);
    }
  }
}

static int CandidateBad(Data *d, int mark, Grid *grid, const char *phase)
{
  int i,j,k,n, bad = 0;
  DOM_LOOP(k,j,i) {
    double *u = d->Uc[k][j][i];
    int cell_bad = !isfinite(u[RHO]) || !isfinite(u[ENG])
                || u[RHO] <= 0.0 || u[ENG] <= 0.0;
    NVAR_LOOP(n) if (!isfinite(u[n])) cell_bad = 1;
    if (cell_bad) {
      bad = 1;
      if (mark) d->flag[k][j][i] |= FLAG_CONS2PRIM_FAIL;
    }
  }
  /* All stored faces are blended when staggered MHD is active. */
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) if (!isfinite(d->Vs[n][k][j][i])) bad = 1;
#endif
  // if (bad) PCPDumpBadCells(d, grid, phase);
  return GlobalBad(bad);
}

int RMHD_PCP_UpdateStage(Data *d, double dt, timeStep *Dts, Grid *grid)
{
  static Data_Arr input, high, low, bin, bhi, blo;
  timeStep saved_dts = *Dts;
  RBox box;
  int i,j,k,n, it;
  double theta;
  if (input == NULL) {
    input = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
    high  = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
    low   = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
#ifdef STAGGERED_MHD
    bin = ARRAY_4D(DIMENSIONS,NX3_TOT,NX2_TOT,NX1_TOT,double);
    bhi = ARRAY_4D(DIMENSIONS,NX3_TOT,NX2_TOT,NX1_TOT,double);
    blo = ARRAY_4D(DIMENSIONS,NX3_TOT,NX2_TOT,NX1_TOT,double);
#endif
  }
  RBoxDefine(IBEG,IEND,JBEG,JEND,KBEG,KEND,CENTER,&box);
  if (CandidateBad(d, 1, grid, "input")) return 1;
  /* Input primitive and ghost states must be usable before any flux call. */
  {
    int bad = !(g_gamma > 1.0 && g_gamma <= 2.0);
    double v[NVAR];
    TOT_LOOP(k,j,i) {
      NVAR_LOOP(n) v[n] = d->Vc[n][k][j][i];
      if (!PhysicalPrimitive(v)) bad = 1;
    }
    if (GlobalBad(bad)) {
      // PCPDumpBadCells(d, grid, "input_primitive");
      return 1;
    }
  }
  TOT_LOOP(k,j,i) NVAR_LOOP(n) input[k][j][i][n] = d->Uc[k][j][i][n];
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) bin[n][k][j][i] = d->Vs[n][k][j][i];
#endif

  UpdateStage(d, d->Uc, d->Vs, NULL, dt, Dts, grid);
#ifdef STAGGERED_MHD
  CT_AverageStaggeredFields(d, 0, &box, grid);
#endif
  if (!CandidateBad(d, 0, grid, "high")) return 0; /* no low-order work on a good stage */


  //more works begin here: 
  //1. low-order candidate, 2. global halving sequence, 3. collective recheck

  TOT_LOOP(k,j,i) NVAR_LOOP(n) {
    high[k][j][i][n] = d->Uc[k][j][i][n];
    d->Uc[k][j][i][n] = input[k][j][i][n];
  }
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) {
    bhi[n][k][j][i] = d->Vs[n][k][j][i];
    d->Vs[n][k][j][i] = bin[n][k][j][i];
  }
#endif
  *Dts = saved_dts;
  g_pcp_low_order = 1;
  UpdateStage(d, d->Uc, d->Vs, NULL, dt, Dts, grid);
  g_pcp_low_order = 0;
#ifdef STAGGERED_MHD
  CT_AverageStaggeredFields(d, 0, &box, grid);
#endif
  if (CandidateBad(d, 1, grid, "low")) {
    printLog("PCP_STAGE reject: low-order candidate, stage=%d dt=%.4g\n",g_intStage,dt);
    /* FailSafe owns the full-step backup, including primitives and flags. */
    return 1;
  }
  TOT_LOOP(k,j,i) NVAR_LOOP(n) low[k][j][i][n] = d->Uc[k][j][i][n];
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) blo[n][k][j][i] = d->Vs[n][k][j][i];
#endif

  /* Fixed GLOBAL halving sequence + collective recheck. No assumption that
   * the state-dependent floating-point margin checker defines a convex set.
   * theta=0 is copied exactly: 0*NaN must never poison the valid anchor. */
  for (it=1; it<=2; it++) {
    theta = ldexp(1.0,-it);
    DOM_LOOP(k,j,i) NVAR_LOOP(n) {
      d->Uc[k][j][i][n] = theta == 0.0 ? low[k][j][i][n]
        : (1.0-theta)*low[k][j][i][n] + theta*high[k][j][i][n];
    }
#ifdef STAGGERED_MHD
    DIM_LOOP(n) TOT_LOOP(k,j,i) {
      d->Vs[n][k][j][i] = theta == 0.0 ? blo[n][k][j][i]
        : (1.0-theta)*blo[n][k][j][i] + theta*bhi[n][k][j][i];
    }
#endif
#ifdef STAGGERED_MHD
    CT_AverageStaggeredFields(d, 0, &box, grid);
#endif
    if (!CandidateBad(d, 0, grid, "blend")) {
      printLog("PCP_STAGE limited: stage=%d theta=%.4g\n",g_intStage,theta);
      return 0;
    }
  }
  /* Supply troubled-cell seeds to the full-step local fallback. */
  CandidateBad(d, 1, grid, "blend_final");
  return 1;
}

#endif
