/* Local a-posteriori RK stage repair. Candidate validation runs AFTER
 * RK averaging, CT correction and primitive recovery. Only flagged shared
 * faces/edges use the existing LF/CT fallback; no global state blending.
 * A failed stage is recomputed from its own entry, not from the full step.
 * The bounded repair is not a proof of a PCP low-order anchor. */
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

/* Replace shared face fluxes, not cell updates. FLAG_LF is the undilated
 * FOFC cell mask; FLAG_FLAT also supports the outer full-step fallback.
 * Rebuild BOTH states after the ordinary solver, which may modify them
 * (e.g. its GLM pre-solver), and refresh LF's CT flux data at the same faces. */
void RMHD_FailSafe_CorrectFluxes(const Sweep *s, int beg, int end,
                               double *cmax, Grid *grid)
{
  int face, first;
  if (!g_rmhd_failsafe_low_order) return;
  for (face = beg; face <= end; ) {
    if (!((s->flag[face] | s->flag[face+1]) & (FLAG_FLAT | FLAG_LF))) {
      face++;
      continue;
    }
    first = face++;
    while (face <= end &&
           ((s->flag[face] | s->flag[face+1]) & (FLAG_FLAT | FLAG_LF))) face++;
    RMHD_FailSafe_StageStates(s,first,face-1);
    LF_Solver(s,first,face-1,cmax,grid);
  }
}

/* Shared acceptance policy for stage input and mapper-modified states. */
int RMHD_PCP_CheckConservative(const double *u)
{
  int n;
  NVAR_LOOP(n) if (!isfinite(u[n])) return 0;
#if RMHD_FAILSAFE_PCP_STRICT == YES
  return RMHD_PCP_Check(u);
#else
  if (RMHD_PCP_Check(u)) return 1;
  double E = u[ENG];
#if RMHD_REDUCED_ENERGY == YES
  E += u[RHO];
#endif
  /* Compatibility mode defers the remaining PCP conditions to recovery. */
  return isfinite(E) && u[RHO] > 0.0 && u[ENG] > 0.0
      && E > hypot(u[RHO],hypot(u[MX1],hypot(u[MX2],u[MX3])));
#endif
}

/* Check the actual RK average after its RK combination and CT averaging,
 * before the mapper can change U. No primitive recovery or state repair is
 * performed here. All ranks must make the same decision before returning
 * from AdvanceStep, since the following C2P contains MPI collectives. */
int RMHD_PCP_CheckStage(Data *d, Grid *grid)
{
  int i,j,k, local_cells = 0, deferred_cells = 0;
  int counts[2];
  (void)grid;
  DOM_LOOP(k,j,i) {
    double *u = d->Uc[k][j][i];
    int cell_bad = !RMHD_PCP_CheckConservative(u);
    if (d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL) cell_bad = 1;
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
    }
#if RMHD_FAILSAFE_PCP_STRICT == NO
    else if (!RMHD_PCP_Check(u)) deferred_cells++;
#endif
  }
  if (local_cells) WARNING(printLog("RMHD_STAGE_REJECT: step=%ld stage=%d rank=%d "
                   "cells=%d dt=%.4g pre_c2p=YES\n",
                   g_stepNumber,g_intStage,prank,local_cells,g_dt);)
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
#if EOS != IDEAL || GEOMETRY != CARTESIAN \
 || ((TIME_STEPPING != RK2) && (TIME_STEPPING != RK3))
#error "Local PCP stages require Cartesian IDEAL RMHD and SSP RK2/RK3"
#endif
#if ENTROPY_SWITCH != NO || PARABOLIC_FLUX != NO || RADIATION != NO \
 || PARTICLES != NO || FORCED_TURB != NO || INCLUDE_LES != NO \
 || DUST_FLUID != NO || RING_AVERAGE > 1 || UPDATE_VECTOR_POTENTIAL == YES
#error "Local PCP stages: rollback for these auxiliary fields is not implemented"
#endif
#if defined(CHOMBO) || defined(FARGO) || defined(SHEARINGBOX) || defined(HIGH_ORDER)
#error "Local PCP stages: this integration path is unsupported"
#endif

static Data_Arr stage_u, stage_v, stage_b;
static uint16_t ***stage_flags, ***stage_mask;
static double *stage_cmax;
static timeStep stage_dts;
static int stage_pass, outer_mode;

/* Called after Boundary/initial P2C, once per RK stage. V must be saved:
 * CT_RMHD_CONSISTENT and the mapper can both modify it on rejected trials. */
void RMHD_PCP_BeginStage(Data *d, timeStep *Dts)
{
  int i,j,k,n;
  if (stage_u == NULL) {
    stage_u = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
    stage_v = ARRAY_4D(NVAR,NX3_TOT,NX2_TOT,NX1_TOT,double);
    stage_flags = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    stage_mask = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    stage_cmax = ARRAY_1D(NMAX_POINT,double);
#ifdef STAGGERED_MHD
    stage_b = ARRAY_4D(DIMENSIONS,NX3_TOT,NX2_TOT,NX1_TOT,double);
#endif
  }
  stage_pass = 0;
  outer_mode = g_rmhd_failsafe_low_order;
  stage_dts = *Dts;
  memcpy(stage_cmax,Dts->cmax,NMAX_POINT*sizeof(double));
  TOT_LOOP(k,j,i) {
    stage_flags[k][j][i] = d->flag[k][j][i];
    stage_mask[k][j][i] = 0;
  }
  for(k=0;k<NX3_TOT;k++) for(j=0;j<NX2_TOT;j++)
    memcpy(stage_u[k][j][0],d->Uc[k][j][0],NX1_TOT*NVAR*sizeof(double));
  NVAR_LOOP(n) for(k=0;k<NX3_TOT;k++) for(j=0;j<NX2_TOT;j++)
    memcpy(stage_v[n][k][j],d->Vc[n][k][j],NX1_TOT*sizeof(double));
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) stage_b[n][k][j][i] = d->Vs[n][k][j][i];
#endif
}

/* Every rank calls this after a collective stage rejection. Return one
 * only when all ranks have restored the same stage and will retry it.
 * FOFC marks only failed cells: the two faces touching a marked cell in
 * each direction get first-order LLF fluxes. Do not dilate into the
 * FLAT/MINMOD rings used by the outer full-step failsafe. */
int RMHD_PCP_RetryStage(Data *d, timeStep *Dts)
{
  int i,j,k,n;
  long local[2]={0,0}, total[2];
  /* Count the cumulative mask with the new seeds before restoration. Both
   * counts share the one collective needed to decide whether to retry. */
  DOM_LOOP(k,j,i) {
    int bad=!!(d->flag[k][j][i]&FLAG_CONS2PRIM_FAIL);
    local[0] += bad;
    local[1] += bad || (stage_mask[k][j][i]&FLAG_LF);
  }
#ifdef PARALLEL
  MPI_Allreduce(local,total,2,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
#else
  memcpy(total,local,sizeof total);
#endif
  if (stage_pass == PCP_STAGE_LOCAL_PASSES || total[0] == 0) {
    /* Retain all observed seeds for the outer full-step fallback, while
     * leaving rejected U/V/B intact for its diagnostics. */
    TOT_LOOP(k,j,i) if (stage_mask[k][j][i]&FLAG_LF)
      d->flag[k][j][i] |= FLAG_CONS2PRIM_FAIL;
    WARNING(printLog("PCP_STAGE exhausted: step=%ld stage=%d passes=%d bad=%ld\n",
             g_stepNumber,g_intStage,stage_pass,total[0]);)
    g_rmhd_failsafe_low_order = outer_mode;
    return 0;
  }
  /* No dilation: install local seeds first, then exchange the final flags
   * once. ArrayLib's directional exchange also propagates CT edge/corners. */
  TOT_LOOP(k,j,i) if (d->flag[k][j][i]&FLAG_CONS2PRIM_FAIL)
    stage_mask[k][j][i] |= FLAG_LF;
  *Dts = stage_dts;
  memcpy(Dts->cmax,stage_cmax,NMAX_POINT*sizeof(double));
  for(k=0;k<NX3_TOT;k++) for(j=0;j<NX2_TOT;j++)
    memcpy(d->Uc[k][j][0],stage_u[k][j][0],NX1_TOT*NVAR*sizeof(double));
  NVAR_LOOP(n) for(k=0;k<NX3_TOT;k++) for(j=0;j<NX2_TOT;j++)
    memcpy(d->Vc[n][k][j],stage_v[n][k][j],NX1_TOT*sizeof(double));
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) d->Vs[n][k][j][i] = stage_b[n][k][j][i];
#endif
  TOT_LOOP(k,j,i) d->flag[k][j][i] = stage_flags[k][j][i] | stage_mask[k][j][i];
#ifdef PARALLEL
  AL_Exchange(d->flag[0][0],SZ_uint16_t);
  /* Track received ghost masks too, so EndStage removes remote seeds.
   * Original LF flags are restored from stage_flags during cleanup. */
  TOT_LOOP(k,j,i) stage_mask[k][j][i] |=
    d->flag[k][j][i] & FLAG_LF & ~stage_flags[k][j][i];
#endif
  stage_pass++;
  g_rmhd_failsafe_low_order = 1;
  WARNING(printLog("PCP_STAGE retry: step=%ld stage=%d pass=%d bad=%ld fofc=%ld\n",
           g_stepNumber,g_intStage,stage_pass,total[0],total[1]);)
  return 1;
}

void RMHD_PCP_EndStage(Data *d)
{
  int i,j,k;
  /* No retry means no temporary mask was installed. Accepted stage checks
   * only add diagnostic bits, so there is no full-grid cleanup to do. */
  if (stage_pass == 0) {
    g_rmhd_failsafe_low_order = outer_mode;
    return;
  }
  /* Repair masks are stage-local; preserve original flags and accepted
   * recovery diagnostics, removing only bits introduced by this repair. */
  TOT_LOOP(k,j,i) d->flag[k][j][i] =
    (d->flag[k][j][i] & ~stage_mask[k][j][i]) | stage_flags[k][j][i];
  if (stage_pass)
    WARNING(printLog("PCP_STAGE accepted: step=%ld stage=%d passes=%d\n",
             g_stepNumber,g_intStage,stage_pass);)
  g_rmhd_failsafe_low_order = outer_mode;
}

/* PressureFix may change energy and need not pass the strict PCP test at
 * floor pressure. Require a finite, usable state instead. Before recovery,
 * reject inputs PressureFix cannot repair without inventing density/B.
 * Reduce on all ranks before deciding whether to enter the C2P collective. */
static int PCPRecoveryInvalid(Data *d, int recovered)
{
  int i,j,k,n,err=0;
  DOM_LOOP(k,j,i) {
    int bad=0;
    NVAR_LOOP(n) if (!isfinite(d->Uc[k][j][i][n])) bad=1;
    if (d->Uc[k][j][i][RHO] <= 0.0) bad=1;
#ifdef STAGGERED_MHD
    DIM_EXPAND(
      if (!isfinite(d->Vs[IDIR][k][j][i]) ||
          !isfinite(d->Vs[IDIR][k][j][i-1])) bad=1;,
      if (!isfinite(d->Vs[JDIR][k][j][i]) ||
          !isfinite(d->Vs[JDIR][k][j-1][i])) bad=1;,
      if (!isfinite(d->Vs[KDIR][k][j][i]) ||
          !isfinite(d->Vs[KDIR][k-1][j][i])) bad=1;
    )
#endif
    if (recovered) {
      double v2=0.0;
      NVAR_LOOP(n) if (!isfinite(d->Vc[n][k][j][i])) bad=1;
      for (n=VX1;n<=VX3;n++) v2+=d->Vc[n][k][j][i]*d->Vc[n][k][j][i];
      if (d->Vc[RHO][k][j][i]<=0.0 || d->Vc[PRS][k][j][i]<=0.0 ||
          !isfinite(v2) || v2>=1.0) bad=1;
      if (d->flag[k][j][i] & (FLAG_CONS2PRIM_FAIL | FLAG_PRESSURE_FIX_FAIL)) bad=1;
    }
    if (bad) {d->flag[k][j][i] |= FLAG_CONS2PRIM_FAIL; err=1;}
  }
#ifdef PARALLEL
  {int local=err; MPI_Allreduce(&local,&err,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);}
#endif
  return err;
}

/* Common exit for every RK stage. The final FOFC candidate (not the stage
 * snapshot) is handed to the configured recovery path after local retries
 * are exhausted.  FAILSAFE owns the failed candidate when enabled; with
 * FAILSAFE disabled, use the native RMHD PressureFix explicitly. */
int RMHD_PCP_FinishStage(Data *d, int err, RBox *box, Grid *grid)
{
  int i,j,k;
  (void)box;
  (void)grid;
  if (err) {
#if FAILSAFE == YES
    /* Leave the final candidate and its diagnostic flags intact.  The caller
     * returns this error to the full-step failsafe, which owns the fallback
     * reconstruction and retry. */
    g_rmhd_failsafe_low_order=outer_mode;
    return err;
#else
    int pressure_fix_err = 0;

    /* Give the ordinary mapper its normal opportunity first.  Its recovery
     * path invokes PressureFix and also refreshes all primitive fields. */
    DOM_LOOP(k,j,i)
      d->flag[k][j][i] &= ~(FLAG_CONS2PRIM_FAIL | FLAG_PRESSURE_FIX_FAIL);
    pressure_fix_err = ConsToPrim3D(d->Uc,d->Vc,d->flag,box,grid);
    pressure_fix_err = PCPRecoveryInvalid(d,1) || pressure_fix_err;
    if (!pressure_fix_err) {
      RMHD_PCP_EndStage(d);
      return 0;
    }

    /* The ordinary mapper can reject a state before reaching its pressure
     * fix (for example after a non-finite intermediate).  Retry the failed
     * cells explicitly before declaring the no-failsafe build fatal. */
    pressure_fix_err = 0;
    DOM_LOOP(k,j,i) {
      if ((d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL) ||
          (stage_mask[k][j][i] & FLAG_LF)) {
        int fix_err;
        d->flag[k][j][i] &= ~(FLAG_CONS2PRIM_FAIL | FLAG_PRESSURE_FIX_FAIL);
        d->flag[k][j][i] |= FLAG_NEGATIVE_PRESSURE;
        fix_err = RMHD_PressureFix(d->Uc[k][j][i], d->Vc[k][j][i]);
        if (fix_err) {
          d->flag[k][j][i] |= FLAG_PRESSURE_FIX_FAIL | FLAG_CONS2PRIM_FAIL;
          pressure_fix_err = 1;
        }
      }
    }
#ifdef PARALLEL
    {int local=pressure_fix_err;
     MPI_Allreduce(&local,&pressure_fix_err,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);}
#endif
    pressure_fix_err = PCPRecoveryInvalid(d,1) || pressure_fix_err;
    printLog("PCP_STAGE recovery: step=%ld stage=%d pressure_fix=%s\n",
             g_stepNumber,g_intStage,pressure_fix_err ? "FAILED" : "accepted");
    if (!pressure_fix_err) {
      RMHD_PCP_EndStage(d);
      return 0;
    }

    g_rmhd_failsafe_low_order=outer_mode;
    printLog("! PCP_STAGE: PressureFix failed with FAILSAFE=NO; "
             "step=%ld stage=%d dt=%.17g\n",g_stepNumber,g_intStage,g_dt);
    LogFileFlush();
    QUIT_PLUTO(1);
    return pressure_fix_err;
#endif
  }
  RMHD_PCP_EndStage(d);
  return 0;
}
#endif
