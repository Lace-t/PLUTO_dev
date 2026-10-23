/* Local a-posteriori RK stage repair. Candidate validation runs AFTER
 * RK averaging, CT correction and primitive recovery. Only flagged shared
 * faces/edges use the existing LF/CT fallback; no global state blending.
 * A failed stage is recomputed from its own entry, not from the full step.
 * The bounded repair is not a proof of a PCP low-order anchor. */
#include "pluto.h"

int g_pcp_low_order = 0;

/* User override point for conservative floors before stage PCP checks. */
#if defined(__GNUC__)
__attribute__((weak))
#endif
void UserDefCFloors(Data *d, Grid *grid)
{
  (void)d;
  (void)grid;
}

static int stage_pass;
#if PCP_FOFC == YES
static Data_Arr stage_u;
static int PCPStageCellAffected(int k, int j, int i);
static Data_Arr stage_candidate_u;
#ifdef STAGGERED_MHD
static Data_Arr stage_candidate_b;
#endif
#endif

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

/* Check the actual RK average after its RK combination and CT averaging,
 * before the mapper can change U. No primitive recovery or state repair is
 * performed here. All ranks must make the same decision before returning
 * from AdvanceStep, since the following C2P contains MPI collectives. */
int RMHD_PCP_CheckStage(Data *d, Grid *grid)
{
  int i,j,k,n, local_cells = 0, rejected_cells;
  (void)grid;
#if PCP_FOFC == YES
  if (stage_u != NULL && stage_pass == 0) {
    DOM_LOOP(k,j,i) {
      memcpy(stage_candidate_u[k][j][i],d->Uc[k][j][i],NVAR*sizeof(double));
#ifdef STAGGERED_MHD
      DIM_LOOP(n) stage_candidate_b[n][k][j][i]=d->Vs[n][k][j][i];
#endif
    }
  } else if (stage_u != NULL) {
#ifndef STAGGERED_MHD
    DOM_LOOP(k,j,i) if (!PCPStageCellAffected(k,j,i))
      memcpy(d->Uc[k][j][i],stage_candidate_u[k][j][i],NVAR*sizeof(double));
#endif
  }
#endif
  DOM_LOOP(k,j,i) {
#if PCP_FOFC == YES
    /* Once FOFC changes local fluxes, only this stencil can have a new
     * conservative candidate. The first candidate still needs a full scan. */
    if (stage_u != NULL && stage_pass > 0 && !PCPStageCellAffected(k,j,i)) {
#ifndef STAGGERED_MHD
      continue;
#else
      /* CT retries rebuild the full stage without a cached edge/flux update. */
#endif
    }
#endif
    double *u = d->Uc[k][j][i];
    /* RMHD_PCP_Check is the sole conservative PCP criterion. */
    int cell_bad = !RMHD_PCP_Check(u);
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
  }
  if (local_cells) { WARNING(printLog("RMHD_STAGE_REJECT: step=%ld stage=%d rank=%d "
                   "cells=%d dt=%.4g pre_c2p=YES\n",
                   g_stepNumber,g_intStage,prank,local_cells,g_dt);) }
  rejected_cells = local_cells;
#ifdef PARALLEL
  {
    MPI_Allreduce(&local_cells,&rejected_cells,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  }
#endif
  if (prank == 0 && rejected_cells){
    WARNING(printLog("RMHD_STAGE_CHECK: step=%ld stage=%d rejected=%d "
             "strict=1\n",
             g_stepNumber,g_intStage,rejected_cells);
            )
    }

  return rejected_cells != 0;
}


#if PCP_FOFC == YES
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

static Data_Arr stage_v, stage_b;
static uint16_t ***stage_flags, ***stage_mask;
static double *stage_cmax;
static timeStep stage_dts;
static int outer_mode;
static Data_Arr stage_flux[DIMENSIONS];
static double ***stage_press[DIMENSIONS];
static double ***stage_speed[DIMENSIONS];

/* Store one face's high-order flux and pressure term by global-in-rank
 * coordinates. The cache is stage-local and reused only on FOFC passes. */
static void PCPStageFluxAddress(int dir, int k, int j, int i,
                                int *fk, int *fj, int *fi)
{
  if (dir == IDIR) { *fk=k; *fj=j; *fi=i; }
  else if (dir == JDIR) { *fk=k; *fj=j; *fi=i; }
  else { *fk=k; *fj=j; *fi=i; }
}

int RMHD_PCP_UseFluxCache(void)
{
#ifdef STAGGERED_MHD
  return 0;
#else
  return stage_pass > 0;
#endif
}

void RMHD_PCP_CacheSweepFlux(int dir, int k, int j, int i,
                             const Sweep *s, const double *cmax, int q)
{
  int fk,fj,fi,n;
  if (stage_pass != 0 || dir < 0 || dir >= DIMENSIONS ||
      stage_flux[dir] == NULL) return;
  PCPStageFluxAddress(dir,k,j,i,&fk,&fj,&fi);
  NVAR_LOOP(n) stage_flux[dir][fk][fj][fi][n]=s->flux[q][n];
  stage_press[dir][fk][fj][fi]=s->press[q];
  stage_speed[dir][fk][fj][fi]=cmax[q];
}

void RMHD_PCP_LoadSweepFlux(int dir, int k, int j, int i,
                            Sweep *s, double *cmax, int q)
{
  int fk,fj,fi,n;
  if (dir < 0 || dir >= DIMENSIONS || stage_flux[dir] == NULL) return;
  PCPStageFluxAddress(dir,k,j,i,&fk,&fj,&fi);
  NVAR_LOOP(n) s->flux[q][n]=stage_flux[dir][fk][fj][fi][n];
  s->press[q]=stage_press[dir][fk][fj][fi];
  /* RetryStage restores the stage-entry CFL storage, which may still be
   * uninitialized. Every replayed face needs its original signal speed. */
  cmax[q]=stage_speed[dir][fk][fj][fi];
}

/* Corrected face fluxes affect face-neighbor cell updates. In CT, the
 * associated edge EMFs also update nearby face fields, so include the full
 * one-cell Chebyshev neighborhood as a conservative dependency stencil. */
static int PCPStageCellAffected(int k, int j, int i)
{
  int di,dj,dk;
#ifdef STAGGERED_MHD
  for (dk=-1; dk<=1; dk++) {
    if (DIMENSIONS < 3 && dk != 0) continue;
    if (k+dk < 0 || k+dk >= NX3_TOT) continue;
    for (dj=-1; dj<=1; dj++) {
      if (DIMENSIONS < 2 && dj != 0) continue;
      if (j+dj < 0 || j+dj >= NX2_TOT) continue;
      for (di=-1; di<=1; di++) {
        if (i+di < 0 || i+di >= NX1_TOT) continue;
        if (stage_mask[k+dk][j+dj][i+di] & FLAG_LF) return 1;
      }
    }
  }
#else
  if (stage_mask[k][j][i] & FLAG_LF) return 1;
  for (di=-1; di<=1; di+=2)
    if (i+di >= 0 && i+di < NX1_TOT &&
        (stage_mask[k][j][i+di] & FLAG_LF)) return 1;
#if DIMENSIONS >= 2
  for (dj=-1; dj<=1; dj+=2)
    if (j+dj >= 0 && j+dj < NX2_TOT &&
        (stage_mask[k][j+dj][i] & FLAG_LF)) return 1;
#endif
#if DIMENSIONS == 3
  for (dk=-1; dk<=1; dk+=2)
    if (k+dk >= 0 && k+dk < NX3_TOT &&
        (stage_mask[k+dk][j][i] & FLAG_LF)) return 1;
#endif
#endif
  return 0;
}

/* Called after Boundary/initial P2C, once per RK stage. V must be saved:
 * CT_RMHD_CONSISTENT and the mapper can both modify it on rejected trials. */
void RMHD_PCP_BeginStage(Data *d, timeStep *Dts)
{
  int i,j,k,n;
  if (stage_u == NULL) {
    stage_u = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
    stage_candidate_u = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
#ifdef STAGGERED_MHD
    stage_candidate_b = ARRAY_4D(DIMENSIONS,NX3_TOT,NX2_TOT,NX1_TOT,double);
#endif
    stage_v = ARRAY_4D(NVAR,NX3_TOT,NX2_TOT,NX1_TOT,double);
    stage_flags = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    stage_mask = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    stage_cmax = ARRAY_1D(NMAX_POINT,double);
    for (n=0; n<DIMENSIONS; n++) {
#ifndef STAGGERED_MHD
      stage_flux[n] = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
      stage_press[n] = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);
      stage_speed[n] = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,double);
#endif
    }
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
  long local[3]={0,0,0}, total[3];
  /* Count bad, cumulative and newly marked cells before restoration. All
   * counts share the one collective needed to decide whether to retry. */
  DOM_LOOP(k,j,i) {
    int bad=!!(d->flag[k][j][i]&FLAG_CONS2PRIM_FAIL);
    local[0] += bad;
    local[1] += bad || (stage_mask[k][j][i]&FLAG_LF);
    local[2] += bad && !(stage_mask[k][j][i]&FLAG_LF);
  }
#ifdef PARALLEL
  MPI_Allreduce(local,total,3,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
#else
  memcpy(total,local,sizeof total);
#endif
  if (stage_pass == PCP_STAGE_LOCAL_PASSES || total[0] == 0 ||
      (stage_pass > 0 && total[2] == 0)) {
    /* The input snapshot, dt and LF operator have not changed. Repeating
     * this candidate without any new global mask cells cannot repair it. */
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

/* Own the FOFC retry loop so the RK driver does not rerun its whole stage.
 * The cache path replays the first high-order fluxes and only invokes LF on
 * interfaces touching the cumulative failure mask. */
int RMHD_PCP_RunStage(Data *d, timeStep *Dts, double dt,
                      Data_Arr Ubase, double ****Bbase, int combination,
                      double w0_stage, double wc_stage, RBox *box, Grid *grid)
{
  int i,j,k,n,err;
  RMHD_PCP_BeginStage(d,Dts);
  do {
    UpdateStage(d,d->Uc,d->Vs,NULL,dt,Dts,grid);
#ifdef STAGGERED_MHD
    /* The predictor also needs CT synchronization and primitive recovery. */
    if (combination == 0) CT_AverageStaggeredFields(d,1,box,grid);
#endif
    if (combination == 1) {
      DOM_LOOP(k,j,i) NVAR_LOOP(n)
        d->Uc[k][j][i][n] = w0_stage*Ubase[k][j][i][n]
                          + wc_stage*d->Uc[k][j][i][n];
#ifdef STAGGERED_MHD
      DIM_LOOP(n) TOT_LOOP(k,j,i)
        d->Vs[n][k][j][i] = w0_stage*Bbase[n][k][j][i]
                          + wc_stage*d->Vs[n][k][j][i];
      CT_AverageStaggeredFields(d,1,box,grid);
#endif
    } else if (combination == 2) {
      DOM_LOOP(k,j,i) NVAR_LOOP(n)
        d->Uc[k][j][i][n] = (Ubase[k][j][i][n]
                           + 2.0*d->Uc[k][j][i][n])/3.0;
#ifdef STAGGERED_MHD
      DIM_LOOP(n) TOT_LOOP(k,j,i)
        d->Vs[n][k][j][i] = (Bbase[n][k][j][i]
                           + 2.0*d->Vs[n][k][j][i])/3.0;
      CT_AverageStaggeredFields(d,1,box,grid);
#endif
    }
    /* The cache contains the first fully combined RK candidate. A retry
     * replays cached high-order fluxes plus local LF replacements; retain
     * the bitwise first candidate outside the flux-dependency stencil. */
    if (stage_pass > 0
#ifndef STAGGERED_MHD
        ) DOM_LOOP(k,j,i)
#else
        && RMHD_PCP_UseFluxCache()) DOM_LOOP(k,j,i)
#endif
      if (!PCPStageCellAffected(k,j,i))
        memcpy(d->Uc[k][j][i],stage_candidate_u[k][j][i],NVAR*sizeof(double));
#ifdef STAGGERED_MHD
    if (stage_pass > 0 && RMHD_PCP_UseFluxCache()) {
      /* CT faces are restored with the matching candidate cell stencil. */
      for (n=0;n<DIMENSIONS;n++) TOT_LOOP(k,j,i)
        if (!PCPStageCellAffected(k,j,i)) d->Vs[n][k][j][i]=stage_candidate_b[n][k][j][i];
    }
#endif
    UserDefCFloors(d,grid);
    err=RMHD_PCP_CheckStage(d,grid);
#if !defined(STAGGERED_MHD) || CT_RMHD_CONSISTENT != YES
    if (!err)
      err=RMHD_PCP_ConsToPrim3D_CheckModified(d->Uc,d->Vc,d->flag,box,grid);
#endif
  } while (err && RMHD_PCP_RetryStage(d,Dts));
  if (err) err=RMHD_PCP_TryPressureFix(d,err,box,grid);
  return RMHD_PCP_FinishStage(d,err,box,grid);
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

#endif /* PCP_FOFC: stage snapshots and local retry machinery */

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

/* Recover failed cells directly at floor pressure; ordinary C2P is only
 * needed for healthy cells whose primitives are still from stage entry. */
static int PCPRecoverAfterFOFC(Data *d)
{
  static double **v;
  int i,j,k,n,first,err=0;
  if (v == NULL) v=ARRAY_2D(NX1_TOT,NVAR,double);
  for(k=KBEG;k<=KEND;k++) for(j=JBEG;j<=JEND;j++) {
    for(i=IBEG;i<=IEND;i++) NVAR_LOOP(n) v[i][n]=d->Vc[n][k][j][i];
    for(i=IBEG;i<=IEND;) {
      if ((d->flag[k][j][i]&FLAG_CONS2PRIM_FAIL) ||
          !RMHD_PCP_Check(d->Uc[k][j][i])) {
        double *u=d->Uc[k][j][i];
        d->flag[k][j][i]&=~(FLAG_CONS2PRIM_FAIL|FLAG_PRESSURE_FIX_FAIL);
        if (RMHD_PressureFix(u,v[i])) {
          d->flag[k][j][i]|=FLAG_CONS2PRIM_FAIL|FLAG_PRESSURE_FIX_FAIL;
          err=1;
        } else {
#if NSCL > 0
          NSCL_LOOP(n) v[i][n]=u[n]/u[RHO];
#endif
#ifdef GLM_MHD
          v[i][PSI_GLM]=u[PSI_GLM];
#endif
        }
        i++;
      } else {
        first=i++;
        while(i<=IEND && !(d->flag[k][j][i]&FLAG_CONS2PRIM_FAIL)
              && RMHD_PCP_Check(d->Uc[k][j][i])) i++;
        g_j=j;g_k=k;
        if (ConsToPrim(d->Uc[k][j],v,first,i-1,d->flag[k][j])) err=1;
      }
    }
    for(i=IBEG;i<=IEND;i++) NVAR_LOOP(n) d->Vc[n][k][j][i]=v[i][n];
  }
#ifdef PARALLEL
  {int local=err;MPI_Allreduce(&local,&err,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);}
#endif
  return err;
}

/* Try PressureFix after local FOFC retries are exhausted.
 * Failure restores U/V so the outer full-step failsafe sees the original
 * rejected candidate and its failure seeds.
 * Recovery at floor pressure uses finite/physical checks, not G_epsilon. */
int RMHD_PCP_TryPressureFix(Data *d, int err, RBox *box, Grid *grid)
{
#if PCP_PRESSURE_FIX == YES
  int i,j,k,n;
  int outer_mode=g_rmhd_failsafe_low_order;
  (void)box;
  if (err) {
    int pressure_fix_err = 0;
    static Data_Arr rejected_u, rejected_v;
    static uint16_t ***rejected_flags;
    long changed=0, global_changed;
    double added=0.0, max_relative=0.0, global_added, global_relative;

    /* Do not let ordinary mapper floors invent density or mask NaNs.
     * All ranks make this decision before entering its collectives. */
    if (PCPRecoveryInvalid(d,0)) {
      g_rmhd_failsafe_low_order=outer_mode;
      return err;
    }
    if (rejected_u == NULL) {
      rejected_u=ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
      rejected_v=ARRAY_4D(NVAR,NX3_TOT,NX2_TOT,NX1_TOT,double);
      rejected_flags=ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    }
    DOM_LOOP(k,j,i) {
      NVAR_LOOP(n) {
        rejected_u[k][j][i][n]=d->Uc[k][j][i][n];
        rejected_v[n][k][j][i]=d->Vc[n][k][j][i];
      }
      rejected_flags[k][j][i]=d->flag[k][j][i];
    }

    pressure_fix_err = PCPRecoverAfterFOFC(d);
    pressure_fix_err = PCPRecoveryInvalid(d,1) || pressure_fix_err;
    if (!pressure_fix_err) DOM_LOOP(k,j,i) {
      double before=rejected_u[k][j][i][ENG];
      double delta=d->Uc[k][j][i][ENG]-before;
      if (delta != 0.0) {
        changed++;
        added+=delta*grid->dV[k][j][i];
        max_relative=fmax(max_relative,fabs(delta)/fmax(fabs(before),1.e-300));
      }
    }
    global_changed=changed; global_added=added; global_relative=max_relative;
#ifdef PARALLEL
    MPI_Allreduce(&changed,&global_changed,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&added,&global_added,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&max_relative,&global_relative,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
    // printLog("PCP_STAGE recovery: step=%ld stage=%d pressure_fix=%s "
    //          "energy_changed=%ld dE_volume=%.17g max_relative=%.17g\n",
    //          g_stepNumber,g_intStage,pressure_fix_err ? "FAILED" : "accepted",
    //          global_changed,global_added,global_relative);
    if (!pressure_fix_err) {
      return 0;
    }

    /* Restore the rejected candidate for the outer fallback/dump. Keep
     * original FOFC seeds as well as newly observed recovery failures. */
    DOM_LOOP(k,j,i) {
      NVAR_LOOP(n) {
        d->Uc[k][j][i][n]=rejected_u[k][j][i][n];
        d->Vc[n][k][j][i]=rejected_v[n][k][j][i];
      }
      d->flag[k][j][i] |= rejected_flags[k][j][i];
    }
    g_rmhd_failsafe_low_order=outer_mode;
    return pressure_fix_err;
  }
  return 0;
#else
  (void)d; (void)box; (void)grid;
  return err;
#endif
}

/* Complete after FOFC and the optional post-FOFC PressureFix attempt. */
int RMHD_PCP_FinishStage(Data *d, int err, RBox *box, Grid *grid)
{
  (void)box; (void)grid;
  if (err) {
#if FAILSAFE == NO
    printLog("! PCP_STAGE: recovery exhausted with FAILSAFE=NO\n");
    LogFileFlush(); QUIT_PLUTO(1);
#endif
    return err;
  }
#if PCP_FOFC == YES
  RMHD_PCP_EndStage(d);
#endif
  return 0;
}
