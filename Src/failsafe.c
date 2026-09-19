#include "pluto.h"

#if (PHYSICS == RMHD)
/* Retry mode enables LOCAL face/edge selection from synchronized flags.
 * It never forces unmarked interfaces or edges to first order. */
int g_rmhd_failsafe_low_order = 0;

/*
 * Write a rank-local snapshot before the recovery path changes the state.
 * The normal PLUTO log is intentionally rank-0 oriented, which hides the
 * cell that caused a globally reduced stage failure.  Keep this diagnostic
 * local and bounded: count all active FLAG_CONS2PRIM_FAIL cells, writing
 * at most 1024 rows per rank with the physical primitive/conservative
 * variables and the positive staggered face fields at each reported cell.
 * Files are written in the run directory and closed before MPI abort.
 */
static void FailSafeDumpBadCells (Data *d, Grid *grid, const char *phase)
{
  char fname[256];
  FILE *fp;
  int i, j, k, n, nbad = 0, nwritten = 0;

  DOM_LOOP(k,j,i) if (d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL) nbad++;
  if (nbad == 0) return; /* Do not create empty files on healthy ranks. */

  snprintf(fname, sizeof(fname),
           "bad_cells_%s_step%08ld_stage%02d_rank%04d.dat",
           phase, g_stepNumber, g_intStage, prank);
  fp = fopen(fname, "w");
  if (fp == NULL) return;

  fprintf(fp, "# rank=%d step=%ld stage=%d time=%.17g dt=%.17g phase=%s\n",
          prank, g_stepNumber, g_intStage, g_time, g_dt, phase);
  fprintf(fp, "# U is the rejected candidate; V may be the previous stage "
              "when the pre-C2P check rejected U.\n");
  fprintf(fp, "# columns: i j k x1 x2 x3 flags v2 "
              "U[rho mx1 mx2 mx3 bx1 bx2 bx3 eng] "
              "V[rho vx1 vx2 vx3 bx1 bx2 bx3 prs]"
#ifdef STAGGERED_MHD
              " Bface[x1 x2 x3]"
#endif
              "\n");

  DOM_LOOP(k,j,i) {
    double v2 = 0.0;
    /* Match the failure flag used by FailSafe, excluding recovered cells. */
    if (!(d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL)) continue;
    if (nwritten >= 1024) continue;
    for (n = VX1; n <= VX3; n++) v2 += d->Vc[n][k][j][i]*d->Vc[n][k][j][i];
    fprintf(fp, "%d %d %d %.17g %.17g %.17g 0x%04x %.17g",
            i, j, k, grid->x[IDIR][i], grid->x[JDIR][j], grid->x[KDIR][k],
            (unsigned)d->flag[k][j][i], v2);
    for (n = RHO; n <= ENG; n++) fprintf(fp, " %.17g", d->Uc[k][j][i][n]);
    for (n = RHO; n <= PRS; n++) fprintf(fp, " %.17g", d->Vc[n][k][j][i]);
#ifdef STAGGERED_MHD
    DIM_LOOP(n) fprintf(fp, " %.17g", d->Vs[n][k][j][i]);
#endif
    fputc('\n', fp);
    nwritten++;
  }
  fprintf(fp, "# bad_cells=%d written_cells=%d%s\n", nbad, nwritten,
          nbad > nwritten ? " (truncated at 1024)" : "");
  fflush(fp);
  fclose(fp);
}
#endif

#if FAILSAFE == YES
#if PHYSICS != RMHD
/* ********************************************************************* */
void FailSafe (Data *d, timeStep *Dts, Grid *grid)
/*!
 *
 * 
 *********************************************************************** */
{
  int  err;
  int  i,j,k,nv;
  int  nretry;

  int    n0  = g_stepNumber; 
  double dt0 = g_dt;
  double t0  = g_time;
 
  static double ***Bss0[3], ***Ess0[3];
  static uint16_t ***flag1;
  static Data_Arr Ucs0, Vcs0;
  
/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (Ucs0 == NULL){
    Ucs0  = ARRAY_4D(NX3_TOT, NX2_TOT, NX1_TOT, NVAR, double);
    Vcs0  = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);
    flag1 = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, uint16_t);
    #ifdef STAGGERED_MHD
    DIM_EXPAND(
      Bss0[IDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      Bss0[JDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      Bss0[KDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    )
    #endif
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_EXPAND(
      Ess0[IDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      Ess0[JDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);  ,
      Ess0[KDIR] = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    )
    #endif
  }

/* --------------------------------------------------------
   1. Backup solution arrays at t = t^n
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i) NVAR_LOOP(nv)   Ucs0[k][j][i][nv] = d->Uc[k][j][i][nv];
  NVAR_LOOP(nv)   TOT_LOOP(k,j,i) Vcs0[nv][k][j][i] = d->Vc[nv][k][j][i];
  #ifdef STAGGERED_MHD
  DIM_LOOP(nv) TOT_LOOP(k,j,i) Bss0[nv][k][j][i] = d->Vs[nv][k][j][i];
  #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
  DIM_LOOP(nv) TOT_LOOP(k,j,i) Ess0[nv][k][j][i] = d->Vs[EX1s+nv][k][j][i];
  #endif
  #endif

/* --------------------------------------------------------
   2. Try step, return to caller if err == 0
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i) d->flag[k][j][i] = 0;
  err = Integrate (d, Dts, grid);
  if (err == 0) return;  /* Step succeeded */

  // #if PHYSICS == RMHD
  // FailSafeDumpBadCells(d, grid, "native_failure");
  // #endif

/* --------------------------------------------------------
   3. Step did not succeed: enable Fail-safe procedure.
      Here 'err' is a globaly reduced variable already, so
      it has the same value for all processors.
   -------------------------------------------------------- */

  if (err > 0) {
    int invalid_zones=0;
    nretry = 1;
    uint16_t filter = 0;

    print("> FailSafe(): err = %d [g_intStage = %d] --> retrying step\n",
               err, g_intStage);
    filter |= FLAG_CONS2PRIM_FAIL;
    if (FlagCheck(d->flag, FLAG_CONS2PRIM_FAIL)) FlagShow(d, filter);

  /* -- 3a. Prepare new flags from the old failures --
   *
   * Build this map in two passes.  Clearing a cell while scanning the
   * failures would erase MINMOD flags written by an earlier failed cell.
   * The failed cell itself gets the TVDLF flux and flat reconstruction;
   * its one-cell Chebyshev neighbourhood gets MINMOD reconstruction.
   */
    TOT_LOOP(k,j,i){

      flag1[k][j][i] = 0;
      if ( (d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL) ) {
        invalid_zones++;
        printLog("  - Flagging zone (%d %d %d)\n",i,j,k);
        printLog("  - U0 = "); ShowState(Ucs0[k][j][i],0);
 
        int ib = i - 2*INCLUDE_IDIR, ie = i + 2*INCLUDE_IDIR;
        int jb = j - 2*INCLUDE_JDIR, je = j + 2*INCLUDE_JDIR;
        int kb = k - 2*INCLUDE_KDIR, ke = k + 2*INCLUDE_KDIR;
	
        ib = MAX(0, ib); ie = MIN(NX1_TOT-1,ie);
	       jb = MAX(0, jb); je = MIN(NX2_TOT-1,je);
	       kb = MAX(0, kb); ke = MIN(NX3_TOT-1,ke);
        int ii, jj, kk;

        flag1[k][j][i] |= FLAG_HLL;
        for (ii = ib; ii <= ie; ii++){
        for (jj = jb; jj <= je; jj++){
        for (kk = kb; kk <= ke; kk++){
          if ( abs(ii-i) <= 1 && abs(jj-j) <= 1 && abs(kk-k) <= 1 ) {
             flag1[k][j][i] |= FLAG_FLAT;
          }
          else flag1[k][j][i] |= FLAG_MINMOD;
       	}}}
      }
    }

  /* -- 3b. Overwrite d->flag with new values  -- */

    TOT_LOOP(k,j,i) d->flag[k][j][i] = flag1[k][j][i];

  /* -- 3c. Exchange data, Reduce and print info -- */

    #ifdef PARALLEL
    AL_Exchange (d->flag[0][0], SZ_uint16_t);
//    MPI_Allreduce (&invalid_zones, &nv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//    invalid_zones = nv;
    #endif
    if (invalid_zones > 0) printLog("  Total number of flagged zones = %d\n",invalid_zones);

  /* -- 3d. Restore solutionn array before repeating step -- */

    g_stepNumber = n0;
    g_time       = t0;
    g_dt         = dt0;

    print ("  --> step:%d; t = %10.4e; dt = %10.4e;\n",  g_stepNumber,
                                                         g_time, g_dt);
    LogFileFlush();

    TOT_LOOP(k,j,i) NVAR_LOOP(nv)   d->Uc[k][j][i][nv] = Ucs0[k][j][i][nv];
    NVAR_LOOP(nv)   TOT_LOOP(k,j,i) d->Vc[nv][k][j][i] = Vcs0[nv][k][j][i];
    #ifdef STAGGERED_MHD
    DIM_LOOP(nv) TOT_LOOP(k,j,i) d->Vs[nv][k][j][i] = Bss0[nv][k][j][i];
    #if (PHYSICS == ResRMHD) && (DIVE_CONTROL == CONSTRAINED_TRANSPORT)
    DIM_LOOP(nv) TOT_LOOP(k,j,i) d->Vs[EX1s+nv][k][j][i] = Ess0[nv][k][j][i];
    #endif
    #endif

  /* -- 3e. Repeat step -- */

    err = Integrate (d, Dts, grid);
    /* err is globally reduced: every rank finishes its dump before abort. */
    if (err) {
      // FailSafeDumpBadCells(d, grid, "failsafe_final");
#ifdef PARALLEL
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    if (err && FlagCheck(d->flag, FLAG_CONS2PRIM_FAIL)) {
      printLog ("  - Final flag:\n");
      FlagShow(d, filter);
      printLog ("  ! Step did not succeed\n");
      QUIT_PLUTO(1);
    }
  }
  
}
#else /* PHYSICS == RMHD */

#if RMHD_FAILSAFE_MAX_RETRIES < 0 || RMHD_FAILSAFE_MAX_RETRIES > 30
#error "RMHD_FAILSAFE_MAX_RETRIES must lie in [0,30]"
#endif
#if SHOCK_FLATTENING != MULTID
#error "RMHD FAILSAFE requires SHOCK_FLATTENING=MULTID for local order reduction"
#endif
#if PARTICLES != NO || FORCED_TURB != NO || defined(FARGO) || defined(HIGH_ORDER)
#error "RMHD FAILSAFE: auxiliary state rollback is not implemented for this configuration"
#endif

/* Accumulate fallback flags, keeping earlier troubled regions on later
 * retries. Exchange failure seeds BEFORE dilation so a failed cell next to
 * an MPI boundary also downgrades active cells on the neighbouring rank.
 * The failed cell selects HLL; radius one uses FLAT, radius two MINMOD.
 * Flags are applied to shared fluxes by the ordinary reconstruction/solver.
 */
void RMHD_LocalFallbackFlag(Data *d, uint16_t ***fallback)
{
  int i,j,k,ii,jj,kk;
#ifdef PARALLEL
  AL_Exchange(d->flag[0][0], SZ_uint16_t);
#endif
  TOT_LOOP(k,j,i) {
    if (!(d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL)) continue;
    fallback[k][j][i] |= FLAG_HLL;
    for (kk = MAX(0,k-2*INCLUDE_KDIR);
         kk <= MIN(NX3_TOT-1,k+2*INCLUDE_KDIR); kk++) {
    for (jj = MAX(0,j-2*INCLUDE_JDIR);
         jj <= MIN(NX2_TOT-1,j+2*INCLUDE_JDIR); jj++) {
    for (ii = MAX(0,i-2*INCLUDE_IDIR);
         ii <= MIN(NX1_TOT-1,i+2*INCLUDE_IDIR); ii++) {
      if (abs(ii-i) <= 1 && abs(jj-j) <= 1 && abs(kk-k) <= 1)
        fallback[kk][jj][ii] |= FLAG_FLAT;
      else
        fallback[kk][jj][ii] |= FLAG_MINMOD;
    }}}
  }
}

/* Normal step -> local fallback at dt0 -> local fallback at dt0/2^n.
 * Every trial starts from the same full-step fluid state. Only successful
 * trials commit U/V/B and Dts; g_dt remains the accepted step for main().
 * Integrate's stage failures must already be collective (as in RK C2P and
 * PCP); the final reduction also handles rank-local errors at its return.
 */
void FailSafe(Data *d, timeStep *Dts, Grid *grid)
{
  int i,j,k,n,err,attempt;
  long int n0 = g_stepNumber;
  int stage0 = g_intStage;
  double t0 = g_time, dt0 = g_dt, trial_dt;
  timeStep dts0 = *Dts;
  static Data_Arr U0, V0;
  static uint16_t ***flag0, ***fallback;
  static double *cmax0;
#ifdef STAGGERED_MHD
  static Data_Arr B0;
#endif
#if UPDATE_VECTOR_POTENTIAL == YES
  static Data_Arr A0;
  double ***A[3] = {d->Ax1, d->Ax2, d->Ax3};
#endif
#ifdef GLM_MHD
  double ch0 = glm_ch;
#endif

  if (U0 == NULL) {
    U0 = ARRAY_4D(NX3_TOT,NX2_TOT,NX1_TOT,NVAR,double);
    V0 = ARRAY_4D(NVAR,NX3_TOT,NX2_TOT,NX1_TOT,double);
    flag0 = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    fallback = ARRAY_3D(NX3_TOT,NX2_TOT,NX1_TOT,uint16_t);
    cmax0 = ARRAY_1D(NMAX_POINT,double);
#ifdef STAGGERED_MHD
    B0 = ARRAY_4D(DIMENSIONS,NX3_TOT,NX2_TOT,NX1_TOT,double);
#endif
#if UPDATE_VECTOR_POTENTIAL == YES
    A0 = ARRAY_4D(3,NX3_TOT,NX2_TOT,NX1_TOT,double);
#endif
  }
  TOT_LOOP(k,j,i) {
    flag0[k][j][i] = d->flag[k][j][i];
    fallback[k][j][i] = 0;
    NVAR_LOOP(n) {
      U0[k][j][i][n] = d->Uc[k][j][i][n];
      V0[n][k][j][i] = d->Vc[n][k][j][i];
    }
  }
  memcpy(cmax0,Dts->cmax,NMAX_POINT*sizeof(double));
#ifdef STAGGERED_MHD
  DIM_LOOP(n) TOT_LOOP(k,j,i) B0[n][k][j][i] = d->Vs[n][k][j][i];
#endif
#if UPDATE_VECTOR_POTENTIAL == YES
  for (n=0; n<3; n++) if (A[n] != NULL)
    TOT_LOOP(k,j,i) A0[n][k][j][i] = A[n][k][j][i];
#endif

  for (attempt=0; attempt<=(RMHD_FAILSAFE_MAX_RETRIES+1)*RMHD_FAILSAFE_LOCAL_PASSES; attempt++) {
    trial_dt = ldexp(dt0,-(MAX(0,attempt-1)/RMHD_FAILSAFE_LOCAL_PASSES));
    if (!(trial_dt > 0.0) || !isfinite(trial_dt) || t0+trial_dt == t0) {
      printLog("RMHD_FAILSAFE stop: time step cannot advance time, dt=%.4g\n",
               trial_dt);
      break;  /* Keep the last failed state for diagnostics. */
    }

    //expand bad region
    g_rmhd_failsafe_low_order = (attempt > 0);
    if (attempt > 0) {
      RMHD_LocalFallbackFlag(d,fallback);
      {
        long counts[2] = {0,0}, totals[2];
        DOM_LOOP(k,j,i) {
          counts[0] += (fallback[k][j][i] & FLAG_FLAT) != 0;
          counts[1]++;
        }
#ifdef PARALLEL
        MPI_Allreduce(counts, totals, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
        totals[0] = counts[0]; totals[1] = counts[1];
#endif
        // printLog("RMHD_FAILSAFE local mask: attempt=%d halvings=%d flat_cells=%ld total_cells=%ld\n",
        //          attempt, (attempt-1)/RMHD_FAILSAFE_LOCAL_PASSES, totals[0], totals[1]);
      }
      g_stepNumber = n0;
      g_time = t0;
      g_intStage = stage0;
      *Dts = dts0;
      memcpy(Dts->cmax,cmax0,NMAX_POINT*sizeof(double));
      TOT_LOOP(k,j,i) {
        d->flag[k][j][i] = flag0[k][j][i] | fallback[k][j][i];
        NVAR_LOOP(n) {
          d->Uc[k][j][i][n] = U0[k][j][i][n];
          d->Vc[n][k][j][i] = V0[n][k][j][i];
        }
      }
#ifdef STAGGERED_MHD
      DIM_LOOP(n) TOT_LOOP(k,j,i) d->Vs[n][k][j][i] = B0[n][k][j][i];
#endif
#if UPDATE_VECTOR_POTENTIAL == YES
      for (n=0; n<3; n++) if (A[n] != NULL)
        TOT_LOOP(k,j,i) A[n][k][j][i] = A0[n][k][j][i];
#endif
#ifdef GLM_MHD
      glm_ch = ch0;
#endif
#ifdef PARALLEL
      AL_Exchange(d->flag[0][0], SZ_uint16_t);
#endif
      // printLog("RMHD_FAILSAFE retry: attempt=%d local_fallback=YES "
      //          "halvings=%d dt=%.17g original_dt=%.17g\n",
      //          attempt,attempt-1,trial_dt,dt0);
    }
    g_dt = trial_dt;
    err = Integrate(d,Dts,grid);
#ifdef PARALLEL
    {
      int local_bad = (err != 0);
      MPI_Allreduce(&local_bad,&err,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    }
#endif
    if (err == 0) {
      if (attempt > 0)
        printLog("RMHD_FAILSAFE accepted: attempt=%d dt=%.4g\n",attempt,g_dt);
      g_rmhd_failsafe_low_order = 0;
      return;
    }
    // printLog("RMHD_FAILSAFE failed: attempt=%d stage=%d dt=%.17g err=%d\n",
    //          attempt,g_intStage,g_dt,err);
    // {
    //   char phase[48];
    //   snprintf(phase,sizeof(phase),"failsafe_attempt%02d",attempt);
    //   FailSafeDumpBadCells(d,grid,phase);
    // }
  }

  /* All ranks dump before any rank aborts, even if the failure has no C2P
   * flag (for example a rejected PCP stage or unusable ghost primitive). */
  //FailSafeDumpBadCells(d,grid,"failsafe_final");
  if (FlagCheck(d->flag,FLAG_CONS2PRIM_FAIL)) FlagShow(d,FLAG_CONS2PRIM_FAIL);
  printLog("! RMHD_FAILSAFE exhausted: step=%ld stage=%d original_dt=%.4g "
           "last_dt=%.4g max_halvings=%d\n",
           n0,g_intStage,dt0,g_dt,RMHD_FAILSAFE_MAX_RETRIES);
  LogFileFlush();
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  QUIT_PLUTO(1);
}
#endif  /* PHYSICS != RMHD */
#endif  /* FAILSAFE == YES */

/* ********************************************************************* */
void FlagShow(Data *d, uint16_t filter)
/*!
 *  Show which bits of the d->flag structure are set.
 *  The variable filter can be used to select only certain bits.
 *
 * \param [in] d        pointer to PLUTO data structure
 * \param [in] filter   a filter-pass variable used to select printing
 *                      of certain bits only.
 * 
 *********************************************************************** */
{
  int i,j,k;
  int nbits = 8*sizeof(uint16_t);
  int bit_position;
  uint16_t ***flag = d->flag;
  int isBitSet, isFilterSet;

  g_dir = IDIR;
  DOM_LOOP(k,j,i){
    g_j = j;
    g_k = k;

  /* ------------------------------------------------------
     1. Check bits to see which of them are on. 
        Enable show_flag only if they correspond to those
        specified by filter.
     ------------------------------------------------------ */
 
    int show_flag = 0;
    for (bit_position = 0; bit_position < nbits; bit_position++){
      isBitSet    = flag[k][j][i] & (1 << bit_position);
      isFilterSet = filter & (1 << bit_position);
      if (isBitSet && isFilterSet) show_flag = 1;
    }

    if (show_flag){
      printLog ("> FlagShow(): (%d, %d, %d)",i,j,k); Where(i,NULL);
/*
      printLog ("              active flags: ");
      if (flag[k][j][i] & FLAG_NEGATIVE_DENSITY) printLog ("/negative density"); 
      if (flag[k][j][i] & FLAG_NEGATIVE_ENERGY)  printLog ("/negative energy"); 
      if (flag[k][j][i] & FLAG_NEGATIVE_PRESSURE) printLog ("/negative pressure"); 
      if (flag[k][j][i] & FLAG_PRESSURE_FIX_FAIL) printLog ("/press. fix failed"); 
      if (flag[k][j][i] & FLAG_CONS2PRIM_FAIL) printLog ("/cons2prim fail"); 
      if (flag[k][j][i] & FLAG_HLL) printLog ("/hll"); 
      if (flag[k][j][i] & FLAG_MINMOD) printLog ("/minmod"); 
      if (flag[k][j][i] & FLAG_FLAT) printLog ("/flat"); 
      if (flag[k][j][i] & FLAG_SPLIT_CELL) printLog ("/splitCell");
      if (flag[k][j][i] & FLAG_ENTROPY) printLog ("/entropy"); 
*/

      printLog ("            ");
      printLog ("| MM  |");
      printLog (" FL  |");
      printLog (" HLL |");
      printLog (" Entr|");
      printLog (" SC  |");
      printLog (" IB  |");
      printLog (" C2P |");
      printLog (" nP  |");
      printLog (" nE  |");
      printLog (" nD  |");
      printLog ("prfix|");
      printLog ("\n            ");
      for (bit_position = 0; bit_position < nbits; bit_position++){
        printLog ("+-----");
      }
      printLog ("+\n            ");

      for (bit_position = 0; bit_position < nbits; bit_position++){
        isBitSet = flag[k][j][i] & (1 << bit_position);
        if (isBitSet) printLog ("|  *  ");
        else          printLog ("|     ");
      }
      printLog ("|\n");
    }
  }
}

/* ********************************************************************* */
int FlagCheck(uint16_t ***flag, int f)
/*
 *
 *
 *********************************************************************** */
{
  int i,j,k;

  DOM_LOOP(k,j,i){
    if (flag[k][j][i] & f) return 1;	   
  }
  return 0;
}
