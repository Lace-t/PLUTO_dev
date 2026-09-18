/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Module header file for relativistic MHD (RMHD).

  Set label, indexes and basic prototyping for the relativistic 
  MHD module.

  \authors A. Mignone (andrea.mignone@unito.it)\n
           C. Zanni   (zanni@oato.inaf.it)\n
           G. Mattia
  \date   Jul 03, 2021
*/
/* ///////////////////////////////////////////////////////////////////// */

/* *********************************************************
    Set flow variable indices.
    Extra vector components, when not needed, point to the
    last element (255) of the array stored by startup.c.  
   ********************************************************* */

#define  RHO 0
#define  MX1 1
#define  MX2 2
#define  MX3 3
#define  BX1 4
#define  BX2 5
#define  BX3 6

#if HAVE_ENERGY
  #define ENG  7
  #define PRS  ENG
#endif

#if DIVB_CONTROL == DIV_CLEANING
    #define PSI_GLM  (ENG + 1)
    #define DIV_COMP 1
#else
  #define DIV_COMP 0
#endif

#define VX1   MX1
#define VX2   MX2
#define VX3   MX3

#if RADIATION
  #define ENR (7 + DIV_COMP + HAVE_ENERGY)
  #define FR1 (ENR + 1)
  #define FR2 (ENR + 2)
  #define FR3 (ENR + 3)

  #include "../Radiation/radiation.h"

  #define NFLX (11 + DIV_COMP + HAVE_ENERGY + IRRADIATION)
#else
  #define NFLX (7 + DIV_COMP + HAVE_ENERGY)
#endif

/* *********************************************************
    Label the different waves in increasing order 
    following the number of vector components.

    IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
               present only in the MHD-GLM formulation.
               We keep them at the END of the enumeration
               so we can skip them in unnecessary loops.
               Please do NOT change them !
   ********************************************************* */

enum KWAVES {
 KFASTM, KFASTP, KENTRP

 #if DIVB_CONTROL != DIV_CLEANING
  , KDIVB
 #endif

  , KSLOWM, KSLOWP, KALFVM, KALFVP

 #if DIVB_CONTROL == DIV_CLEANING  
  , KPSI_GLMM, KPSI_GLMP 
 #endif
};

/* ******************************************************
     Vector potential: these labels are and MUST only
     be used in the STARTUP / INIT  functions;
     they're convenient in obtaining a discretization 
     that preserve divB since the beginning.   
   ****************************************************** */

#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)

#define AX  AX1  /* backward compatibility */
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL
                               
 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3
                               
 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBZ    BX2
 #define iBPHI  BX3
                               
#endif

#if GEOMETRY == POLAR
                               
 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3
                               
 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3
                               
 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3
                               
#endif

#if GEOMETRY == SPHERICAL
                    
 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3
                             
 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3
                               
 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3

#endif

/* ******************************************************************
    Module-specific symbolic constants (switches)
   ****************************************************************** */

#ifndef COUNT_FAILURES  
  #define COUNT_FAILURES  NO /**< When set to YES, count number of failures and 
          write the count to "hlld_fails.dat" (works in parallel as well).
          The number of failures is normalized to the total number of zones: 
          \c nf/nz where \c nz is incremented during the main loop and \c nf
          is incremented when a failure occur. */
#endif

#ifndef RMHD_FAST_EIGENVALUES
  #define RMHD_FAST_EIGENVALUES  NO /**< If set to YES, use approximate (and
                                         faster) expressions when computing the
                                         fast magnetosonic speed, see Sect. 3.3
                                         of Del Zanna,  A&A (2006), 473.
                                         Solutions of quartic equation is avoided
                                         and replace with upper bounds provided by
                                         quadratic equation. */
#endif  

#ifndef RMHD_REDUCED_ENERGY
  #define RMHD_REDUCED_ENERGY    YES  /**< By turning RMHD_REDUCED_ENERGY to YES, 
                                            we let PLUTO evolve the total energy 
                                            minus the mass density contribution. */
#endif 

/* Enable the optional admissibility-preserving reconstruction and the
   physical-constraint-preserving primitive recovery.  The default remains
   the historical PLUTO RMHD pipeline. */
#ifndef PCP_NR
  #define PCP_NR  NO
#endif

/* Opt-in, conservative three-dimensional a-posteriori stage limiter.
   This is not an assertion that the native CT update satisfies Wu's
   a-priori LxF theorem. Unsupported physics is rejected in pcp_stage.c. */
#ifndef PCP_STAGE_LIMITING
  #define PCP_STAGE_LIMITING NO
#endif

#if (PCP_NR == YES)
	#ifndef RMHD_FAILSAFE_PCP_CHECK
	  /* Check RK cell averages before C2P can repair them.
	     Independent of PCP_NR and of the optional stage blending wrapper. */
	  #define RMHD_FAILSAFE_PCP_CHECK YES
	#endif
#endif

#ifndef RMHD_FAILSAFE_PCP_STRICT
  /* NO retains the existing pressure-recovery policy for states passing
     the hard D/energy/q checks but failing the full PCP check. YES rejects
     every failed PCP check; cold production restarts may require a more
     substantial change of the spatial update to run in this mode. */
  #define RMHD_FAILSAFE_PCP_STRICT NO
#endif
#ifndef RMHD_FAILSAFE_LOCAL_PASSES
  /* Accumulate newly discovered bad regions at fixed dt before reducing it. */
  #define RMHD_FAILSAFE_LOCAL_PASSES 2
#endif
#if RMHD_FAILSAFE_LOCAL_PASSES < 1
  #error "RMHD_FAILSAFE_LOCAL_PASSES must be positive"
#endif
#ifndef RMHD_FAILSAFE_MAX_RETRIES
  /* Number of dt halvings AFTER local mask passes at the original dt.
     Independent of PCP_STAGE_LIMITING; zero disables dt reduction.
     Accept explicitly configured legacy retry counts for compatibility. */
  #ifdef PCP_STAGE_MAX_RETRIES
    #define RMHD_FAILSAFE_MAX_RETRIES PCP_STAGE_MAX_RETRIES
  #else
    #define RMHD_FAILSAFE_MAX_RETRIES 1
  #endif
#endif

#if PCP_NR == YES
  /* The first PCP implementation is proved and exercised only for the SSP
     convex-combination stages of RK2 and RK3. */
  #if (TIME_STEPPING != RK2) && (TIME_STEPPING != RK3)
    #error "PCP_NR requires TIME_STEPPING = RK2 or RK3"
  #endif
#endif

                                      

/* ---- Function prototyping ----  */

int  ApproximateFastWaves  (double *, double, double, double *);
int  ConsToPrim   (double **, double **, int, int, uint16_t *);
void ConvertTo4vel (double **, int, int);
void ConvertTo3vel (double **, int, int);
void PrimEigenvectors (double *, double, double, double *, double **, double **);
 
void Flux      (const State *, int, int);
void HLL_Speed (const State *, const State *, double *, double *, int, int);
int  MaxSignalSpeed (const State *, double *, double *, int, int);
void RMHD_LF_Flux (const Sweep *, int, double *, Grid *);

void PrimToCons   (double **, double **, int, int);
void VelocityLimiter (double *, double *, double *);

int  Magnetosonic (double *, double, double, double *);

void RiemannCheck (double *, double *s);

Riemann_Solver LF_Solver, HLL_Solver, HLLC_MB_Solver, HLLC_KB_Solver;
Riemann_Solver HLLD_Solver, HLLEM_Solver;
Riemann_Solver HLL_Linde_Solver, GMUSTA1_Solver; 
Riemann_Solver GFORCE_Solver;

int  RMHD_EnergySolve  (double *, double *);
int  RMHD_EntropySolve (double *, double *);
int  RMHD_PressureFix  (double *, double *);

/* Wu admissibility/limiter helpers also serve PCP_STAGE_LIMITING when the
   optional Cai PCP-NR primitive recovery is disabled. */
int  RMHD_PCP_Check (const double *);
int  RMHD_PCP_Limit (const double *, double *);
int  RMHD_PCP_NR    (const double *, double *);
void RMHD_PCP_SetFallback (const double *, double *);
void RMHD_PCP_LimitStates (const Sweep *, int, int);
extern int g_pcp_low_order;
/* Set only while the RMHD full-step failsafe is rebuilding a failed step.
 * Unlike the local HLL/FLAT flags, this selects a complete first-order
 * FV+CT operator: piecewise-constant shared-normal-field states, LF fluxes,
 * and first-order edge-EMF reconstruction. */
extern int g_rmhd_failsafe_low_order;
void RMHD_FailSafe_StageStates (const Sweep *, int, int);
int RMHD_PCP_UpdateStage (Data *, double, timeStep *, Grid *);
int RMHD_PCP_CheckStage (Data *, Grid *);
void RMHD_PCP_StageStates (const Sweep *, int, int);

#if DIVB_CONTROL == EIGHT_WAVES
 void POWELL_DIVB_SOURCE(const Sweep *, int, int, Grid *);
 void HLL_DIVB_SOURCE (const Sweep *, double **, int, int, Grid *);
#elif DIVB_CONTROL == DIV_CLEANING
 #include "MHD/GLM/glm.h"
#elif DIVB_CONTROL == CONSTRAINED_TRANSPORT
  #include "MHD/CT/ct.h"
#endif

#ifndef NEW_RMHD_FLUXES
  #define NEW_RMHD_FLUXES  NO
#endif
