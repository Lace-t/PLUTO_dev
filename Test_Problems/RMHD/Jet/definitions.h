#define  PHYSICS                        RMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define  MACH                           0
#define  LORENTZ                        1
#define  RHOJ                           2
#define  RHOA                           3
#define  SIGMA_POL                      4
#define  SIGMA_TOR                      5

/* [Beg] user-defined constants (do not change this line) */

#define  INCLUDE_JDIR                   NO
#define  SHOCK_FLATTENING               MULTID
#define  INTERNAL_BOUNDARY              NO
#define  FAILSAFE                       NO
#define  CT_EN_CORRECTION               YES
#define  INITIAL_SMOOTHING              NO
#define  WARNING_MESSAGES               NO
#define  PCP_NR                         NO

/* [End] user-defined constants (do not change this line) */
