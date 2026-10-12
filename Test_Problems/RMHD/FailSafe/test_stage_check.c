/* Real admissibility checker + real stage gate. In MPI only rank 0 owns
 * the invalid state; every rank must reject without changing U or V. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>

#ifdef COUNT_PCP_CHECKS
static int pcp_checks;
int __real_RMHD_PCP_Check(const double *);
int __wrap_RMHD_PCP_Check(const double *u)
{ pcp_checks++; return __real_RMHD_PCP_Check(u); }
#endif

void printLog(const char *fmt, ...)
{
  va_list args;
  va_start(args,fmt); vfprintf(stdout,fmt,args); va_end(args);
}
void PlutoError(int bad, char *msg) { if (bad) { puts(msg); abort(); } }
void LogFileFlush(void) { fflush(stdout); }

static void SetEnergy(double *u, double E)
{
#if RMHD_REDUCED_ENERGY == YES
  u[ENG] = E-u[RHO];
#else
  u[ENG] = E;
#endif
}

int main(int argc, char **argv)
{
  Data d = {0};
  Grid grid = {0};
  int i,j,k,n,test;
  double saved[NVAR];
#ifdef PARALLEL
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&prank);
#else
  (void)argc; (void)argv;
#endif
  NX1_TOT = NX2_TOT = NX3_TOT = NMAX_POINT = 8;
  IBEG = JBEG = KBEG = 2; IEND = JEND = KEND = 5;
  d.Uc = ARRAY_4D(8,8,8,NVAR,double);
  d.Vc = ARRAY_4D(NVAR,8,8,8,double);
  d.Vs = ARRAY_4D(DIMENSIONS,8,8,8,double);
  d.flag = ARRAY_3D(8,8,8,uint16_t);
  for (test=0;test<6;test++) {
    TOT_LOOP(k,j,i) {
      NVAR_LOOP(n) { d.Uc[k][j][i][n]=0.0; d.Vc[n][k][j][i]=123.0+n; }
      d.Uc[k][j][i][RHO]=1.0;
      SetEnergy(d.Uc[k][j][i],3.0);
      DIM_LOOP(n) d.Vs[n][k][j][i]=0.0;
      d.flag[k][j][i]=0;
    }
    if (prank==0) {
      double *u=d.Uc[3][3][3];
      if (test==1) u[RHO]=-1.0;
      if (test==2) u[MX1]=4.0; /* q<0 despite positive D and E */
      if (test==3) u[BX1]=4.0; /* q>0 but energy below magnetic energy */
      if (test==4) u[MX2]=NAN;
      if (test==5) d.Vs[IDIR][3][3][2]=NAN; /* negative bounding face */
    }
    memcpy(saved,d.Uc[3][3][3],sizeof(saved));
    {
      int rejected = test!=0;
#ifdef COUNT_PCP_CHECKS
      pcp_checks=0;
#endif
      assert(RMHD_PCP_CheckStage(&d,&grid)==rejected);
#ifdef COUNT_PCP_CHECKS
      assert(pcp_checks==64-(test==4 && prank==0));
#endif
    assert(memcmp(saved,d.Uc[3][3][3],sizeof(saved))==0);
    TOT_LOOP(k,j,i) NVAR_LOOP(n) assert(d.Vc[n][k][j][i]==123.0+n);
    if (rejected && prank==0) assert(d.flag[3][3][3]&FLAG_CONS2PRIM_FAIL);
    else DOM_LOOP(k,j,i) assert(d.flag[k][j][i]==0);
    }
  }
  puts("PASS stage gate: valid, negative D, negative q, negative Psi, NaN, CT face; U/V unchanged");
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
}
