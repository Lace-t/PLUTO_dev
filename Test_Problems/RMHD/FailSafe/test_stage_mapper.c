/* Real 3D mapper and PCP policy; inject successful C2P repairs and failures. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>

static int mode;
void printLog(const char *fmt, ...)
{
  va_list args;
  va_start(args,fmt); vprintf(fmt,args); va_end(args);
}
void PlutoError(int bad, char *msg) { if (bad) { puts(msg); abort(); } }
void LogFileFlush(void) { fflush(stdout); }

int ConsToPrim(double **u, double **v, int beg, int end, uint16_t *flag)
{
  int i,n;
  for (i=beg; i<=end; i++) NVAR_LOOP(n) v[i][n]=123.0+n;
  if (prank != 0 || g_j != 3 || g_k != 3) return 0;
  if (mode==1) u[3][ENG]+=1.0;          /* valid repair */
  if (mode==2) u[3][ENG]=-1.0;         /* invalid repair, reported success */
  if (mode==3) u[3][BX1]=4.0;          /* Psi failure only */
  if (mode==4) u[3][MX2]=NAN;          /* nonfinite repair */
  if (mode==5) flag[3]|=FLAG_CONS2PRIM_FAIL; /* flag without return code */
  return mode==6;                     /* return code without flag */
}

int main(int argc, char **argv)
{
  Data d={0}; Grid grid={0}; RBox box;
  int i,j,k,n,err;
#ifdef PARALLEL
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&prank);
#else
  (void)argc; (void)argv;
#endif
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2; IEND=JEND=KEND=5;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);
  d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.Vs=ARRAY_4D(DIMENSIONS,8,8,8,double);
  d.flag=ARRAY_3D(8,8,8,uint16_t);
  RBoxDefine(IBEG,IEND,JBEG,JEND,KBEG,KEND,CENTER,&box);
  for (mode=0; mode<=6; mode++) {
    TOT_LOOP(k,j,i) {
      NVAR_LOOP(n) d.Uc[k][j][i][n]=0.0;
      d.Uc[k][j][i][RHO]=1.0;
      d.Uc[k][j][i][ENG]=3.0;
#if RMHD_REDUCED_ENERGY == YES
      d.Uc[k][j][i][ENG]-=1.0;
#endif
      DIM_LOOP(n) d.Vs[n][k][j][i]=0.0;
      d.flag[k][j][i]=0;
    }
    assert(RMHD_PCP_CheckStage(&d,&grid)==0);
    g_dir=JDIR;
    err=RMHD_PCP_ConsToPrim3D(d.Uc,d.Vc,d.flag,&box,&grid);
    int expected=mode>=2;
#if RMHD_FAILSAFE_PCP_STRICT == NO
    if (mode==3) expected=0;
#endif
    assert(err==expected);
    assert(g_dir==JDIR);
    DOM_LOOP(k,j,i) NVAR_LOOP(n) assert(d.Vc[n][k][j][i]==123.0+n);
    if (prank==0 && expected && mode!=6)
      assert(d.flag[3][3][3]&FLAG_CONS2PRIM_FAIL);
  }
  /* The ordinary mapper retains its original return-code policy. */
  mode=2;
  TOT_LOOP(k,j,i) d.flag[k][j][i]=0;
  assert(ConsToPrim3D(d.Uc,d.Vc,d.flag,&box,&grid)==0);
  puts("PASS stage mapper: unchanged, valid/invalid repair, Psi, NaN, flags, errors");
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
}
