/* Exercise production FailSafe with an integration stub that corrupts every
 * saved field on failure. No changes to the production integrator are needed.
 */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>

static int calls, success_call, no_flags;
static double initial_time;
static Data *test_data;

void printLog(const char *fmt, ...)
{
  va_list args;
  va_start(args,fmt); vfprintf(stdout,fmt,args); va_end(args);
}
void print(const char *fmt, ...) { (void)fmt; }
void LogFileFlush(void) { fflush(stdout); }
void Where(int i, Grid *grid) { (void)i; (void)grid; }
void PlutoError(int bad, char *msg) { if (bad) { puts(msg); abort(); } }

#ifdef PARALLEL
/* Two x-decomposed ranks, three ghost cells, six active cells. Transport
 * real MPI data, emulating the flag-only ArrayLib exchange used by FailSafe.
 */
int AL_Exchange(void *buffer, int size)
{
  uint16_t remote[12*12*12];
  uint16_t *local = buffer;
  int i,j,k;
  (void)size;
  MPI_Sendrecv(local,12*12*12,MPI_UNSIGNED_SHORT,1-prank,0,
               remote,12*12*12,MPI_UNSIGNED_SHORT,1-prank,0,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  for (k=0; k<12; k++) for (j=0; j<12; j++) for (i=0; i<3; i++) {
    if (prank == 0) local[(k*12+j)*12+9+i] = remote[(k*12+j)*12+3+i];
    else           local[(k*12+j)*12+i] = remote[(k*12+j)*12+6+i];
  }
  return 0;
}
#endif

int Integrate(Data *d, timeStep *Dts, Grid *grid)
{
  int i,j,k,n;
  double expected_dt = ldexp(1.0,-(MAX(0,calls-1)/RMHD_FAILSAFE_LOCAL_PASSES));
  (void)grid;
  assert(g_dt == expected_dt && g_time == initial_time);
  assert(g_rmhd_failsafe_low_order == (calls > 0));
  assert(g_stepNumber == 12345678901L && g_intStage == 7);
  assert(Dts->invDt_hyp == 12.0 && Dts->dt_cool == 13.0);
  for (i=0; i<NMAX_POINT; i++) assert(Dts->cmax[i] == 14.0+i);
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) {
      assert(d->Uc[k][j][i][n] == 1.0+n+0.001*(i+12*j+144*k));
      assert(d->Vc[n][k][j][i] == (n >= VX1 && n <= VX3 ? 0.0 : 2.0+n+0.002*(i+12*j+144*k)));
    }
    DIM_LOOP(n) assert(d->Vs[n][k][j][i] == 3.0+n);
    assert(d->Ax1[k][j][i] == 4.0 && d->Ax2[k][j][i] == 5.0);
    assert(d->Ax3[k][j][i] == 6.0);
    assert(!(d->flag[k][j][i] & FLAG_CONS2PRIM_FAIL));
    assert(d->flag[k][j][i] & FLAG_ENTROPY);
  }
  if (calls > 0 && !no_flags) {
    if (prank == 0) {
      assert(d->flag[5][5][8] & FLAG_HLL);
      assert(d->flag[6][6][7] & FLAG_FLAT);
      assert(d->flag[7][7][6] & FLAG_MINMOD);
      if (calls == 1)
        assert(!(d->flag[5][5][3] & (FLAG_HLL|FLAG_FLAT|FLAG_MINMOD)));
      if (calls >= 2) assert(d->flag[3][3][3] & FLAG_HLL);
    }
#ifdef PARALLEL
    else {
      assert(d->flag[5][5][3] & FLAG_FLAT);
      assert(d->flag[5][5][4] & FLAG_MINMOD);
      assert(!(d->flag[5][5][5] & (FLAG_HLL|FLAG_FLAT|FLAG_MINMOD)));
    }
#endif
  }
  /* Mutation happens on success too: successful values must be committed. */
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) { d->Uc[k][j][i][n] = -9.0; d->Vc[n][k][j][i] = -8.0; }
    DIM_LOOP(n) d->Vs[n][k][j][i] = -7.0;
    d->Ax1[k][j][i] = d->Ax2[k][j][i] = d->Ax3[k][j][i] = -6.0;
    d->flag[k][j][i] = 0; /* fallback must survive Integrate clearing flags */
  }
  for (i=0; i<NMAX_POINT; i++) Dts->cmax[i] = -5.0;
  Dts->invDt_hyp = -4.0; Dts->dt_cool = -3.0;
  g_intStage = 2;
  if (!no_flags && prank == 0 && calls != success_call) {
    if (calls == 0) d->flag[5][5][8] = FLAG_CONS2PRIM_FAIL;
    else           d->flag[3][3][3] = FLAG_CONS2PRIM_FAIL;
  }
  printLog("TEST_CALL %d dt=%.17g rank=%d\n",calls,g_dt,prank);
  calls++;
  /* In MPI only rank zero reports failure, testing collective retry/abort. */
  return prank == 0 && calls-1 != success_call;
}

static void CheckExit(void)
{
  if (success_call >= 0) return;
  assert(calls == (initial_time == 0.0 ? (RMHD_FAILSAFE_MAX_RETRIES+1)*RMHD_FAILSAFE_LOCAL_PASSES+1 : 0));
  if (calls) assert(test_data->Uc[5][5][5][RHO] == -9.0);
}

int main(int argc, char **argv)
{
  Data d = {0};
  Grid grid = {0};
  timeStep dts = {0};
  int i,j,k,n;
#ifdef PARALLEL
  int ranks;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&prank);
  MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  assert(ranks == 2);
#endif
  assert(argc == 2);
  success_call = atoi(argv[1]);
  no_flags = success_call == -2;
  initial_time = success_call == -3 ? 1.e30 : 0.0;
  test_data = &d;
  atexit(CheckExit);
  NX1_TOT = NX2_TOT = NX3_TOT = NMAX_POINT = 12;
  IBEG = JBEG = KBEG = 3; IEND = JEND = KEND = 8;
  d.Uc = ARRAY_4D(12,12,12,NVAR,double);
  d.Vc = ARRAY_4D(NVAR,12,12,12,double);
  d.Vs = ARRAY_4D(DIMENSIONS,12,12,12,double);
  d.flag = ARRAY_3D(12,12,12,uint16_t);
  d.Ax1 = ARRAY_3D(12,12,12,double);
  d.Ax2 = ARRAY_3D(12,12,12,double);
  d.Ax3 = ARRAY_3D(12,12,12,double);
  dts.cmax = ARRAY_1D(12,double);
  for (n=0; n<3; n++) {
    grid.x[n] = ARRAY_1D(12,double);
    for (i=0; i<12; i++) grid.x[n][i] = i;
  }
  for (i=0; i<12; i++) dts.cmax[i] = 14.0+i;
  dts.invDt_hyp = 12.0; dts.dt_cool = 13.0;
  TOT_LOOP(k,j,i) {
    NVAR_LOOP(n) {
      d.Uc[k][j][i][n] = 1.0+n+0.001*(i+12*j+144*k);
      d.Vc[n][k][j][i] = n >= VX1 && n <= VX3 ? 0.0 : 2.0+n+0.002*(i+12*j+144*k);
    }
    DIM_LOOP(n) d.Vs[n][k][j][i] = 3.0+n;
    d.Ax1[k][j][i] = 4.0; d.Ax2[k][j][i] = 5.0; d.Ax3[k][j][i] = 6.0;
    d.flag[k][j][i] = FLAG_ENTROPY;
  }
  g_dt = 1.0; g_time = initial_time; g_intStage = 7; g_stepNumber = 12345678901L;
  FailSafe(&d,&dts,&grid);
  assert(g_rmhd_failsafe_low_order == 0);
  assert(success_call >= 0 && calls == success_call+1);
  assert(g_dt == ldexp(1.0,-(MAX(0,success_call-1)/RMHD_FAILSAFE_LOCAL_PASSES)));
  assert(g_time == initial_time && g_stepNumber == 12345678901L);
  assert(d.Uc[5][5][5][RHO] == -9.0 && dts.invDt_hyp == -4.0);
  puts("PASS retry order, full rollback, accumulated flags, accepted state");
#ifdef PARALLEL
  MPI_Finalize();
#endif
  return 0;
}
