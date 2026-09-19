/* Real RK driver and stage rollback, with controlled flux/CT/C2P failures.
 * The separate CT test exercises actual shared-edge locality. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>
static int calls[4], mode, target, maps;
void printLog(const char *fmt, ...) { va_list a; va_start(a,fmt); vprintf(fmt,a); va_end(a); }
void print(const char *fmt, ...) {(void)fmt;}
void PlutoError(int b,char *s) {if(b){puts(s);abort();}}
void LogFileFlush(void){}
void Boundary(const Data *d,int side,Grid *g){(void)d;(void)side;(void)g;}
void FlagShock(const Data *d,Grid *g){(void)d;(void)g;}
void PrimToCons3D(Data_Arr v,Data_Arr u,RBox *b,Grid *g){(void)v;(void)u;(void)b;(void)g;}
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

void UpdateStage(Data *d,Data_Arr u,Data_Arr b,double **f,double dt,timeStep *Dts,Grid *g)
{
 int i,j,k,n,pass=calls[g_intStage]++;
 (void)f;(void)g;
 assert(dt==0.25);
 assert(Dts->cmax[0]==(g_intStage==1?10:20));
 if(g_intStage==target && pass) {
   assert(g_rmhd_failsafe_low_order);
   if(prank==0) {
     assert(d->flag[5][5][8]&FLAG_HLL);
     assert(d->flag[5][5][7]&FLAG_FLAT);
     assert(!(d->flag[5][5][3]&FLAG_FLAT));
   }
#ifdef PARALLEL
   else { assert(d->flag[5][5][3]&FLAG_FLAT); assert(!(d->flag[5][5][5]&FLAG_FLAT)); }
#endif
 }
 /* Restore must recover entry primitives, not failed CT or mapper output. */
 assert(d->Vc[RHO][5][5][5]==(g_intStage==1?1:2));
 DOM_LOOP(k,j,i) {u[k][j][i][RHO]+=0.01; DIM_LOOP(n)b[n][k][j][i]+=0.125;}
 Dts->cmax[0]=20; Dts->invDt_hyp=30;
 if(g_intStage==target && prank==0 && (mode==2 || pass==0)) {
   u[5][5][8][MX1]=100; /* D/E positive: full candidate check must reject. */
 }
}
void CT_AverageStaggeredFields(const Data *d,int correct,RBox *b,Grid *g)
{
 (void)correct;(void)b;(void)g;
 /* Simulate CT internal mapper touching V even on a rejected candidate. */
 d->Vc[RHO][5][5][5]=-99;
}
int RMHD_PCP_ConsToPrim3D(Data_Arr u,Data_Arr v,uint16_t ***flags,RBox *b,Grid *g)
{
 int i,j,k,n,bad=0;
 (void)u;(void)b;(void)g;
 maps++;
 DOM_LOOP(k,j,i) NVAR_LOOP(n) v[n][k][j][i]=n==RHO?2:0;
 if(g_intStage==target && calls[g_intStage]==2 && mode==1 && prank==0) {
   flags[5][5][8]|=FLAG_CONS2PRIM_FAIL;bad=1;
 }
#ifdef PARALLEL
 int global;MPI_Allreduce(&bad,&global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);bad=global;
#endif
 return bad;
}
int main(int argc,char **argv)
{
 Data d={0};Grid grid={0};timeStep ts={0};int i,j,k,n;
#ifdef PARALLEL
 int size;MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&prank);MPI_Comm_size(MPI_COMM_WORLD,&size);assert(size==2);
#endif
 mode=argc>1?atoi(argv[1]):1;target=argc>2?atoi(argv[2]):2;
 NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=12;
 IBEG=JBEG=KBEG=3;IEND=JEND=KEND=8;g_dt=.25;g_gamma=5.0/3;
 d.Uc=ARRAY_4D(12,12,12,NVAR,double);d.Vc=ARRAY_4D(NVAR,12,12,12,double);
 d.Vs=ARRAY_4D(3,12,12,12,double);d.flag=ARRAY_3D(12,12,12,uint16_t);
 ts.cmax=ARRAY_1D(12,double);for(i=0;i<12;i++)ts.cmax[i]=10;
 TOT_LOOP(k,j,i){
   NVAR_LOOP(n){d.Uc[k][j][i][n]=n==RHO?1:n==ENG?3:0;d.Vc[n][k][j][i]=n==RHO?1:0;}
   DIM_LOOP(n)d.Vs[n][k][j][i]=0;d.flag[k][j][i]=0;
 }
 int err=AdvanceStep(&d,&ts,&grid);
 if(mode==2){assert(err);assert(calls[target]==PCP_STAGE_LOCAL_PASSES+1);}
 else {
   assert(!err);assert(calls[target]==(mode==1?3:2));
   assert(calls[target==1?2:1]==1);
   assert(fabs(d.Uc[5][5][5][RHO]-1.01)<1e-14);
   DOM_LOOP(k,j,i)assert(!(d.flag[k][j][i]&(FLAG_FLAT|FLAG_HLL|FLAG_MINMOD|FLAG_CONS2PRIM_FAIL)));
 }
 assert(g_rmhd_failsafe_low_order==0);
 printf("PASS local stage integrator=%d mode=%d target=%d calls=%d,%d,%d rank=%d\n",TIME_STEPPING,mode,target,calls[1],calls[2],calls[3],prank);
#ifdef PARALLEL
 MPI_Finalize();
#endif
 return 0;
}
