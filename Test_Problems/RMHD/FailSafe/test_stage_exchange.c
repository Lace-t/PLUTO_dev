/* Production AL_Exchange on a periodic 2x2x2 rank grid. Check face/edge/
 * corner masks, accumulation, cleanup and one exchange/reduction per retry. */
#include "pluto.h"
#include "globals.h"
/* PLUTO has a global int SZ; ArrayLib uses SZ as a private typedef. */
#define SZ ArrayDescriptor
#include "al_hidden.h"
#undef SZ
#include <assert.h>
#include <stdarg.h>
ArrayDescriptor *sz_stack[AL_MAX_ARRAYS];
int stack_ptr[AL_MAX_ARRAYS];
static int exchanges, reductions;
int __real_AL_Exchange(void *,int);
int __wrap_AL_Exchange(void *p,int s)
{exchanges++;return __real_AL_Exchange(p,s);}
int __wrap_MPI_Allreduce(const void *s,void *r,int n,MPI_Datatype t,MPI_Op op,MPI_Comm c)
{reductions++;return PMPI_Allreduce(s,r,n,t,op,c);}
void printLog(const char *fmt, ...)
{va_list a;va_start(a,fmt);if(prank==0)vprintf(fmt,a);va_end(a);}
void PlutoError(int bad,char *msg){if(bad)MPI_Abort(MPI_COMM_WORLD,2);}
void LogFileFlush(void){fflush(stdout);}
static int same(int x,int y,int z,int p){return x==p && y==p && z==p;}
int main(int argc,char **argv)
{
  Data d={0};timeStep ts={0};ArrayDescriptor descriptor={0};
  int size,coords[3],dims[3]={2,2,2},periods[3]={1,1,1};
  int i,j,k,n,pass;MPI_Comm cart;
  MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&prank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);assert(size==8);
  MPI_Cart_create(MPI_COMM_WORLD,3,dims,periods,0,&cart);
  MPI_Cart_coords(cart,prank,3,coords);
  descriptor.comm=cart;descriptor.ndim=3;descriptor.rank=prank;descriptor.size=size;
  SZ_uint16_t=0;sz_stack[0]=&descriptor;stack_ptr[0]=AL_STACK_USED;
  for(n=0;n<3;n++) {
    int sizes[3]={8,8,8},sub[3]={8,8,8},start[3]={0,0,0};
    sub[2-n]=2;
    MPI_Type_create_subarray(3,sizes,sub,start,MPI_ORDER_C,MPI_UNSIGNED_SHORT,&descriptor.type_rl[n]);
    MPI_Type_commit(&descriptor.type_rl[n]);descriptor.type_lr[n]=descriptor.type_rl[n];
    MPI_Cart_shift(cart,n,1,&descriptor.left[n],&descriptor.right[n]);
    descriptor.bg[n]=2;descriptor.tag1[n]=2*n;descriptor.tag2[n]=2*n+1;
    int stride=n==0?1:n==1?8:64;
    descriptor.sendb1[n]=2*stride*sizeof(uint16_t);
    descriptor.recvb1[n]=6*stride*sizeof(uint16_t);
    descriptor.sendb2[n]=4*stride*sizeof(uint16_t);descriptor.recvb2[n]=0;
  }
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=8;
  IBEG=JBEG=KBEG=2;IEND=JEND=KEND=5;
  d.Uc=ARRAY_4D(8,8,8,NVAR,double);d.Vc=ARRAY_4D(NVAR,8,8,8,double);
  d.Vs=ARRAY_4D(3,8,8,8,double);d.flag=ARRAY_3D(8,8,8,uint16_t);
  ts.cmax=ARRAY_1D(8,double);
  TOT_LOOP(k,j,i) {
    int x=(coords[0]*4+i+6)%8,y=(coords[1]*4+j+6)%8,z=(coords[2]*4+k+6)%8;
    d.flag[k][j][i]=FLAG_MINMOD|(same(x,y,z,0)?FLAG_LF:0);
  }
  g_rmhd_failsafe_low_order=1;
  RMHD_PCP_BeginStage(&d,&ts);
  /* No-retry cleanup preserves initial flags and newly added diagnostics. */
  d.flag[3][3][3]|=FLAG_NEGATIVE_PRESSURE;
  RMHD_PCP_EndStage(&d);
  assert(d.flag[3][3][3]&(FLAG_NEGATIVE_PRESSURE));
  assert(g_rmhd_failsafe_low_order==1);
  RMHD_PCP_BeginStage(&d,&ts);
  for(pass=0;pass<2;pass++) {
    if(pass==0 && prank==0)d.flag[5][5][5]|=FLAG_CONS2PRIM_FAIL;
    if(pass==1 && prank==7)d.flag[2][2][2]|=FLAG_CONS2PRIM_FAIL;
    exchanges=reductions=0;
    assert(RMHD_PCP_RetryStage(&d,&ts));
    assert(exchanges==1 && reductions==1);
    TOT_LOOP(k,j,i) {
      int x=(coords[0]*4+i+6)%8,y=(coords[1]*4+j+6)%8,z=(coords[2]*4+k+6)%8;
      assert(!!(d.flag[k][j][i]&FLAG_LF)==
             (same(x,y,z,0)||same(x,y,z,3)||(pass==1 && same(x,y,z,4))));
      assert(!(d.flag[k][j][i]&FLAG_CONS2PRIM_FAIL));
    }
  }
  RMHD_PCP_EndStage(&d);
  TOT_LOOP(k,j,i) {
    int x=(coords[0]*4+i+6)%8,y=(coords[1]*4+j+6)%8,z=(coords[2]*4+k+6)%8;
    assert(!!(d.flag[k][j][i]&FLAG_LF)==same(x,y,z,0));
    assert(d.flag[k][j][i]&FLAG_MINMOD);
  }
  assert(g_rmhd_failsafe_low_order==1);
  if(prank==0)puts("PASS 8-rank production halo: faces/edges/corners, cumulative mask, ghost cleanup, 1 exchange + 1 reduction/retry");
  for(n=0;n<3;n++)MPI_Type_free(&descriptor.type_rl[n]);
  MPI_Comm_free(&cart);MPI_Finalize();return 0;
}
