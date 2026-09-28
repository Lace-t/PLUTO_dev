/* Differential oracle: unchanged legacy mapper vs split mapper, including
 * per-cell recovery order, exact U/V/flags, and the 3D changed-state gate. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include <stdarg.h>

#define CELLS 128
static double **tracked;
static int order[CELLS], counts[3], force_errors;
void printLog(const char *fmt, ...) { (void)fmt; }
void print(const char *fmt, ...) { (void)fmt; }
void Where(int i,Grid *g) { (void)i; (void)g; }
void LogFileFlush(void) {}
void PlutoError(int bad,char *msg) { if(bad) { puts(msg); abort(); } }
static int Record(const double *u,int method)
{
  int i;
  counts[method-1]++;
  if(tracked) for(i=0;i<CELLS;i++) if(u==tracked[i]) {
    order[i]=10*order[i]+method;
    return i;
  }
  return -1;
}
int __real_RMHD_EnergySolve(double *,double *);
int __wrap_RMHD_EnergySolve(double *u,double *v)
{
  int i=Record(u,1);
  if(force_errors && i>=0 && i%4==0) return 3;
  return __real_RMHD_EnergySolve(u,v);
}
int __real_RMHD_PCP_NR(const double *,double *);
int __wrap_RMHD_PCP_NR(const double *u,double *v)
{
  int i=Record(u,2);
  if(force_errors && i>=0 && i%8==0) return 1;
  return __real_RMHD_PCP_NR(u,v);
}
int __real_RMHD_PressureFix(double *,double *);
int __wrap_RMHD_PressureFix(double *u,double *v)
{
  int i=Record(u,3);
  if(force_errors && i>=0 && i%16==0) return 1;
  return __real_RMHD_PressureFix(u,v);
}

static void Fixture(double *u,int cell)
{
  double v[NVAR]={0}; int kind=cell%16;
  v[RHO]=0.5+cell*0.01; v[PRS]=0.1+cell*0.001;
  v[VX1]=0.5*sin(cell); v[VX2]=0.2; v[VX3]=-0.1;
  v[BX1]=0.3+cell*0.02; v[BX2]=-0.2; v[BX3]=0.5;
#if NTRACER > 0
  v[TRC]=0.37;
#endif
#ifdef GLM_MHD
  v[PSI_GLM]=0.01*cell;
#endif
  /* These positive sub-floor states retain the historical native policy. */
  if(kind==13) { v[PRS]=0.1*g_smallPressure; v[BX1]=0; }
  if(kind==14) v[RHO]=0.1*g_smallDensity;
  PrimToConsLoc(v,u);
  if(kind==8) u[ENG]=-1;
  if(kind==9) u[RHO]=-1;
  if(kind==10) u[MX1]=1.e8;
  if(kind==11) u[ENG]=NAN;
  if(kind==12) u[RHO]=0;
  if(kind==15) u[MX2]=NAN;
}

int main(void)
{
  double **a,**b,**va,**vb,initial[CELLS][NVAR];
  uint16_t fa[CELLS],fb[CELLS]; int changed[CELLS],nchanged,seen[CELLS];
  int i,n,mode,expected_order[CELLS],expected_counts[3],ea,eb;
  NX1_TOT=NX2_TOT=NX3_TOT=8; NMAX_POINT=CELLS;
  IBEG=JBEG=KBEG=2; IEND=JEND=KEND=5;
  g_gamma=5.0/3; g_smallDensity=1.e-10; g_smallPressure=1.e-8;
  a=ARRAY_2D(CELLS,NVAR,double); b=ARRAY_2D(CELLS,NVAR,double);
  va=ARRAY_2D(CELLS,NVAR,double); vb=ARRAY_2D(CELLS,NVAR,double);
  for(mode=0;mode<2;mode++) {
    force_errors=mode;
    for(i=0;i<CELLS;i++) {
      Fixture(initial[i],i);
      memcpy(a[i],initial[i],sizeof initial[i]); memcpy(b[i],initial[i],sizeof initial[i]);
      for(n=0;n<NVAR;n++)va[i][n]=vb[i][n]=0.001*(n+1);
      fa[i]=fb[i]=FLAG_MINMOD | (i%7==0 ? FLAG_PCP_NR_FAIL:0)
                            | (i%11==0 ? FLAG_CONS2PRIM_FAIL:0);
    }
    tracked=a; memset(order,0,sizeof order); memset(counts,0,sizeof counts);
    ea=ConsToPrim(a,va,0,CELLS-1,fa);
    memcpy(expected_order,order,sizeof order); memcpy(expected_counts,counts,sizeof counts);
    tracked=b; memset(order,0,sizeof order); memset(counts,0,sizeof counts);
    eb=RMHD_PCP_ConsToPrimRow(b,vb,0,CELLS-1,fb,changed,&nchanged);
    assert(ea==eb && !memcmp(order,expected_order,sizeof order));
    assert(!memcmp(counts,expected_counts,sizeof counts));
    assert(!memcmp(fa,fb,sizeof fa));
    memset(seen,0,sizeof seen);
    for(i=0;i<nchanged;i++) { assert(!seen[changed[i]]); seen[changed[i]]=1; }
    for(i=0;i<CELLS;i++) {
      assert(!memcmp(a[i],b[i],NVAR*sizeof(double)));
      assert(!memcmp(va[i],vb[i],NVAR*sizeof(double)));
      assert(seen[i]==(memcmp(initial[i],b[i],NVAR*sizeof(double))!=0));
    }
  }
  /* All-normal row: no repair or changed-state work. */
  force_errors=0; tracked=b; memset(counts,0,sizeof counts);
  for(i=0;i<CELLS;i++) { Fixture(b[i],0); fb[i]=0; }
  assert(!RMHD_PCP_ConsToPrimRow(b,vb,0,CELLS-1,fb,changed,&nchanged));
  assert(nchanged==0 && counts[0]==CELLS && counts[1]==0 && counts[2]==0);

  /* Real 3D wrapper, independent legacy row + post-modification PCP oracle. */
  {
    Data d={0}; Grid g={0}; RBox box; int j,k,row=0,expected_error=0;
    d.Uc=ARRAY_4D(8,8,8,NVAR,double); d.Vc=ARRAY_4D(NVAR,8,8,8,double);
    d.flag=ARRAY_3D(8,8,8,uint16_t);
    tracked=NULL; force_errors=0;
    RBoxDefine(2,5,2,5,2,5,CENTER,&box);
    /* Prime the mapper's scratch V with valid values before failure tests. */
    DOM_LOOP(k,j,i) { Fixture(d.Uc[k][j][i],0); d.flag[k][j][i]=0; }
    assert(!RMHD_PCP_ConsToPrim3D_CheckModified(d.Uc,d.Vc,d.flag,&box,&g));
    for(i=0;i<CELLS;i++)Fixture(b[i],0);
    ConsToPrim(b,vb,0,CELLS-1,fb);
    DOM_LOOP(k,j,i) {
      int p=row++;
      Fixture(a[p],p); memcpy(initial[p],a[p],NVAR*sizeof(double));
      memcpy(d.Uc[k][j][i],a[p],NVAR*sizeof(double));
      /* Row scratch slots carry previous successful recovery, as in production. */
      memcpy(va[p],vb[0],NVAR*sizeof(double));
      fa[p]=d.flag[k][j][i]=0;
    }
    /* Invalid failures can leave partial V, which is intentionally unspecified
     * across scratch reuse. Compare flags/U on all cells and V on successes. */
    expected_error=ConsToPrim(a,va,0,row-1,fa);
    for(i=0;i<row;i++) {
      if(memcmp(initial[i],a[i],NVAR*sizeof(double)) && !RMHD_PCP_CheckConservative(a[i])) {
        fa[i]|=FLAG_CONS2PRIM_FAIL;
        if(a[i][RHO]<=0)fa[i]|=FLAG_NEGATIVE_DENSITY;
        if(a[i][ENG]<=0)fa[i]|=FLAG_NEGATIVE_ENERGY;
      }
      if(fa[i]&FLAG_CONS2PRIM_FAIL)expected_error=1;
    }
    g_dir=JDIR;
    assert(RMHD_PCP_ConsToPrim3D_CheckModified(d.Uc,d.Vc,d.flag,&box,&g)==expected_error);
    assert(g_dir==JDIR); row=0;
    DOM_LOOP(k,j,i) {
      int p=row++;
      assert(d.flag[k][j][i]==fa[p]);
      assert(!memcmp(a[p],d.Uc[k][j][i],NVAR*sizeof(double)));
      if(!(fa[p]&FLAG_CONS2PRIM_FAIL)) for(n=0;n<NVAR;n++)
        assert(!memcmp(&va[p][n],&d.Vc[n][k][j][i],sizeof(double)));
    }
  }
  puts("PASS split recovery: exact row U/V/flags, recovery order/counts, sparse modifications, 3D PCP oracle");
  return 0;
}
