/* Exercise production face selection, constant states and LLF formula.
 * A linear test flux isolates locality/conservation from RMHD EOS recovery. */
#include "pluto.h"
#include "globals.h"
#include <assert.h>
int g_rmhd_failsafe_low_order;
static int ct_visits[12];
void PlutoError(int bad, char *msg) { if (bad) abort(); }
void PrimToConsLoc(double *v, double *u)
{ int n; NVAR_LOOP(n) u[n]=v[n]; }
void SoundSpeed2(const State *s,int beg,int end,int pos,Grid *g) {}
int MaxSignalSpeed(const State *s,double *lo,double *hi,int beg,int end)
{ int i; for(i=beg;i<=end;i++){lo[i]=-2;hi[i]=2;} return 0; }
void Flux(const State *s,int beg,int end)
{ int i,n; for(i=beg;i<=end;i++){NFLX_LOOP(n)s->flux[i][n]=s->u[i][n];s->prs[i]=s->v[i][PRS];} }
void CT_Flux(const Sweep *s,int beg,int end,Grid *g)
{ int i; for(i=beg;i<=end;i++)ct_visits[i]++; }
int main(void)
{
  Sweep s={0}; Grid grid={0}; double cmax[12];
  int i,n,seed,dir,enabled;
  NMAX_POINT=12;
  s.stateC.v=ARRAY_2D(12,NVAR,double);
  s.stateL.v=ARRAY_2D(12,NVAR,double);s.stateR.v=ARRAY_2D(12,NVAR,double);
  s.stateL.u=ARRAY_2D(12,NVAR,double);s.stateR.u=ARRAY_2D(12,NVAR,double);
  s.stateL.flux=ARRAY_2D(12,NVAR,double);s.stateR.flux=ARRAY_2D(12,NVAR,double);
  s.stateL.prs=ARRAY_1D(12,double);s.stateR.prs=ARRAY_1D(12,double);
  s.flux=ARRAY_2D(12,NVAR,double);s.press=ARRAY_1D(12,double);
  s.SL=ARRAY_1D(12,double);s.SR=ARRAY_1D(12,double);
  s.Bn=ARRAY_1D(12,double);s.flag=ARRAY_1D(12,uint16_t);
  for(dir=0;dir<3;dir++) for(seed=2;seed<=9;seed++) for(enabled=0;enabled<2;enabled++) {
    g_dir=dir;SetVectorIndices(dir);g_rmhd_failsafe_low_order=enabled;
    for(i=0;i<12;i++) {
      s.flag[i]=i==seed?FLAG_LF:0;ct_visits[i]=0;cmax[i]=99;s.Bn[i]=100+i;
      NVAR_LOOP(n) {
        s.stateC.v[i][n]=10*i+n;
        s.stateL.v[i][n]=s.stateR.v[i][n]=-999;
        s.flux[i][n]=1000+i+n;
      }
    }
    RMHD_FailSafe_CorrectFluxes(&s,2,8,cmax,&grid);
    for(i=0;i<12;i++) {
      int selected=enabled && i>=2 && i<=8 && (i==seed-1 || i==seed);
      assert(ct_visits[i]==selected);
      assert(cmax[i]==(selected?2:99));
      NVAR_LOOP(n) {
        if(selected) {
          double l=n==BXn?s.Bn[i]:s.stateC.v[i][n];
          double r=n==BXn?s.Bn[i]:s.stateC.v[i+1][n];
          assert(s.stateL.v[i][n]==l && s.stateR.v[i][n]==r);
          assert(s.flux[i][n]==0.5*(l+r-2*(r-l)));
        } else {
          assert(s.stateL.v[i][n]==-999 && s.stateR.v[i][n]==-999);
          assert(s.flux[i][n]==1000+i+n);
        }
      }
    }
    /* Shared face flux differences telescope, including corrected faces. */
    double change=0;
    for(i=3;i<=8;i++)change+=s.flux[i-1][RHO]-s.flux[i][RHO];
    assert(change==s.flux[2][RHO]-s.flux[8][RHO]);
  }
  puts("PASS FOFC flux: all directions, boundary seeds, both constant states, LLF, CT refresh, locality, conservation");
  return 0;
}
