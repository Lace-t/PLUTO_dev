#include "pluto.h"
#include "globals.h"
#include <assert.h>
#include "../../../Src/MHD/CT/ct_emf_average.c"
int g_rmhd_failsafe_low_order;
void printLog(const char *fmt, ...) {}
void PlutoError(int bad, char *msg) { if (bad) abort(); }
int main(void)
{
  int i,j,k,dir,n,mask;
  double left[10],right[10];
  uint16_t flag[10];
  double ***q;
  NX1_TOT=NX2_TOT=NX3_TOT=NMAX_POINT=10;
  q=ARRAY_3D(10,10,10,double);
  for(k=0;k<10;k++) for(j=0;j<10;j++) for(i=0;i<10;i++)
    q[k][j][i]=100*k+10*j+i;
  for(mask=0;mask<2;mask++) for(dir=0;dir<3;dir++) {
    for(n=0;n<10;n++) {flag[n]=0;left[n]=12345+n;right[n]=23456+n;}
    flag[4]=mask ? FLAG_LF : FLAG_FLAT;
    g_rmhd_failsafe_low_order=0;
    CT_FailSafeEdgeStates(q,flag,2,2,2,dir,left,right);
    for(n=0;n<10;n++) {assert(left[n]==12345+n);assert(right[n]==23456+n);}
    g_rmhd_failsafe_low_order=1;
    CT_FailSafeEdgeStates(q,flag,2,2,2,dir,left,right);
    for(n=0;n<10;n++) {
      if(n==3 || n==4) {
        double a=dir==0?q[2][2][n]:(dir==1?q[2][n][2]:q[n][2][2]);
        double b=dir==0?q[2][2][n+1]:(dir==1?q[2][n+1][2]:q[n+1][2][2]);
        assert(left[n]==a);assert(right[n]==b);
      } else {assert(left[n]==12345+n);assert(right[n]==23456+n);}
    }
  }
  puts("Local CT: all directions, both edge sides, unmarked states unchanged PASS");
  return 0;
}
