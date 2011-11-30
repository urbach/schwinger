#include <stdlib.h>
#include "lattice.h"

/* index arrays for neighbours */
int * left1;
int * right1;
int * left2;
int * right2;

/* routine to initialise the geometry index routines */
/* takes care of the periodic boundary conditions    */
int init_lattice(const int NX, const int NY) {
  int i;
  div_t dg;
  int N = NX*NY;
  left1  = (int *)malloc(N*sizeof(int));
  left2  = (int *)malloc(N*sizeof(int));
  right1 = (int *)malloc(N*sizeof(int));
  right2 = (int *)malloc(N*sizeof(int));

  left1[0]=NX-1;
  right1[0]=1;
  left2[0]=N-1;
  right2[0]=NX;

  for(i=0; i<N; i++) {
    /* div from stdlib for integer division                                  */
    /* div_t has two members, quot for the quotient and rem for the reminder */
    dg=div(i, NX);
    if(dg.rem==0) {
      left1[i]=i+NX-1;
      right1[i]=i+1;
    }
    else {
      if(dg.rem==(NX-1)) {
	left1[i]=i-1;
	right1[i]=i-NX+1;
      }
      else {
	left1[i]=i-1;
	right1[i]=i+1;
      }
    } 

    left2[i]=i-NX;
    if(left2[i]<0) {
      left2[i] += N;
    }

    right2[i]=i+NX;
    if(right2[i]>=N) {
      right2[i] -= N;
    }
  }
  return(0);
}

int idx(const int x1, const int x2, const int NX) {
  return x1 + NX*x2;    
}

