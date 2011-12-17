#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "complex/complex.h"
#include "fields.h"
#include "lattice.h"
#include "linalg.h"
#include "hmc.h"
#include "dirac.h"
#include "measurements.h"
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

double mean_plaquette()
{
 int i;
 double mp = 0.0;
 for(i=0; i<GRIDPOINTS; i++)
 {
  mp += cos(gauge1[i] + gauge2[right1[i]] - gauge1[right2[i]] - gauge2[i]);
 }
 return mp/(double)GRIDPOINTS;
}

double polyakov_loop()
{
 //Here we calculate the Polyakov loop in X2 direction and average it over all x1
  int x1, x2;
  complex double pl;
  double apl = 0.0;
  calculatelinkvars();      
  for(x1=0; x1<X1; x1++) {
    pl = 1.0 + I*0.0;
    for(x2=0; x2<X2; x2++) {
      pl *= link2[idx(x1, x2, X1)];
    }
    apl += creal(pl); 
  };
  return apl/(double)X1;
}

spinor S0[GRIDPOINTS], S[GRIDPOINTS], R1[GRIDPOINTS], R2[GRIDPOINTS];
double chiral_condensate()
{
 complex double q1, q2;
 int i0 = 0;
 q1 = 0.0 + I*0.0;
 q2 = 0.0 + I*0.0;

 //To calculate D^{-1}(x,x), we invert solve the equation D R  = S
 //Where the source S is only nonzero at x and for different spinor components
 //The required number of sources is the number of spinor components
 // Source 1
 set_zero(S0);
 S0[i0].s1 = 1.0 + I*0.0;
 gam5D_wilson(S, S0);
 cg(R1, S, ITER_MAX, DELTACG, &gam5D_SQR_wilson); //Inverting the Dirac operator on source 1
 q1 += R1[i0].s1;
 // Source 2
 set_zero(S0);
 S0[i0].s2 = 1.0 + I*0.0;
 gam5D_wilson(S, S0);
 cg(R2, S, ITER_MAX, DELTACG, &gam5D_SQR_wilson); //Inverting the Dirac operator on source 2
 q2 += R2[i0].s2;
 
 if(fabs(cimag(q1 - q2))>sqrt(DELTACG))
 {
  printf("\n Imaginary part of chiral condensate detected!!! \n"); 
 };
 //q1 and q2 are the diagonal components (11 and 22) of the propagator
 //q1 - q2 is tr(gamma_5 D^-1)
 return creal(q1 - q2);
}

double pion_correlation_function(int t)
{
  double C = 0;
  for (int x = 0; x < X1; x ++)
  {
    int index = idx(x, t, X1);
    complex double a = R1[index].s1;
    complex double b = R1[index].s2;
    complex double c = R2[index].s1;
    complex double d = R2[index].s2;
    C += fabs(a * cconj(a) + b * cconj(b) + c * cconj(c) + d * cconj(d));
  }
  return C;
}

double pcac_correlation_function(int t)
{
  complex double C = 0;
  for (int x = 0; x < X1; x ++)
  {
    int index = idx(x, t, X1);
    complex double a = R1[index].s1;
    complex double b = R1[index].s2;
    complex double c = R2[index].s1;
    complex double d = R2[index].s2;
    //C += I * (a * cconj(c) - c * cconj(a) + b * cconj(d) - d * cconj(b));
    C += -cimag(a * cconj(c)) - cimag(b * cconj(d));
  }
  C *= 2;
  if (fabs(cimag(C)) > 1e-6)
    printf("\nImaginary part of PCAC correlation function detected!\n");
  return fabs(C);
}

double topological_charge()
{
  double tmp = 0;
  for (int i = 0; i < GRIDPOINTS; i ++)
    tmp += gauge1[i] + gauge2[right1[i]] - gauge1[right2[i]] - gauge2[i];
  return 0.5 * tmp / M_PI;
}
