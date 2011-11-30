#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif
#include <complex.h>
#include "lattice.h"
#include "linalg.h"
#include "rand/ranlxd.h"
#include "fields.h"

spinor g_fermion[GRIDPOINTS], g_fermion2[GRIDPOINTS], g_X[GRIDPOINTS], g_gam5DX[GRIDPOINTS], g_temp[GRIDPOINTS];
double gauge1[GRIDPOINTS], gauge2[GRIDPOINTS], gauge1_old[GRIDPOINTS], gauge2_old[GRIDPOINTS];
complex double link1[GRIDPOINTS], link2[GRIDPOINTS];
double gp1[GRIDPOINTS], gp2[GRIDPOINTS];
double s_g, s_g_old;

// S_G = e^-2 * (1 - 0.5 * (U_P + U_P^dag)),
// U_P = exp(gauge1[i] + gauge2[right1[i]] - gauge1[right2[i]] - gauge2[i]) ?
double S_G(int i)
{
 return (-beta*cos(gauge1[i] + gauge2[right1[i]] - gauge1[right2[i]] - gauge2[i]));
}

double DS_G1(int i)
{
 return (beta*(-sin(gauge1[left2[i]]+gauge2[right1[left2[i]]]-gauge1[i]-gauge2[left2[i]])
		+sin(gauge1[i]+gauge2[right1[i]]-gauge1[right2[i]]-gauge2[i])));
}

double DS_G2(int i)
{
  return (beta*(sin(gauge1[left1[i]]+gauge2[i]-gauge1[left1[right2[i]]]-gauge2[left1[i]])
		-sin(gauge1[i]+gauge2[right1[i]]-gauge1[right2[i]]-gauge2[i])));
}

void coldstart()
{
 int i;
 for(i=0; i<GRIDPOINTS; i++)
 {
  gauge1[i]     = 0.0;
  gauge1_old[i] = 0.0;
  gauge2[i]     = 0.0;
  gauge2_old[i] = 0.0;
 };
 calculatelinkvars();
 s_g=0;
 for(i=0; i<GRIDPOINTS; i++)
  s_g += S_G(i);
 s_g_old = s_g;
}

void hotstart()
{
 int i;
 double r[GRIDPOINTS*2];
 ranlxd(r, GRIDPOINTS*2);
 for(i=0; i<GRIDPOINTS; i++)
 {
  gauge1[i]=2*M_PI*r[i]-M_PI;
  gauge2[i]=2*M_PI*r[i+GRIDPOINTS]-M_PI;
  gauge1_old[i]=gauge1[i];
  gauge2_old[i]=gauge2[i];
 }
 calculatelinkvars();
 s_g=0;
 for(i=0; i<GRIDPOINTS; i++)
  s_g += S_G(i);
 s_g_old=s_g;
}

int calculatelinkvars()
{
 int i;
 for(i=0; i<GRIDPOINTS; i++)
 {
  link1[i] = cos(gauge1[i]) + I*sin(gauge1[i]);
  link2[i] = cos(gauge2[i]) + I*sin(gauge2[i]);
 };
 return(0);
}
