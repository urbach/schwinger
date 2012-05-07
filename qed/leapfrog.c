#include <stdlib.h>
#include "hmc.h"
#include "leapfrog.h"
#include "dirac.h"
#include "fields.h"

/*  leap frog */
void leapfrog(const int nsteps, const double dtau) {
  int l;

  /* first phase: \Delta\Tau / 2 step for p */
  update_momenta(0.5*dtau); 

  /*  second phase: iterate with steps of \Delta\Tau */
  for(l = 0; l < nsteps-1; l++) {
    update_gauge(dtau);
    update_momenta(dtau);
  }
  /* a last one for the fields (because N steps for fields, */
  /*      and N-1 steps for impulses) */
  update_gauge(dtau);

  /*  last phase: \Delta\Tau / 2 step for p */
  update_momenta(dtau*0.5);
}

void update_momenta(const double dtau) 
{
  int i;
  g_cgiterations1 += cg(g_X, g_fermion, ITER_MAX, DELTACG, &gam5D_SQR_wilson);
  gam5D_wilson(g_gam5DX, g_X);
#ifdef OMP
#pragma omp parallel for
#endif
  for(i = 0; i < GRIDPOINTS; i++) {
    gp1[i] = gp1[i] - dtau*(DS_G1(i) - trX_dQ_wilson_dalpha1_X(i));
    gp2[i] = gp2[i] - dtau*(DS_G2(i) - trX_dQ_wilson_dalpha2_X(i));
  }
  return;
}

void update_gauge(const double dtau) {
  int i;
#ifdef OMP
#pragma omp parallel for
#endif
  for(i = 0; i < GRIDPOINTS; i++) {
    gauge1[i] = gauge1[i] + dtau*gp1[i];
    gauge2[i] = gauge2[i] + dtau*gp2[i];
  }
  calculatelinkvars();
  return;
}
