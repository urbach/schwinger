/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Nils Christian,
 *               Carsten Urbach
 *
 * This file is part of a Schwinger code for the Helmholtz summer school
 * 2011 in Dubna
 *
 * this is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * this software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this code.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hmc.h"
#include "2MN_integrator.h"
#include "leapfrog.h"
#include "fields.h"
#include "dirac.h"

double gauge_force = 0;
double PF1_force = 0;
double PF2_force = 0;

//#define _DEBUG_

/* inner update defined below */
void MN2_inner_update(const int n2, const double tau, const double lambda);

/* Second order Minimal Norm Integration scheme     */
/* for two timescales, non-recursive implementation */

void MN2_integrator(const int n1, const int n2, const double tau, const double lambda) {
  int i;
  double dtau1 = tau/((double)n1);
  double oneminus2lambda = (1.-2.*lambda);

  /* here we need a halfstep, so 2*0.5*lambda*dtau */
  update_momenta_fermion(lambda*dtau1);

  for(i = 1; i < n1; i++){
    MN2_inner_update(n2, 0.5*dtau1, lambda);
    update_momenta_fermion(oneminus2lambda*dtau1);
    MN2_inner_update(n2, 0.5*dtau1, lambda);
    update_momenta_fermion(2*lambda*dtau1);
  }
  MN2_inner_update(n2, 0.5*dtau1, lambda);
  update_momenta_fermion(oneminus2lambda*dtau1);
  MN2_inner_update(n2, 0.5*dtau1, lambda);

  /* here we need a halfstep, so 2*0.5*lambda*dtau */
  update_momenta_fermion(lambda*dtau1);

  return;
}


void MN2_inner_update(const int n2, const double tau, const double lambda) {
  int j;
  double oneminus2lambda = (1.-2.*lambda);
  double dtau = tau/((double)n2);
  
  update_momenta_gauge(lambda*dtau);

  for(j = 1; j < n2; j++) {
    update_gauge(0.5*dtau);
    update_momenta_gauge(oneminus2lambda*dtau);
    update_gauge(0.5*dtau);
    update_momenta_gauge(2.*lambda*dtau);
  }
  update_gauge(0.5*dtau);
  update_momenta_gauge(oneminus2lambda*dtau);
  update_gauge(0.5*dtau);

  update_momenta_gauge(2*lambda*dtau);

  return;
}

void update_momenta_fermion(const double dtau) {
  int i;
  double f1=0., f2=0., sqsum = 0.;
  g_cgiterations1 += cg(g_X, g_fermion, ITER_MAX, DELTACG, &gam5D_SQR_wilson);
  gam5D_wilson(g_gam5DX, g_X);
#pragma omp parallel for private(f1,f2)
  for(i = 0; i < GRIDPOINTS; i++) {
    f1 = trX_dQ_wilson_dalpha1_X(i);
    f2 = trX_dQ_wilson_dalpha2_X(i);
#ifdef _DEBUG_
    sqsum = f1*f1 + f2*f2;
#endif
    gp1[i] = gp1[i] - dtau*(- f1);
    gp2[i] = gp2[i] - dtau*(- f2);
  }
#ifdef _DEBUG_
  printf("fermion momenta\t ||f|| = %e,\t dtau*||f|| = %e\n",
 	 sqrt(sqsum)/((double)GRIDPOINTS), dtau*sqrt(sqsum)/((double)GRIDPOINTS));
#endif
  return;
}

void update_momenta_PF2(const double dtau) {
  int i;
  double f1=0., f2=0.;
  double sqsum = 0.;
  g_cgiterations2 += cg(g_X, g_fermion2, ITER_MAX, DELTACG, &gam5D_SQR_wilson);
  gam5D_wilson(g_gam5DX, g_X);
#pragma omp parallel for private(f1,f2)
  for(i = 0; i < GRIDPOINTS; i++) {
    f1 = g_musqr*trX_dQ_wilson_dalpha1_X(i);
    f2 = g_musqr*trX_dQ_wilson_dalpha2_X(i);
#ifdef _DEBUG_
    sqsum = f1*f1 + f2*f2;
#endif
    gp1[i] = gp1[i] - dtau*(- f1);
    gp2[i] = gp2[i] - dtau*(- f2);
  }
  PF2_force += dtau * sqrt(sqsum);
#ifdef _DEBUG_
  printf("force PF2 \t ||f|| = %e,\t dtau*||f|| = %e\n", 
	  sqrt(sqsum)/((double)GRIDPOINTS), dtau*sqrt(sqsum)/((double)GRIDPOINTS));
#endif
  return;
}

void update_momenta_PF1(const double dtau) {
  int i;
  double f1=0., f2=0.;
  double sqsum = 0.;
  g_cgiterations1 += cg(g_X, g_fermion, ITER_MAX, DELTACG, &gam5D_SQR_musqr_wilson);
  gam5D_wilson(g_gam5DX, g_X);
#pragma omp parallel for private(f1,f2)
  for(i = 0; i < GRIDPOINTS; i++) {
    f1 = trX_dQ_wilson_dalpha1_X(i);
    f2 = trX_dQ_wilson_dalpha2_X(i);
#ifdef _DEBUG_
    sqsum = f1*f1 + f2*f2;
#endif
    gp1[i] = gp1[i] - dtau*(- f1);
    gp2[i] = gp2[i] - dtau*(- f2);
  }
  PF1_force += dtau * sqrt(sqsum);
#ifdef _DEBUG_
  printf("force PF1 \t ||f|| = %e,\t dtau*||f|| = %e\n", 
	  sqrt(sqsum)/((double)GRIDPOINTS), dtau*sqrt(sqsum)/((double)GRIDPOINTS));
#endif
  return;
}

void update_momenta_gauge(const double dtau) {
  int i;
  double f1=0., f2=0., sqsum = 0.;
#pragma omp parallel for private (f1,f2)
  for(i = 0; i < GRIDPOINTS; i++) {
    f1 = DS_G1(i);
    f2 = DS_G2(i);
    sqsum = f1*f1 + f2*f2;
    gp1[i] = gp1[i] - dtau*f1;
    gp2[i] = gp2[i] - dtau*f2;
  }
  gauge_force += dtau * sqrt(sqsum);
#ifdef _DEBUG_
  printf("gauge momenta\t ||f|| = %e,\t dtau*||f|| = %e\n", 
	 sqrt(sqsum)/((double)GRIDPOINTS), dtau*sqrt(sqsum)/((double)GRIDPOINTS));
#endif
  return;
}
