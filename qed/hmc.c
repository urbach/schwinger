/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2011 Nils Christian,
 *               Pavel Buividovic, Carsten Urbach
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
#include <complex.h>
#include "rand/ranlxd.h"
#include "rand/gauss.h"
#include "linalg.h"
#include "fields.h"
#include "lattice.h"
#include "dirac.h"
#include "2MN_integrator.h"
#include "leapfrog.h"
#include "leapfrog2.h"
#include "rec_lf_integrator.h"
#include "hmc.h"


int R;
int g_cgiterations1;
int g_cgiterations2;

int update() //Basic HMC update step
{
  double squnrm;
  int i, acc;
  double exphdiff;
  
  /* the new impulses and the 'generator' of the arbitrary pseudofield */
  /* calculate the hamiltonian of this state: new impulses + action */
  /* g_X is ab-used a bit - here it is \xi = (gamma5 D)^{-1} \phi */
  
  ham_old = s_g_old;
  for(i=0; i<GRIDPOINTS; i++) {
    gp1[i] = gauss();
    gp2[i] = gauss();
    ham_old += 0.5*(gp1[i]*gp1[i] + gp2[i]*gp2[i]);
  }
  
  /* Now create the field and calculate its contributions to the action (end of the 'misuse') */
  /* squnrm is the fermion part of the action : */
  /*   S = R^dagger * R  =  g_fermion^dag * D^{-1 dag} * D^{-1} * g_fermion = g_fermion Q^-1 g_fermion */

  /* PF1 det(1/(Q^2 + mu^2)) */
  for(i=0; i<GRIDPOINTS; i++) {
    g_X[i].s1 = (gauss() + I*gauss())/sqrt(2); //Gaussian fields R
    g_X[i].s2 = (gauss() + I*gauss())/sqrt(2);
  }
  squnrm = square_norm(g_X);
  
  // step iv): g_fermion = \phi = K^dag * g_X = K^dag * \xi
  gam5D_wilson(g_fermion, g_X);
  assign_diff_mul(g_fermion, g_X, 0.+I*sqrt(g_musqr));
  ham_old += squnrm;

  /* PF2 det((Q^2 + mu^2)/Q^2) */
  if(no_timescales > 2) {
    for(i=0; i<GRIDPOINTS; i++) {
      g_X[i].s1 = (gauss() + I*gauss())/sqrt(2); //Gaussian fields R
      g_X[i].s2 = (gauss() + I*gauss())/sqrt(2);
    }
    squnrm = square_norm(g_X);

    cg(g_fermion2, g_X, ITER_MAX, DELTACG, &gam5D_SQR_musqr_wilson);    
    gam5D_wilson(g_gam5DX, g_fermion2);
    assign_add_mul(g_gam5DX, g_fermion2, 0.+I*sqrt(g_musqr));
    gam5D_wilson(g_fermion2, g_gam5DX);
    ham_old += squnrm;
  }
  // Add the part for the fermion fields

  // Do the molecular dynamic chain
  /* the simple LF scheme */

  /* the second order minimal norm multi-timescale integrator*/
  /* MN2_integrator(g_steps, 2, g_steps*g_stepsize, 0.2); */

  /* This is the recursive implementation */
  /* in can be found in rec_lf_integrator.c|h */
  if (no_timescales == 1)
    leapfrog(n_steps[0], tau/n_steps[0]);
  else
    integrate_leap_frog(tau/n_steps[no_timescales-1], no_timescales-1, no_timescales, n_steps, 1, up_momenta);
  
  // Calculate the new action and hamiltonian
  ham = 0;
  s_g = 0;
  for (i=0; i<GRIDPOINTS; i++) {
    s_g += S_G(i);
    ham += 0.5*(gp1[i]*gp1[i] + gp2[i]*gp2[i]);
  }
  /* Sum_ij [(g_fermion^*)_i (Q^-1)_ij (g_fermion)_j]  =  Sum_ij [(g_fermion^*)_i (g_X)_i] */
  ham += s_g;
  // add in the part for the fermion fields.
  cg(g_X, g_fermion, ITER_MAX, DELTACG, &gam5D_SQR_musqr_wilson);
  ham += scalar_prod_r(g_fermion, g_X);
  
  if(no_timescales > 2) {
    cg(g_gam5DX, g_fermion2, ITER_MAX, DELTACG, &gam5D_SQR_wilson);
    gam5D_SQR_musqr_wilson(g_X, g_temp, g_gam5DX);
    ham += scalar_prod_r(g_fermion2, g_X);
  }

  exphdiff = exp(ham_old-ham);
  acc = accept(exphdiff);
 
  for(i=0; i<GRIDPOINTS; i++) {
    gauge1_old[i]=gauge1[i];
    gauge2_old[i]=gauge2[i];
  }
 
  s_g_old = s_g;
  return(acc);
}

int accept(const double exphdiff)
{
  int acc=0, i;
  double r[1];

  // the acceptance step
  if(exphdiff>=1) {
    acc = 1; 
    R += 1;
  }
  else {
    ranlxd(r,1);
    if(r[0]<exphdiff) {
      acc = 1;
      R += 1;
    }
    else {
      // get the old values for phi, cause the configuration was not accepted
      for (i=0; i<GRIDPOINTS; i++)
	{
	  gauge1[i]=gauge1_old[i];
	  gauge2[i]=gauge2_old[i];
	};
      calculatelinkvars();
      s_g = s_g_old;
    }
  }
  return acc;
}







