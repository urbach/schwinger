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
#include "hmc.h"
#include "2MN_integrator.h"
#include "leapfrog.h"
#include "leapfrog2.h"
#include "fields.h"
#include "dirac.h"

/* inner update defined below */
void LF_inner_update(const int n2, const double tau);

/* Second order Minimal Norm Integration scheme     */
/* for two timescales, non-recursive implementation */

void LF_integrator(const int n1, const int n2, const double tau) {
  int i;
  double dtau1 = tau/((double)n1);

  update_momenta_fermion(0.5*dtau1);

  for(i = 1; i < n1; i++){
    LF_inner_update(n2, dtau1);
    update_momenta_fermion(dtau1);
  }
  LF_inner_update(n2, dtau1);

  update_momenta_fermion(0.5*dtau1);

  return;
}


void LF_inner_update(const int n2, const double tau) {
  int j;
  double dtau = tau/((double)n2);
  
  update_momenta_gauge(0.5*dtau);

  for(j = 1; j < n2; j++) {
    update_gauge(dtau);
    update_momenta_gauge(dtau);
  }
  update_gauge(dtau);

  update_momenta_gauge(0.5*dtau);
  return;
}

