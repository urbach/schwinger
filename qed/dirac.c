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
#include <math.h>
#ifndef M_PI
# define M_PI    3.14159265358979323846f
#endif
#include "fields.h"
#include <complex.h>
#include "lattice.h"
#include "linalg.h"
#include "rand/ranlxd.h"
#include "dirac.h"

// out = K^dag * in
// K_nm = (M + 4r) \delta_nm - 0.5 * \sum_\mu  ((r - \gamma_\mu) \delta_m,n+\mu + (r + \gamma_\mu) \delta_m,n-\mu)?
// Rothe, p. 56
// gamma_1/2 -> sigma_1/2, gamma_5 -> sigma_3, gamma_5 vorne dran
void gam5D_wilson(spinor *out, spinor *in) {
  int i;
  double factor = (2*g_R + g_mass);

#ifdef OMP
#pragma omp parallel for 
#endif
  for(i=0; i<GRIDPOINTS; i++) {
    complex double link1_i = link1[i];
    spinor in_i = in[i];
    int left1_i = left1[i];
    int left2_i = left2[i];
    int right1_i = right1[i];
    int right2_i = right2[i];
    spinor in_right1_i = in[right1_i];
    spinor in_right2_i = in[right2_i];
    spinor in_left1_i = in[left1_i];
    spinor in_left2_i = in[left2_i];
    complex double cconj_link1_left1_i = cconj(link1[left1_i]);
    complex double cconj_link2_left2_i = cconj(link2[left2_i]);
    out[i].s1 = factor * in_i.s1 - 
    0.5*(link1_i*(g_R*in_right1_i.s1 - in_right1_i.s2) +
         cconj_link1_left1_i * ( g_R*in_left1_i.s1  +   in_left1_i.s2)  +
         link2[i] * ( g_R*in_right2_i.s1 + I*in_right2_i.s2) +
         cconj_link2_left2_i * ( g_R*in_left2_i.s1  - I*in_left2_i.s2));
    
    out[i].s2 = - factor * in_i.s2 -
    0.5*(link1_i * (   in_right1_i.s1 - g_R*in_right1_i.s2) -
         cconj_link1_left1_i * (in_left1_i.s1  + g_R*in_left1_i.s2)  +
         link2[i] * ( I*in_right2_i.s1 - g_R*in_right2_i.s2) -
         cconj_link2_left2_i * (I*in_left2_i.s1  + g_R*in_left2_i.s2));
  }
  return;
}

/* Q^2 */
void gam5D_SQR_wilson(spinor *out, spinor *temp, spinor *in) {
  gam5D_wilson(temp, in);
  gam5D_wilson(out, temp);
}

/* Q^2 + mu^2 */
/* should have better condition than gam5D_SQR_wilson() */
void gam5D_SQR_musqr_wilson(spinor *out, spinor *temp, spinor *in) {
  gam5D_wilson(temp, in);
  gam5D_wilson(out, temp);
  assign_add_mul(out, in, g_musqr);
}



// assumes that in g_gam5DX is gam5D * X 
// (X^*)(dM/da)(M)(X) + (X^*)(M)(dM/da)(X) = 2Re{(X^*)(dM/da)(M)(X)}
double trX_dQ_wilson_dalpha1_X(int j)
{
  return(creal(I*((cconj(link1[j]) 
		   * (cconj(g_X[right1[j]].s1)*(g_R*g_gam5DX[j].s1 +  g_gam5DX[j].s2) -
		      cconj(g_X[right1[j]].s2)*(  g_gam5DX[j].s1   + g_R*g_gam5DX[j].s2))) 
		  -
		  (link1[j] 
		   * (cconj(g_X[j].s1) * (g_R*g_gam5DX[right1[j]].s1-  g_gam5DX[right1[j]].s2) + 
		      cconj(g_X[j].s2) * (  g_gam5DX[right1[j]].s1-g_R*g_gam5DX[right1[j]].s2)))
		  )
	       )
	 );
}

double trX_dQ_wilson_dalpha2_X(int j)
{
  return(creal(I*((cconj(link2[j]) 
		   * (cconj(g_X[right2[j]].s1)*(g_R*g_gam5DX[j].s1  - I*g_gam5DX[j].s2) -
		      cconj(g_X[right2[j]].s2)*(I*g_gam5DX[j].s1    + g_R*g_gam5DX[j].s2))) 
		  -
		  (link2[j] 
		   * (cconj(g_X[j].s1) * (g_R*g_gam5DX[right2[j]].s1+I*g_gam5DX[right2[j]].s2) +
		      cconj(g_X[j].s2) * (I*g_gam5DX[right2[j]].s1-g_R*g_gam5DX[right2[j]].s2)))
		  )
	       )
	 );
}
