/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2011 Carsten Urbach
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
#include "leapfrog.h"
#include "rec_lf_integrator.h"

void integrate_leap_frog(const double tau, const int S, const int no_timescales, 
			 int n[], const int halfstep,  up_m * up_momenta) {
  int i;
  double eps;

  if(S == no_timescales-1) {
    // coarsest level: perform one step on each scale first
    eps = tau/((double)n[S]);

    for(i = S; i > 0; i--) {
      up_momenta[i](0.5*eps);
      eps /= ((double)n[i-1]);
    }
    up_momenta[0](0.5*eps);
  }

  // current stepsize
  eps = tau/((double)n[S]);
  if(S == 0) {
    // finest level, no further recursion
    for(i = 1; i < n[0]; i++) {
      update_gauge(eps);
      up_momenta[0](eps);
    }
    update_gauge(eps);
    if(halfstep != 1) {
      up_momenta[0](eps);
    }
  }
  else {
    for(i = 1; i < n[S]; i++){
      integrate_leap_frog(eps, S-1, no_timescales, n, 0, up_momenta);
      up_momenta[S](eps);
    }
    if(S == no_timescales-1) {
      integrate_leap_frog(eps, S-1, no_timescales, n, 1, up_momenta);
    }
    else integrate_leap_frog(eps, S-1, no_timescales, n, halfstep, up_momenta);
    if(halfstep != 1 && S != no_timescales-1) {
      up_momenta[S](eps);
    }
  }

  if(S == no_timescales-1) {
    for(i = S; i > 0; i--) {
      up_momenta[i](0.5*eps);
      eps /= ((double)n[i-1]);
    }
    up_momenta[0](0.5*eps);
  }
}
