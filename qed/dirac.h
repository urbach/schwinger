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


#ifndef _DIRAC_H
#define _DIRAC_H

/**********************************************************************************************/
/*** Implementation of the Dirac operator, its square and derivatives over the gauge fields ***/
/**********************************************************************************************/

#include "linalg.h"

void   gam5D_wilson(spinor *out, spinor *in);                   //out = gamma_5*D*in, gamma_5*D has the advantage of being hermitian
void   gam5D_SQR_wilson(spinor *out, spinor *temp, spinor *in); //out = (gamma_5*D)^2 in, temp is used for internal calculations
void gam5D_SQR_musqr_wilson(spinor *out, spinor *temp, spinor *in);
double trX_dQ_wilson_dalpha1_X(int j);                          //Derivative of \phi (gamma_5*D)^(-2) \phi over the 1st component of the gauge field at site i, phi is the pseudofermion
double trX_dQ_wilson_dalpha2_X(int j);                          //Derivative of \phi (gamma_5*D)^(-2) \phi over the 2nd component of the gauge field at site i

//At present the boundary conditions for fermions are periodic in all directions!!!

#endif
