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

#ifndef _2MN_INTEGRATOR_H
#define _2MN_INTEGRATOR_H

extern double gauge_force;
extern double PF1_force;
extern double PF2_force;

void MN2_integrator(const int n1, const int n2, const double tau, const double lambda);
void update_momenta_fermion(const double dtau);
void update_momenta_gauge(const double dtau);

void update_momenta_PF2(const double dtau);
void update_momenta_PF1(const double dtau);
#endif
