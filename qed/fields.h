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

#ifndef _FIELDS_H
#define _FIELDS_H

/*****************************************************************************/
/*** This unit defines the field variables used by the HMC procedure,      ***/
/*** procedures for their initialization and the action of the gauge field ***/
/*****************************************************************************/

#include "lattice.h" /* for GRIDPOINTS */
#include "linalg.h"  /* for spinor     */
#include "rand/ranlxd.h"

/* the Wilson-parameter (usually set to 1) */
#define g_R 1.0

/* declared in qed.c */
extern double g_mass;
extern double g_musqr;
extern double beta;

extern spinor g_fermion[GRIDPOINTS]; //Pseudofermion field
extern spinor g_fermion2[GRIDPOINTS]; //Pseudofermion field
extern spinor g_X[GRIDPOINTS];       //g_X = (gamma_5 D)^{-2} g_fermion
extern spinor g_gam5DX[GRIDPOINTS];  //gam5D_wilson * g_X
extern spinor g_temp[GRIDPOINTS];

//Here the names like gauge1, gauge2 etc. mean the first and the second component of the gauge field
extern double gauge1[GRIDPOINTS], gauge2[GRIDPOINTS];         //Non-compact real-valued gauge fields
extern double gauge1_old[GRIDPOINTS], gauge2_old[GRIDPOINTS]; //Old values of the gauge field
extern complex double link1[GRIDPOINTS], link2[GRIDPOINTS];   //Compact lattice gauge fields: link = exp(I*gauge)
extern double gp1[GRIDPOINTS], gp2[GRIDPOINTS];               //Momenta corresponding to the gauge fields

extern double s_g, s_g_old; //Action for the gauge fields in the beginning and in the end of the MD trajectory

void coldstart(); //Cold start for the gauge fields
void hotstart();  //Hot start for the gauge fields

int calculatelinkvars(); //Sets the values of the compact gauge fields from noncompact ones

double S_G(int i);   //Action of the gauge fields
double DS_G1(int i); //Derivative of the action over the 1st component of the noncompact gauge field
double DS_G2(int i); //Derivative of the action over the 2nd component of the noncompact gauge field
#endif
