#ifndef _HMC_H
#define _HMC_H

#include "rec_lf_integrator.h"

/***********************************************************************************/
/*** This unit implement the basic HMC update step and the necessary procedures ****/
/***********************************************************************************/

/* Maximal number of iterations for CG method */
#define ITER_MAX 5000
/* Tolerance for CG method */
#define DELTACG 1.e-16
/* Step size and the number of steps for the leapfrog */
/* declared in qed.c */

extern double ham, ham_old;

/* The following we need for the recursive integration scheme */
/* the integration step numbers for the different time scales */
extern int n_steps[3];
/* list of function pointers to the momentum update functions */
extern up_m up_momenta[3];
extern int no_timescales;
extern double tau;

extern int R;  // Counter of all accepted configurations
extern int g_cgiterations1, g_cgiterations2;

int  update(); //Basic HMC update step
int  accept(const double exphdiff); //Accept or reject the trajectory depending on exphdiff

#endif
