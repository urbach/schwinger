#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "lattice.h"
#include "hmc.h"
#include "fields.h"
#include "measurements.h"
#include "rand/gauss.h"
#include "dirac.h"

/* global variables */
int g_thermalize = 2;      //Number of MC updates for thermalization
int g_measurements = 100;  //Number of measurements (statistically independent configurations)
int g_intermediate = 0;    //Number of MC updates between the measurements
/* extern in fields.h   */
double g_mass = 0.1;       //Fermion mass
double beta   = 1.0;       //Coupling constant for the gauge field
/* extern in hmc.h      */
int    g_steps    = 10;    //Number of steps in the molecular dynamics chain
double g_stepsize = 0.1;   //Size of each step

void rand_spinor(spinor *P)
{
 int i;
 for(i=0; i<GRIDPOINTS; i++)
 {
  P[i].s1 = (gauss() + I*gauss())/sqrt(2);
  P[i].s2 = (gauss() + I*gauss())/sqrt(2);
 }    
}

int main(int argc, char **argv) 
{
 /* Initialize the random number generator */
 rlxd_init(2, 12345); 
 /* Initialize the lattice geometry */
 init_lattice(X1, X2);
  
 //Testing the hermiticity of gam5D_wilson
 //Operator is Hermitian if for any vectors x and y (y, Ax) = (Ay, x)
 //1. Write the procedure which fills a spinor field with gaussian-distributed complex random numbers such that <z z*> = 1
 //2. Initialize two random spinor fields X and Y (arrays of type spinor and size GRIDPOINTS) 
 //3. Check that (y, Ax) = (Ay, x)
 
 spinor X[GRIDPOINTS], Y[GRIDPOINTS], tmp[GRIDPOINTS];
 rand_spinor(X);
 rand_spinor(Y);
 //Calculating p1 = (y, Ax)
 gam5D_wilson(tmp, X);
 complex double p1 = scalar_prod(Y, tmp);
 //Calculating p2 = (Ax, y)
 gam5D_wilson(tmp, Y);
 complex double p2 = scalar_prod(tmp, X);
 //Check that p1 and p2 are equal
 double epsilon = cabs(p1 - p2);
 printf("\n\n\t Hermiticity check: |p1 - p2| = %2.2E - %s \n\n", epsilon, (epsilon < 1E-14)? "OK":"No");

 //Testing the performance of the Conjugate Gradient Solver
 //1. In order to monitor the progress of the solver, #define MONITOR_CG_PROGRESS in linalg.h
 //   This will print the squared norm of the residue at each iteration to the screen
 //2. Initialize the gauge fields with the coldstart() procedure
 //3. Generate a random spinor field  Y and use cg(X, Y, ITER_MAX, DELTACG, &gam5D_SQR_wilson)
 //   to solve the equation (gamma_5 D)^2 X = Y
 //   (gamma_5 D)^2 is implemented as gam5D_SQR_wilson(out, temp, in), see Dirac.h
 //4. Now apply gam5D_SQR_wilson to X, subtract Y and calculate the norm of the result
 //5. Do the same with hotstart() and see how the number of CG iterations changes

 coldstart();
 spinor Y1[GRIDPOINTS];
 rand_spinor(Y);
 cg(X, Y, ITER_MAX, DELTACG, &gam5D_SQR_wilson);
 gam5D_SQR_wilson(Y1, tmp, X);
 diff(tmp, Y1, Y);
 epsilon = sqrt(square_norm(tmp));
 printf("\n\n\t Test of CG inverter: |Q X - Y| = %2.2E - %s \n\n", epsilon, (epsilon < 1E-6)? "OK":"No");

 free(left1);
 free(left2);
 free(right1);
 free(right2);
  
 system("PAUSE");
 return 0;
}
