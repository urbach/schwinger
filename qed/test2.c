#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "lattice.h"
#include "hmc.h"
#include "fields.h"
#include "rand/gauss.h"
#include "measurements.h"
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

double hamiltonian()
{
 double h = 0.0;
 spinor tmp[GRIDPOINTS];
 int i;
 for(i=0; i<GRIDPOINTS; i++)
 {
  h += 0.5*(gp1[i]*gp1[i] + gp2[i]*gp2[i]); //Calculate the gauge-field part of the hamiltonian
  h += S_G(i);
 };
 //Calculate the fermionic contribution 0.5*(g_fermion, (gamma_5 D)^{-2} g_fermion)
 i = cg(tmp, g_fermion, ITER_MAX, DELTACG, &gam5D_SQR_wilson);
 h += 0.5*scalar_prod_r(g_fermion, tmp);
 return h;
}

int main(int argc, char **argv) 
{
  /* Initialize the random number generator */
  rlxd_init(2, 12345); 
  /* Initialize the lattice geometry */
  init_lattice(X1, X2);
  /* Initialize the gauge fields */
  coldstart();

  //Test the reversibility of the leapfrog procedure
  //1. Initialize the momenta conjugate to the gauge fields (gp1 and gp2 arrays)
  // and the pseudofermion fields. Fill gp1 and gp2 with random gaussian-distributed numbers.
  // In order to initialize the pseudofermion field, first generate the random 
  // spinor field \xi distributed according to P(\xi) \sim \exp{-1/2 \xi\xi*}
  // and then calculate the pseudofermion field as \phi = gam5D_wilson \xi
  // Use gauss() and rand_spinor() procedures. Store the pseudofermion field in the g_fermion array. 
  //2. Implement the procedure for calculating the hamiltonian of the fields gauge1,2, gp1,2 and g_fermion
  //   Calculate the hamiltonian for the initial field configuration
  //3. Do the molecular dynamics chain with, say, 10 steps of size 0.1 and calculate the hamiltonian of the new state
  //4. Reverse the momenta gp1 and gp2 - just multiply them by -1
  //5. Again do the molecular dynamics chain with the same number of steps of the same size
  //6. Again calculate the hamiltonian of the obtained state and compare it 
  //   with the hamiltonian of the initial state
  //7. Repeat everything, but initialize the gauge fields with the hotstart() and see how the values of 
  //   the hamiltonian before and after the MD chain
  int i;
  rand_spinor(g_X); //Initialize the \xi = (gamma5_D)^{-1} \phi
  gam5D_wilson(g_fermion, g_X); //Calculate the  pseudofermion field
  for(i=0; i<GRIDPOINTS; i++)
  {
   gp1[i] = gauss(); //Initialize the momenta for the gauge fields
   gp2[i] = gauss(); 
  };
  ham_old = hamiltonian();
  printf("Hamiltonian for the initial state:\t %2.4lf \n", ham_old);
  
  // Do the molecular dynamic chain with the step size g_stepsize
  leapfrog(g_steps, g_stepsize);
  printf("Hamiltonian after %i leapfrog steps of size %2.4lf:\t %2.4lf \n", g_steps, g_stepsize, hamiltonian());
  
  //Reverse the momenta
  for(i=0; i<GRIDPOINTS; i++)
  {
   gp1[i] *= -1.0;
   gp2[i] *= -1.0;
  };
  printf("Now reversing the momenta...\n");
  // Do the molecular dynamic chain with the step size g_stepsize
  leapfrog(g_steps, g_stepsize);
  printf("Hamiltonian after %i leapfrog steps of size %2.4lf:\t %2.4lf \n", g_steps, g_stepsize, hamiltonian());
  

  free(left1);
  free(left2);
  free(right1);
  free(right2);
  
  system("PAUSE");
  return 0;
}

