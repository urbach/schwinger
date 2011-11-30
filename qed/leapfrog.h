#ifndef _LEAPFROG_H
#define _LEAPFROG_H

//This is the implementation of the Leapfrog integrator ... 

void leapfrog(const int nsteps, const double dtau);

void update_momenta(const double dtau);
void update_gauge(const double dtau);

#endif
