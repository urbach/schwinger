#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "lattice.h"
#include "linalg.h"

void add(spinor *Q, spinor *R, spinor *S)
{
  int ix;
  spinor *q,*r,*s;
  
#ifdef OMP
#pragma omp parallel for private(q,r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++)
  {
    q = (spinor *) Q + ix;
    r = (spinor *) R + ix;
    s = (spinor *) S + ix;
    
    (*q).s1=(*r).s1+(*s).s1;
    (*q).s2=(*r).s2+(*s).s2;
  }
}

double square_norm(spinor *S)
{
  int ix;
  static double ds;
  spinor *s;

  ds=0.0;

  for (ix=0;ix<GRIDPOINTS;ix++)
  {
    s=(spinor *)S + ix;
    ds+=creal(conj((*s).s1)*(*s).s1+conj((*s).s2)*(*s).s2);
  }

  return ds;
}

void assign(spinor *R, spinor *S)
{
  int ix;
  spinor *r,*s;
  
#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s1=(*s).s1;
    
    (*r).s2=(*s).s2;
  }
}

void assign_add_mul(spinor *P, spinor *Q, complex double c)
{
  int ix;
  spinor *r,*s;

  
#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    r=(spinor *) P + ix;
    s=(spinor *) Q + ix;
    
    (*r).s1 = (*r).s1 + c*(*s).s1;
    (*r).s2 = (*r).s2 + c*(*s).s2;
  }
}

void assign_add_mul_r(spinor *P, spinor *Q, double c)
{
  int ix;
  static double fact;

  fact=c;
#ifdef OMP
#pragma omp parallel for
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){

    P[ix].s1+=fact*Q[ix].s1;
    P[ix].s2+=fact*Q[ix].s2;
  }
}

void assign_diff_mul(spinor *R, spinor *S, complex double c){
  int ix;
  spinor *r, *s;

#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++)
  {
    s = (spinor *) S + ix;
    r = (spinor *) R + ix;

    (*r).s1=(*r).s1-c*(*s).s1;
    (*r).s2=(*r).s2-c*(*s).s2;
  }
}

void assign_mul_add_r(spinor *R, spinor *S, double c)
{
  int ix;
  static double fact;
  spinor *r,*s;
  
  fact=c;
  
#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s1=fact*(*r).s1+(*s).s1;
    (*r).s2=fact*(*r).s2+(*s).s2;
  }
}

void diff(spinor *Q, spinor *R, spinor *S)
{
  int ix;
  spinor *q,*r,*s;
  
#ifdef OMP
#pragma omp parallel for private(q,r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    q=(spinor *) Q + ix;
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*q).s1=(*r).s1-(*s).s1;
    (*q).s2=(*r).s2-(*s).s2;
  }
}

void mul_r(spinor *R, double c, spinor *S)
{
  int ix;
  spinor *r,*s;

#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s1=c*(*s).s1; 
    (*r).s2=c*(*s).s2;
  }
}

void mul_c(spinor *R, complex double c, spinor *S)
{
  int ix;
  spinor *r,*s;

#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    r=(spinor *) R + ix;
    s=(spinor *) S + ix;
    
    (*r).s1=c*(*s).s1; 
    (*r).s2=c*(*s).s2;
  }
}

complex double scalar_prod(spinor *S, spinor *R)
{
  int ix;
  static complex double ds;
  spinor *s,*r;
  
  /* Real Part */

  ds=0.0 + I*0.0;
  
#ifdef OMP
#pragma omp parallel for private(r,s)
#endif
  for (ix=0;ix<GRIDPOINTS;ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    ds+=conj((*s).s1)*(*r).s1+conj((*s).s2)*(*r).s2;
  }

  return(ds);
}

double scalar_prod_r(spinor *S, spinor *R)
{
  int ix;
  static double ds;
  spinor *s,*r;
  
  /* Real Part */

  ds=0.0;
  
  for (ix=0;ix<GRIDPOINTS;ix++){
    s=(spinor *) S + ix;
    r=(spinor *) R + ix;
    
    //ds+=conj((*s).s1)*(*r).s1+conj((*s).s2)*(*r).s2;
    ds+=creal(conj((*s).s1)*(*r).s1+conj((*s).s2)*(*r).s2);
  }

  return(ds);
}

/************ Conjugate gradient ****************/
/***    Solves the equation f*P = Q           ***/
/************************************************/

int cg(spinor *P, spinor *Q, int max_iter, double eps_sq, matrix_mult f)
{
 double normsq, pro, err, alpha_cg, beta_cg;
 int iteration;
 spinor r[GRIDPOINTS], p[GRIDPOINTS];
 spinor x[GRIDPOINTS], q2p[GRIDPOINTS];
 spinor tmp1[GRIDPOINTS];
 
 /* Initial guess for the solution - zero works well here */ 
 set_zero(x);
  
 /* initialize residue r and search vector p */
 assign(r, Q); /* r = Q - f*x, x=0 */
 assign(p, r);
 normsq = square_norm(r);
  
 double max_rel_err = eps_sq * square_norm(Q);
  
 /* main loop */
#ifdef MONITOR_CG_PROGRESS
  printf("\n\n Starting CG iterations...\n");
#endif
 for(iteration=0; iteration<max_iter; iteration++)
 {
  f(q2p, tmp1, p);
  pro = scalar_prod_r(p, q2p);
  /*  Compute alpha_cg(i+1)   */
  alpha_cg = normsq/pro;
  /*  Compute x_(i+1) = x_i + alpha_cg(i+1) p_i    */
  assign_add_mul_r(x, p,  alpha_cg);
  /*  Compute r_(i+1) = r_i - alpha_cg(i+1) Q2p_i   */
  assign_add_mul_r(r, q2p, -alpha_cg);
  /* Check whether the precision is reached ... */
  err=square_norm(r);
#ifdef MONITOR_CG_PROGRESS
  printf("\t CG iteration %i, |r|^2 = %2.2E\n", iteration, err);
#endif
  if(err<=max_rel_err)
  {
#ifdef MONITOR_CG_PROGRESS
   printf("Required precision reached, stopping CG iterations...\n\n");
#endif
   assign(P, x);
   return(iteration);
  };
     
  /* Compute beta_cg(i+1) */
  beta_cg = err/normsq;
  /* Compute p_(i+1) = r_i+1 + beta_(i+1) p_i     */
  assign_mul_add_r(p, r, beta_cg);
  normsq = err;
 }
 fprintf(stderr, "WARNING: CG didn't converge after %d iterations!\n", max_iter);
 return (-1);
}


spinor gamma5_i(spinor s)
{
 spinor ret;
 ret.s1 = s.s1;
 ret.s2 = -1.0*s.s2;
 return ret;
}

void gamma5(spinor *out, spinor *in)
{
 int i;
 for(i=0; i<GRIDPOINTS; i++)
 {
  out[i]=gamma5_i(in[i]);
 }
}

void set_zero(spinor *P)
{
 int i;
 for(i=0; i<GRIDPOINTS; i++)
 {
  P[i].s1 = 0.0 + I*0.0;
  P[i].s2 = 0.0 + I*0.0;
 };
}

