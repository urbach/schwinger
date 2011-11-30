#ifndef _LINALG_H
#define _LINALG_H

#include "complex/complex.h"

//If this is #defined, CG will print the norm of each residue
#undef MONITOR_CG_PROGRESS

/***********************************************************/
/**** Implementation of linear algebra on spinor fields ****/
/***********************************************************/

//Basic structure for the fermion field - two-component spinor in two dimensions
typedef struct spinor {
 complex double s1, s2;
} spinor;

//Function type for linear operators acting on spinor space
typedef void (*matrix_mult) (spinor * const, spinor * const, spinor * const);

void set_zero(spinor *P);                                     //Set the field P to zero

void assign_add_mul(spinor *P, spinor *Q, complex double c);  // P = P + c Q
void assign(spinor *R, spinor *S);                            // R = S
double scalar_prod_r(spinor *S, spinor *R);                   // Re(S*, R)
complex double scalar_prod(spinor *S, spinor *R);             // (S*, R)
void assign_mul_add_r(spinor *R, spinor *S, double c);        // R = c R + S, c is real
void assign_add_mul_r(spinor *P, spinor *Q, double c);        // P = P + c Q, c is real
void assign_add_mul(spinor *P, spinor *Q, complex double c);  // P = P + c Q
void assign_diff_mul(spinor *R, spinor *S, complex double c); // R = R - c S
void diff(spinor *Q, spinor *R, spinor *S);                   // Q = R - S
void mul_r(spinor *R, double c, spinor *S);                   // R = c S, c is real
void mul_c(spinor *R, complex double c, spinor *S);           // R = c S, c is complex
double square_norm(spinor *P);                                // (P, P*)
void add(spinor *Q, spinor *R, spinor *S);                    // Q = R + S

spinor gamma5_i(spinor s);                                    //Multiplies s by gamma_5
void gamma5(spinor *out, spinor *in);                         //out = gamma_5 in

//Conjugate gradient method...
int cg(spinor *P, spinor *Q, int max_iter, double eps_sq, matrix_mult f);

#endif
