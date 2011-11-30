#ifndef _LATTICE_H
#define _LATTICE_H

/***********************************************/
/***** This unit defines lattice geometry  *****/
/***********************************************/    

#define X1 (12)                                       //Lattice size in direction 1
#define X2 (X1)                                      //Lattice size in direction 2
#define GRIDPOINTS (X1*X2)                           //Total number of lattice sites

//These are arrays of neighbors to the left and to the right of a lattice point with index i in direction 1 or 2, respectively
//E.g. left1[i] is the left neighboring site of the site i in the direction 1
extern int *left1;     
extern int *right1;
extern int *left2;
extern int *right2;

int  init_lattice(const int NX, const int NY);      //This procedure initializes the above left/right arrays
int  idx(const int x1, const int x2, const int NX); //Converts Cartesian coordinates x1 and x2 into a single integer site index

#endif
