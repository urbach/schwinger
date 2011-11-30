#include <stdlib.h>
#include <math.h>
#include "gauss.h"
#include "../rand/ranlxd.h"
#include <math.h>
#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

double gauss()
{
 double r[2];
 ranlxd(r, 2);
 return (sqrt(-2.0*log(r[0]))*cos(2.0*M_PI*r[1]));
}
