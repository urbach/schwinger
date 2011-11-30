#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./complex/complex.h"

int main(int argc, char *argv[])
{
    double complex a, b, c;
    
    a = 0.0 + I*1.0;
    b = 1.0 + I*0.0;
    
    c = a + b;
    
    printf("c = %2.2lf + I*%2.2lf, |c| = %2.2lf\n", creal(c), cimag(c), cabs(c));
    printf("conj(c) = %2.2lf + I*%2.2lf\n", creal(cconj(c)), cimag(cconj(c)));
    printf("cos(c) = %2.2lf + I*%2.2lf\n", creal(ccos(c)), cimag(ccos(c)));
    
    system("PAUSE");
    return EXIT_SUCCESS;
}
