#if !defined GAULEG_H
#define GAULEG_H

#include <stdio.h>
#include <math.h>

void gauleg(double x1, double x2, double x[], double w[], int n);
/*gauleg.c
  Given the lower and upper limits of integration x1 and x2, and given n,
  this routine returns arrays x[n] and w[n] of length n, containing 
  the abscissas and weights of the Guass-Legendre n-point quadrature formula.
*/

#endif
