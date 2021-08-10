#if !defined RAB_H
#define RAB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "const.h"
#include "gauleg.h"

typedef struct rab_data{
  double ki;             // wave number 
  double E0;             // power coefficient
  double complex er,ea;  // normalized polarization coefficient for electric field     
  double complex hr,ha;  // normalized polarization coefficient for magnetic field
  double cos_t,sin_t;    // cos(theta),sin(theta)
  double cos_p,sin_p;    // cos(phi),sin(phi), theta,phi:rotation parameter
  int nt;                // sampling number of gau-leg guadrature
  // gau-leg data
  double *wt;
  double *ct;
  double *st;  
  // trapezoidal data
  double *cp; 
  double *sp;  
}RAbD;

typedef struct rab{
  double lambda0;          // wavelength of plane wave in vacuum [m]
  double ni;               // reflactive index of surrounding
  double NA;               // Numerical Aperture
  double power;            // incident beam power passing through pupil plane [W]
  double complex e0r;      // radial-component of polarization coefficient include phase
  double complex e0a;      // azimuthal-component of polarization coefficient include phase
  double fx,fy,fz;         // translation vector (fx,fy,fz)
  double theta,phi;        // rotation parameter
  int    nn;               // gauss-legendre quadrature sampling number
  RAbD data;
}RAb;

void read_data_rab(char *rfile,RAb *rab);
void print_data_rab(RAb *rab);
void setup_rab(RAb *rab);  // calculate coefficient and memory allocation 
void  free_rab(RAb *rab);  // memory free 

void calc_rab_EH(double _Complex *e,double _Complex *h,double *x,RAb *rab);

#endif
