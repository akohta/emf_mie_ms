#if !defined FPW_H
#define FPW_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "const.h"
#include "gauleg.h"


typedef struct fpw_data{
  double ki;            // wave number 
  double E0;            // power coefficient
  double complex ex,ey; // normalized polarization coefficient for electric field 
  double complex hx,hy; // normalized polarization coefficient for magnetic field
  double cos_t,sin_t;   // cos(theta),sin(theta)
  double cos_p,sin_p;   // cos(phi),sin(phi), theta,phi:rotation parameter
  int nt;               // sampling number for gauss-legendre quadrature
  
  // data for gauss-legendre 
  double *wt,*ct,*st;  
  // data for trapezoidial  
  double *cp,*sp;  

}FpwD;

typedef struct fpw{
  double lambda0;          // wavelength of plane wave in vacuum [m]
  double ni;               // reflactive index of surrounding
  double NA;               // Numerical Aperture
  double power;            // incident beam power passing through pupil plane [W]
  double complex e0x;      // x-component of polarization coefficient include phase
  double complex e0y;      // y-component of polarization ceofficient include phase
  double fx,fy,fz;         // translation vector (fx,fy,fz) 
  double theta,phi;        // rotation parameter
  int    nn;               // gauss-legendre integration sampling number 
  FpwD data;
}Fpw;

void read_data_fpw(char *rfile,Fpw *fpw);
void print_data_fpw(Fpw *fpw);
void setup_Fpw(Fpw *fpw);  // calculate parameter and memory allocation
void free_Fpw(Fpw *fpw);   // memory free 

void calc_fpw_EH(double complex *e,double complex *h,double *x,Fpw *fpw);

#endif
