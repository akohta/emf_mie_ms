#if !defined RAB_H
#define RAB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "mfb_const.h"
#include "gauleg.h"
#include "osu_mksa.h"
#include "my_utils.h"

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

void print_data_rab(RAb *rab);
void print_data_rab_mksa(RAb *rab);
void setup_rab(RAb *rab);  // calculate coefficient and memory allocation 
void  free_rab(RAb *rab);  // memory free 

void calc_rab_EH(double _Complex *e,double _Complex *h,double *x,RAb *rab);
void calc_rab_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,RAb *rab);
// directional devirative ( define : df/dv=df/dx v_x + df/dy v_y + df/fz v_z )
// dedv[0]=dE_x/dv, dedv[1]=dE_y/dv, dedv[2]=dE_z/dv, dhdv[0]=dH_x/dv, dhdv[1]=dH_y/dv, dhdv[2]=dH_z/dv

#endif
