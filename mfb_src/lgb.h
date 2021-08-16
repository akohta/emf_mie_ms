#if !defined LGB_H
#define LGB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "mfb_const.h"
#include "gauleg.h"
#include "osu_mksa.h"
#include "my_utils.h"

typedef struct lgb_data{
  double ki;            // wave number 
  double E0;            // power coefficient
  double complex ex,ey; // normalized polarization coefficient for electric field
  double complex hx,hy; // normalized polarization coefficient for magnetic field
  double cos_t,sin_t;   // cos(theta),sin(theta)
  double cos_p,sin_p;   // cos(phi),sin(phi), theta,phi:rotation parameter
  int nt;               // sampling number for gau-leg quadrature
  
  // data for gau-leg
  double *wt,*ct,*st;
  // data for trapezoidal
  double *cp,*sp; 
  // phase data
  double complex *ph_p;
}LGbD;

typedef struct lgb{
  double lambda0;          // wavelength of plane wave in vacuum [m]
  double ni;               // reflactive index of surrounding
  double NA;               // Numerical Aperture
  double power;            // incident beam power passig through pupil plane [W]
  double complex e0x;      // x-component of polarization coefficient include phase
  double complex e0y;      // y-component of polarization coefficient include phase
  double fx,fy,fz;         // translation vector (fx,fy,fz)
  double theta,phi;        // rotation parameter
  int    lg_m;             // mode number for spiral phase modulation
  int    nn;               // sampling number for gau-leg quadrature
  LGbD data;
}LGb;

void print_data_lgb(LGb *lgb);
void print_data_lgb_mksa(LGb *lgb);
void setup_LGb(LGb *lgb);  // calculate coefficient and memory allocation
void free_LGb(LGb *fpw);   // memory free 

void calc_lgb_EH(double complex *e,double complex *h,double *x,LGb *lgb);
void calc_lgb_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,LGb *lgb);
// directional devirative ( define : df/dv=df/dx v_x + df/dy v_y + df/fz v_z )
// dedv[0]=dE_x/dv, dedv[1]=dE_y/dv, dedv[2]=dE_z/dv, dhdv[0]=dH_x/dv, dhdv[1]=dH_y/dv, dhdv[2]=dH_z/dv


#endif
