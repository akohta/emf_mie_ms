#if !defined BSB_H
#define BSB_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "const.h"

typedef struct bsb_data{
  double ki;             // wave number 
  double E0;             // power coefficient
  double complex ex,ey;  // normalized polarization coefficient include initial phase 
  double complex hx,hy;  // normalized polarization coefficient include initial phase
  double cos_t,sin_t;    // cos(theta),sin(theta)
  double cos_p,sin_p;    // cos(phi),sin(phi),theta,phi:rotation parameter
  double ct,st;          // cos(theta_d),sin(theta_d)
  // trapezoidal data
  double *cp; 
  double *sp;
}BsbD;

typedef struct bsb{
  double lambda0;          // wavelength of plane wave in vacuum [m]
  double ni;               // reflactive index of surrounding
  double d_angle;          // deflection angle [rad] (for axicon:=arcsin(n_axi/ni*sin(angle_axi))-angle_axi)
  double power;            // beam power per unit square metre [W/m^2]
  double complex e0x;      // x-component of polarization coefficient include phase
  double complex e0y;      // y-component of polarization coefficient include phase
  double fx,fy,fz;         // translation vector (fx,fy,fz)
  double theta,phi;        // rotation parameter 
  BsbD data;
}Bsb;

void read_data_bsb(char *rfile,Bsb *bsb);
void print_data_bsb(Bsb *bsb);
void setup_Bsb(Bsb *bsb);  // calculate coefficient and memory allocation 
void free_Bsb(Bsb *bsb);   // memory free

void calc_bsb_EH(double complex *e,double complex *h,double *x,Bsb *bsb);

#endif
