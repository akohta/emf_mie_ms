#if !defined IPW_H
#define IPW_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "const.h"

typedef struct incident_planewave_data{
  double E0;           // power coefficient
  double ki;           // wave number 
  double cos_t,sin_t;  // cos(theta),sin(theta)
  double cos_p,sin_p;  // sin(phi),cos(phi)
  
  double complex ex,ey,ez;  // normalized electric field vector 
  double complex hx,hy,hz;  // normalized magnetic field vector
}IpwD;

typedef struct incident_planewave{
  double lambda0;       // wavelength of plane wave in vacuum [m]
  double ni;            // reflactive index of surroundings
  double power;         // planewave power [W/m^2]
  double complex e0x;   // x component of polarization coefficient include phase
  double complex e0y;   // y component of polarization coefficient include phase
  double fx,fy,fz;      // translation vector (fx,fy,fz)
  double theta;         // theta,phi are rotation parameter.  
  double phi;           // ex. unit wave vector k : kx=sin(theta)*cos(phi),ky=sin(theta)*sin(phi),kz=cos(theta)
  IpwD data;            // data
}Ipw;

void read_data_ipw(char *rfile,Ipw *ipw);
void print_data_ipw(Ipw *ipw);
void setup_ipw(Ipw *ipw);

void calc_ipw_EH(double complex *e,double complex *h,double *x,Ipw *ipw);

#endif
