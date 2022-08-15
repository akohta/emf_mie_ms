#if !defined FGB_H
#define FGB_H

#include <complex.h>

typedef struct focused_gaussian_beam_data_t{
  double E0;           // power coefficient
  double ki;           // wave number 
  double complex ex,ey;// field coefficients
  double cos_t,sin_t;  // cos(theta),sin(theta)
  double cos_p,sin_p;  // sin(phi),cos(phi)
  
  double l,i_l;        // l=ki*omega0*omega0 (confocal parameter). i_l=1.0/l;
  double i_w0;         // i_w0=1.0/omega0;
}FgbD;

typedef struct focused_gaussian_beam_t{
  double lambda0;       // wavelength in vacuum 
  double ni;            // refractive index of surroundings
  double omega0;        // beam waist radius
  double power;         // total power transmitted by the beam
  double complex e0x;   // x component of polarization coefficient include initial phase
  double complex e0y;   // y component of polarization coefficient include initial phase
  double fx,fy,fz;      // focal position (fx,fy,fz)
  double theta,phi;     // theta,phi are rotation parameter.  

  FgbD data;            // data
}Fgb;

void print_data_fgb(Fgb *fgb);
void print_data_fgb_mksa(Fgb *fgb);
void setup_fgb(Fgb *fgb);

void calc_fgb_EH(double complex *e,double complex *h,double *x,Fgb *fgb);
void calc_fgb_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,Fgb *fgb);
// directional devirative ( define : df/dv=df/dx v_x + df/dy v_y + df/fz v_z )
// dedv[0]=dE_x/dv, dedv[1]=dE_y/dv, dedv[2]=dE_z/dv, dhdv[0]=dH_x/dv, dhdv[1]=dH_y/dv, dhdv[2]=dH_z/dv

#endif
