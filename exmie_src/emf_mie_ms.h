#if !defined EMF_MIE_MS_H_
#define EMF_MIE_MS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "multi_fbeam.h"
#include "rctbess.h"
#include "gsl/gsl_specfunc.h"

// default sphere datafile name
#define fn_sphr "msphr.txt"

// for iterative operation 
#define ito_max 256
#define ito_eps 1.0e-12
#define ito_breakcount 2

// struct definition
typedef struct discrete_data{
  int    nt,np;                  // sampling point number
  double *xt,*wt,*xp,*wp;        // gauss-legender weight and point data
  double complex *eri,*hri;      // incident field on sphere surface. Er(r,theta,phi), Hr(r,theta,phi) 
  double complex *ers,*hrs;      // scattered field on sphere surface. for iterative solutions.
  int l_max;                     // limit of Ricatti-Bessel function order 
  double *cab;                   // coefficient for A_lm B_lm
  double complex *ca,*cb,*cc,*cd;// coefficient for a_lm, b_lm, c_lm, d_lm
  double complex *Alm,*Blm;      // coefficient A_lm, B_lm
}DDT;

typedef struct sphere_data{
  double a;             // radius of sphere [m]
  double complex ns;    // reflactive index of sphere
  double xs,ys,zs;      // sphere center (xs,ys,zs)
  int bsn;              // basic sampling number for Gauss-Legendre quadrature
  int bdv;              // division number of sphere surface ( per M_PI )
  int l_limit;          // limit of order number 
  
  DDT ddt;              // discrete data
}SPD;

typedef struct sphere_objects{
  int  n_sphr;          // number of spheres
  SPD *sp;              // sphere data 
  Bobj bm;              // multi_fbeam data
}MSPD;

// ---- emf_mie_ms.c ----
void  read_data_ms(MSPD *msp);    // seach the sphere datafile (defined fn_sphr) and load 
void print_data_ms(MSPD *msp);    // print loaded data 
void setup_ms(MSPD *msp);         // allocate memory and setup coefficients
void  free_ms(MSPD *msp);         // free allocated memory
void iterative_ops_ms(MSPD *msp); // solve multiple scattering

// ---- emf_mie_ms_field.c ----
void  incident_EH_ms(double complex *e,double complex *h,double *r,MSPD *msp); // calculate incident field 
void  internal_EH_ms(double complex *e,double complex *h,double *r,MSPD *msp); // calculate internal field 
void scattered_EH_ms(double complex *e,double complex *h,double *r,MSPD *msp); // calculate scattered field
void     total_EH_ms(double complex *e,double complex *h,double *r,MSPD *msp); // calculate total field
// e[0]=Ex,e[1]=Ey,e[2]=Ez,h[0]=Hx,h[1]=Hy,h[2]=Hz,r[0]=x,r[1]=y,r[2]=z

// ---- emf_mie_ms_force.c ----
void  force_ms(int id,double *vf,MSPD *msp);                  // calculate radiation force of the sphere 
void torque_ms(int id,double *vn,MSPD *msp);                  // calculate radiation torque of the sphere 
void force_torque_ms(int id,double *vf,double *vn,MSPD *msp); // calculate radiation force and torque of the sphere
// id:sphereid, vf[0]=Fx,vf[1]=Fy,vf[2]=Fz,vn[0]=Nx,vn[1]=Ny,vn[2]=Nz

// ---- emf_mie_ms_dat.c ----
void write_dat_ms(char *fn,MSPD *msp); // write datafile
void  read_dat_ms(char *fn,MSPD *msp); // read datafile and allocate memory

#endif