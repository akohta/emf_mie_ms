#if !defined MULTI_FBEAM_H
#define MULTI_FBEAM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include "const.h"
#include "gauleg.h"
#include "ipw.h"  
#include "fpw.h"  
#include "lgb.h"  
#include "bsb.h"  
#include "bslgb.h" 
#include "rab.h"   

// default datafile name. 
#define fn_ipw "ipw.txt" // plane wave
#define fn_fpw "fpw.txt" // focused plane wave
#define fn_lgb "lgb.txt" // focused plane wave with spiral phase modulation
#define fn_bsb "bsb.txt" // bessel beam
#define fn_blg "blg.txt" // bessel beam with spiral phase modulation
#define fn_rab "rab.txt" // focused radial-azimuthal polarized beam

typedef struct beam_data{
  Ipw *ipw;
  Fpw *fpw;
  LGb *lgb;
  Bsb *bsb;
  BsLGb *blg;
  RAb *rab;
}Bdata;

typedef struct beam_object{
  int n_ipw;  char fname_ipw[28]; // defined beam number, readed beam datafile name
  int n_fpw;  char fname_fpw[28]; 
  int n_lgb;  char fname_lgb[28];
  int n_bsb;  char fname_bsb[28];
  int n_blg;  char fname_blg[28];
  int n_rab;  char fname_rab[28];

  double n_0; // refractive index of surroundings 
  double lambda_0; // wavelength in vacuum
  double omega;  // angular frequency, omega=2 pi f=2 pi c0/lambda_0. c0 : speed of light in vacuum

  Bdata bd;
}Bobj;

void init_mfb(Bobj *obj);       // initalize 
void  read_data_mfb(Bobj *obj); // beam datafile is automatically searched and readed. search path is current directry.
void print_data_mfb(Bobj *obj); // print defined beam data
void setup_mfb(Bobj *obj);      // memory allocation and calculation of coefficients.
void  free_mfb(Bobj *obj);      // free allocated memory

void calc_mfb_EH(double complex *e,double complex *h,double *x,Bobj *obj);
// e[0]=Ex,e[1]=Ey,e[2]=Ez, h[0]=Hx,h[1]=Hy,h[2]=Hz, x[0]=x,x[1]=y,x[2]=z
// for real electromagneticfield  Re(e*exp(-i omega t)), Re(h*exp(-i omega t)). Re(z) : real part of z

// manually set beam datafile. select proper function. returnd value is defined beam numbers.
int read_data_mfb_ipw(char *fname,Bobj *obj); // for plane wave
int read_data_mfb_fpw(char *fname,Bobj *obj); // for focused plane wave
int read_data_mfb_lgb(char *fname,Bobj *obj); // for focused plane wave with spiral phase modulation
int read_data_mfb_bsb(char *fname,Bobj *obj); // for bessel beam
int read_data_mfb_blg(char *fname,Bobj *obj); // for bessel beam with spiral phase modulation
int read_data_mfb_rab(char *fname,Bobj *obj); // for focused radial-azimuthal polarization beam

#endif
