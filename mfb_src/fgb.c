#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fgb.h"
#include "mfb_const.h"
#include "osu_mksa.h"

void print_data_fgb(Fgb *fgb)
{
  printf("-- first-order focused Gaussian beam --\n");
  printf("wave length in vacuum                       : %15.14g\n",fgb->lambda0);
  printf("refractive index of surroundings            : %15.14g\n",fgb->ni);
  printf("beam waist radius                           : %15.14g\n",fgb->omega0);
  printf("incident beam power                         : %15.14g\n",fgb->power);
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fgb->e0x),cimag(fgb->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fgb->e0y),cimag(fgb->e0y));
  printf("x-component of translation vector           : %15.14g\n",fgb->fx);
  printf("y-component of translation vector           : %15.14g\n",fgb->fy);
  printf("z-component of translation vector           : %15.14g\n",fgb->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",fgb->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",fgb->phi);
  
  if(fgb->lambda0/fgb->ni > fgb->omega0){
    printf("Note that the analysis results will be a good approximation as long as the beam waist radius is greater than the wavelength.\n");
  }
}

void print_data_fgb_mksa(Fgb *fgb)
{
  printf("-- first-order focused Gaussian beam, MKSA system --\n");
  printf("wave length in vacuum                    [m]: %15.14g\n",OSUtoMKSA_length(fgb->lambda0));
  printf("refractive index of surroundings            : %15.14g\n",fgb->ni);
  printf("beam waist radius                        [m]: %15.14g\n",OSUtoMKSA_length(fgb->omega0));
  printf("incident beam power                      [W]: %15.14g\n",OSUtoMKSA_power(fgb->power));
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fgb->e0x), cimag(fgb->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fgb->e0y), cimag(fgb->e0y));
  printf("x-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(fgb->fx));
  printf("y-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(fgb->fy));
  printf("z-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(fgb->fz));
  printf("rotation parameter theta               [rad]: %15.14g\n",fgb->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",fgb->phi); 

  if(fgb->lambda0/fgb->ni > fgb->omega0){
    printf("Note that the analysis results will be a good approximation as long as the beam waist radius is greater than the wavelength.\n");
  }
}

void setup_fgb(Fgb *fgb)
{
  double ce;
  
  fgb->data.E0=sqrt(4.0*fgb->power/(M_PI*fgb->omega0*fgb->omega0*fgb->ni));
  fgb->data.ki=2.0*M_PI*fgb->ni/fgb->lambda0;
  fgb->data.l=fgb->data.ki*fgb->omega0*fgb->omega0;
  fgb->data.i_l=1.0/fgb->data.l;
  fgb->data.i_w0=1.0/fgb->omega0;
  
  ce=1.0/sqrt(pow(cabs(fgb->e0x),2)+pow(cabs(fgb->e0y),2));
  fgb->data.ex=fgb->data.E0*ce*fgb->e0x;
  fgb->data.ey=fgb->data.E0*ce*fgb->e0y;

  fgb->data.sin_t=sin(fgb->theta);
  fgb->data.cos_t=cos(fgb->theta);
  fgb->data.sin_p=sin(fgb->phi);
  fgb->data.cos_p=cos(fgb->phi);
}

void calc_fgb_EH(double complex *e,double complex *h,double *x,Fgb *fgb)
{
  void calc_fgb_eh(double complex *e,double complex *h,double *x,Fgb *fgb);
  
  double xn[3],xb,yb,zb;
  double cos_p,sin_p,cos_t,sin_t;
  double complex te[3],th[3];
  
  xb=x[0]-fgb->fx;  yb=x[1]-fgb->fy;  zb=x[2]-fgb->fz;
  cos_p=fgb->data.cos_p;    sin_p=fgb->data.sin_p;
  cos_t=fgb->data.cos_t;    sin_t=fgb->data.sin_t;
  
  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  calc_fgb_eh(te,th,xn,fgb);

  e[0]= te[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+te[1]*sin_p*cos_p*(cos_t-1.0)        +te[2]*sin_t*cos_p;
  e[1]= te[0]*sin_p*cos_p*(cos_t-1.0)        +te[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+te[2]*sin_t*sin_p;
  e[2]=-te[0]*sin_t*cos_p                    -te[1]*sin_t*sin_p                    +te[2]*cos_t;
  h[0]= th[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+th[1]*sin_p*cos_p*(cos_t-1.0)        +th[2]*sin_t*cos_p;
  h[1]= th[0]*sin_p*cos_p*(cos_t-1.0)        +th[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+th[2]*sin_t*sin_p;
  h[2]=-th[0]*sin_t*cos_p                    -th[1]*sin_t*sin_p                    +th[2]*cos_t;
}

void calc_fgb_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,Fgb *fgb)
{
  void calc_fgb_eh_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,Fgb *fgb);
  
  double xn[3],xb,yb,zb,nn[3];
  double cos_p,sin_p,cos_t,sin_t;
  double complex te[3],th[3],tde[3],tdh[3];
  xb=x[0]-fgb->fx;  yb=x[1]-fgb->fy;  zb=x[2]-fgb->fz;

  cos_p=fgb->data.cos_p;    sin_p=fgb->data.sin_p;
  cos_t=fgb->data.cos_t;    sin_t=fgb->data.sin_t;

  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  nn[0]= v[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+v[1]*sin_p*cos_p*(cos_t-1.0)        -v[2]*sin_t*cos_p;
  nn[1]= v[0]*sin_p*cos_p*(cos_t-1.0)        +v[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)-v[2]*sin_t*sin_p;
  nn[2]= v[0]*sin_t*cos_p                    +v[1]*sin_t*sin_p                    +v[2]*cos_t;

  calc_fgb_eh_dv(te,th,tde,tdh,xn,nn,fgb);

  e[0]= te[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+te[1]*sin_p*cos_p*(cos_t-1.0)        +te[2]*sin_t*cos_p;
  e[1]= te[0]*sin_p*cos_p*(cos_t-1.0)        +te[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+te[2]*sin_t*sin_p;
  e[2]=-te[0]*sin_t*cos_p                    -te[1]*sin_t*sin_p                    +te[2]*cos_t;
  h[0]= th[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+th[1]*sin_p*cos_p*(cos_t-1.0)        +th[2]*sin_t*cos_p;
  h[1]= th[0]*sin_p*cos_p*(cos_t-1.0)        +th[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+th[2]*sin_t*sin_p;
  h[2]=-th[0]*sin_t*cos_p                    -th[1]*sin_t*sin_p                    +th[2]*cos_t;

  dedv[0]= tde[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+tde[1]*sin_p*cos_p*(cos_t-1.0)        +tde[2]*sin_t*cos_p;
  dedv[1]= tde[0]*sin_p*cos_p*(cos_t-1.0)        +tde[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+tde[2]*sin_t*sin_p;
  dedv[2]=-tde[0]*sin_t*cos_p                    -tde[1]*sin_t*sin_p                    +tde[2]*cos_t;
  dhdv[0]= tdh[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+tdh[1]*sin_p*cos_p*(cos_t-1.0)        +tdh[2]*sin_t*cos_p;
  dhdv[1]= tdh[0]*sin_p*cos_p*(cos_t-1.0)        +tdh[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+tdh[2]*sin_t*sin_p;
  dhdv[2]=-tdh[0]*sin_t*cos_p                    -tdh[1]*sin_t*sin_p                    +tdh[2]*cos_t;

}

/////////////////////////////////////////////////////////////////////////
void calc_fgb_eh(double complex *e,double complex *h,double *x,Fgb *fgb)
{
  double complex qc,p0,cez,tex,tey,thx,thy;
  double rho2,zt0,zt1,zt2,kz;
  
  rho2=(x[0]*x[0]+x[1]*x[1])*fgb->data.i_w0*fgb->data.i_w0;
  zt0=x[0]*fgb->data.i_l;
  zt1=x[1]*fgb->data.i_l;
  zt2=x[2]*fgb->data.i_l;
  qc=1.0/(I+2.0*zt2);
  p0=I*qc*cexp(-I*qc*rho2);
  kz=fgb->data.ki*x[2];
  cez=cos(kz)+I*sin(kz);

  tex=fgb->data.ex*conj(p0)*cez;
  tey=fgb->data.ey*conj(p0)*cez;
  thx=fgb->ni*tex;
  thy=fgb->ni*tey;

  e[0]= tex;
  e[1]= tey;
  e[2]=-2.0*conj(qc)*(zt0*tex-zt1*tey);
  h[0]=-thy;
  h[1]= thx;
  h[2]=-2.0*conj(qc)*(zt0*thy+zt1*thx);
}

void calc_fgb_eh_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,Fgb *fgb)
{
  double complex qc,cqc,p0,cez,tex,tey,thx,thy;
  double complex c0,c1,ce0,ddx,ddy,ddz,ddn,tec,thc;
  double rho2,zt0,zt1,zt2,kz,i_w02;
  
  i_w02=fgb->data.i_w0*fgb->data.i_w0;
  rho2=(x[0]*x[0]+x[1]*x[1])*i_w02;
  zt0=x[0]*fgb->data.i_l;
  zt1=x[1]*fgb->data.i_l;
  zt2=x[2]*fgb->data.i_l;
  qc=1.0/(I+2.0*zt2);
  cqc=conj(qc);
  p0=I*qc*cexp(-I*qc*rho2);
  kz=fgb->data.ki*x[2];
  cez=cos(kz)+I*sin(kz);

  tex=fgb->data.ex*conj(p0)*cez;
  tey=fgb->data.ey*conj(p0)*cez;
  thx=fgb->ni*tex;
  thy=fgb->ni*tey;

  e[0]= tex;
  e[1]= tey;
  e[2]=-2.0*cqc*(zt0*tex-zt1*tey);
  h[0]=-thy;
  h[1]= thx;
  h[2]=-2.0*cqc*(zt0*thy+zt1*thx);
  
  c0=I*cqc*rho2;
  c1=cqc*cqc;
  ce0=cexp(c0);
  ddx=2.0*x[0]*i_w02*c1*ce0*cez;
  ddy=2.0*x[1]*i_w02*c1*ce0*cez;
  ddz=I*cqc*(2.0*fgb->data.i_l*cqc*(1.0+c0)-I*fgb->data.ki)*ce0*cez;
  ddn=ddx*v[0]+ddy*v[1]+ddz*v[2];
  tec= zt0*fgb->data.ex-zt1*fgb->data.ey;
  thc=(zt0*fgb->data.ey+zt1*fgb->data.ex)*fgb->ni;
  
  dedv[0]= fgb->data.ex*ddn;
  dedv[1]= fgb->data.ey*ddn;
  dedv[2]=-2.0*cqc*( tex*fgb->data.i_l+ddx*tec)*v[0]
          -2.0*cqc*(-tey*fgb->data.i_l+ddy*tec)*v[1]
          -2.0*cqc*(e[2]*fgb->data.i_l+ddz*tec)*v[2];

  dhdv[0]=-fgb->ni*dedv[1];
  dhdv[1]= fgb->ni*dedv[0];
  dhdv[2]=-2.0*cqc*( thy*fgb->data.i_l+ddx*thc)*v[0]
          -2.0*cqc*( thx*fgb->data.i_l+ddy*thc)*v[1]
          -2.0*cqc*(h[2]*fgb->data.i_l+ddz*thc)*v[2];
}
