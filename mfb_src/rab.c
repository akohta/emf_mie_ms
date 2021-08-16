#include "rab.h"

void calc_rab_EH(double _Complex *e,double _Complex *h,double *x,RAb *rab)
{
  void calc_rab_eh(double complex *e,double complex *h,double *x,RAb *rab);
  
  int i;
  double xn[3],xb,yb,zb;
  double cos_p,sin_p,cos_t,sin_t;
  double _Complex te[3],th[3];
  
  xb=x[0]-rab->fx;  yb=x[1]-rab->fy;  zb=x[2]-rab->fz;
  cos_p=rab->data.cos_p;    sin_p=rab->data.sin_p;
  cos_t=rab->data.cos_t;    sin_t=rab->data.sin_t;
  
  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;
  
  calc_rab_eh(te,th,xn,rab);
  
  e[0]= te[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+te[1]*sin_p*cos_p*(cos_t-1.0)        +te[2]*sin_t*cos_p;
  e[1]= te[0]*sin_p*cos_p*(cos_t-1.0)        +te[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+te[2]*sin_t*sin_p;
  e[2]=-te[0]*sin_t*cos_p                    -te[1]*sin_t*sin_p                    +te[2]*cos_t;
  h[0]= th[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+th[1]*sin_p*cos_p*(cos_t-1.0)        +th[2]*sin_t*cos_p;
  h[1]= th[0]*sin_p*cos_p*(cos_t-1.0)        +th[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+th[2]*sin_t*sin_p;
  h[2]=-th[0]*sin_t*cos_p                    -th[1]*sin_t*sin_p                    +th[2]*cos_t;
  
  for(i=0;i<3;i++){
    e[i]*=rab->data.E0;    h[i]*=rab->data.E0;
  }
}

void calc_rab_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,RAb *rab)
{
  void calc_rab_eh_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,RAb *rab);
  
  int i;
  double xn[3],xb,yb,zb,nn[3];
  double cos_p,sin_p,cos_t,sin_t;
  double complex te[3],th[3],tde[3],tdh[3];
  xb=x[0]-rab->fx;  yb=x[1]-rab->fy;  zb=x[2]-rab->fz;

  cos_p=rab->data.cos_p;    sin_p=rab->data.sin_p;
  cos_t=rab->data.cos_t;    sin_t=rab->data.sin_t;

  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  nn[0]= v[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+v[1]*sin_p*cos_p*(cos_t-1.0)        -v[2]*sin_t*cos_p;
  nn[1]= v[0]*sin_p*cos_p*(cos_t-1.0)        +v[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)-v[2]*sin_t*sin_p;
  nn[2]= v[0]*sin_t*cos_p                    +v[1]*sin_t*sin_p                    +v[2]*cos_t;

  calc_rab_eh_dv(te,th,tde,tdh,xn,nn,rab);

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

  for(i=0;i<3;i++){
    e[i]*=rab->data.E0;    h[i]*=rab->data.E0;
    dedv[i]*=rab->data.E0;    dhdv[i]*=rab->data.E0;
  }
}

void print_data_rab(RAb *rab)
{
  printf("-- focused radially and azimuthally polarized beam --\n");
  printf("wave length of incident beam in vacuum      : %15.14g\n",rab->lambda0);
  printf("refractive index of surrounding             : %15.14g\n",rab->ni);
  printf("numerical aperture                          : %15.14g\n",rab->NA);
  printf("incident beam power                         : %15.14g\n",rab->power);
  printf("radial-component of polarization coefficient:%7.6g+%7.6gi\n",creal(rab->e0r),cimag(rab->e0r));
  printf("azimuthal-component of polarization coef    :%7.6g+%7.6gi\n",creal(rab->e0a),cimag(rab->e0a));
  printf("x-component of translation vector           : %15.14g\n",rab->fx);
  printf("y-component of translation vector           : %15.14g\n",rab->fy);
  printf("z-component of translation vector           : %15.14g\n",rab->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",rab->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",rab->phi);
  printf("gauss-legendre integration sampling number  : %15d\n",rab->nn);
}

void print_data_rab_mksa(RAb *rab)
{
  printf("-- focused radially and azimuthally polarized beam, MKSA system --\n");
  printf("wave length of incident beam in vacuum   [m]: %15.14g\n",OSUtoMKSA_length(rab->lambda0));
  printf("refractive index of surrounding             : %15.14g\n",rab->ni);
  printf("numerical aperture                          : %15.14g\n",rab->NA);
  printf("incident beam power                      [W]: %15.14g\n",OSUtoMKSA_power(rab->power));
  printf("radial-component of polarization coefficient:%7.6g+%7.6gi\n",creal(rab->e0r),cimag(rab->e0r));
  printf("azimuthal-component of polarization coef    :%7.6g+%7.6gi\n",creal(rab->e0a),cimag(rab->e0a));
  printf("x-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(rab->fx));
  printf("y-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(rab->fy));
  printf("z-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(rab->fz));
  printf("rotation parameter theta               [rad]: %15.14g\n",rab->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",rab->phi);
  printf("gauss-legendre integration sampling number  : %15d\n",rab->nn);
}

void setup_rab(RAb *rab)
{
  double kr,ce,*xt,h,apd;
  int i,j,s,st,nnc;

  kr=rab->NA/rab->ni;
  rab->data.ki=2.0*M_PI*rab->ni/rab->lambda0;
  ce=sqrt(1.0/M_PI)/(kr*sqrt(pow(cabs(rab->e0r),2)+pow(cabs(rab->e0a),2)));
  rab->data.er= rab->e0r*ce;
  rab->data.ea= rab->e0a*ce;
  rab->data.hr=-rab->data.ea*rab->ni;
  rab->data.ha= rab->data.er*rab->ni;
  rab->data.E0=sqrt(2.0*rab->power*rab->ni)/rab->lambda0;

  rab->data.cos_t=cos(rab->theta);
  rab->data.sin_t=sin(rab->theta);
  rab->data.cos_p=cos(rab->phi);
  rab->data.sin_p=sin(rab->phi);
  rab->data.nt=rab->nn;
  // memory allocation 
  rab->data.ct=(double *)m_alloc2(rab->data.nt,sizeof(double),"rab.c,setup_rab(),rab->data.ct");
  rab->data.st=(double *)m_alloc2(rab->data.nt,sizeof(double),"rab.c,setup_rab(),rab->data.st");
  rab->data.wt=(double *)m_alloc2(rab->data.nt,sizeof(double),"rab.c,setup_rab(),rab->data.wt");
  // gauleg data
  xt=(double *)m_alloc2(rab->data.nt,sizeof(double),"rab.c,setup_rab(),xt");
  gauleg( 0.0 ,kr   ,xt,rab->data.wt,rab->data.nt);
  for(i=0;i<rab->data.nt;i++){    rab->data.st[i]=xt[i];         rab->data.ct[i]=sqrt(1.0-xt[i]*xt[i]);    rab->data.wt[i]*=xt[i];  } 
  free(xt);                                                
  // trapezoidral data
  nnc=(2<<(TRAP_MNUM))+1;
  rab->data.cp=(double *)m_alloc2(nnc,sizeof(double),"rab.c,setup_rab(),rab->data.cp");
  rab->data.sp=(double *)m_alloc2(nnc,sizeof(double),"rab.c,setup_rab(),rab->data.sp");
  h=2.0*M_PI;
  rab->data.cp[0]=-1.0;
  rab->data.sp[0]= 0.0;
  j=1;
  st=1;
  for(i=0;i<=TRAP_MNUM;i++){
    h*=0.5;
    for(s=1;s<=st;s++){
      apd=-M_PI+(double)(2*s-1)*h;
      rab->data.cp[j]=cos(apd);
      rab->data.sp[j]=sin(apd);
      j++;
    }
    st*=2;
  }
}

void free_rab(RAb *rab)
{
  free(rab->data.ct);
  free(rab->data.st);
  free(rab->data.wt);
  free(rab->data.cp);
  free(rab->data.sp);
}

/////////////////////////////////////////////////////////////////////////
void calc_rab_eh(double complex *e,double complex *h,double *x,RAb *rab)
{
  void int_phi_rab(double complex *eh,double ct,double st,double *x,RAb *rab);

  double complex eh[6];
  int j,n;
  
  for(n=0;n<3;n++){
    e[n]=0.0;    h[n]=0.0;
  }
  
  for(j=0;j<rab->data.nt;j++){
  int_phi_rab(eh,rab->data.ct[j],rab->data.st[j],x,rab);
  
  e[0]+=eh[0]*rab->data.wt[j];    e[1]+=eh[1]*rab->data.wt[j];    e[2]+=eh[2]*rab->data.wt[j];
    h[0]+=eh[3]*rab->data.wt[j];    h[1]+=eh[4]*rab->data.wt[j];    h[2]+=eh[5]*rab->data.wt[j];
  }
}

void int_phi_rab(double complex *eh,double ct,double st,double *x,RAb *rab)
{
  void ac_eh_rab(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,RAb *rab);

  double h,abst,abss;
  double complex Ia[6],It[6],Is[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=rab->ni*rab->ni;
  double cmu=1.0;

  cc=0;
  h=2.0*M_PI;
  ac_eh_rab(Ia,ct,st,rab->data.cp[0],rab->data.sp[0],x,rab);
  abst=0.0;
  for(i=0;i<6;i++){
    It[i]=h*Ia[i];
    if(i<3) abst+=cep*creal(It[i]*conj(It[i]));
    else    abst+=cmu*creal(It[i]*conj(It[i]));
  }
  f_order=abst;
  j=1;
  sc=1;
  while(cc<=TRAP_MNUM){
    h*=0.5;
    for(i=0;i<6;i++) Is[i]=0.0;
    for(s=1;s<=sc;s++){
      ac_eh_rab(Ia,ct,st,rab->data.cp[j],rab->data.sp[j],x,rab);
      j++;
      for(i=0;i<6;i++) Is[i]+=Ia[i];
    }
    abss=0.0;
    for(i=0;i<6;i++){
      It[i]=0.5*It[i]+Is[i]*h;
      if(i<3) abss+=cep*creal(It[i]*conj(It[i]));
      else    abss+=cmu*creal(It[i]*conj(It[i]));
    }
    cr_coef=sqrt((double)sc)*f_order;
    if((cc>TRAP_LNUM && fabs(abst-abss) < TRAP_EPS*cr_coef) || (abss==0.0 && abst==0.0) ){
      break;
    }
    abst=abss;
    if(abst>f_order) f_order=abst;
    sc*=2;
    cc++;
  }
  if(cc>=TRAP_MNUM) {
    printf("phi integration limit over! Exit...\n");
    exit(1);
  }
  for(i=0;i<6;i++) eh[i]=It[i];
}

void ac_eh_rab(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,RAb *rab)
{
  double tmp=rab->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=(cos(tmp)+I*sin(tmp));

  eh[0]=( rab->data.er*cos_t*cos_p-rab->data.ea*sin_p)*ce;
  eh[1]=( rab->data.er*cos_t*sin_p+rab->data.ea*cos_p)*ce;
  eh[2]=(-rab->data.er*sin_t)*ce;

  eh[3]=( rab->data.hr*cos_t*cos_p-rab->data.ha*sin_p)*ce;
  eh[4]=( rab->data.hr*cos_t*sin_p+rab->data.ha*cos_p)*ce;
  eh[5]=(-rab->data.hr*sin_t)*ce;
}

void calc_rab_eh_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,RAb *rab)
{
  void int_phi_rab_dv(double complex *eh,double complex *deh,double ct,double st,double *x,double *v,RAb *rab);
  
  double complex Et[6],dEt[6];
  int i,j;

  for(i=0;i<3;i++){
    e[i]=0.0;    h[i]=0.0;
    dedv[i]=0.0;    dhdv[i]=0.0;
  }
  for(j=0;j<rab->data.nt;j++){
    int_phi_rab_dv(Et,dEt,rab->data.ct[j],rab->data.st[j],x,v,rab);

    e[0]+=Et[0]*rab->data.wt[j];
    e[1]+=Et[1]*rab->data.wt[j];
    e[2]+=Et[2]*rab->data.wt[j];
    h[0]+=Et[3]*rab->data.wt[j];
    h[1]+=Et[4]*rab->data.wt[j];
    h[2]+=Et[5]*rab->data.wt[j];

    dedv[0]+=dEt[0]*rab->data.wt[j];
    dedv[1]+=dEt[1]*rab->data.wt[j];
    dedv[2]+=dEt[2]*rab->data.wt[j];
    dhdv[0]+=dEt[3]*rab->data.wt[j];
    dhdv[1]+=dEt[4]*rab->data.wt[j];
    dhdv[2]+=dEt[5]*rab->data.wt[j];
  }
}

void int_phi_rab_dv(double complex *eh,double complex *deh,double ct,double st,double *x,double *v,RAb *rab)
{
  void ac_eh_rab_dv(double complex *eh,double complex *deh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,double *v,RAb *rab);
  
  double h,abst,abss;
  double complex Ia[6],It[6],Is[6],dIa[6],dIt[6],dIs[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=rab->ni*rab->ni;
  double cmu=1.0;

  cc=0;
  h=2.0*M_PI;
  ac_eh_rab_dv(Ia,dIa,ct,st,rab->data.cp[0],rab->data.sp[0],x,v,rab);
  abst=0.0;
  for(i=0;i<6;i++){
    It[i]=h*Ia[i];
    dIt[i]=h*dIa[i];
    if(i<3) abst+=cep*creal(It[i]*conj(It[i]));
    else    abst+=cmu*creal(It[i]*conj(It[i]));
  }
  f_order=abst;
  j=1;
  sc=1;
  while(cc<=TRAP_MNUM){
    h*=0.5;
    for(i=0;i<6;i++){
      Is[i]=0.0;
      dIs[i]=0.0;
    }
    for(s=1;s<=sc;s++){
      ac_eh_rab_dv(Ia,dIa,ct,st,rab->data.cp[j],rab->data.sp[j],x,v,rab);
      j++;
      for(i=0;i<6;i++){
        Is[i]+=Ia[i];
        dIs[i]+=dIa[i];
      }
    }
    abss=0.0;
    for(i=0;i<6;i++){
      It[i]=0.5*It[i]+Is[i]*h;
      dIt[i]=0.5*dIt[i]+dIs[i]*h;
      if(i<3) abss+=cep*creal(It[i]*conj(It[i]));
      else    abss+=cmu*creal(It[i]*conj(It[i]));
    }
    cr_coef=sqrt((double)sc)*f_order;
    if((cc>TRAP_LNUM && fabs(abst-abss) < TRAP_EPS*cr_coef) || (abss==0.0 && abst==0.0) ){
      break;
    }
    abst=abss;
    if(abst>f_order) f_order=abst;
    sc*=2;
    cc++;
  }
  if(cc>=TRAP_MNUM) {
    printf("phi integration limit over! Exit...\n");
    exit(1);
  }
  for(i=0;i<6;i++){
    eh[i]=It[i];
    deh[i]=dIt[i];
  }
}

void ac_eh_rab_dv(double complex *eh,double complex *deh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,double *v,RAb *rab)
{
  double tmp=rab->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=(cos(tmp)+I*sin(tmp));
  double complex cdx,cdy,cdz,gfc;

  cdx=I*rab->data.ki*sin_t*cos_p;
  cdy=I*rab->data.ki*sin_t*sin_p;
  cdz=I*rab->data.ki*cos_t;
  gfc=cdx*v[0]+cdy*v[1]+cdz*v[2];

  eh[0]=( rab->data.er*cos_t*cos_p-rab->data.ea*sin_p)*ce;
  eh[1]=( rab->data.er*cos_t*sin_p+rab->data.ea*cos_p)*ce;
  eh[2]=(-rab->data.er*sin_t)*ce;

  eh[3]=( rab->data.hr*cos_t*cos_p-rab->data.ha*sin_p)*ce;
  eh[4]=( rab->data.hr*cos_t*sin_p+rab->data.ha*cos_p)*ce;
  eh[5]=(-rab->data.hr*sin_t)*ce;

  deh[0]=eh[0]*gfc;
  deh[1]=eh[1]*gfc;
  deh[2]=eh[2]*gfc;
  deh[3]=eh[3]*gfc;
  deh[4]=eh[4]*gfc;
  deh[5]=eh[5]*gfc;
}
