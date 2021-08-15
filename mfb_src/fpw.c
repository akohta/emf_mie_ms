#include "fpw.h"


void calc_fpw_EH(double complex *e,double complex *h,double *x,Fpw *fpw)
{
  void calc_fpw_eh(double complex *e,double complex *h,double *x,Fpw *fpw);
  
  int i;
  double xn[3],xb,yb,zb;
  double cos_p,sin_p,cos_t,sin_t;
  double complex te[3],th[3];
  xb=x[0]-fpw->fx;  yb=x[1]-fpw->fy;  zb=x[2]-fpw->fz;
  
  cos_p=fpw->data.cos_p;    sin_p=fpw->data.sin_p;
  cos_t=fpw->data.cos_t;    sin_t=fpw->data.sin_t;
  
  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  calc_fpw_eh(te,th,xn,fpw);

  e[0]= te[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+te[1]*sin_p*cos_p*(cos_t-1.0)        +te[2]*sin_t*cos_p;
  e[1]= te[0]*sin_p*cos_p*(cos_t-1.0)        +te[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+te[2]*sin_t*sin_p;
  e[2]=-te[0]*sin_t*cos_p                    -te[1]*sin_t*sin_p                    +te[2]*cos_t;
  h[0]= th[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+th[1]*sin_p*cos_p*(cos_t-1.0)        +th[2]*sin_t*cos_p;
  h[1]= th[0]*sin_p*cos_p*(cos_t-1.0)        +th[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+th[2]*sin_t*sin_p;
  h[2]=-th[0]*sin_t*cos_p                    -th[1]*sin_t*sin_p                    +th[2]*cos_t;
  
  for(i=0;i<3;i++){
    e[i]*=fpw->data.E0;    h[i]*=fpw->data.E0;
  }
}

void calc_fpw_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,Fpw *fpw)
{
  void calc_fpw_eh_dv(double complex *e,double complex *h,double complex *de,double complex *dh,double *x,double *v,Fpw *fpw);
          
  int i;
  double xn[3],xb,yb,zb,nn[3];
  double cos_p,sin_p,cos_t,sin_t;
  double complex te[3],th[3],tde[3],tdh[3];
  xb=x[0]-fpw->fx;  yb=x[1]-fpw->fy;  zb=x[2]-fpw->fz;

  cos_p=fpw->data.cos_p;    sin_p=fpw->data.sin_p;
  cos_t=fpw->data.cos_t;    sin_t=fpw->data.sin_t;

  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  nn[0]= v[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+v[1]*sin_p*cos_p*(cos_t-1.0)        -v[2]*sin_t*cos_p;
  nn[1]= v[0]*sin_p*cos_p*(cos_t-1.0)        +v[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)-v[2]*sin_t*sin_p;
  nn[2]= v[0]*sin_t*cos_p                    +v[1]*sin_t*sin_p                    +v[2]*cos_t;

  calc_fpw_eh_dv(te,th,tde,tdh,xn,nn,fpw);

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
    e[i]*=fpw->data.E0;      h[i]*=fpw->data.E0;
    dedv[i]*=fpw->data.E0;   dhdv[i]*=fpw->data.E0;
  }
}

void print_data_fpw(Fpw *fpw)
{
  printf("-- focused plane wave --\n");
  printf("wave length of incident beam in vacuum      : %15.14g\n",fpw->lambda0);
  printf("refractive index of surrounding             : %15.14g\n",fpw->ni);
  printf("numerical aperture                          : %15.14g\n",fpw->NA);
  printf("incident beam power                         : %15.14g\n",fpw->power);
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fpw->e0x),cimag(fpw->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fpw->e0y),cimag(fpw->e0y));
  printf("x-component of translation vector           : %15.14g\n",fpw->fx);
  printf("y-component of translation vector           : %15.14g\n",fpw->fy);
  printf("z-component of translation vector           : %15.14g\n",fpw->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",fpw->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",fpw->phi);
  printf("gauss-legendre integration sampling number  : %15d\n",fpw->nn);
}

void print_data_fpw_mksa(Fpw *fpw)
{
  printf("-- focused plane wave, MKSA system --\n");
  printf("wave length of incident beam in vacuum   [m]: %15.14g\n",OSUtoMKSA_length(fpw->lambda0));
  printf("refractive index of surrounding             : %15.14g\n",fpw->ni);
  printf("numerical aperture                          : %15.14g\n",fpw->NA);
  printf("incident beam power                      [W]: %15.14g\n",OSUtoMKSA_power(fpw->power));
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fpw->e0x),cimag(fpw->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(fpw->e0y),cimag(fpw->e0y));
  printf("x-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(fpw->fx));
  printf("y-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(fpw->fy));
  printf("z-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(fpw->fz));
  printf("rotation parameter theta               [rad]: %15.14g\n",fpw->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",fpw->phi);
  printf("gauss-legendre integration sampling number  : %15d\n",fpw->nn);  
}

void setup_Fpw(Fpw *fpw)
{
  double kr,ce,*xt,h,apd;
  int i,j,s,st,nnc;
  
  kr=fpw->NA/fpw->ni;
  fpw->data.ki=2.0*M_PI*fpw->ni/fpw->lambda0;
  ce=sqrt(1.0/M_PI)/(kr*sqrt(pow(cabs(fpw->e0x),2)+pow(cabs(fpw->e0y),2)));
  fpw->data.ex=fpw->e0x*ce;
  fpw->data.ey=fpw->e0y*ce;
  fpw->data.hx=-fpw->data.ey*fpw->ni;
  fpw->data.hy= fpw->data.ex*fpw->ni;
  fpw->data.E0=sqrt(2.0*fpw->power*fpw->ni)/fpw->lambda0;

  fpw->data.cos_t=cos(fpw->theta);
  fpw->data.sin_t=sin(fpw->theta);
  fpw->data.cos_p=cos(fpw->phi);
  fpw->data.sin_p=sin(fpw->phi);
  fpw->data.nt=fpw->nn;
  // gau-leg 
  fpw->data.ct=(double *)m_alloc2(fpw->data.nt,sizeof(double),"fpw.c,setup_Fpw(),fpw->data.ct");
  fpw->data.st=(double *)m_alloc2(fpw->data.nt,sizeof(double),"fpw.c,setup_Fpw(),fpw->data.st");
  fpw->data.wt=(double *)m_alloc2(fpw->data.nt,sizeof(double),"fpw.c,setup_Fpw(),fpw->data.wt");
  xt=(double *)m_alloc2(fpw->data.nt,sizeof(double),"fpw.c,setup_Fpw(),fpw->data.xt");
  gauleg( 0.0 ,kr,xt,fpw->data.wt,fpw->data.nt);
  for(i=0;i<fpw->data.nt;i++){    fpw->data.st[i]=xt[i];         fpw->data.ct[i]=sqrt(1.0-xt[i]*xt[i]);  fpw->data.wt[i]*=xt[i];  } 
  free(xt);                                         // free
  // trapezoid
  nnc=(2<<(TRAP_MNUM))+1;
  fpw->data.cp=(double *)m_alloc2(nnc,sizeof(double),"fpw.c,setup_Fpw(),fpw->data.cp");
  fpw->data.sp=(double *)m_alloc2(nnc,sizeof(double),"fpw.c,setup_Fpw(),fpw->data.sp");
  h=2.0*M_PI;
  fpw->data.cp[0]=-1.0;
  fpw->data.sp[0]= 0.0;
  j=1;
  st=1;
  for(i=0;i<=TRAP_MNUM;i++){
    h*=0.5;
    for(s=1;s<=st;s++){
      apd=-M_PI+(double)(2*s-1)*h;
      fpw->data.cp[j]=cos(apd);
      fpw->data.sp[j]=sin(apd);
      j++;
    }
    st*=2;
  }
}

void free_Fpw(Fpw *fpw)
{
  free(fpw->data.ct);
  free(fpw->data.st);
  free(fpw->data.wt);
  free(fpw->data.cp);
  free(fpw->data.sp);
}


////////////////////////////////////////////////////////////////////
void calc_fpw_eh(double complex *e,double complex *h,double *x,Fpw *fpw)
{
  void int_phi_fpw(double complex *eh,double ct,double st,double *x,Fpw *fpw);
  
  double complex eh[6];
  int j;
  
  e[0]=0.0;  e[1]=0.0;  e[2]=0.0;
  h[0]=0.0;  h[1]=0.0;  h[2]=0.0;
  for(j=0;j<fpw->data.nt;j++){
    int_phi_fpw(eh,fpw->data.ct[j],fpw->data.st[j],x,fpw);
    
    e[0]+=eh[0]*fpw->data.wt[j];    e[1]+=eh[1]*fpw->data.wt[j];    e[2]+=eh[2]*fpw->data.wt[j];
    h[0]+=eh[3]*fpw->data.wt[j];    h[1]+=eh[4]*fpw->data.wt[j];    h[2]+=eh[5]*fpw->data.wt[j];
  }
}

void int_phi_fpw(double complex *eh,double ct,double st,double *x,Fpw *fpw)
{
  void ac_eh_fpw(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,Fpw *fpw);
  
  double h,cpa,spa,abst,abss;
  double complex Ia[6],It[6],Is[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=fpw->ni*fpw->ni;
  double cmu=1.0;

  cc=0;
  h=2.0*M_PI;
  cpa=fpw->data.cp[0];  spa=fpw->data.sp[0];
  ac_eh_fpw(Ia,ct,st,cpa,spa,x,fpw);
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
      cpa=fpw->data.cp[j];  spa=fpw->data.sp[j];
      j++;
      ac_eh_fpw(Ia,ct,st,cpa,spa,x,fpw);
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

void ac_eh_fpw(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,Fpw *fpw)
{
  double tmp=fpw->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=(cos(tmp)+I*sin(tmp));

  eh[0]=( fpw->data.ex*(sin_p*sin_p+cos_t*cos_p*cos_p)+fpw->data.ey*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[1]=( fpw->data.ex*sin_p*cos_p*(cos_t-1.0)        +fpw->data.ey*(cos_t*sin_p*sin_p+cos_p*cos_p))*ce;
  eh[2]=(-fpw->data.ex*sin_t*cos_p                    -fpw->data.ey*sin_t*sin_p)*ce;
  
  eh[3]=( fpw->data.hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+fpw->data.hy*sin_p*cos_p*(cos_t-1.0))*ce;                                            
  eh[4]=( fpw->data.hx*sin_p*cos_p*(cos_t-1.0)        +fpw->data.hy*(sin_p*sin_p*cos_t+cos_p*cos_p))*ce;                                    
  eh[5]=(-fpw->data.hx*sin_t*cos_p                    -fpw->data.hy*sin_t*sin_p)*ce;                                              
}

void calc_fpw_eh_dv(double complex *e,double complex *h,double complex *de,double complex *dh,double *x,double *v,Fpw *fpw)
{
  void int_phi_fpw_dv(double complex *eh,double complex *deh,double ct,double st,double *x,double *v,Fpw *fpw);
  
  double complex Et[6],dEt[6];
  double cos_t,sin_t;
  int i,j;

  for(i=0;i<3;i++){
    e[i]=0.0;    h[i]=0.0;
    de[i]=0.0;    dh[i]=0.0;
  }
  for(j=0;j<fpw->data.nt;j++){
    cos_t=fpw->data.ct[j];
    sin_t=fpw->data.st[j];
    
    int_phi_fpw_dv(Et,dEt,cos_t,sin_t,x,v,fpw);

    e[0]+=Et[0]*fpw->data.wt[j];
    e[1]+=Et[1]*fpw->data.wt[j];
    e[2]+=Et[2]*fpw->data.wt[j];
    h[0]+=Et[3]*fpw->data.wt[j];
    h[1]+=Et[4]*fpw->data.wt[j];
    h[2]+=Et[5]*fpw->data.wt[j];

    de[0]+=dEt[0]*fpw->data.wt[j];
    de[1]+=dEt[1]*fpw->data.wt[j];
    de[2]+=dEt[2]*fpw->data.wt[j];
    dh[0]+=dEt[3]*fpw->data.wt[j];
    dh[1]+=dEt[4]*fpw->data.wt[j];
    dh[2]+=dEt[5]*fpw->data.wt[j];
  } 
  
}

void int_phi_fpw_dv(double complex *eh,double complex *deh,double ct,double st,double *x,double *v,Fpw *fpw)
{
  void ac_eh_fpw_dv(double complex *eh,double complex *deh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,double *v,Fpw *fpw);
  
  double h,cpa,spa,abst,abss;
  double complex Ia[6],It[6],Is[6],dIa[6],dIt[6],dIs[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=fpw->ni*fpw->ni;
  double cmu=1.0;

  cc=0;
  h=2.0*M_PI;
  cpa=fpw->data.cp[0];  spa=fpw->data.sp[0];
  ac_eh_fpw_dv(Ia,dIa,ct,st,cpa,spa,x,v,fpw);
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
      cpa=fpw->data.cp[j];  spa=fpw->data.sp[j];
      j++;
      ac_eh_fpw_dv(Ia,dIa,ct,st,cpa,spa,x,v,fpw);
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

void ac_eh_fpw_dv(double complex *eh,double complex *deh,double cos_t,double sin_t,double cos_p,double sin_p,double *x,double *v,Fpw *fpw)
{
  double tmp=fpw->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=(cos(tmp)+I*sin(tmp));
  double complex cdx,cdy,cdz,gfc;

  cdx=I*fpw->data.ki*sin_t*cos_p;
  cdy=I*fpw->data.ki*sin_t*sin_p;
  cdz=I*fpw->data.ki*cos_t;
  gfc=cdx*v[0]+cdy*v[1]+cdz*v[2];

  eh[0]=( fpw->data.ex*(sin_p*sin_p+cos_t*cos_p*cos_p)+fpw->data.ey*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[1]=( fpw->data.ex*sin_p*cos_p*(cos_t-1.0)        +fpw->data.ey*(cos_t*sin_p*sin_p+cos_p*cos_p))*ce;
  eh[2]=(-fpw->data.ex*sin_t*cos_p                    -fpw->data.ey*sin_t*sin_p)*ce;

  eh[3]=( fpw->data.hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+fpw->data.hy*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[4]=( fpw->data.hx*sin_p*cos_p*(cos_t-1.0)        +fpw->data.hy*(sin_p*sin_p*cos_t+cos_p*cos_p))*ce;
  eh[5]=(-fpw->data.hx*sin_t*cos_p                    -fpw->data.hy*sin_t*sin_p)*ce;

  deh[0]=eh[0]*gfc;
  deh[1]=eh[1]*gfc;
  deh[2]=eh[2]*gfc;
  deh[3]=eh[3]*gfc;
  deh[4]=eh[4]*gfc;
  deh[5]=eh[5]*gfc;
}
