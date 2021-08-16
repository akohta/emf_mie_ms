#include "lgb.h"

void calc_lgb_EH(double _Complex *e,double _Complex *h,double *x,LGb *lgb)
{
  void calc_lgb_eh(double complex *e,double complex *h,double *x,LGb *lgb);
  
  int i;
  double xn[3],xb,yb,zb;
  double cos_p,sin_p,cos_t,sin_t;
  double _Complex te[3],th[3];
  xb=x[0]-lgb->fx;  yb=x[1]-lgb->fy;  zb=x[2]-lgb->fz;

  cos_p=lgb->data.cos_p;    sin_p=lgb->data.sin_p;
  cos_t=lgb->data.cos_t;    sin_t=lgb->data.sin_t;
  
  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  calc_lgb_eh(te,th,xn,lgb);

  e[0]= te[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+te[1]*sin_p*cos_p*(cos_t-1.0)        +te[2]*sin_t*cos_p;
  e[1]= te[0]*sin_p*cos_p*(cos_t-1.0)        +te[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+te[2]*sin_t*sin_p;
  e[2]=-te[0]*sin_t*cos_p                    -te[1]*sin_t*sin_p                    +te[2]*cos_t;
  h[0]= th[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+th[1]*sin_p*cos_p*(cos_t-1.0)        +th[2]*sin_t*cos_p;
  h[1]= th[0]*sin_p*cos_p*(cos_t-1.0)        +th[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+th[2]*sin_t*sin_p;
  h[2]=-th[0]*sin_t*cos_p                    -th[1]*sin_t*sin_p                    +th[2]*cos_t;

  for(i=0;i<3;i++){
    e[i]*=lgb->data.E0;
    h[i]*=lgb->data.E0;
  }
}

void calc_lgb_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,LGb *lgb)
{
  void calc_lgb_eh_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,LGb *lgb);
  
  int i;
  double xn[3],xb,yb,zb,nn[3];
  double cos_p,sin_p,cos_t,sin_t;
  double complex te[3],th[3],tde[3],tdh[3];
  xb=x[0]-lgb->fx;  yb=x[1]-lgb->fy;  zb=x[2]-lgb->fz;

  cos_p=lgb->data.cos_p;    sin_p=lgb->data.sin_p;
  cos_t=lgb->data.cos_t;    sin_t=lgb->data.sin_t;

  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  nn[0]= v[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+v[1]*sin_p*cos_p*(cos_t-1.0)        -v[2]*sin_t*cos_p;
  nn[1]= v[0]*sin_p*cos_p*(cos_t-1.0)        +v[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)-v[2]*sin_t*sin_p;
  nn[2]= v[0]*sin_t*cos_p                    +v[1]*sin_t*sin_p                    +v[2]*cos_t;

  calc_lgb_eh_dv(te,th,tde,tdh,xn,nn,lgb);

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
    e[i]*=lgb->data.E0;
    h[i]*=lgb->data.E0;
    dedv[i]*=lgb->data.E0;
    dhdv[i]*=lgb->data.E0;
  }
}

void print_data_lgb(LGb *lgb)
{
  printf("-- focused plane wave with spiral phase modulation --\n");
  printf("wave length of incident beam in vacuum      : %15.14g\n",lgb->lambda0);
  printf("refractive index of surrounding             : %15.14g\n",lgb->ni);
  printf("numerical aperture                          : %15.14g\n",lgb->NA);
  printf("incident beam power                         : %15.14g\n",lgb->power);
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(lgb->e0x),cimag(lgb->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(lgb->e0y),cimag(lgb->e0y));
  printf("x-component of translation vector           : %15.14g\n",lgb->fx);
  printf("y-component of translation vector           : %15.14g\n",lgb->fy);
  printf("z-component of translation vector           : %15.14g\n",lgb->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",lgb->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",lgb->phi);
  printf("mode number for spiral phase modulation     : %15d\n",lgb->lg_m);
  printf("gauss-legendre integration sampling number  : %15d\n",lgb->nn);
}

void print_data_lgb_mksa(LGb *lgb)
{
  printf("-- focused plane wave with spiral phase modulation, MKSA system --\n");
  printf("wave length of incident beam in vacuum   [m]: %15.14g\n",OSUtoMKSA_length(lgb->lambda0));
  printf("refractive index of surrounding             : %15.14g\n",lgb->ni);
  printf("numerical aperture                          : %15.14g\n",lgb->NA);
  printf("incident beam power                      [W]: %15.14g\n",OSUtoMKSA_power(lgb->power));
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(lgb->e0x),cimag(lgb->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(lgb->e0y),cimag(lgb->e0y));
  printf("x-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(lgb->fx));
  printf("y-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(lgb->fy));
  printf("z-component of translation vector        [m]: %15.14g\n",OSUtoMKSA_length(lgb->fz));
  printf("rotation parameter theta               [rad]: %15.14g\n",lgb->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",lgb->phi);
  printf("mode number for spiral phase modulation     : %15d\n",lgb->lg_m);
  printf("gauss-legendre integration sampling number  : %15d\n",lgb->nn);
}

void setup_LGb(LGb *lgb)
{
  double complex phase_lgb(int mode,double phi);
  
  double kr,ce,*xt,h,apd;
  int i,j,s,st,nnc;
  
  kr=lgb->NA/lgb->ni;
  lgb->data.ki=2.0*M_PI*lgb->ni/lgb->lambda0;
  ce=sqrt(1.0/M_PI)/(kr*sqrt(pow(cabs(lgb->e0x),2)+pow(cabs(lgb->e0y),2)));
  lgb->data.ex= lgb->e0x*ce;
  lgb->data.ey= lgb->e0y*ce;
  lgb->data.hx=-lgb->data.ey*lgb->ni;
  lgb->data.hy= lgb->data.ex*lgb->ni;
  lgb->data.E0=sqrt(2.0*lgb->power*lgb->ni)/lgb->lambda0;
  lgb->data.cos_t=cos(lgb->theta);  lgb->data.sin_t=sin(lgb->theta);
  lgb->data.cos_p=cos(lgb->phi);    lgb->data.sin_p=sin(lgb->phi);
  lgb->data.nt=lgb->nn;
  // gau-leg data 
  lgb->data.ct=(double *)m_alloc2(lgb->data.nt,sizeof(double),"lgb.c,setup_LGb(),lgb->data.ct");
  lgb->data.st=(double *)m_alloc2(lgb->data.nt,sizeof(double),"lgb.c,setup_LGb(),lgb->data.st");
  lgb->data.wt=(double *)m_alloc2(lgb->data.nt,sizeof(double),"lgb.c,setup_LGb(),lgb->data.wt");
  xt=(double *)m_alloc2(lgb->data.nt,sizeof(double),"lgb.c,setup_LGb(),lgb->data.xt");
  gauleg(  0.0,kr,xt,lgb->data.wt,lgb->data.nt);
  for(i=0;i<lgb->data.nt;i++){    lgb->data.st[i]=xt[i];         lgb->data.ct[i]=sqrt(1.0-xt[i]*xt[i]);  lgb->data.wt[i]*=xt[i];  } 
  free(xt);                                                
  // trapezoid data
  nnc=(2<<(TRAP_MNUM))+1;
  lgb->data.cp=(double *)m_alloc2(nnc,sizeof(double),"lgb.c,setup_LGb(),lgb->data.cp");
  lgb->data.sp=(double *)m_alloc2(nnc,sizeof(double),"lgb.c,setup_LGb(),lgb->data.sp");
  lgb->data.ph_p=(double complex *)m_alloc2(nnc,sizeof(double complex),"lgb.c,setup_LGb(),lgb->data.ph_p");
  h=2.0*M_PI;
  lgb->data.cp[0]=-1.0;
  lgb->data.sp[0]= 0.0;
  lgb->data.ph_p[0]=phase_lgb(lgb->lg_m,-M_PI);
  j=1;
  st=1;
  for(i=0;i<=TRAP_MNUM;i++){
    h*=0.5;
    for(s=1;s<=st;s++){
      apd=-M_PI+(double)(2*s-1)*h;
      lgb->data.cp[j]=cos(apd);
      lgb->data.sp[j]=sin(apd);
      lgb->data.ph_p[j]=phase_lgb(lgb->lg_m,apd);
      j++;
    }
    st*=2;
  }
}

void free_LGb(LGb *lgb)
{
  free(lgb->data.ct);
  free(lgb->data.st);
  free(lgb->data.wt);
  free(lgb->data.cp);
  free(lgb->data.sp);
  free(lgb->data.ph_p);
}

////////////////////////////////////////////////////////////////////////////
double complex phase_lgb(int mode,double phi)
{
  double tmp,pp;
  if(mode==0) return 1.0;
  else {
    if(phi<0.0) pp=phi+2.0*M_PI;
    else        pp=phi;
    tmp=fmod((double)mode*pp,2.0*M_PI);
    return cos(tmp)+I*sin(tmp);
  }
}

void calc_lgb_eh(double complex *e,double complex *h,double *x,LGb *lgb)
{
  void int_phi_lgb(double _Complex *eh,double ct,double st,double *x,LGb *lgb);
  
  double complex eh[6];
  int j;
  
  e[0]=0.0;  e[1]=0.0;  e[2]=0.0;
  h[0]=0.0;  h[1]=0.0;  h[2]=0.0;

  for(j=0;j<lgb->data.nt;j++){
  int_phi_lgb(eh,lgb->data.ct[j],lgb->data.st[j],x,lgb);
  
    e[0]+=eh[0]*lgb->data.wt[j];    e[1]+=eh[1]*lgb->data.wt[j];    e[2]+=eh[2]*lgb->data.wt[j];
    h[0]+=eh[3]*lgb->data.wt[j];    h[1]+=eh[4]*lgb->data.wt[j];    h[2]+=eh[5]*lgb->data.wt[j];
  }
}

void int_phi_lgb(double _Complex *eh,double ct,double st,double *x,LGb *lgb)
{
  void ac_eh_lgb(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double complex phase,double *x,LGb *lgb);
  
  double h,abst,abss;
  double complex Ia[6],It[6],Is[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=lgb->ni*lgb->ni;
  double cmu=1.0;

  cc=0;
  h=2.0*M_PI;
  ac_eh_lgb(Ia,ct,st,lgb->data.cp[0],lgb->data.sp[0],lgb->data.ph_p[0],x,lgb);
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
      ac_eh_lgb(Ia,ct,st,lgb->data.cp[j],lgb->data.sp[j],lgb->data.ph_p[j],x,lgb);
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

void ac_eh_lgb(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double complex phase,double *x,LGb *lgb)
{
  double tmp=lgb->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=phase*(cos(tmp)+I*sin(tmp));

  eh[0]=( lgb->data.ex*(sin_p*sin_p+cos_t*cos_p*cos_p)+lgb->data.ey*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[1]=( lgb->data.ex*sin_p*cos_p*(cos_t-1.0)        +lgb->data.ey*(cos_t*sin_p*sin_p+cos_p*cos_p))*ce;
  eh[2]=(-lgb->data.ex*sin_t*cos_p                    -lgb->data.ey*sin_t*sin_p)*ce;
  
  eh[3]=( lgb->data.hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+lgb->data.hy*sin_p*cos_p*(cos_t-1.0))*ce;                                                                          
  eh[4]=( lgb->data.hx*sin_p*cos_p*(cos_t-1.0)        +lgb->data.hy*(sin_p*sin_p*cos_t+cos_p*cos_p))*ce;                                                                          
  eh[5]=(-lgb->data.hx*sin_t*cos_p                    -lgb->data.hy*sin_t*sin_p)*ce;                                              
}

void calc_lgb_eh_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,LGb *lgb)
{
  void int_phi_lgb_dv(double complex *eh,double complex *deh,double ct,double st,double *x,double *v,LGb *lgb);

  double complex Et[6],dEt[6];
  double cos_t,sin_t;
  int i,j;

  for(i=0;i<3;i++){
    e[i]=0.0;    h[i]=0.0;
    dedv[i]=0.0;    dhdv[i]=0.0;
  }
  for(j=0;j<lgb->data.nt;j++){
    cos_t=lgb->data.ct[j];
    sin_t=lgb->data.st[j];

    int_phi_lgb_dv(Et,dEt,cos_t,sin_t,x,v,lgb);

    e[0]+=Et[0]*lgb->data.wt[j];
    e[1]+=Et[1]*lgb->data.wt[j];
    e[2]+=Et[2]*lgb->data.wt[j];
    h[0]+=Et[3]*lgb->data.wt[j];
    h[1]+=Et[4]*lgb->data.wt[j];
    h[2]+=Et[5]*lgb->data.wt[j];

    dedv[0]+=dEt[0]*lgb->data.wt[j];
    dedv[1]+=dEt[1]*lgb->data.wt[j];
    dedv[2]+=dEt[2]*lgb->data.wt[j];
    dhdv[0]+=dEt[3]*lgb->data.wt[j];
    dhdv[1]+=dEt[4]*lgb->data.wt[j];
    dhdv[2]+=dEt[5]*lgb->data.wt[j];
  }
}

void int_phi_lgb_dv(double complex *eh,double complex *deh,double ct,double st,double *x,double *v,LGb *lgb)
{
  void ac_eh_lgb_dv(double complex *eh,double complex *deh,double cos_t,double sin_t,double cos_p,double sin_p,double complex phase,double *x,double *v,LGb *lgb);
  
  double h,abst,abss;
  double complex Ia[6],It[6],Is[6],dIa[6],dIt[6],dIs[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=lgb->ni*lgb->ni;
  double cmu=1.0;

  cc=0;
  h=2.0*M_PI;
  ac_eh_lgb_dv(Ia,dIa,ct,st,lgb->data.cp[0],lgb->data.sp[0],lgb->data.ph_p[0],x,v,lgb);
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
      ac_eh_lgb_dv(Ia,dIa,ct,st,lgb->data.cp[j],lgb->data.sp[j],lgb->data.ph_p[j],x,v,lgb);
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

void ac_eh_lgb_dv(double complex *eh,double complex *deh,double cos_t,double sin_t,double cos_p,double sin_p,double complex phase,double *x,double *v,LGb *lgb)
{
  double tmp=lgb->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=phase*(cos(tmp)+I*sin(tmp));
  double complex cdx,cdy,cdz,gfc;

  cdx=I*lgb->data.ki*sin_t*cos_p;
  cdy=I*lgb->data.ki*sin_t*sin_p;
  cdz=I*lgb->data.ki*cos_t;
  gfc=cdx*v[0]+cdy*v[1]+cdz*v[2];

  eh[0]=( lgb->data.ex*(sin_p*sin_p+cos_t*cos_p*cos_p)+lgb->data.ey*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[1]=( lgb->data.ex*sin_p*cos_p*(cos_t-1.0)        +lgb->data.ey*(cos_t*sin_p*sin_p+cos_p*cos_p))*ce;
  eh[2]=(-lgb->data.ex*sin_t*cos_p                    -lgb->data.ey*sin_t*sin_p)*ce;

  eh[3]=( lgb->data.hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+lgb->data.hy*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[4]=( lgb->data.hx*sin_p*cos_p*(cos_t-1.0)        +lgb->data.hy*(sin_p*sin_p*cos_t+cos_p*cos_p))*ce;
  eh[5]=(-lgb->data.hx*sin_t*cos_p                    -lgb->data.hy*sin_t*sin_p)*ce;

  deh[0]=eh[0]*gfc;
  deh[1]=eh[1]*gfc;
  deh[2]=eh[2]*gfc;
  deh[3]=eh[3]*gfc;
  deh[4]=eh[4]*gfc;
  deh[5]=eh[5]*gfc;
}
