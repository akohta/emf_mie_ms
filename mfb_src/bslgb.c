#include "bslgb.h"

void calc_bslgb_EH(double _Complex *e,double _Complex *h,double *x,BsLGb *bsb)
{
  void calc_bslgb_eh(double complex *e,double complex *h,double *x,BsLGb *bsb);
  
  int i;
  double xn[3],xb,yb,zb;
  double cos_p,sin_p,cos_t,sin_t;
  double _Complex te[3],th[3];
  xb=x[0]-bsb->fx;  yb=x[1]-bsb->fy;  zb=x[2]-bsb->fz;

  cos_p=bsb->data.cos_p;    sin_p=bsb->data.sin_p;
  cos_t=bsb->data.cos_t;    sin_t=bsb->data.sin_t;
  
  xn[0]= xb*(cos_t*cos_p*cos_p+sin_p*sin_p)+yb*sin_p*cos_p*(cos_t-1.0)        -zb*sin_t*cos_p;
  xn[1]= xb*sin_p*cos_p*(cos_t-1.0)        +yb*(cos_t*sin_p*sin_p+cos_p*cos_p)-zb*sin_t*sin_p;
  xn[2]= xb*sin_t*cos_p                    +yb*sin_t*sin_p                    +zb*cos_t;

  calc_bslgb_eh(te,th,xn,bsb);
  
  e[0]= te[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+te[1]*sin_p*cos_p*(cos_t-1.0)        +te[2]*sin_t*cos_p;
  e[1]= te[0]*sin_p*cos_p*(cos_t-1.0)        +te[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+te[2]*sin_t*sin_p;
  e[2]=-te[0]*sin_t*cos_p                    -te[1]*sin_t*sin_p                    +te[2]*cos_t;
  h[0]= th[0]*(cos_t*cos_p*cos_p+sin_p*sin_p)+th[1]*sin_p*cos_p*(cos_t-1.0)        +th[2]*sin_t*cos_p;
  h[1]= th[0]*sin_p*cos_p*(cos_t-1.0)        +th[1]*(cos_t*sin_p*sin_p+cos_p*cos_p)+th[2]*sin_t*sin_p;
  h[2]=-th[0]*sin_t*cos_p                    +th[1]*sin_t*sin_p                    +th[2]*cos_t;

  for(i=0;i<3;i++){
    e[i]*=bsb->data.E0;
    h[i]*=bsb->data.E0;
  }
}

void read_data_bslgb(char *rfile,BsLGb *bsb)
{
  FILE *fp;
  if((fp=fopen(rfile,"rt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  char buf[256]="";  int tmpi;  double tmpd,tmpd2;
  fgets(buf,256,fp);  fgets(buf,256,fp);

  printf("-- beam parameter --\n");
  fscanf(fp,"%d",&tmpi);
  if(tmpi!=4){ printf("beam type %d is not supported\n",tmpi); exit(1);}
  else                                                printf("incident beam type                          : Spiral phase modulation Bessel beam \n");
  fscanf(fp,"%lf",&tmpd);  bsb->lambda0=tmpd;         printf("wave length of incident beam in vacuum   [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->ni     =tmpd;         printf("refractive index of surrounding             : %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->d_angle=tmpd;         printf("deflection angle                       [rad]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->power  =tmpd;         printf("incident beam power                  [W/m^2]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);
  fscanf(fp,"%lf",&tmpd2); bsb->e0x    =tmpd+I*tmpd2; printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",tmpd,tmpd2);
  fscanf(fp,"%lf",&tmpd);
  fscanf(fp,"%lf",&tmpd2); bsb->e0y    =tmpd+I*tmpd2; printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",tmpd,tmpd2);
  fscanf(fp,"%lf",&tmpd);  bsb->fx     =tmpd;         printf("x-component of translation vector        [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->fy     =tmpd;         printf("y-component of translation vector        [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->fz     =tmpd;         printf("z-component of translation vector        [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->theta  =tmpd;         printf("rotation parameter theta               [rad]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  bsb->phi    =tmpd;         printf("rotation parameter phi                 [rad]: %15.14g\n",tmpd);
  fscanf(fp,"%d" ,&tmpi);  bsb->lg_m   =tmpi;         printf("mode number for spiral phase modulation     : %15d\n",tmpi);

  fclose(fp);
}

void print_data_bslgb(BsLGb *bsb)
{
  printf("-- Bessel beam with spiral phase modulation --\n");
  printf("wave length of incident beam in vacuum   [m]: %15.14g\n",bsb->lambda0);
  printf("refractive index of surrounding             : %15.14g\n",bsb->ni);
  printf("deflection angle                       [rad]: %15.14g\n",bsb->d_angle);
  printf("incident beam power                  [W/m^2]: %15.14g\n",bsb->power);
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(bsb->e0x),cimag(bsb->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(bsb->e0y),cimag(bsb->e0y));
  printf("x-component of translation vector        [m]: %15.14g\n",bsb->fx);
  printf("y-component of translation vector        [m]: %15.14g\n",bsb->fy);
  printf("z-component of translation vector        [m]: %15.14g\n",bsb->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",bsb->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",bsb->phi);
  printf("mode number for spiral phase modulation     : %15d\n",bsb->lg_m);

}

void setup_BsLGb(BsLGb *bsb)
{
  double complex phase_bslgb(int mode,double phi);
  
  double ce,h,apd;
  int nnc,i,j,s,st;
  
  bsb->data.st=sin(bsb->d_angle);
  bsb->data.ct=cos(bsb->d_angle);
  bsb->data.ki=2.0*M_PI*bsb->ni/bsb->lambda0;

  ce=1.0/sqrt(pow(cabs(bsb->e0x),2)+pow(cabs(bsb->e0y),2));
  bsb->data.ex= bsb->e0x*ce;
  bsb->data.ey= bsb->e0y*ce;
  bsb->data.hx=-bsb->data.ey*bsb->ni/Z0;
  bsb->data.hy= bsb->data.ex*bsb->ni/Z0;
  bsb->data.E0=sqrt(2.0*Z0*bsb->power/bsb->ni);
  bsb->data.cos_t=cos(bsb->theta);  bsb->data.sin_t=sin(bsb->theta);
  bsb->data.cos_p=cos(bsb->phi);    bsb->data.sin_p=sin(bsb->phi);
  // trapezoidal data
  nnc=(2<<(TRAP_MNUM))+1;
  bsb->data.cp=(double *)malloc(sizeof(double)*nnc);
  bsb->data.sp=(double *)malloc(sizeof(double)*nnc);
  bsb->data.ph_p=(double complex *)malloc(sizeof(double complex)*nnc);
  h=2.0*M_PI;
  bsb->data.cp[0]=-1.0;
  bsb->data.sp[0]= 0.0;
  bsb->data.ph_p[0]=phase_bslgb(bsb->lg_m,-M_PI);
  j=1;
  st=1;
  for(i=0;i<=TRAP_MNUM;i++){
    h*=0.5;
    for(s=1;s<=st;s++){
      apd=-M_PI+(double)(2*s-1)*h;
      bsb->data.cp[j]=cos(apd);
      bsb->data.sp[j]=sin(apd);
      bsb->data.ph_p[j]=phase_bslgb(bsb->lg_m,apd);
      j++;
    }
    st*=2;
  }
}

void free_BsLGb(BsLGb *bsb)
{
  free(bsb->data.cp);
  free(bsb->data.sp);
  free(bsb->data.ph_p);
}

///////////////////////////////////////////////////////////////////////
double complex phase_bslgb(int mode,double phi)
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

void calc_bslgb_eh(double complex *e,double complex *h,double *x,BsLGb *bsb)
{
  void int_phi_bslgb(double complex *eh,double ct,double st,double *x,BsLGb *bsb);
  
  double complex eh[6];
  double cos_t,sin_t;
  
  cos_t=bsb->data.ct;  sin_t=bsb->data.st;
  int_phi_bslgb(eh,cos_t,sin_t,x,bsb);
  e[0]=eh[0]*sin_t;  e[1]=eh[1]*sin_t;  e[2]=eh[2]*sin_t;
  h[0]=eh[3]*sin_t;  h[1]=eh[4]*sin_t;  h[2]=eh[5]*sin_t;
}

void int_phi_bslgb(double complex *eh,double ct,double st,double *x,BsLGb *bsb)
{
  void ac_eh_bslgb(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double complex phase,double *x,BsLGb *bsb);

  double h,abst,abss;
  double complex Ia[6],It[6],Is[6];
  int i,j,s,sc,cc;
  double f_order,cr_coef;
  double cep=bsb->ni*bsb->ni*epsilon0;
  double cmu=mu0;
  
  cc=0;
  h=2.0*M_PI;
  ac_eh_bslgb(Ia,ct,st,bsb->data.cp[0],bsb->data.sp[0],bsb->data.ph_p[0],x,bsb);
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
      ac_eh_bslgb(Ia,ct,st,bsb->data.cp[j],bsb->data.sp[j],bsb->data.ph_p[j],x,bsb);
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
    if((cc>1 && fabs(abst-abss) < EPS*cr_coef) || (abss==0.0 && abst==0.0) ){
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

void ac_eh_bslgb(double complex *eh,double cos_t,double sin_t,double cos_p,double sin_p,double complex phase,double *x,BsLGb *bsb)
{
  double tmp=bsb->data.ki*(sin_t*cos_p*x[0]+sin_t*sin_p*x[1]+cos_t*x[2]);
  double complex ce=phase*(cos(tmp)+I*sin(tmp));

  eh[0]=( bsb->data.ex*(sin_p*sin_p+cos_t*cos_p*cos_p)+bsb->data.ey*sin_p*cos_p*(cos_t-1.0))*ce;
  eh[1]=( bsb->data.ex*sin_p*cos_p*(cos_t-1.0)        +bsb->data.ey*(cos_t*sin_p*sin_p+cos_p*cos_p))*ce;
  eh[2]=(-bsb->data.ex*sin_t*cos_p                    -bsb->data.ey*sin_t*sin_p)*ce;
  
  eh[3]=( bsb->data.hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+bsb->data.hy*sin_p*cos_p*(cos_t-1.0))*ce;                                                                          
  eh[4]=( bsb->data.hx*sin_p*cos_p*(cos_t-1.0)        +bsb->data.hy*(sin_p*sin_p*cos_t+cos_p*cos_p))*ce;                                                                          
  eh[5]=(-bsb->data.hx*sin_t*cos_p                    -bsb->data.hy*sin_t*sin_p)*ce;                                              
}
