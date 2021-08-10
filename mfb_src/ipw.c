#include "ipw.h"

void calc_ipw_EH(double complex *e,double complex *h,double *x,Ipw *ipw)
{
  double xc=x[0]-ipw->fx;
  double yc=x[1]-ipw->fy;
  double zc=x[2]-ipw->fz;
  double te=ipw->data.ki*(xc*ipw->data.sin_t*ipw->data.cos_p+yc*ipw->data.sin_t*ipw->data.sin_p+zc*ipw->data.cos_t);
  double complex Ee=(cos(te)+I*sin(te));
  
  e[0]=ipw->data.ex*Ee*ipw->data.E0;  e[1]=ipw->data.ey*Ee*ipw->data.E0;  e[2]=ipw->data.ez*Ee*ipw->data.E0;
  h[0]=ipw->data.hx*Ee*ipw->data.E0;  h[1]=ipw->data.hy*Ee*ipw->data.E0;  h[2]=ipw->data.hz*Ee*ipw->data.E0;
}

void read_data_ipw(char *rfile,Ipw *ipw)
{
  FILE *fp;
  if((fp=fopen(rfile,"rt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  char buf[256]="";  int tmpi;  double tmpd,tmpd2;
  fgets(buf,256,fp);  fgets(buf,256,fp);

  printf("-- wave parameter --\n");
  fscanf(fp,"%d",&tmpi); 
  if(tmpi!=0){ printf("beam type %d is not supported\n",tmpi); exit(1);}
  else                                                printf("incident beam type                          : plane wave\n");
  fscanf(fp,"%lf",&tmpd);  ipw->lambda0=tmpd;         printf("wave length of incident beam in vacuum   [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  ipw->ni     =tmpd;         printf("refractive index of surrounding             : %15.14g\n",tmpd); 
  fscanf(fp,"%lf",&tmpd);  ipw->power  =tmpd;         printf("incident beam power                  [W/m^2]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd); 
  fscanf(fp,"%lf",&tmpd2); ipw->e0x    =tmpd+I*tmpd2; printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",tmpd,tmpd2);
  fscanf(fp,"%lf",&tmpd);
  fscanf(fp,"%lf",&tmpd2); ipw->e0y    =tmpd+I*tmpd2; printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",tmpd,tmpd2);
  fscanf(fp,"%lf",&tmpd);  ipw->fx     =tmpd;         printf("x-component of translation vector        [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  ipw->fy     =tmpd;         printf("y-component of translation vector        [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  ipw->fz     =tmpd;         printf("z-component of translation vector        [m]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  ipw->theta  =tmpd;         printf("rotation parameter theta               [rad]: %15.14g\n",tmpd);
  fscanf(fp,"%lf",&tmpd);  ipw->phi    =tmpd;         printf("rotation parameter phi                 [rad]: %15.14g\n",tmpd);

  fclose(fp);
}

void print_data_ipw(Ipw *ipw)
{
  printf("-- plane wave --\n");
  printf("wave length of incident beam in vacuum   [m]: %15.14g\n",ipw->lambda0);
  printf("refractive index of surrounding             : %15.14g\n",ipw->ni);
  printf("incident beam power                  [W/m^2]: %15.14g\n",ipw->power);
  printf("x-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(ipw->e0x),cimag(ipw->e0x));
  printf("y-component of polarization coefficient     :%7.6g+%7.6gi\n",creal(ipw->e0y),cimag(ipw->e0y));
  printf("x-component of translation vector        [m]: %15.14g\n",ipw->fx);
  printf("y-component of translation vector        [m]: %15.14g\n",ipw->fy);
  printf("z-component of translation vector        [m]: %15.14g\n",ipw->fz);
  printf("rotation parameter theta               [rad]: %15.14g\n",ipw->theta);
  printf("rotation parameter phi                 [rad]: %15.14g\n",ipw->phi);
}

void setup_ipw(Ipw *ipw)
{
  double cs=1.0/sqrt(pow(cabs(ipw->e0x),2)+pow(cabs(ipw->e0y),2));
  double complex ex,ey,hx,hy;
  double sin_t,cos_t,sin_p,cos_p;
  ex= ipw->e0x*cs;
  ey= ipw->e0y*cs;
  hx=-ey*ipw->ni/Z0;
  hy= ex*ipw->ni/Z0;

  ipw->data.E0=sqrt(ipw->power*2.0*Z0/ipw->ni);
  ipw->data.ki=2.0*M_PI*ipw->ni/ipw->lambda0;

  sin_t=sin(ipw->theta);  cos_t=cos(ipw->theta);
  sin_p=sin(ipw->phi  );  cos_p=cos(ipw->phi  );
  ipw->data.sin_t=sin_t;  ipw->data.cos_t=cos_t;
  ipw->data.sin_p=sin_p;  ipw->data.cos_p=cos_p;
  
  ipw->data.ex=ex*(sin_p*sin_p+cos_p*cos_p*cos_t)+ey*(sin_p*cos_p*cos_t-sin_p*cos_p);
  ipw->data.ey=ex*(sin_p*cos_p*cos_t-sin_p*cos_p)+ey*(sin_p*sin_p*cos_t+cos_p*cos_p);
  ipw->data.ez=-(ex*cos_p+ey*sin_p)*sin_t;
  ipw->data.hx=hx*(sin_p*sin_p+cos_p*cos_p*cos_t)+hy*(sin_p*cos_p*cos_t-sin_p*cos_p);
  ipw->data.hy=hx*(sin_p*cos_p*cos_t-sin_p*cos_p)+hy*(sin_p*sin_p*cos_t+cos_p*cos_p);
  ipw->data.hz=-(hx*cos_p+hy*sin_p)*sin_t;
}
