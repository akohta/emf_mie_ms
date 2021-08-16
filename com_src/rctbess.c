#include "rctbess.h"

int msta1(double x,int mp)
{
  double a0,f0,f1,f;
  int i,n0,n1,nn;

  a0 = fabs(x);
  n0 = (int)(1.1*a0)+1;
  f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-mp;
  n1 = n0+5;
  f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-mp;
  for (i=0;i<20;i++) {
    nn = n1-(n1-n0)/(1.0-f0/f1);
    f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-mp;
    if (abs(nn-n1) < 1) break;
    n0 = n1;
    f0 = f1;
    n1 = nn;
    f1 = f;
  }
  return nn;
}

int msta2(double x,int n,int mp)
{
  double a0,ejn,hmp,f0,f1,f,obj;
  int i,n0,n1,nn;

  a0 = fabs(x);
  hmp = 0.5*mp;
  ejn = 0.5*log10(6.28*n)-n*log10(1.36*a0/n);
  if (ejn <= hmp) {
    obj = mp;
    n0 = (int)(1.1*a0);
    if (n0 < 1) n0 = 1;
  }
  else {
    obj = hmp+ejn;
    n0 = n;
  }
  f0 = 0.5*log10(6.28*n0)-n0*log10(1.36*a0/n0)-obj;
  n1 = n0+5;
  f1 = 0.5*log10(6.28*n1)-n1*log10(1.36*a0/n1)-obj;
  for (i=0;i<20;i++) {
    nn = n1-(n1-n0)/(1.0-f0/f1);
    f = 0.5*log10(6.28*nn)-nn*log10(1.36*a0/nn)-obj;
    if (abs(nn-n1) < 1) break;
    n0 = n1;
    f0 = f1;
    n1 = nn;
    f1 = f;
  }
  return nn+10;
}

void rctjd(int n,double x,int *mn,double *rj,double *dj)
{
  double rj0,rj1,f0,f1,f,cs;
  int i,m;
  for(i=0;i<=n;i++){
    rj[i]=0.0;    dj[i]=0.0;
  }
  *mn=n;
  if(x<0.0){
    printf("rctjd(x) negative argument not supported\n");
    exit(1);
  }  
  else if(fabs(x)<1.0e-100){
    rj[0]=sin(x);
    dj[0]=cos(x);
    for(i=1;i<n;i++){
      rj[i]=0.0;      dj[i]=0.0;
    }
    *mn=0;
  }
  else {
    rj[0]=sin(x);
    rj[1]=rj[0]/x-cos(x);
    rj0=rj[0];
    rj1=rj[1];
    if(n>=2){
      m=msta1(x,300);
      if(m<n)      *mn=m;
      else m=msta2(x,n,15);
      f=0.0;
      f0=0.0;
      f1=1.0e-100;
      for(i=m;i>=0;i--){
        f=(double)(2*i+3)*f1/x-f0;
        if(i<=*mn) rj[i]=f;
        f0=f1;
        f1=f;
      }
      if(fabs(rj0) > fabs(rj1)) cs=rj0/f;
      else cs=rj1/f0;
      for(i=0;i<=*mn;i++)       rj[i]=cs*rj[i];
    }
    dj[0]=cos(x);
    for(i=1;i<=*mn;i++)      dj[i]=-(double)i*rj[i]/x+rj[i-1];
  }
}

void rctyd(int n,double x,int *mn,double *ry,double *dy)
{
  double rf0,rf1,rf2;
  int i;
  for(i=0;i<=n;i++){
    ry[i]=0.0;    dy[i]=0.0;
  }
  *mn=n;
  if(x<0.0){ 
    printf("rctyd(x) negative argument of x not supported\n");
    exit(1);
  }
  else if(x<1.0e-60){
    ry[0]= cos(x);
    dy[0]= sin(x);
    for(i=1;i<=n;i++){
      ry[i]= 1.0e300;
      dy[i]= 1.0e300;
    }
    *mn=0;
  }
  else {
    ry[0]=cos(x);
    ry[1]=ry[0]/x+sin(x);
    rf0=ry[0];
    rf1=ry[1];
    for(i=2;i<=n;i++){
      rf2=(double)(2*i-1)*rf1/x-rf0;
      if(fabs(rf2)>1.0e300) {
	*mn=i-1;
	break;
      }
      ry[i]=rf2;
      rf0=rf1;
      rf1=rf2;
    }
    dy[0]=sin(x);
    for(i=1;i<=*mn;i++){
      dy[i]=-(double)i*ry[i]/x+ry[i-1];
    }
  }
}

void rctjc(int n,double _Complex z,int *mn,double _Complex *rj,double _Complex *dj)
{
  double _Complex rj0,rj1,f0,f1,f,cs;
  double a0;
  int i,m;
  for(i=0;i<=n;i++){
    rj[i]=0.0;    dj[i]=0.0;
  }
  a0=cabs(z);
  *mn=n;
  rj[0]=csin(z);
  rj[1]=rj[0]/z-ccos(z);
  rj0=rj[0];
  rj1=rj[1];
  if(n>=2){
    m=msta1(a0,300);
    if(m<n)      *mn=m;
    else m=msta2(a0,n,15);
    f=0.0;
    f0=0.0;
    f1=1.0e-100;
    for(i=m;i>=0;i--){
      f=(double)(2*i+3)*f1/z-f0;
      if(i<=*mn) rj[i]=f;
      f0=f1;
      f1=f;
    }
    if(cabs(rj0) > cabs(rj1)) cs=rj0/f;
    else cs=rj1/f0;
    for(i=0;i<=*mn;i++)       rj[i]=cs*rj[i];
  }
  dj[0]=ccos(z);
  for(i=1;i<=*mn;i++)      dj[i]=-(double)i*rj[i]/z+rj[i-1];
}

void rctyc(int n,double _Complex z,int *mn,double _Complex *ry,double _Complex *dy)
{
  double _Complex rf0,rf1,rf2;
  int i;
  for(i=0;i<=n;i++){
    ry[i]=0.0;    dy[i]=0.0;
  }
  *mn=n;
  ry[0]=ccos(z);
  ry[1]=ry[0]/z+csin(z);
  rf0=ry[0];
  rf1=ry[1];
  for(i=2;i<=n;i++){
    rf2=(double)(2*i-1)*rf1/z-rf0;
    if(cabs(rf2)>1.0e300) {
      *mn=i-1;
      break;
    }
    ry[i]=rf2;
    rf0=rf1;
    rf1=rf2;
  }
  dy[0]=csin(z);
  for(i=1;i<=*mn;i++){
    dy[i]=-(double)i*ry[i]/z+ry[i-1];
  }
}

void rcth1d(int n,double x,int *mn,double _Complex *ch,double _Complex *dch)
{
  double *psi =(double *)malloc(sizeof(double)*(n+1)); if( psi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1d(),psi. Exit...\n"); exit(1);}
  double *dpsi=(double *)malloc(sizeof(double)*(n+1)); if(dpsi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1d(),dpsi. Exit...\n"); exit(1);}
  double *chi =(double *)malloc(sizeof(double)*(n+1)); if( chi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1d(),chi. Exit...\n"); exit(1);}
  double *dchi=(double *)malloc(sizeof(double)*(n+1)); if(dchi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1d(),dchi. Exit...\n"); exit(1);}
  int nj,ny,i;
  for(i=0;i<=n;i++){
    ch[i]=0.0;    dch[i]=0.0;
  }
  *mn=n;
  rctjd( n,x,&nj,psi,dpsi);  if(nj<n ) *mn=nj;
  rctyd(nj,x,&ny,chi,dchi);  if(ny<nj) *mn=ny;
  for(i=0;i<=*mn;i++){
    ch[i] = psi[i]-I* chi[i];
    dch[i]=dpsi[i]-I*dchi[i];
  }
  free( psi);  free( chi);
  free(dpsi);  free(dchi);
}

void rcth2d(int n,double x,int *mn,double _Complex *ch,double _Complex *dch)
{
  double *psi =(double *)malloc(sizeof(double)*(n+1)); if( psi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2d(),psi. Exit...\n"); exit(1);}
  double *dpsi=(double *)malloc(sizeof(double)*(n+1)); if(dpsi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2d(),dpsi. Exit...\n"); exit(1);}
  double *chi =(double *)malloc(sizeof(double)*(n+1)); if( chi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2d(),chi. Exit...\n"); exit(1);}
  double *dchi=(double *)malloc(sizeof(double)*(n+1)); if(dchi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2d(),dchi. Exit...\n"); exit(1);}
  int nj,ny,i;
  for(i=0;i<=n;i++){
    ch[i]=0.0;    dch[i]=0.0;
  }
  *mn=n;
  rctjd( n,x,&nj,psi,dpsi);  if(nj<n ) *mn=nj;
  rctyd(nj,x,&ny,chi,dchi);  if(ny<nj) *mn=ny;
  for(i=0;i<=*mn;i++){
    ch[i] = psi[i]+I* chi[i];
    dch[i]=dpsi[i]+I*dchi[i];
  }
  free( psi);  free( chi);
  free(dpsi);  free(dchi);
}

void rcth1c(int n,double _Complex z,int *mn,double _Complex *ch,double _Complex *dch)
{
  double _Complex *psi =(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if( psi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1c(),psi. Exit...\n"); exit(1);}
  double _Complex *dpsi=(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if(dpsi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1c(),dpsi. Exit...\n"); exit(1);}
  double _Complex *chi =(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if( chi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1c(),chi. Exit...\n"); exit(1);}
  double _Complex *dchi=(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if(dchi==NULL){ printf("failed to allocate memory. rctbess.c,rcth1c(),dchi. Exit...\n"); exit(1);}
  int nj,ny,i;
  for(i=0;i<=n;i++){
    ch[i]=0.0;    dch[i]=0.0;
  }
  *mn=n;
  rctjc( n,z,&nj,psi,dpsi);  if(nj<n ) *mn=nj;
  rctyc(nj,z,&ny,chi,dchi);  if(ny<nj) *mn=ny;
  for(i=0;i<=*mn;i++){
    ch[i] = psi[i]-I* chi[i];
    dch[i]=dpsi[i]-I*dchi[i];
  }
  free( psi);  free( chi);
  free(dpsi);  free(dchi);
}

void rcth2c(int n,double _Complex z,int *mn,double _Complex *ch,double _Complex *dch)
{
  double _Complex *psi =(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if( psi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2c(),psi. Exit...\n"); exit(1);}
  double _Complex *dpsi=(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if(dpsi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2c(),dpsi. Exit...\n"); exit(1);}
  double _Complex *chi =(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if( chi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2c(),chi. Exit...\n"); exit(1);}
  double _Complex *dchi=(double _Complex *)malloc(sizeof(double _Complex)*(n+1)); if(dchi==NULL){ printf("failed to allocate memory. rctbess.c,rcth2c(),dchi. Exit...\n"); exit(1);}
  int nj,ny,i;
  for(i=0;i<=n;i++){
    ch[i]=0.0;    dch[i]=0.0;
  }
  *mn=n;
  rctjc( n,z,&nj,psi,dpsi);  if(nj<n ) *mn=nj;
  rctyc(nj,z,&ny,chi,dchi);  if(ny<nj) *mn=ny;
  for(i=0;i<=*mn;i++){
    ch[i] = psi[i]+I* chi[i];
    dch[i]=dpsi[i]+I*dchi[i];
  }
  free( psi);  free( chi);
  free(dpsi);  free(dchi);
}
