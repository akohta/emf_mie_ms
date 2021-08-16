#include "emf_mie_ms.h"

void  read_data_ms(MSPD *msp)
{
  FILE *fp;
  if((fp=fopen(fn_sphr,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",fn_sphr);    exit(1);  }
  char buf[256]="";  int tmpi;  double tmpd,tmpd2;
  fgets(buf,256,fp);  fgets(buf,256,fp);

  int num,nc;
  fscanf(fp,"%d\n",&tmpi);   num=tmpi;
  fgets(buf,256,fp);
  if(num==0){    printf("No sphere defined. Exit...\n"); exit(1);  }
  
  msp->n_sphr=num;
  msp->sp=(SPD *)m_alloc2(num,sizeof(SPD),"read_data_ms(),msp->sp");
  
  for(nc=0;nc<num;nc++){
    fscanf(fp,"%lf",&tmpd);  msp->sp[nc].a      =tmpd;
    fscanf(fp,"%lf",&tmpd); 
    fscanf(fp,"%lf",&tmpd2); msp->sp[nc].ns     =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);  msp->sp[nc].xs     =tmpd;
    fscanf(fp,"%lf",&tmpd);  msp->sp[nc].ys     =tmpd;
    fscanf(fp,"%lf",&tmpd);  msp->sp[nc].zs     =tmpd; 
    fscanf(fp,"%d",&tmpi);   msp->sp[nc].bsn    =tmpi; 
    fscanf(fp,"%d",&tmpi);   msp->sp[nc].bdv    =tmpi; 
    fscanf(fp,"%d",&tmpi);   msp->sp[nc].l_limit=tmpi; 
  }
  fclose(fp);
  
  // multi fbeam
  init_mfb(&(msp->bm));       // initialize
  read_data_mfb(&(msp->bm));  // search and read beam datafile
}

void print_data_ms(MSPD *msp)
{
  int nc;
  
  print_data_mfb(&(msp->bm)); // print beam data

  // print sphere data  
  printf("---- sphere data ( %s ) ----\n",fn_sphr);
  printf("number of spheres                             : %16d\n",msp->n_sphr);
  for(nc=0;nc<msp->n_sphr;nc++){
    printf("  Sphere ID %d\n",nc);
    printf("radius of sphere                              : %16.15g\n",msp->sp[nc].a);
    printf("refractive index of sphere                    : %7.6g+%7.6gi\n",creal(msp->sp[nc].ns),cimag(msp->sp[nc].ns));
    printf("x-coordinate of sphere center                 : %16.15g\n",msp->sp[nc].xs);
    printf("y-coordinate of sphere center                 : %16.15g\n",msp->sp[nc].ys);
    printf("z-coordinate of sphere center                 : %16.15g\n",msp->sp[nc].zs);
    printf("basic sampling number on sphere surface       : %16d\n",msp->sp[nc].bsn);
    printf("division number for sphere surface    (per PI): %16d\n",msp->sp[nc].bdv);
    printf("limit of order number l                       : %16d\n",msp->sp[nc].l_limit);
  }
  printf("\n");
  
  //printf("continue? (y/n) : ");  if(getchar()!='y'){ printf("Exit\n");  exit(0);}
}

void print_data_ms_mksa(MSPD *msp)
{
  int nc;
  
  print_data_mfb_mksa(&(msp->bm)); // print beam data

  // print sphere data  
  printf("---- sphere data ( %s ), MKSA system ----\n",fn_sphr);
  printf("number of spheres                             : %16d\n",msp->n_sphr);
  for(nc=0;nc<msp->n_sphr;nc++){
    printf("  Sphere ID %d\n",nc);
    printf("radius of sphere                           [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[nc].a));
    printf("refractive index of sphere                    : %7.6g+%7.6gi\n",creal(msp->sp[nc].ns),cimag(msp->sp[nc].ns));
    printf("x-coordinate of sphere center              [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[nc].xs));
    printf("y-coordinate of sphere center              [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[nc].ys));
    printf("z-coordinate of sphere center              [m]: %16.15g\n",OSUtoMKSA_length(msp->sp[nc].zs));
    printf("basic sampling number on sphere surface       : %16d\n",msp->sp[nc].bsn);
    printf("division number for sphere surface    (per PI): %16d\n",msp->sp[nc].bdv);
    printf("limit of order number l                       : %16d\n",msp->sp[nc].l_limit);
  }
  printf("\n");
  
  //printf("continue? (y/n) : ");  if(getchar()!='y'){ printf("Exit\n");  exit(0);}
}

void setup_ms(MSPD *msp)
{
  void check_position(MSPD *msp);
  void setup_sp(SPD *sp);
  void setup_coefficient_dispd(SPD *sp,Bobj *bm);
  void initialize_eh_r(SPD *sp,Bobj *bm);
  void coefficient(SPD *sp);
  void check_l_limit_ms(MSPD *msp);
  
  int i;
  
  check_position(msp);
  
  // multi_fbeam
  setup_mfb(&(msp->bm));
  
  // spheres
  for(i=0;i<msp->n_sphr;i++){
    setup_sp(&(msp->sp[i]));
	setup_coefficient_dispd(&(msp->sp[i]),&(msp->bm));
	initialize_eh_r(&(msp->sp[i]),&(msp->bm)); 
	coefficient(&(msp->sp[i]));
  }
  
  check_l_limit_ms(msp);

}

void  free_ms(MSPD *msp)
{
  void free_sp(SPD *sp);

  int i;
  
  // spheres
  for(i=0;i<msp->n_sphr;i++){
    free_sp(&(msp->sp[i]));
  }
  free(msp->sp);  msp->n_sphr=0;
  
  // multi_fbeam
  free_mfb(&(msp->bm));
}

void iterative_ops_ms(MSPD *msp)
{
  void field_s_ehr(int src,int obj,MSPD *msp); 
  void coefficient(SPD *sp); 
  
  int i,j,t,nn,sbc,num=msp->n_sphr;
  double vf[3],f1,*f0;
  int *bc;
  
  if(num<2) return;
  
  bc=(int *)m_alloc2(num,sizeof(int),"iterative_ops_ms(),*mc");
  f0=(double *)m_alloc2(num,sizeof(double),"iterative_ops_ms(),*f0");

  for(t=0;t<num;t++){
    force_ms(t,vf,msp);
    f0[t]=vf[0]*vf[0]+vf[1]*vf[1]+vf[2]*vf[2];
    bc[t]=ito_breakcount;
  }

  printf("iterative operation start (convergence criterion : cv < %g)\n",ito_eps);
  for(nn=0;nn<ito_max;nn++){
    for(i=0;i<num;i++)
      for(j=0;j<num;j++) if(i!=j) field_s_ehr(i,j,msp);
    for(i=0;i<num;i++) coefficient(&(msp->sp[i]));
    printf("%3d, cv : ",nn);
    for(t=0;t<num;t++){
      force_ms(t,vf,msp);
      f1=vf[0]*vf[0]+vf[1]*vf[1]+vf[2]*vf[2];
      if(fabs(f1/f0[t]-1.0)<ito_eps)  bc[t]--;
      printf("%g\t",fabs(f1/f0[t]-1.0));
      f0[t]=f1;
    }
    printf("\n");
    sbc=0;
    for(t=0;t<num;t++) if(bc[t]<=0) sbc++;
    if(sbc==num) break;
  }

  free(bc);  free(f0);
}

///////////////////////////////////////////////////////////////////////
void check_position(MSPD *msp)
{
  double r,rs;
  int i,j;
  
  for(i=0;i<msp->n_sphr;i++){
    for(j=i+1;j<msp->n_sphr;j++){
      r =msp->sp[i].a+msp->sp[j].a;
      rs=sqrt(pow(msp->sp[j].xs-msp->sp[i].xs,2)+pow(msp->sp[j].ys-msp->sp[i].ys,2)+pow(msp->sp[j].zs-msp->sp[i].zs,2));
      if(rs<r){
        printf("Sphere Position Check Error! sphere_id[%d] and sphere_id[%d] is overlaped. Exit...\n",i,j);
        exit(1);
      }
    }
  }
}

void setup_sp(SPD *sp)
{
  void gauleg_dv(double a,double b,double *x,double *w,int bn,int dv);

  int mn;
  int nt=  sp->bsn*sp->bdv;
  int np=2*sp->bsn*sp->bdv;
  
  sp->ddt.nt=nt;
  sp->ddt.np=np;
  sp->ddt.xp=(double *)m_alloc2(np,sizeof(double),"setup_sp(),sp->ddt.xp");
  sp->ddt.wp=(double *)m_alloc2(np,sizeof(double),"setup_sp(),sp->ddt.wp");
  sp->ddt.xt=(double *)m_alloc2(nt,sizeof(double),"setup_sp(),sp->ddt.xt");
  sp->ddt.wt=(double *)m_alloc2(nt,sizeof(double),"setup_sp(),sp->ddt.wt");
  gauleg_dv(0.0,    M_PI,sp->ddt.xt,sp->ddt.wt,sp->bsn,  sp->bdv);
  gauleg_dv(0.0,2.0*M_PI,sp->ddt.xp,sp->ddt.wp,sp->bsn,2*sp->bdv);
  
  sp->ddt.eri=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.eri");
  sp->ddt.hri=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.hri");
  sp->ddt.ers=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.ers");
  sp->ddt.hrs=(double complex *)m_alloc2(np*nt,sizeof(double complex),"setup_sp(),sp->ddt.hrs");
  
  mn=sp->l_limit;
  sp->ddt.cab=(double *)m_alloc2(mn+1,sizeof(double),"setup_sp(),sp->ddt.cab");
  sp->ddt.ca=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.ca");
  sp->ddt.cb=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.cb");
  sp->ddt.cc=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.cc");
  sp->ddt.cd=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_sp(),sp->ddt.cd");

  sp->ddt.Alm=(double complex *)m_alloc2(mn*(mn+2),sizeof(double complex),"setup_sp(),sp->ddt.Alm");
  sp->ddt.Blm=(double complex *)m_alloc2(mn*(mn+2),sizeof(double complex),"setup_sp(),sp->ddt.Blm");
}

void gauleg_dv(double a,double b,double *x,double *w,int bn,int dv)
{
  double xt[bn],wt[bn];
  gauleg(-1.0, 1.0,xt,wt,bn);
  
  double h,dh,x0,x1,cx,cc;
  int d,i,j;
  h=b-a;
  dh=h/(double)dv;
  x1=a;
  j=0;
  for(d=0;d<dv;d++){
    x0=x1;
    x1=x0+dh;
    
    cx=0.5*(x1-x0);
    cc=0.5*(x1+x0);
    for(i=0;i<bn;i++){
      x[j]= cx*xt[i]+cc;
      w[j]= cx*wt[i];
      j++;
    }
  }
}

void free_sp(SPD *sp)
{
  free(sp->ddt.xp);  free(sp->ddt.wp);
  free(sp->ddt.xt);  free(sp->ddt.wt);
  
  free(sp->ddt.eri);  free(sp->ddt.hri);
  free(sp->ddt.ers);  free(sp->ddt.hrs);

  free(sp->ddt.cab);
  free(sp->ddt.ca );  free(sp->ddt.cb );
  free(sp->ddt.cc );  free(sp->ddt.cd );

  free(sp->ddt.Alm);  free(sp->ddt.Blm);
}

void setup_coefficient_dispd(SPD *sp,Bobj *bm)
{
  double complex z,nr,b_ac,b_bd,b_t;
  double x,a2;
  int mn,nn,i;
  double complex *xi,*dxi;
  double complex *psic,*dpsic;
  
  mn=sp->l_limit;
  x=2.0*M_PI*bm->n_0*sp->a/bm->lambda_0;
  xi =(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_coefficient_dispd(),xi");
  dxi=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_coefficient_dispd(),dxi");
  rcth1d(mn,x,&nn,xi,dxi);
  sp->ddt.l_max=nn;
  z=2.0*M_PI*sp->ns*sp->a/bm->lambda_0;
  psic =(double complex *)m_alloc2(mn+1,sizeof(double complex),"setup_coefficient_dispd(),psic");
  dpsic=(double complex *)m_alloc2(mn+1,sizeof(double complex),"setpp_coefficient?dispd(),dpsic");
  rctjc(mn,z,&nn,psic,dpsic);
  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  
  a2=pow(sp->a,2);
  for(i=1;i<=mn;i++){
    sp->ddt.cab[i]=a2/((double)(i*(i+1))*creal(xi[i]));
  }
  
  nr=sp->ns/bm->n_0;
  for(i=1;i<=mn;i++){
    b_ac=1.0/(nr*psic[i]*dxi[i]-dpsic[i]*xi[i]);
    b_bd=1.0/(psic[i]*dxi[i]-nr*dpsic[i]*xi[i]);
    b_t =dxi[i]*creal(xi[i])-xi[i]*creal(dxi[i]);
    sp->ddt.ca[i]=(dpsic[i]*creal(xi[i])-nr*psic[i]*creal(dxi[i]))*b_ac;
    sp->ddt.cb[i]=(nr*dpsic[i]*creal(xi[i])-psic[i]*creal(dxi[i]))*b_bd;
    sp->ddt.cc[i]=b_t*b_ac/nr;
    sp->ddt.cd[i]=b_t*b_bd;
  }
  free(xi);   free(dxi);
  free(psic);  free(dpsic);
}

void initialize_eh_r(SPD *sp,Bobj *bm)
{
  double complex e[3],h[3];
  double r,theta,phi,x[3],sin_t,cos_t,sin_p,cos_p;
  int i,j,nt,np;
  
  nt=sp->ddt.nt;
  np=sp->ddt.np;
  r=sp->a;
  for(i=0;i<nt;i++){
    theta=sp->ddt.xt[i];
    sin_t=sin(theta);    cos_t=cos(theta);
    #pragma omp parallel for schedule(dynamic) private(phi,sin_p,cos_p,x,e,h) // OpenMP parallel for
    for(j=0;j<np;j++){
      phi=sp->ddt.xp[j];
      sin_p=sin(phi);      cos_p=cos(phi);
      x[0]=r*sin_t*cos_p+sp->xs;
      x[1]=r*sin_t*sin_p+sp->ys;
      x[2]=r*cos_t      +sp->zs;
      calc_mfb_EH(e,h,x,bm);
      sp->ddt.eri[i*np+j]=e[0]*sin_t*cos_p+e[1]*sin_t*sin_p+e[2]*cos_t;
      sp->ddt.hri[i*np+j]=h[0]*sin_t*cos_p+h[1]*sin_t*sin_p+h[2]*cos_t;
    }
  }
}


void coefficient(SPD *sp)
{
  int ti,lm,np,nt,tt,l,m,i,j,t;
  size_t ms,lmax;
  
  lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  np=sp->ddt.np;
  nt=sp->ddt.nt;

  for(ti=0;ti<lm*(lm+2);ti++){
    sp->ddt.Alm[ti]=0.0;     sp->ddt.Blm[ti]=0.0;
  }
  
  #pragma omp parallel private(i,j,l,m,t,tt,ti) // pragma parallel 
  {
  double complex Yp,Ym;
  double theta,sin_t,cos_t,flag;
  
  double *sphPlm=(double *)m_alloc2(ms,sizeof(double),"coefficient(),*sphPlm");
  double complex *e_phim=(double complex *)m_alloc2(np,sizeof(double complex),"coefficient(),*e_phim");
  double complex *tmpAlm=(double complex *)m_alloc2(lm*(lm+2),sizeof(double complex),"coefficient(),*tmpAlm");
  double complex *tmpBlm=(double complex *)m_alloc2(lm*(lm+2),sizeof(double complex),"coefficient(),*tmpBlm");
  
  #pragma omp for schedule(dynamic) // parallel for
  for(i=0;i<nt;i++){
    theta=sp->ddt.xt[i];
    sin_t=sin(theta);    cos_t=cos(theta);
    tt=0;
    m=0;
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm);
    for(l=1;l<=lm;l++){
      Yp=sphPlm[gsl_sf_legendre_array_index(l,m)];
      for(j=0;j<np;j++){
        tmpAlm[tt]+=(sp->ddt.eri[i*np+j]+sp->ddt.ers[i*np+j])*conj(Yp)*sp->ddt.wp[j];
        tmpBlm[tt]+=(sp->ddt.hri[i*np+j]+sp->ddt.hrs[i*np+j])*conj(Yp)*sp->ddt.wp[j];
      }
      tt++;
    }
    for(m=1;m<=lm;m++){
      for(t=0;t<np;t++) e_phim[t]=cos((double)m*sp->ddt.xp[t])+I*sin((double)m*sp->ddt.xp[t]);
      if(m%2==0) flag=1.0;
      else flag=-1.0;
      for(l=m;l<=lm;l++){
        for(j=0;j<np;j++){
          Yp=sphPlm[gsl_sf_legendre_array_index(l,m)]*e_phim[j];
          Ym=flag*conj(Yp);
          tmpAlm[tt  ]+=(sp->ddt.eri[i*np+j]+sp->ddt.ers[i*np+j])*conj(Yp)*sp->ddt.wp[j];
          tmpBlm[tt  ]+=(sp->ddt.hri[i*np+j]+sp->ddt.hrs[i*np+j])*conj(Yp)*sp->ddt.wp[j];
          tmpAlm[tt+1]+=(sp->ddt.eri[i*np+j]+sp->ddt.ers[i*np+j])*conj(Ym)*sp->ddt.wp[j];
          tmpBlm[tt+1]+=(sp->ddt.hri[i*np+j]+sp->ddt.hrs[i*np+j])*conj(Ym)*sp->ddt.wp[j];
        }
        tt+=2;
      }
    }
    for(ti=0;ti<lm*(lm+2);ti++){
      #pragma omp critical // pragma omp critical
      sp->ddt.Alm[ti]+=tmpAlm[ti]*sin_t*sp->ddt.wt[i];
      #pragma omp critical // pragma omp critical
      sp->ddt.Blm[ti]+=tmpBlm[ti]*sin_t*sp->ddt.wt[i];
      tmpAlm[ti]=0.0;      tmpBlm[ti]=0.0;
    }
  }
  free(sphPlm);
  free(e_phim);
  free(tmpAlm);  free(tmpBlm);
  } // pragma parallel end
  
  tt=0;  m=0;
  for(l=1;l<=lm;l++){
    sp->ddt.Alm[tt]*=sp->ddt.cab[l]; 
    sp->ddt.Blm[tt]*=sp->ddt.cab[l];
    tt++;
  }
  for(m=1;m<=lm;m++){
    for(l=m;l<=lm;l++){
      sp->ddt.Alm[tt  ]*=sp->ddt.cab[l];      sp->ddt.Alm[tt+1]*=sp->ddt.cab[l];
      sp->ddt.Blm[tt  ]*=sp->ddt.cab[l];      sp->ddt.Blm[tt+1]*=sp->ddt.cab[l];
      tt+=2;
    }
  }
  for(i=0;i<nt;i++){
    for(j=0;j<np;j++){
      sp->ddt.ers[i*np+j]=0.0;
      sp->ddt.hrs[i*np+j]=0.0;
    }
  }
 
}

void check_l_limit_ms(MSPD *msp)
{
  int i;
  for(i=0;i<msp->n_sphr;i++){
    if(msp->sp[i].l_limit>msp->sp[i].ddt.l_max){
      printf("Overflow and underflow problem of Riccati-Bessel function occurred sphere id %d. Check the data precision.\n",i);
      printf("Available order number is less than %d.\n",msp->sp[i].ddt.l_max); 
    }   
  }
}

void field_s_ehr(int src,int obj,MSPD *msp)
{
  void scattered_EH(double complex *e,double complex *h,double *xb,SPD *sp,Bobj *bm);
  
  double complex es[3],hs[3];
  double cos_t,sin_t,cos_p,sin_p;
  double r[3];
  double a=msp->sp[obj].a;
  int np=msp->sp[obj].ddt.np;
  int nt=msp->sp[obj].ddt.nt;
  int i,j;

  #pragma omp parallel for schedule(dynamic) private(j,cos_t,sin_t,cos_p,sin_p,r,es,hs)
  for(i=0;i<nt;i++){
    cos_t=cos(msp->sp[obj].ddt.xt[i]);    sin_t=sin(msp->sp[obj].ddt.xt[i]);
    for(j=0;j<np;j++){
      cos_p=cos(msp->sp[obj].ddt.xp[j]);      sin_p=sin(msp->sp[obj].ddt.xp[j]);
      r[0]=a*sin_t*cos_p+msp->sp[obj].xs;
      r[1]=a*sin_t*sin_p+msp->sp[obj].ys;
      r[2]=a*cos_t      +msp->sp[obj].zs;
      scattered_EH(es,hs,r,&(msp->sp[src]),&(msp->bm));
      msp->sp[obj].ddt.ers[i*np+j]+=es[0]*sin_t*cos_p+es[1]*sin_t*sin_p+es[2]*cos_t;
      msp->sp[obj].ddt.hrs[i*np+j]+=hs[0]*sin_t*cos_p+hs[1]*sin_t*sin_p+hs[2]*cos_t;
    }
  }
}

void scattered_EH(double complex *e,double complex *h,double *xb,SPD *sp,Bobj *bm)
{
  double complex er,et,ep,hr,ht,hp,Yp,Ym,dYp,dYm,dep,expi;
  double r,rxy,r2,cos_t,sin_t,cos_p,sin_p,ker,ke,flag,i_ne,ne,x,y,z,i_sin_t;
  int l,m,tt,nn,lm,ai;
  size_t ms,lmax=(size_t)sp->l_limit;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"scattered_EH(),*sphPlm");
  double *dsphPlm=(double *)m_alloc2(ms,sizeof(double),"scattered_EH(),*dsphPlm");
  double complex *xi =(double complex *)m_alloc2(lm+1,sizeof(double complex),"scattered_EH(),*xi");
  double complex *dxi=(double complex *)m_alloc2(lm+1,sizeof(double complex),"scattered_EH(),*dxi");

  x=xb[0]-sp->xs;  y=xb[1]-sp->ys;  z=xb[2]-sp->zs;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  rxy=sqrt(x*x+y*y);
  if(rxy==0.0){ // x==0,y==0
    x=z*0.7e-7;
    y=z*0.7e-7;
    r2=x*x+y*y+z*z;    r=sqrt(r2);
    rxy=sqrt(x*x+y*y); 
  }
  cos_t=z/r;    sin_t=rxy/r;   i_sin_t=r/rxy;
  cos_p=x/rxy;  sin_p=y/rxy;
  ke =2.0*M_PI*bm->n_0/bm->lambda_0;
  ker=ke*r;
  ne=bm->n_0;
  i_ne=1.0/(bm->n_0);
  
  rcth1d(lm,ker,&nn,xi,dxi);
  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  dep=cos_p+I*sin_p; expi=1.0;
  er=0.0;  et=0.0;  ep=0.0;
  hr=0.0;  ht=0.0;  hp=0.0;
  tt=0;  m=0;
  gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm,dsphPlm);
  for(l=1;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    dYp=dsphPlm[ai];
    er+=(double)(l*(l+1))*sp->ddt.ca[l]*sp->ddt.Alm[tt]*xi[l]*Yp;
    hr+=(double)(l*(l+1))*sp->ddt.cb[l]*sp->ddt.Blm[tt]*xi[l]*Yp;
    et+=sp->ddt.ca[l]*sp->ddt.Alm[tt]*dxi[l]*dYp-(double)m*i_ne*sp->ddt.cb[l]*sp->ddt.Blm[tt]*xi[l]*Yp*i_sin_t;
    ht+=sp->ddt.cb[l]*sp->ddt.Blm[tt]*dxi[l]*dYp+(double)m*  ne*sp->ddt.ca[l]*sp->ddt.Alm[tt]*xi[l]*Yp*i_sin_t;
    ep+=(double)m*sp->ddt.ca[l]*sp->ddt.Alm[tt]*dxi[l]*Yp*i_sin_t-i_ne*sp->ddt.cb[l]*sp->ddt.Blm[tt]*xi[l]*dYp;
    hp+=(double)m*sp->ddt.cb[l]*sp->ddt.Blm[tt]*dxi[l]*Yp*i_sin_t+  ne*sp->ddt.ca[l]*sp->ddt.Alm[tt]*xi[l]*dYp;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    if(m%2==0)flag= 1.0;
    else      flag=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      dYp=dsphPlm[ai]*expi;
      Ym =flag*conj( Yp);
      dYm=flag*conj(dYp);
      er+=(double)(l*(l+1))*sp->ddt.ca[l]*xi[l]*(sp->ddt.Alm[tt  ]*Yp+sp->ddt.Alm[tt+1]*Ym);
      et+=dxi[l]*sp->ddt.ca[l]*(sp->ddt.Alm[tt  ]*dYp+sp->ddt.Alm[tt+1]*dYm)
        -(double)m*i_ne*sp->ddt.cb[l]*xi[l]*i_sin_t*(sp->ddt.Blm[tt  ]*Yp-sp->ddt.Blm[tt+1]*Ym);
      ep+=(double)m*sp->ddt.ca[l]*dxi[l]*i_sin_t*(sp->ddt.Alm[tt  ]*Yp-sp->ddt.Alm[tt+1]*Ym)
        -i_ne*sp->ddt.cb[l]*xi[l]*(sp->ddt.Blm[tt  ]*dYp+sp->ddt.Blm[tt+1]*dYm);
      hr+=(double)(l*(l+1))*sp->ddt.cb[l]*xi[l]*(sp->ddt.Blm[tt  ]*Yp+sp->ddt.Blm[tt+1]*Ym); 
      ht+=sp->ddt.cb[l]*dxi[l]*(sp->ddt.Blm[tt  ]*dYp+sp->ddt.Blm[tt+1]*dYm)
        +(double)m*ne*sp->ddt.ca[l]*xi[l]*i_sin_t*(sp->ddt.Alm[tt  ]*Yp-sp->ddt.Alm[tt+1]*Ym);
      hp+=(double)m*sp->ddt.cb[l]*dxi[l]*i_sin_t*(sp->ddt.Blm[tt  ]*Yp-sp->ddt.Blm[tt+1]*Ym)
        +ne*sp->ddt.ca[l]*xi[l]*(sp->ddt.Alm[tt  ]*dYp+sp->ddt.Alm[tt+1]*dYm);
      tt+=2;
    }
  }
  er/=r2;
  et*=ke/r;
  ep*=I*ke/r;
  hr/=r2;
  ht*=ke/r;
  hp*=I*ke/r;

  e[0]=er*sin_t*cos_p+et*cos_t*cos_p-ep*sin_p;
  e[1]=er*sin_t*sin_p+et*cos_t*sin_p+ep*cos_p;
  e[2]=er*cos_t-et*sin_t;
  h[0]=hr*sin_t*cos_p+ht*cos_t*cos_p-hp*sin_p;
  h[1]=hr*sin_t*sin_p+ht*cos_t*sin_p+hp*cos_p;
  h[2]=hr*cos_t-ht*sin_t;

  free(sphPlm);  free(dsphPlm);
  free(xi);      free(dxi);
}

void internal_EH(double complex *e,double complex *h,double *xb,SPD *sp,Bobj *bm)
{
  void internal_EH_r0(double complex *e,double complex *h,SPD *sp,Bobj *bm);
  
  double complex er,et,ep,hr,ht,hp,Yp,Ym,dYp,dYm,dep,expi,ke,ker,i_ne,ne;
  double r,rxy,r2,cos_t,sin_t,cos_p,sin_p,flag,x,y,z,i_sin_t;
  int l,m,tt,nn,lm,ai;
  size_t ms,lmax=(size_t)sp->ddt.l_max;
  ms=gsl_sf_legendre_array_n(lmax);
  lm=(int)lmax;
  
  double *sphPlm =(double *)m_alloc2(ms,sizeof(double),"internal_EH(),*sphPlm");
  double *dsphPlm=(double *)m_alloc2(ms,sizeof(double),"internal_EH(),*dsphPlm");
  double complex *psi =(double complex *)m_alloc2(lm+1,sizeof(double complex),"internal_EH(),*psi");
  double complex *dpsi=(double complex *)m_alloc2(lm+1,sizeof(double complex),"internal_EH(),*dpsi");
  x=xb[0]-sp->xs;  y=xb[1]-sp->ys;  z=xb[2]-sp->zs;
  r2=x*x+y*y+z*z;    r=sqrt(r2);
  rxy=sqrt(x*x+y*y); 
  if(rxy==0.0){ // x==0,y==0
	if(z==0.0){
      internal_EH_r0(e,h,sp,bm);
      return;
    } 
    x=z*0.7e-7;
    y=z*0.7e-7;
    r2=x*x+y*y+z*z;    r=sqrt(r2);
    rxy=sqrt(x*x+y*y); 
  }
  cos_t=z/r;    sin_t=rxy/r;   i_sin_t=r/rxy;
  cos_p=x/rxy;  sin_p=y/rxy;

  ke =2.0*M_PI*sp->ns/bm->lambda_0;
  ker=ke*r;
  ne=sp->ns;
  i_ne=1.0/(sp->ns);

  rctjc(lm,ker,&nn,psi,dpsi);
  if(nn<sp->ddt.l_max) sp->ddt.l_max=nn;
  dep=cos_p+I*sin_p; expi=1.0;
  er=0.0;  et=0.0;  ep=0.0;
  hr=0.0;  ht=0.0;  hp=0.0;
  tt=0;  m=0;
  gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,lmax,cos_t,-1,sphPlm,dsphPlm);
  for(l=1;l<=lm;l++){
    ai=gsl_sf_legendre_array_index(l,m);
    Yp = sphPlm[ai];
    dYp=dsphPlm[ai];
    er+=(double)(l*(l+1))*sp->ddt.cc[l]*sp->ddt.Alm[tt]*psi[l]*Yp;
    et+=sp->ddt.cc[l]*sp->ddt.Alm[tt]*dpsi[l]*dYp-(double)m*i_ne*sp->ddt.cd[l]*sp->ddt.Blm[tt]*psi[l]*Yp*i_sin_t;
    ep+=(double)m*sp->ddt.cc[l]*sp->ddt.Alm[tt]*dpsi[l]*Yp*i_sin_t-i_ne*sp->ddt.cd[l]*sp->ddt.Blm[tt]*psi[l]*dYp;
    hr+=(double)(l*(l+1))*sp->ddt.cd[l]*sp->ddt.Blm[tt]*psi[l]*Yp;
    ht+=sp->ddt.cd[l]*sp->ddt.Blm[tt]*dpsi[l]*dYp+(double)m*ne*sp->ddt.cc[l]*sp->ddt.Alm[tt]*psi[l]*Yp*i_sin_t;
    hp+=(double)m*sp->ddt.cd[l]*sp->ddt.Blm[tt]*dpsi[l]*Yp*i_sin_t+ne*sp->ddt.cc[l]*sp->ddt.Alm[tt]*psi[l]*dYp;
    tt++;
  }
  for(m=1;m<=lm;m++){
    expi*=dep;
    if(m%2==0)flag= 1.0;
    else      flag=-1.0;
    for(l=m;l<=lm;l++){
      ai=gsl_sf_legendre_array_index(l,m);
      Yp = sphPlm[ai]*expi;
      dYp=dsphPlm[ai]*expi;
      Ym =flag*conj( Yp);
      dYm=flag*conj(dYp);
      er+=(double)(l*(l+1))*sp->ddt.cc[l]*psi[l]*(sp->ddt.Alm[tt  ]*Yp +sp->ddt.Alm[tt+1]*Ym); 
      et+=sp->ddt.cc[l]*dpsi[l]*(sp->ddt.Alm[tt  ]*dYp +sp->ddt.Alm[tt+1]*dYm)
        -(double)m*i_ne*sp->ddt.cd[l]*psi[l]*i_sin_t*(sp->ddt.Blm[tt  ]*Yp-sp->ddt.Blm[tt+1]*Ym);
      ep+=(double)m*sp->ddt.cc[l]*dpsi[l]*i_sin_t*(sp->ddt.Alm[tt  ]*Yp-sp->ddt.Alm[tt+1]*Ym)
        -i_ne*sp->ddt.cd[l]*psi[l]*(sp->ddt.Blm[tt  ]*dYp+sp->ddt.Blm[tt+1]*dYm);
      hr+=(double)(l*(l+1))*sp->ddt.cd[l]*psi[l]*(sp->ddt.Blm[tt  ]*Yp+sp->ddt.Blm[tt+1]*Ym);
      ht+=sp->ddt.cd[l]*dpsi[l]*(sp->ddt.Blm[tt  ]*dYp+sp->ddt.Blm[tt+1]*dYm)
        +(double)m*ne*sp->ddt.cc[l]*psi[l]*i_sin_t*(sp->ddt.Alm[tt  ]*Yp-sp->ddt.Alm[tt+1]*Ym);
      hp+=(double)m*sp->ddt.cd[l]*dpsi[l]*i_sin_t*(sp->ddt.Blm[tt  ]*Yp-sp->ddt.Blm[tt+1]*Ym)
        +ne*sp->ddt.cc[l]*psi[l]*(sp->ddt.Alm[tt  ]*dYp+sp->ddt.Alm[tt+1]*dYm);
      tt+=2;
    }
  }
  er/=r2;
  et*=ke/r;
  ep*=I*ke/r;
  hr/=r2;
  ht*=ke/r;
  hp*=I*ke/r;

  e[0]=er*sin_t*cos_p+et*cos_t*cos_p-ep*sin_p;
  e[1]=er*sin_t*sin_p+et*cos_t*sin_p+ep*cos_p;
  e[2]=er*cos_t-et*sin_t;
  h[0]=hr*sin_t*cos_p+ht*cos_t*cos_p-hp*sin_p;
  h[1]=hr*sin_t*sin_p+ht*cos_t*sin_p+hp*cos_p;
  h[2]=hr*cos_t-ht*sin_t;

  free(sphPlm);  free(dsphPlm);
  free(psi);      free(dpsi);
}

// verification function for sphere center field 
void internal_EH_r0(double complex *e,double complex *h,SPD *sp,Bobj *bm)
{
  double complex pickup_Alm(int l,int m,SPD *sp);
  double complex pickup_Blm(int l,int m,SPD *sp);

  double complex ke,er,et,ep,hr,ht,hp,c1m1,c1p1,c10,d1m1,d1p1,d10,Y10,Y1p1,Y1m1,dY10,dY1p1,dY1m1;  
  double cos_t,sin_t,i_sin_t,cos_p,sin_p;
  
  ke =2.0*M_PI*sp->ns/bm->lambda_0;
  
  // assume r=(x,0,0) then x to 0
  cos_t=0.0;    sin_t=1.0;   i_sin_t=1.0;
  cos_p=1.0;    sin_p=0.0;

  c10 =sp->ddt.cc[1]*pickup_Alm(1,0,sp);
  c1p1=sp->ddt.cc[1]*pickup_Alm(1,1,sp);
  c1m1=sp->ddt.cc[1]*pickup_Alm(1,-1,sp);
  d10 =sp->ddt.cd[1]*pickup_Blm(1,0,sp);
  d1p1=sp->ddt.cd[1]*pickup_Blm(1,1,sp);
  d1m1=sp->ddt.cd[1]*pickup_Blm(1,-1,sp);

  Y10 = 0.5*sqrt(3.0/M_PI)*cos_t;
  Y1p1=-0.5*sqrt(3.0/(2.0*M_PI))*sin_t*(cos_p+I*sin_p);
  Y1m1= 0.5*sqrt(3.0/(2.0*M_PI))*sin_t/(cos_p+I*sin_p);
  
  dY10 =-0.5*sqrt(3.0/M_PI)*sin_t;
  dY1p1=-0.5*sqrt(3.0/(2.0*M_PI))*cos_t*(cos_p+I*sin_p);
  dY1m1= 0.5*sqrt(3.0/(2.0*M_PI))*cos_t/(cos_p+I*sin_p);
  
  er=2.0*ke*ke/(3.0)*(c10*Y10 +c1p1*Y1p1 +c1m1*Y1m1);
  et=2.0*ke*ke/(3.0)*(c10*dY10+c1p1*dY1p1+c1m1*dY1m1);
  ep=2.0*I*ke*ke*i_sin_t/(3.0)*(c1p1*Y1p1-c1m1*Y1m1);
  hr=2.0*ke*ke/(3.0)*(d10*Y10 +d1p1*Y1p1 +d1m1*Y1m1);
  ht=2.0*ke*ke/(3.0)*(d10*dY10+d1p1*dY1p1+d1m1*dY1m1);
  hp=2.0*I*ke*ke*i_sin_t/(3.0)*(d1p1*Y1p1-d1m1*Y1m1);
  
  e[0]=er*sin_t*cos_p+et*cos_t*cos_p-ep*sin_p;
  e[1]=er*sin_t*sin_p+et*cos_t*sin_p+ep*cos_p;
  e[2]=er*cos_t-et*sin_t;
  h[0]=hr*sin_t*cos_p+ht*cos_t*cos_p-hp*sin_p;
  h[1]=hr*sin_t*sin_p+ht*cos_t*sin_p+hp*cos_p;
  h[2]=hr*cos_t-ht*sin_t;

}

double complex pickup_Alm(int l,int m,SPD *sp)
{
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;
  
  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.Alm[l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.Alm[ll  ];
      else     return sp->ddt.Alm[ll+1];
    }
  }
}

double complex pickup_Blm(int l,int m,SPD *sp)
{
  int lm=sp->l_limit;
  int am=abs(m);
  int ll;

  if(l<am)  return 0.0;
  else {
    if(m==0)
      return sp->ddt.Blm[l-1];
    else {
      ll=2*lm*am-am*am+am-2-lm+2*l;
      if (m>0) return sp->ddt.Blm[ll  ];
      else     return sp->ddt.Blm[ll+1];
    }
  }
}


