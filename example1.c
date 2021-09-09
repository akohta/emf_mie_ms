// calculation example using a datafile outputed by 'mie_ms_solver'.
// verification of radiation force and torque.
#include "emf_mie_ms.h"

// radiation force and torque are calculated by surface integral of maxwell stress tensor
int force_torque_integral(int i,double *vf,double *vn,MSPD *msp);

int main(int argc,char *argv[]) 
{
  MSPD msp;
  double complex e[3],h[3],cet,fv;
  double vf[3],vn[3],r[3],t;
  int i;
  
  read_dat_ms(argv[1],&msp); // read data file
  print_data_ms(&msp);       // print data
  print_data_ms_mksa(&msp);  // print data in MKSA system of units
  
  r[0]= 0.0; // set x-coordinate 
  r[1]= 0.0; // set y-coordinate
  r[2]=-1.5; // set z-coordinate
  total_EH_ms(e,h,r,&msp); // calclation of total field ( add incident field to scattered field )
  printf("Electromagnetic field at r=( % g,% g,% g )\n",r[0],r[1],r[2]);
  fv=OSUtoMKSA_ElectricField(e[0]);
  printf("Ex = % 15.14e %+15.14e I (=% 15.14e %+15.14e I [V/m](MKSA))\n",creal(e[0]),cimag(e[0]),creal(fv),cimag(fv));
  fv=OSUtoMKSA_ElectricField(e[1]);
  printf("Ey = % 15.14e %+15.14e I (=% 15.14e %+15.14e I [V/m](MKSA))\n",creal(e[1]),cimag(e[1]),creal(fv),cimag(fv));
  fv=OSUtoMKSA_ElectricField(e[2]);
  printf("Ez = % 15.14e %+15.14e I (=% 15.14e %+15.14e I [V/m](MKSA))\n",creal(e[2]),cimag(e[2]),creal(fv),cimag(fv));
  fv=OSUtoMKSA_MagneticField(h[0]);
  printf("Hx = % 15.14e %+15.14e I (=% 15.14e %+15.14e I [A/m](MKSA))\n",creal(h[0]),cimag(h[0]),creal(fv),cimag(fv));
  fv=OSUtoMKSA_MagneticField(h[1]);
  printf("Hy = % 15.14e %+15.14e I (=% 15.14e %+15.14e I [A/m](MKSA))\n",creal(h[1]),cimag(h[1]),creal(fv),cimag(fv));
  fv=OSUtoMKSA_MagneticField(h[2]);
  printf("Hz = % 15.14e %+15.14e I (=% 15.14e %+15.14e I [A/m](MKSA))\n",creal(h[2]),cimag(h[2]),creal(fv),cimag(fv));
  
  t=1.0e-3; // set time 
  cet=cexp(-I*msp.bm.omega*t);
  printf("Real electromagnetic field at t=%g\n",t);
  fv=OSUtoMKSA_ElectricField(e[0]*cet);
  printf("Ex = % 15.14e (=% 15.14e [V/m](MKSA))\n",creal(e[0]*cet),creal(fv));
  fv=OSUtoMKSA_ElectricField(e[1]*cet);
  printf("Ey = % 15.14e (=% 15.14e [V/m](MKSA))\n",creal(e[1]*cet),creal(fv));
  fv=OSUtoMKSA_ElectricField(e[2]*cet);
  printf("Ez = % 15.14e (=% 15.14e [V/m](MKSA))\n",creal(e[2]*cet),creal(fv));
  fv=OSUtoMKSA_MagneticField(h[0]*cet);
  printf("Hx = % 15.14e (=% 15.14e [A/m](MKSA))\n",creal(h[0]*cet),creal(fv));
  fv=OSUtoMKSA_MagneticField(h[1]*cet);
  printf("Hy = % 15.14e (=% 15.14e [A/m](MKSA))\n",creal(h[1]*cet),creal(fv));
  fv=OSUtoMKSA_MagneticField(h[2]*cet);
  printf("Hz = % 15.14e (=% 15.14e [A/m](MSKA))\n",creal(h[2]*cet),creal(fv));
  printf("\n");
  
  printf("Radiation force and torque\n");
  for(i=0;i<msp.n_sphr;i++){
    force_torque_ms(i,vf,vn,&msp);
    printf("Mie coefficient\n");
    printf("sphere id %2d, F=( % 15.14g,% 15.14g,% 15.14g )\n",i,vf[0],vf[1],vf[2]);
    printf("          %2d, N=( % 15.14g,% 15.14g,% 15.14g )\n",i,vn[0],vn[1],vn[2]);
    printf("          %2d, F=( % 15.14g,% 15.14g,% 15.14g ) [ N ](MKSA)\n",i,OSUtoMKSA_Force(vf[0]),OSUtoMKSA_Force(vf[1]),OSUtoMKSA_Force(vf[2]));
    printf("          %2d, N=( % 15.14g,% 15.14g,% 15.14g ) [N m](MKSA)\n",i,OSUtoMKSA_Torque(vn[0]),OSUtoMKSA_Torque(vn[1]),OSUtoMKSA_Torque(vn[2]));
    
    if(force_torque_integral(i,vf,vn,&msp)){ // for verification
      printf("Surface integral\n");
      printf("sphere id %2d, F=( % 15.14g,% 15.14g,% 15.14g ) [ N ]\n",i,vf[0],vf[1],vf[2]);
      printf("          %2d, N=( % 15.14g,% 15.14g,% 15.14g ) [N m]\n",i,vn[0],vn[1],vn[2]);
      printf("          %2d, F=( % 15.14g,% 15.14g,% 15.14g ) [ N ](MKSA)\n",i,OSUtoMKSA_Force(vf[0]),OSUtoMKSA_Force(vf[1]),OSUtoMKSA_Force(vf[2]));
      printf("          %2d, N=( % 15.14g,% 15.14g,% 15.14g ) [N m](MKSA)\n",i,OSUtoMKSA_Torque(vn[0]),OSUtoMKSA_Torque(vn[1]),OSUtoMKSA_Torque(vn[2]));
    }
  }
  
  free_ms(&msp); // free allocated memory
  return 0;
}

int force_torque_integral(int id,double *vf,double *vn,MSPD *msp)
{
  const int nc=80;
  double xt[nc],wt[nc],xp[2*nc],wp[2*nc];
  double complex e[3],h[3];
  double rc,r[3],Tx[3],Ty[3],Tz[3],sin_t,cos_t,sin_p,cos_p,eps,mu,aex2,aey2,aez2,ahx2,ahy2,ahz2,ne2,nh2;
  double tfx,tfy,tfz,tnx,tny,tnz;
  int i,j;
  
  vf[0]=0.0;  vf[1]=0.0;  vf[2]=0.0;
  vn[0]=0.0;  vn[1]=0.0;  vn[2]=0.0;
  
  if(msp->n_sphr!=1){
    //printf("this code can analize single sphere only. Return...\n");
    return 0;
  }
    
  eps=msp->bm.n_0*msp->bm.n_0;
  mu=1.0;
  
  gauleg(0.0,M_PI,xt,wt,nc);
  gauleg(0.0,2*M_PI,xp,wp,nc*2);
  
  rc=msp->sp[id].a*2.0;
  
  for(i=0;i<2*nc;i++){ // phi 0 to 2pi
    sin_p=sin(xp[i]);    cos_p=cos(xp[i]);
    tfx=0.0;    tfy=0.0;    tfz=0.0;
    tnx=0.0;    tny=0.0;    tnz=0.0;
    for(j=0;j<nc;j++){ // theta 0 to pi
      sin_t=sin(xt[j]);      cos_t=cos(xt[j]);
      
      r[0]=rc*sin_t*cos_p;
      r[1]=rc*sin_t*sin_p;
      r[2]=rc*cos_t;
      total_EH_ms(e,h,r,msp);
      aex2=creal(e[0]*conj(e[0]));      aey2=creal(e[1]*conj(e[1]));      aez2=creal(e[2]*conj(e[2]));
      ahx2=creal(h[0]*conj(h[0]));      ahy2=creal(h[1]*conj(h[1]));      ahz2=creal(h[2]*conj(h[2]));
      ne2=aex2+aey2+aez2;      nh2=ahx2+ahy2+ahz2;
      // maxwell stress tensor
      Tx[0]=0.5*(eps*aex2+mu*ahx2)-0.25*(eps*ne2+mu*nh2);
      Tx[1]=0.5*(eps*creal(e[0]*conj(e[1]))+mu*creal(h[0]*conj(h[1])));
      Tx[2]=0.5*(eps*creal(e[0]*conj(e[2]))+mu*creal(h[0]*conj(h[2])));
      Ty[0]=0.5*(eps*creal(e[1]*conj(e[0]))+mu*creal(h[1]*conj(h[0])));
      Ty[1]=0.5*(eps*aey2+mu*ahy2)-0.25*(eps*ne2+mu*nh2);
      Ty[2]=0.5*(eps*creal(e[1]*conj(e[2]))+mu*creal(h[1]*conj(h[2])));
      Tz[0]=0.5*(eps*creal(e[2]*conj(e[0]))+mu*creal(h[2]*conj(h[0])));
      Tz[1]=0.5*(eps*creal(e[2]*conj(e[1]))+mu*creal(h[2]*conj(h[1])));
      Tz[2]=0.5*(eps*aez2+mu*ahz2)-0.25*(eps*ne2+mu*nh2);
      
      tfx+=(Tx[0]*sin_t*cos_p+Tx[1]*sin_t*sin_p+Tx[2]*cos_t)*rc*rc*sin_t*wt[j];
      tfy+=(Ty[0]*sin_t*cos_p+Ty[1]*sin_t*sin_p+Ty[2]*cos_t)*rc*rc*sin_t*wt[j];
      tfz+=(Tz[0]*sin_t*cos_p+Tz[1]*sin_t*sin_p+Tz[2]*cos_t)*rc*rc*sin_t*wt[j];
      
      tnx+=( (Tz[0]*sin_t*cos_p+Tz[1]*sin_t*sin_p+Tz[2]*cos_t)*sin_t*sin_p
            -(Ty[0]*sin_t*cos_p+Ty[1]*sin_t*sin_p+Ty[2]*cos_t)*cos_t)*rc*rc*rc*sin_t*wt[j];
      tny+=( (Tx[0]*sin_t*cos_p+Tx[1]*sin_t*sin_p+Tx[2]*cos_t)*cos_t
            -(Tz[0]*sin_t*cos_p+Tz[1]*sin_t*sin_p+Tz[2]*cos_t)*sin_t*cos_p)*rc*rc*rc*sin_t*wt[j];
      tnz+=( (Ty[0]*sin_t*cos_p+Ty[1]*sin_t*sin_p+Ty[2]*cos_t)*sin_t*cos_p
            -(Tx[0]*sin_t*cos_p+Tx[1]*sin_t*sin_p+Tx[2]*cos_t)*sin_t*sin_p)*rc*rc*rc*sin_t*wt[j];
    }
    vf[0]+=tfx*wp[i];
    vf[1]+=tfy*wp[i];
    vf[2]+=tfz*wp[i];
    vn[0]+=tnx*wp[i];
    vn[1]+=tny*wp[i];
    vn[2]+=tnz*wp[i];
  }
  
  return 1;
}
