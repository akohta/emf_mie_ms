#include "emf_mie_ms.h"

void absorbed_energy_ms(int id,double *p,MSPD *msp)
{
  void absorb_eng(double *p,SPD *sp,Bobj *bm);
  
  if(id>=msp->n_sphr || id<0){
    printf("absorbed_energy_ms(), Sphere ID error! ID = %d. Exit...\n",id);
    exit(1);
  }
  absorb_eng(p,&(msp->sp[id]),&(msp->bm));
}

///////////////////////////////////////////////////////////////////////
void absorb_eng(double *p,SPD *sp,Bobj *bm)
{
  double complex pickup_Alm(int l,int m,SPD *sp); // emf_mie_ms.c
  double complex pickup_Blm(int l,int m,SPD *sp); // emf_mie_ms.c

  double complex Alm,Blm,alm,blm;
  double k0=2.0*M_PI/bm->lambda_0;
  double ne=bm->n_0;
  int lm=sp->l_limit;
  int l,m;
  
  *p=0.0;
  
  m=0;
  for(l=1;l<lm;l++){
    Alm=pickup_Alm(l,m,sp);          alm=sp->ddt.ca[l]*Alm;
    Blm=pickup_Blm(l,m,sp);          blm=sp->ddt.cb[l]*Blm;
    *p+=(double)(l*(l+1))*(creal(ne*ne*alm*conj(alm))+creal(blm*conj(blm))+creal(ne*ne*Alm*conj(alm)+Blm*conj(blm)));
  }

  for(m=1;m<lm;m++){
    for(l=m;l<lm;l++){
      Alm=pickup_Alm(l,m,sp);          alm=sp->ddt.ca[l]*Alm;
      Blm=pickup_Blm(l,m,sp);          blm=sp->ddt.cb[l]*Blm;
      *p+=(double)(l*(l+1))*(creal(ne*ne*alm*conj(alm))+creal(blm*conj(blm))+creal(ne*ne*Alm*conj(alm)+Blm*conj(blm)));

      Alm=pickup_Alm(l,-m,sp);         alm=sp->ddt.ca[l]*Alm;
      Blm=pickup_Blm(l,-m,sp);         blm=sp->ddt.cb[l]*Blm;
      *p+=(double)(l*(l+1))*(creal(ne*ne*alm*conj(alm))+creal(blm*conj(blm))+creal(ne*ne*Alm*conj(alm)+Blm*conj(blm)));
    }
  }
  
  *p*=-0.5*k0*k0*ne;
}
