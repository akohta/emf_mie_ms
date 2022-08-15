#include <stdio.h>
#include <stdlib.h>
#include "ipw.h"
#include "fpw.h"
#include "lgb.h"
#include "bsb.h"
#include "bslgb.h"
#include "rab.h"
#include "fgb.h"
#include "my_utils.h"
#include "emf_mie_ms.h"

void write_dat_ms(char *fn,MSPD *msp)
{
  FILE *fp;
  int i,np,nt,mn;
  
  if((fp=fopen(fn,"wb"))==NULL){    printf("write_dat_ms(), Failed to open the file %s. Exit...\n",fn);    exit(1);  }
  
  if(fwrite(msp,sizeof(MSPD),1,fp)!=1){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write msp. exit...\n");
    exit(1);
  }
  // beam data
  if(fwrite(msp->bm.bd.ipw,sizeof(Ipw),msp->bm.n_ipw,fp)!=msp->bm.n_ipw){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write ipw. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.fpw,sizeof(Fpw),msp->bm.n_fpw,fp)!=msp->bm.n_fpw){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write fpw. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.lgb,sizeof(LGb),msp->bm.n_lgb,fp)!=msp->bm.n_lgb){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write lgb. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.bsb,sizeof(Bsb),msp->bm.n_bsb,fp)!=msp->bm.n_bsb){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write bsb. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.blg,sizeof(BsLGb),msp->bm.n_blg,fp)!=msp->bm.n_blg){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write blg. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.rab,sizeof(RAb),msp->bm.n_rab,fp)!=msp->bm.n_rab){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write rab. exit...\n");
    exit(1);
  }
  if(fwrite(msp->bm.bd.fgb,sizeof(Fgb),msp->bm.n_fgb,fp)!=msp->bm.n_fgb){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write fgb. exit...\n");
    exit(1);
  }
  // sphere data  
  if(fwrite(msp->sp,sizeof(SPD),msp->n_sphr,fp)!=msp->n_sphr){
    printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write sp. exit...\n");
    exit(1);
  }
  for(i=0;i<msp->n_sphr;i++){
    np=msp->sp[i].ddt.np;
    nt=msp->sp[i].ddt.nt;
    mn=msp->sp[i].l_limit;
    if(fwrite(msp->sp[i].ddt.eri,sizeof(double complex),np*nt,fp)!=np*nt){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write eri. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.hri,sizeof(double complex),np*nt,fp)!=np*nt){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write hri. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cab,sizeof(double),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write cab. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.ca,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write ca. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cb,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write cb. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cc,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write cc. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.cd,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write cd. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.Alm,sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write Alm. exit...\n");
      exit(1);
    }
    if(fwrite(msp->sp[i].ddt.Blm,sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
      printf("emf_mie_ms_dat.c, write_dat_ms(), failed to write Blm. exit...\n");
      exit(1);
    }
  }
  fclose(fp);
}

void  read_dat_ms(char *fn,MSPD *msp)
{
  void setup_sp(SPD *sp); // emf_mie_ms.c

  FILE *fp;
  int mn,i,nt,np;
  if((fp=fopen(fn,"rb"))==NULL){    printf("read_dat_ms(), Failed to open the %s. Exit...\n",fn);    exit(1);  }
  
  if(fread(msp,sizeof(MSPD),1,fp)!=1){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read msp. exit...\n");
    exit(1);
  }
    ;
  // beam data
  msp->bm.bd.ipw=(Ipw *)m_alloc2(msp->bm.n_ipw,sizeof(Ipw),"read_dat_ms(),msp->bm.bd.ipw");
  if(fread(msp->bm.bd.ipw,sizeof(Ipw),msp->bm.n_ipw,fp)!=msp->bm.n_ipw){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read ipw. exit...\n");
    exit(1);
  }
  msp->bm.bd.fpw=(Fpw *)m_alloc2(msp->bm.n_fpw,sizeof(Fpw),"read_dat_ms(),msp->bm.bd.fpw");
  if(fread(msp->bm.bd.fpw,sizeof(Fpw),msp->bm.n_fpw,fp)!=msp->bm.n_fpw){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read fpw. exit...\n");
    exit(1);
  }
  msp->bm.bd.lgb=(LGb *)m_alloc2(msp->bm.n_lgb,sizeof(LGb),"read_dat_ms(),msp->bm.bd.lgb");
  if(fread(msp->bm.bd.lgb,sizeof(LGb),msp->bm.n_lgb,fp)!=msp->bm.n_lgb){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read lgb. exit...\n");
    exit(1);
  }
  msp->bm.bd.bsb=(Bsb *)m_alloc2(msp->bm.n_bsb,sizeof(Bsb),"read_dat_ms(),msp->bm.bd.bsb");
  if(fread(msp->bm.bd.bsb,sizeof(Bsb),msp->bm.n_bsb,fp)!=msp->bm.n_bsb){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read bsb. exit...\n");
    exit(1);
  }
  msp->bm.bd.blg=(BsLGb *)m_alloc2(msp->bm.n_blg,sizeof(BsLGb),"read_dat_ms(),msp->bm.bd.blg");
  if(fread(msp->bm.bd.blg,sizeof(BsLGb),msp->bm.n_blg,fp)!=msp->bm.n_blg){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read blg. exit...\n");
    exit(1);
  }
  msp->bm.bd.rab=(RAb *)m_alloc2(msp->bm.n_rab,sizeof(RAb),"read_dat_ms(),msp->bm.bd.rab");
  if(fread(msp->bm.bd.rab,sizeof(RAb),msp->bm.n_rab,fp)!=msp->bm.n_rab){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read rab. exit...\n");
    exit(1);
  }
  msp->bm.bd.fgb=(Fgb *)m_alloc2(msp->bm.n_fgb,sizeof(Fgb),"read_dat_ms(),msp->bm.bd.fgb");
  if(fread(msp->bm.bd.fgb,sizeof(Fgb),msp->bm.n_fgb,fp)!=msp->bm.n_fgb){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read fgb. exit...\n");
    exit(1);
  } 
  setup_mfb(&(msp->bm));
  // sphere data
  msp->sp=(SPD *)m_alloc2(msp->n_sphr,sizeof(SPD),"read_dat_ms(),msp->sp");
  if(fread(msp->sp,sizeof(SPD),msp->n_sphr,fp)!=msp->n_sphr){
    printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read sp. exit...\n");
    exit(1);
  }
  for(i=0;i<msp->n_sphr;i++){
    setup_sp(&(msp->sp[i]));
    nt=msp->sp[i].ddt.nt;
    np=msp->sp[i].ddt.np;
    mn=msp->sp[i].l_limit;
    if(fread(msp->sp[i].ddt.eri,sizeof(double complex),nt*np,fp)!=nt*np){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read eri. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.hri,sizeof(double complex),nt*np,fp)!=nt*np){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read hri. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cab,sizeof(double),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read cab. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.ca,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read ca. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cb,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read cb. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cc,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read cc. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.cd,sizeof(double complex),mn+1,fp)!=mn+1){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read cd. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.Alm,sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read Alm. exit...\n");
      exit(1);
    }
    if(fread(msp->sp[i].ddt.Blm,sizeof(double complex),mn*(mn+2),fp)!=mn*(mn+2)){
      printf("emf_mie_ms_dat.c, read_dat_ms(), failed to read Blm. exit...\n");
      exit(1);
    }
  }
  fclose(fp); 
}
