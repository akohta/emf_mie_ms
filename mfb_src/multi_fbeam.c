#include "multi_fbeam.h"

void init_mfb(Bobj *obj)
{
  obj->n_ipw=0;  sprintf(obj->fname_ipw,"%s","null");
  obj->bd.ipw=(Ipw *)m_alloc2(0,sizeof(Ipw),"multi_fbeam.c,init_mfb(),obj->bd.ipw");
  obj->n_fpw=0;  sprintf(obj->fname_fpw,"%s","null");
  obj->bd.fpw=(Fpw *)m_alloc2(0,sizeof(Fpw),"multi_fbeam.c,init_mfb(),obj->bd.fpw");
  obj->n_lgb=0;  sprintf(obj->fname_lgb,"%s","null");
  obj->bd.lgb=(LGb *)m_alloc2(0,sizeof(LGb),"multi_fbeam.c,init_mfb(),obj->bd.lgb");
  obj->n_bsb=0;  sprintf(obj->fname_bsb,"%s","null");
  obj->bd.bsb=(Bsb *)m_alloc2(0,sizeof(Bsb),"multi_fbeam.c,init_mfb(),obj->bd.bsb");
  obj->n_blg=0;  sprintf(obj->fname_blg,"%s","null");
  obj->bd.blg=(BsLGb *)m_alloc2(0,sizeof(BsLGb),"multi_fbeam.c,init_mfb(),obj->bd.blg");
  obj->n_rab=0;  sprintf(obj->fname_rab,"%s","null");  
  obj->bd.rab=(RAb *)m_alloc2(0,sizeof(RAb),"multi_fbeam.c,init_mfb(),obj->bd.rab"); 

}

void read_data_mfb(Bobj *obj)
{
  int wn=0;
  
  wn+=read_data_mfb_ipw(fn_ipw,obj); //
  wn+=read_data_mfb_fpw(fn_fpw,obj); //
  wn+=read_data_mfb_lgb(fn_lgb,obj); //
  wn+=read_data_mfb_bsb(fn_bsb,obj); //
  wn+=read_data_mfb_blg(fn_blg,obj); //
  wn+=read_data_mfb_rab(fn_rab,obj); //
  
  if(wn==0){   
    printf("read_data_mfb() No beam defined, check beam data file : %s, %s, %s, %s, %s, %s.\n",
                fn_ipw,fn_fpw,fn_lgb,fn_bsb,fn_blg,fn_rab);
    exit(1);
  }
}

void print_data_mfb(Bobj *obj)
{
  int i;
  // ipw
  for(i=0;i<obj->n_ipw;i++){
    printf(" \"%s\" No.%02d ",obj->fname_ipw,i);    print_data_ipw(&(obj->bd.ipw[i]));
  }
  // fpw
  for(i=0;i<obj->n_fpw;i++){
    printf(" \"%s\" No.%02d ",obj->fname_fpw,i);    print_data_fpw(&(obj->bd.fpw[i]));
  }
  // lgb
  for(i=0;i<obj->n_lgb;i++){
    printf(" \"%s\" No.%02d ",obj->fname_lgb,i);    print_data_lgb(&(obj->bd.lgb[i]));
  }
  // bsb
  for(i=0;i<obj->n_bsb;i++){
    printf(" \"%s\" No.%02d ",obj->fname_bsb,i);    print_data_bsb(&(obj->bd.bsb[i]));
  }
  // blg
  for(i=0;i<obj->n_blg;i++){
    printf(" \"%s\" No.%02d ",obj->fname_blg,i);    print_data_bslgb(&(obj->bd.blg[i]));
  }
  // rab
  for(i=0;i<obj->n_rab;i++){
    printf(" \"%s\" No.%02d ",obj->fname_rab,i);    print_data_rab(&(obj->bd.rab[i]));
  }

  printf("\n");
}

void print_data_mfb_mksa(Bobj *obj)
{
  int i;
  // ipw
  for(i=0;i<obj->n_ipw;i++){
    printf(" \"%s\" No.%02d ",obj->fname_ipw,i);    print_data_ipw_mksa(&(obj->bd.ipw[i]));
  }
  // fpw
  for(i=0;i<obj->n_fpw;i++){
    printf(" \"%s\" No.%02d ",obj->fname_fpw,i);    print_data_fpw_mksa(&(obj->bd.fpw[i]));
  }
  // lgb
  for(i=0;i<obj->n_lgb;i++){
    printf(" \"%s\" No.%02d ",obj->fname_lgb,i);    print_data_lgb_mksa(&(obj->bd.lgb[i]));
  }
  // bsb
  for(i=0;i<obj->n_bsb;i++){
    printf(" \"%s\" No.%02d ",obj->fname_bsb,i);    print_data_bsb_mksa(&(obj->bd.bsb[i]));
  }
  // blg
  for(i=0;i<obj->n_blg;i++){
    printf(" \"%s\" No.%02d ",obj->fname_blg,i);    print_data_bslgb_mksa(&(obj->bd.blg[i]));
  }
  // rab
  for(i=0;i<obj->n_rab;i++){
    printf(" \"%s\" No.%02d ",obj->fname_rab,i);    print_data_rab_mksa(&(obj->bd.rab[i]));
  }

  printf("\n");
} 

void setup_mfb(Bobj *obj)
{
  void check_data_mfb(Bobj *obj);
  
  int i;
  
  check_data_mfb(obj);
  // ipw
  for(i=0;i<obj->n_ipw;i++)      setup_ipw(&(obj->bd.ipw[i])); 
  // fpw
  for(i=0;i<obj->n_fpw;i++)      setup_Fpw(&(obj->bd.fpw[i]));
  // lgb
  for(i=0;i<obj->n_lgb;i++)      setup_LGb(&(obj->bd.lgb[i]));
  // bsb
  for(i=0;i<obj->n_bsb;i++)      setup_Bsb(&(obj->bd.bsb[i]));
  // blg
  for(i=0;i<obj->n_blg;i++)      setup_BsLGb(&(obj->bd.blg[i]));
  // rab
  for(i=0;i<obj->n_rab;i++)      setup_rab(&(obj->bd.rab[i]));
  
}

void free_mfb(Bobj *obj)
{
  int i;
  // ipw
  free(obj->bd.ipw);  obj->n_ipw=0; 
  // fpw
  for(i=0;i<obj->n_fpw;i++)  free_Fpw(&(obj->bd.fpw[i]));
  free(obj->bd.fpw);     obj->n_fpw=0; 
  // lgb
  for(i=0;i<obj->n_lgb;i++)  free_LGb(&(obj->bd.lgb[i]));
  free(obj->bd.lgb);     obj->n_lgb=0;
  // bsb
  for(i=0;i<obj->n_bsb;i++)  free_Bsb(&(obj->bd.bsb[i]));
  free(obj->bd.bsb);     obj->n_bsb=0;
  // blg
  for(i=0;i<obj->n_blg;i++)  free_BsLGb(&(obj->bd.blg[i]));
  free(obj->bd.blg);     obj->n_blg=0;
  // rab
  for(i=0;i<obj->n_rab;i++)  free_rab(&(obj->bd.rab[i]));
  free(obj->bd.rab);     obj->n_rab=0;

}

int read_data_mfb_ipw(char *fname,Bobj *obj)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){ 
    //printf("failed to read %s. no beams defined.\n",fname);
    return 0;
  }  
  char buf[256]="";    int i,nn;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the line. exit...\n");
    exit(1);
  }
  
  sprintf(obj->fname_ipw,"%s",fname); 
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the nn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;
  
  free(obj->bd.ipw);
  obj->n_ipw=nn;     obj->bd.ipw=(Ipw *)m_alloc2(nn,sizeof(Ipw),"multi_fbeam.c,read_data_mfb_ipw(),obj->bd.ipw");

  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the lambda0. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].lambda0=tmpd;      
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the ni. exit...\n");  
      exit(1);
    }
    obj->bd.ipw[i].ni     =tmpd;   
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].power  =tmpd;  
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the real(e0x). exit...\n");
      exit(1);
    } 
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the imag(e0x). exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].e0x    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the read(e0y). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the imag(e0y). exit...\n");
      exit(1); 
    }
    obj->bd.ipw[i].e0y    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the fx. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].fx     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].fy     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].fz     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].theta  =tmpd;
    if(fscanf(fp,"%lf\n",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_ipw(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.ipw[i].phi    =tmpd;       
  }
  fclose(fp);
  return nn;
}

int read_data_mfb_fpw(char *fname,Bobj *obj)
{
  FILE *fp;  
  if((fp=fopen(fname,"rt"))==NULL){ 
  //printf("failed to read %s. no beams defined.\n",fname);
  return 0; 
  }
  char buf[256]="";  int i,tmpi,nn;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the line. exit...\n");
    exit(1);
  }
  
  sprintf(obj->fname_fpw,"%s",fname);
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the nn. exit...\n");
    exit(1); 
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;
  
  free(obj->bd.fpw);
  obj->n_fpw=nn;     obj->bd.fpw=(Fpw *)m_alloc2(nn,sizeof(Fpw),"multi_fbeam.c,read_data_mfb_fpw(),obj->bd.fpw");

  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the lambda0. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].lambda0=tmpd;   
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the ni. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].ni     =tmpd;   
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the NA. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].NA     =tmpd;   
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].power  =tmpd;   
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the real(e0x). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the imag(e0x). exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].e0x    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the real(e0y). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the imag(e0y). exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].e0y    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the fx. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].fx     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].fy     =tmpd;
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].fz     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].theta  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].phi    =tmpd;         
    if(fscanf(fp,"%d\n",&tmpi)!=1){
      printf("multi_fbeam.c, read_data_mfb_fpw(), failed to read the nn. exit...\n");
      exit(1);
    }
    obj->bd.fpw[i].nn     =tmpi;         
  }
  fclose(fp);
  return nn;
}

int read_data_mfb_lgb(char *fname,Bobj *obj)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){
    //printf("failed to read %s. no beams defined.\n",fname);
    return 0;
  }
  char buf[256]="";  int i,nn,tmpi;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the line. exit...\n");
    exit(1);
  }
  
  sprintf(obj->fname_lgb,"%s",fname);
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the bn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;
  
  free(obj->bd.lgb);
  obj->n_lgb=nn;     obj->bd.lgb=(LGb *)m_alloc2(nn,sizeof(LGb),"multi_fbeam.c,read_data_mfb_lgb(),obj->bd.lgb");

  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the lambda0. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].lambda0=tmpd;      
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the ni. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].ni     =tmpd;       
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the NA. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].NA     =tmpd;      
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].power  =tmpd;      
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the real(e0x). exit...\n");
      exit(1);
    } 
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the imag(e0x). exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].e0x    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the real(e0y). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the imag(e0y). exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].e0y    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the fx. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].fx     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].fy     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].fz     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].theta  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].phi    =tmpd;         
    if(fscanf(fp,"%d" ,&tmpi)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the lg_m. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].lg_m   =tmpi;         
    if(fscanf(fp,"%d\n" ,&tmpi)!=1){
      printf("multi_fbeam.c, read_data_mfb_lgb(), failed to read the nn. exit...\n");
      exit(1);
    }
    obj->bd.lgb[i].nn     =tmpi;         
  }
  fclose(fp);
  return nn;
}

int read_data_mfb_bsb(char *fname,Bobj *obj)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){ 
    //printf("failed to read %s. no beams defined.\n",fname);
    return 0;
  }
  char buf[256]="";  int i,nn;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the line. exit...\n");
    exit(1);
  }
  sprintf(obj->fname_bsb,"%s",fname);
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the bn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;
  
  free(obj->bd.bsb);
  obj->n_bsb=nn;     obj->bd.bsb=(Bsb *)m_alloc2(nn,sizeof(Bsb),"multi_fbeam.c,read_data_mfb_bsb(),obj->bd.bsb");
 
  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the lambda0. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].lambda0=tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the ni. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].ni     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the d_angle. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].d_angle=tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].power  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the real(e0x). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the imag(e0x). exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].e0x    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the real(e0y). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the imag(e0y). exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].e0y    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the fx. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].fx     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].fy     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].fz     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].theta  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_bsb(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.bsb[i].phi    =tmpd;         
  }
  fclose(fp);
  return nn;
}

int read_data_mfb_blg(char *fname,Bobj *obj)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){
    //printf("failed to read %s. no beams defined.\n",fname);
    return 0;
  }
  char buf[256]="";  int i,nn,tmpi;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the line. exit...\n");
    exit(1);
  }
  
  sprintf(obj->fname_blg,"%s",fname);
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the bn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;
  
  free(obj->bd.blg);
  obj->n_blg=nn;  obj->bd.blg=(BsLGb *)m_alloc2(nn,sizeof(BsLGb),"multi_fbeam.c,read_data_mfb_blg(),obj->bd.blg");

  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the lambda0. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].lambda0=tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the ni. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].ni     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the d_angle. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].d_angle=tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].power  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the real(e0x). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the imag(e0x). exit...\n");
      exit(1);
    }
    obj->bd.blg[i].e0x    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the real(e0y). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the imag(e0y). exit...\n");
      exit(1);
    }
    obj->bd.blg[i].e0y    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the fx. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].fx     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].fy     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].fz     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].theta  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].phi    =tmpd;         
    if(fscanf(fp,"%d" ,&tmpi)!=1){
      printf("multi_fbeam.c, read_data_mfb_blg(), failed to read the lg_m. exit...\n");
      exit(1);
    }
    obj->bd.blg[i].lg_m   =tmpi;         
  }
  fclose(fp);
  return nn;
}

int read_data_mfb_rab(char *fname,Bobj *obj)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){
    //printf("failed to read %s. no beams defined.\n",fname);
    return 0;
  }
  char buf[256]="";  int i,nn,tmpi;  double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the line. exit...\n");
    exit(1);
  }
  
  sprintf(obj->fname_rab,"%s",fname);
  if(fscanf(fp,"%d\n",&nn)!=1){
    printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the bn. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the line. exit...\n");
    exit(1);
  }
  if(nn==0) return 0;
  
  free(obj->bd.rab);
  obj->n_rab=nn;  obj->bd.rab=(RAb *)m_alloc2(nn,sizeof(RAb),"multi_fbeam.c,read_data_mfb_rab(),obj->bd.rab");

  for(i=0;i<nn;i++){
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the labmda0. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].lambda0=tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the ni. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].ni     =tmpd;          
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the NA. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].NA     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the power. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].power  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the real(e0r). exit...\n");
      exit(1);
    } 
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the imag(e0r). exit...\n");
      exit(1);
    }
    obj->bd.rab[i].e0r    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the real(e0a). exit...\n");
      exit(1);
    }
    if(fscanf(fp,"%lf",&tmpd2)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the imag(e0a). exit...\n");
      exit(1);
    }
    obj->bd.rab[i].e0a    =tmpd+I*tmpd2; 
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the fx. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].fx     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the fy. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].fy     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the fz. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].fz     =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the theta. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].theta  =tmpd;         
    if(fscanf(fp,"%lf",&tmpd)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the phi. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].phi    =tmpd;         
    if(fscanf(fp,"%d\n",&tmpi)!=1){
      printf("multi_fbeam.c, read_data_mfb_rab(), failed to read the nn. exit...\n");
      exit(1);
    }
    obj->bd.rab[i].nn     =tmpi;         
  }
  fclose(fp);
  return nn;
}

void calc_mfb_EH(double complex *e,double complex *h,double *x,Bobj *obj)
{
  double complex te[3],th[3];
  int i,j;
  for(i=0;i<3;i++){
    e[i]=0.0;    h[i]=0.0;
  }
  // ipw
  for(i=0;i<obj->n_ipw;i++){
    calc_ipw_EH(te,th,x,&(obj->bd.ipw[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  // fpw
  for(i=0;i<obj->n_fpw;i++){
    calc_fpw_EH(te,th,x,&(obj->bd.fpw[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  // lgb
  for(i=0;i<obj->n_lgb;i++){
    calc_lgb_EH(te,th,x,&(obj->bd.lgb[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  // bsb
  for(i=0;i<obj->n_bsb;i++){
    calc_bsb_EH(te,th,x,&(obj->bd.bsb[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  
  // blg
  for(i=0;i<obj->n_blg;i++){
    calc_bslgb_EH(te,th,x,&(obj->bd.blg[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  // rab
  for(i=0;i<obj->n_rab;i++){
    calc_rab_EH(te,th,x,&(obj->bd.rab[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];      h[j]+=th[j];
    }
  }
  
}

void calc_mfb_EH_dv(double complex *e,double complex *h,double complex *dedv,double complex *dhdv,double *x,double *v,Bobj *obj)
{
  double complex te[3], th[3],tde[3],tdh[3];
  int i, j;
  
  for(i=0;i<3;i++){
    e[i]=0.0;    h[i]=0.0;
    dedv[i]=0.0;    dhdv[i]=0.0;
  }
  // ipw
  for(i=0;i<obj->n_ipw;i++){
    calc_ipw_EH_dv(te,th,tde,tdh,x,v,&(obj->bd.ipw[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j]+=th[j];
      dedv[j]+=tde[j];    dhdv[j]+=tdh[j];
    }
  }
  // fpw
  for(i=0;i<obj->n_fpw;i++){
    calc_fpw_EH_dv(te,th,tde,tdh,x,v,&(obj->bd.fpw[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j] += th[j];
      dedv[j]+=tde[j];    dhdv[j]+=tdh[j];
    }
  }
  // lgb
  for(i=0;i<obj->n_lgb;i++){
    calc_lgb_EH_dv(te,th,tde,tdh,x,v,&(obj->bd.lgb[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j]+=th[j];
      dedv[j]+=tde[j];    dhdv[j]+=tdh[j];
    }
  }
  // bsb
  for(i=0;i<obj->n_bsb;i++){
    calc_bsb_EH_dv(te,th,tde,tdh,x,v,&(obj->bd.bsb[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j]+=th[j];
      dedv[j]+=tde[j];    dhdv[j]+=tdh[j];
    }
  }
  // blg
  for(i=0;i<obj->n_blg;i++){
    calc_bslgb_EH_dv(te,th,tde,tdh,x,v,&(obj->bd.blg[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j]+=th[j];
      dedv[j]+=tde[j];    dhdv[j]+=tdh[j];
    }
  }
  // rab
  for(i=0;i<obj->n_rab;i++){
    calc_rab_EH_dv(te,th,tde,tdh,x,v,&(obj->bd.rab[i]));
    for(j=0;j<3;j++){
      e[j]+=te[j];    h[j]+=th[j];
      dedv[j]+=tde[j];    dhdv[j]+=tdh[j];
    }
  }

}

////////////////////////////////////////////////////////////
void check_data_mfb(Bobj *obj)
{
  int i;
  int st=0;
  obj->n_0     =0.0;
  obj->lambda_0=0.0;
  //ipw
  for(i=0;i<obj->n_ipw;i++){
    if(st==0){
      obj->n_0     =obj->bd.ipw[i].ni;
      obj->lambda_0=obj->bd.ipw[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.ipw[i].ni     ){ printf("data error ipw ni, exit\n");     exit(1);} 
      if(obj->lambda_0!=obj->bd.ipw[i].lambda0){ printf("data error ipw lambda0, exit\n"); exit(1);}
    }
  }  
  //fpw
  for(i=0;i<obj->n_fpw;i++){
    if(st==0){
      obj->n_0     =obj->bd.fpw[i].ni;
      obj->lambda_0=obj->bd.fpw[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.fpw[i].ni     ){ printf("data error fpw ni, exit\n");     exit(1);}
      if(obj->lambda_0!=obj->bd.fpw[i].lambda0){ printf("data error fpw lambda0, exit\n"); exit(1);}
    }
  }
  // lgb
  for(i=0;i<obj->n_lgb;i++){
    if(st==0){
      obj->n_0     =obj->bd.lgb[i].ni;
      obj->lambda_0=obj->bd.lgb[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.lgb[i].ni     ){ printf("data error lgb ni, exit\n");     exit(1);}
      if(obj->lambda_0!=obj->bd.lgb[i].lambda0){ printf("data error lgb lambda0, exit\n"); exit(1);}
    }
  }
  // bsb
  for(i=0;i<obj->n_bsb;i++){
    if(st==0){
      obj->n_0     =obj->bd.bsb[i].ni;
      obj->lambda_0=obj->bd.bsb[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.bsb[i].ni     ){ printf("data error bsb ni, exit\n");     exit(1);}
      if(obj->lambda_0!=obj->bd.bsb[i].lambda0){ printf("data error bsb lambda0, exit\n"); exit(1);}
    }
  }
  // blg
  for(i=0;i<obj->n_blg;i++){
    if(st==0){
      obj->n_0     =obj->bd.blg[i].ni;
      obj->lambda_0=obj->bd.blg[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.blg[i].ni     ){ printf("data error blg ni, exit\n");     exit(1);}
      if(obj->lambda_0!=obj->bd.blg[i].lambda0){ printf("data error blg lambda0, exit\n"); exit(1);}
    }
  }
  // rab
  for(i=0;i<obj->n_rab;i++){
    if(st==0){
      obj->n_0     =obj->bd.rab[i].ni;
      obj->lambda_0=obj->bd.rab[i].lambda0;
      st++;
    }
    else {
      if(obj->n_0     !=obj->bd.rab[i].ni     ){ printf("data error rab ni, exit\n");     exit(1);}
      if(obj->lambda_0!=obj->bd.rab[i].lambda0){ printf("data error rab lambda0, exit\n"); exit(1);}
    }
  }
  
  if(st==0){
    printf("no beams defined. check datafile name. Exit...\n");
    exit(0);
  }
  obj->omega=2.0*M_PI/obj->lambda_0;
}
