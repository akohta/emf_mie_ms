#include "multi_fbeam.h"

void init_mfb(Bobj *obj)
{
  obj->n_ipw=0;  sprintf(obj->fname_ipw,"%s","null");  obj->bd.ipw=(Ipw *)malloc(sizeof(Ipw)*0);
  obj->n_fpw=0;  sprintf(obj->fname_fpw,"%s","null");  obj->bd.fpw=(Fpw *)malloc(sizeof(Fpw)*0);
  obj->n_lgb=0;  sprintf(obj->fname_lgb,"%s","null");  obj->bd.lgb=(LGb *)malloc(sizeof(LGb)*0);
  obj->n_bsb=0;  sprintf(obj->fname_bsb,"%s","null");  obj->bd.bsb=(Bsb *)malloc(sizeof(Bsb)*0);
  obj->n_blg=0;  sprintf(obj->fname_blg,"%s","null");  obj->bd.blg=(BsLGb *)malloc(sizeof(BsLGb)*0);
  obj->n_rab=0;  sprintf(obj->fname_rab,"%s","null");  obj->bd.rab=(RAb *)malloc(sizeof(RAb)*0); 

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
  for(i=0;i<obj->n_fpw;i++)  free_Fpw(&(obj->bd.fpw[i]));     free(obj->bd.fpw);     obj->n_fpw=0; 
  // lgb
  for(i=0;i<obj->n_lgb;i++)  free_LGb(&(obj->bd.lgb[i]));     free(obj->bd.lgb);     obj->n_lgb=0;
  // bsb
  for(i=0;i<obj->n_bsb;i++)  free_Bsb(&(obj->bd.bsb[i]));     free(obj->bd.bsb);     obj->n_bsb=0;
  // blg
  for(i=0;i<obj->n_blg;i++)  free_BsLGb(&(obj->bd.blg[i]));   free(obj->bd.blg);     obj->n_blg=0;
  // rab
  for(i=0;i<obj->n_rab;i++)  free_rab(&(obj->bd.rab[i]));     free(obj->bd.rab);     obj->n_rab=0;

}

int read_data_mfb_ipw(char *fname,Bobj *obj)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){ 
    //printf("failed to read %s. no beams defined.\n",fname);
    return 0;
  }  
  char buf[256]="";    int i,nn;  double tmpd,tmpd2;
  fgets(buf,256,fp);  fgets(buf,256,fp);
  
  sprintf(obj->fname_ipw,"%s",fname); 
  fscanf(fp,"%d\n",&nn);  fgets(buf,256,fp);
  if(nn==0) return 0;
  
  free(obj->bd.ipw);
  obj->n_ipw=nn;     obj->bd.ipw=(Ipw *)malloc(sizeof(Ipw)*nn);

  for(i=0;i<nn;i++){
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].lambda0=tmpd;      
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].ni     =tmpd;   
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].power  =tmpd;  
    fscanf(fp,"%lf",&tmpd); 
    fscanf(fp,"%lf",&tmpd2); obj->bd.ipw[i].e0x    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.ipw[i].e0y    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].fx     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].fy     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].fz     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.ipw[i].theta  =tmpd;         
    fscanf(fp,"%lf\n",&tmpd);  obj->bd.ipw[i].phi    =tmpd;       
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
  fgets(buf,256,fp);  fgets(buf,256,fp);
  
  sprintf(obj->fname_fpw,"%s",fname);
  fscanf(fp,"%d\n",&nn);  fgets(buf,256,fp);
  if(nn==0) return 0;
  
  free(obj->bd.fpw);
  obj->n_fpw=nn;     obj->bd.fpw=(Fpw *)malloc(sizeof(Fpw)*nn);

  for(i=0;i<nn;i++){
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].lambda0=tmpd;   
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].ni     =tmpd;   
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].NA     =tmpd;   
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].power  =tmpd;   
    fscanf(fp,"%lf",&tmpd); 
    fscanf(fp,"%lf",&tmpd2); obj->bd.fpw[i].e0x    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.fpw[i].e0y    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].fx     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].fy     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].fz     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].theta  =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.fpw[i].phi    =tmpd;         
    fscanf(fp,"%d\n",&tmpi); obj->bd.fpw[i].nn     =tmpi;         
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
  fgets(buf,256,fp);  fgets(buf,256,fp);
  
  sprintf(obj->fname_lgb,"%s",fname);
  fscanf(fp,"%d\n",&nn);  fgets(buf,256,fp);
  if(nn==0) return 0;
  
  free(obj->bd.lgb);
  obj->n_lgb=nn;     obj->bd.lgb=(LGb *)malloc(sizeof(LGb)*nn);

  for(i=0;i<nn;i++){
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].lambda0=tmpd;      
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].ni     =tmpd;       
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].NA     =tmpd;      
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].power  =tmpd;      
    fscanf(fp,"%lf",&tmpd); 
    fscanf(fp,"%lf",&tmpd2);  obj->bd.lgb[i].e0x    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2);  obj->bd.lgb[i].e0y    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].fx     =tmpd;         
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].fy     =tmpd;         
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].fz     =tmpd;         
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].theta  =tmpd;         
    fscanf(fp,"%lf",&tmpd);   obj->bd.lgb[i].phi    =tmpd;         
    fscanf(fp,"%d" ,&tmpi);   obj->bd.lgb[i].lg_m   =tmpi;         
    fscanf(fp,"%d\n" ,&tmpi); obj->bd.lgb[i].nn     =tmpi;         
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
  fgets(buf,256,fp);  fgets(buf,256,fp);
  sprintf(obj->fname_bsb,"%s",fname);
  fscanf(fp,"%d\n",&nn);  fgets(buf,256,fp);
  if(nn==0) return 0;
  
  free(obj->bd.bsb);
  obj->n_bsb=nn;     obj->bd.bsb=(Bsb *)malloc(sizeof(Bsb)*nn);
 
  for(i=0;i<nn;i++){
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].lambda0=tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].ni     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].d_angle=tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].power  =tmpd;         
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.bsb[i].e0x    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.bsb[i].e0y    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].fx     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].fy     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].fz     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].theta  =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.bsb[i].phi    =tmpd;         
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
  fgets(buf,256,fp);  fgets(buf,256,fp);
  sprintf(obj->fname_blg,"%s",fname);
  fscanf(fp,"%d\n",&nn);   fgets(buf,256,fp);
  if(nn==0) return 0;
  
  free(obj->bd.blg);
  obj->n_blg=nn;  obj->bd.blg=(BsLGb *)malloc(sizeof(BsLGb)*nn);

  for(i=0;i<nn;i++){
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].lambda0=tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].ni     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].d_angle=tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].power  =tmpd;         
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.blg[i].e0x    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.blg[i].e0y    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].fx     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].fy     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].fz     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].theta  =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.blg[i].phi    =tmpd;         
    fscanf(fp,"%d" ,&tmpi);  obj->bd.blg[i].lg_m   =tmpi;         
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
  fgets(buf,256,fp);  fgets(buf,256,fp);
  sprintf(obj->fname_rab,"%s",fname);
  fscanf(fp,"%d\n",&nn);   fgets(buf,256,fp);
  if(nn==0) return 0;
  
  free(obj->bd.rab);
  obj->n_rab=nn;  obj->bd.rab=(RAb *)malloc(sizeof(RAb)*nn);

  for(i=0;i<nn;i++){
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].lambda0=tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].ni     =tmpd;          
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].NA     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].power  =tmpd;         
    fscanf(fp,"%lf",&tmpd); 
    fscanf(fp,"%lf",&tmpd2); obj->bd.rab[i].e0r    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);
    fscanf(fp,"%lf",&tmpd2); obj->bd.rab[i].e0a    =tmpd+I*tmpd2; 
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].fx     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].fy     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].fz     =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].theta  =tmpd;         
    fscanf(fp,"%lf",&tmpd);  obj->bd.rab[i].phi    =tmpd;         
    fscanf(fp,"%d\n",&tmpi); obj->bd.rab[i].nn     =tmpi;         
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
  obj->omega=2.0*M_PI*c0/obj->lambda_0;
}
