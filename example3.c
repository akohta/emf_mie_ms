#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "emf_mie_ms.h"
#include "my_utils.h"
#include <sys/stat.h>
#include <errno.h>  
#include <png.h>

typedef struct image_data{
  char dir_name[64];      // directory name to output image
  int scale;              // number for enlarge the output image
  
  int m;                  // sampling number 
  double rang;            // range of sampling

  int ts;                 // time step per cycle
  
  double complex *ve,*vh; // electromagnetic field data
  double me[3],mh[3];     // maximum amplitude of each field component
}IMD;

void directory_name(char *src,char *nn);
void make_directory(char *dir_name);

void eh_field_x(IMD *id,MSPD *sp);
void eh_field_y(IMD *id,MSPD *sp);
void eh_field_z(IMD *id,MSPD *sp);
void output_field(char *pl,IMD *id,MSPD *sp);

// color table
png_byte ct1[9][3]={{0x00,0x00,0x90},{0x00,0x0f,0xff},{0x00,0x90,0xff},{0x0f,0xff,0xee},
                    {0xff,0xff,0xff},{0xff,0xee,0x00},{0xff,0x70,0x00},{0xee,0x00,0x00},{0x7f,0x00,0x00}};
/*                    
png_byte ct1[9][3]={{0x00,0x00,0x90},{0x00,0x0f,0xff},{0x00,0x90,0xff},{0x0f,0xff,0xee},
                    {0x90,0xff,0x70},{0xff,0xee,0x00},{0xff,0x70,0x00},{0xee,0x00,0x00},{0x7f,0x00,0x00}};  
*/

int main(int argc,char *argv[])
{
  MSPD msp;
  IMD id;
  double mf;
  int sn;
  
  if(argc!=2 && argc!=4){
    printf("Usage : %s datafile_name [sampling_number multplier_factor](optional)\n",argv[0]);
    printf("default sampling number 200, multiplier factor 4 (range is -4*lambda0 to 4*lambda0)\n");
    exit(0);
  }
  else if(argc==4){
    sn=atoi(argv[2]);
    mf=atof(argv[3]);
  }
  else{
    sn=200;
    mf=4.0;
  }

  read_dat_ms(argv[1],&msp); // read data file
  print_data_ms(&msp);       // print data
  
  directory_name(argv[1],id.dir_name); // remove file-extension from argv[1] and add "_images"
  id.scale=1;                          // number for enlarge the output image
  id.m=sn;                             // sampling number 
  id.rang=mf*msp.bm.lambda_0;          // range of sampling
  id.ts=40;                            // time step per cycle
  
  make_directory(id.dir_name);

  id.ve=(double complex *)m_alloc2(id.m*id.m*3,sizeof(double complex),"example3.c, ve");
  id.vh=(double complex *)m_alloc2(id.m*id.m*3,sizeof(double complex),"example3.c, vh");

  // x=0 plane
  eh_field_x(&id,&msp);
  output_field("yz",&id,&msp);
  // y=0 plane
  eh_field_y(&id,&msp);
  output_field("xz",&id,&msp);
  // z=0 plane
  eh_field_z(&id,&msp);
  output_field("xy",&id,&msp);

  free(id.ve);
  free(id.vh);
  
  free_ms(&msp);
  return 0;
}

void directory_name(char *src,char *nn)
{
  int s1,s2;
  char *sd,fo[64]={},buf[54]={};
  
  s1=strlen(src);
  if(s1>54){
    printf("example3.c, directory_name(), directory name is too long. exit...\n");
    exit(1);
  }
  sprintf(fo,"%s",src);
  sd=strrchr(fo,'.');
  if(sd!=NULL){
    s2=strlen(sd);
    strncpy(buf,src,s1-s2);
    sprintf(fo,"%s_images",buf);
  }
  sprintf(nn,"%s",fo);
}

void make_directory(char *dir_name)
{
  int ret;
  
  ret=mkdir(dir_name,S_IRWXU|S_IRWXG);
  if(ret!=0 && errno!=EEXIST){
    printf("failed to make directory. Exit..");
    exit(1);
  }
}

void eh_field_x(IMD *id,MSPD *sp)
{
  double complex e[3],h[3];
  double x[3],dr;
  int i,j,d;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  for(i=0;i<3;i++){
    id->me[i]=0.0;
    id->mh[i]=0.0;
  }
  
  // x=0 plane  
  x[0]=0.0;
  #pragma omp parallel for schedule(dynamic) firstprivate(x) private(j,d,e,h) 
  for(i=0;i<id->m;i++){
    x[2]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[1]=-id->rang+(double)j*dr;
      total_EH_ms(e,h,x,sp); // total field
      
      #pragma omp critical
      for(d=0;d<3;d++){
        if(cabs(e[d])>id->me[d]) id->me[d]=cabs(e[d]);
        if(cabs(h[d])>id->mh[d]) id->mh[d]=cabs(h[d]);
      }
      
      for(d=0;d<3;d++){
        id->ve[i*id->m*3+j*3+d]=e[d];
        id->vh[i*id->m*3+j*3+d]=h[d];
      }
    }
  }
}

void eh_field_y(IMD *id,MSPD *sp)
{
  double complex e[3],h[3];
  double x[3],dr;
  int i,j,d;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  for(i=0;i<3;i++){
    id->me[i]=0.0;
    id->mh[i]=0.0;
  }
  
  // y=0 plane  
  x[1]=0.0;
  #pragma omp parallel for schedule(dynamic) firstprivate(x) private(j,d,e,h) 
  for(i=0;i<id->m;i++){
    x[2]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[0]=-id->rang+(double)j*dr;
      total_EH_ms(e,h,x,sp); // total field
      
      #pragma omp critical
      for(d=0;d<3;d++){
        if(cabs(e[d])>id->me[d]) id->me[d]=cabs(e[d]);
        if(cabs(h[d])>id->mh[d]) id->mh[d]=cabs(h[d]);
      }
      
      for(d=0;d<3;d++){
        id->ve[i*id->m*3+j*3+d]=e[d];
        id->vh[i*id->m*3+j*3+d]=h[d];
      }
    }
  }
}

void eh_field_z(IMD *id,MSPD *sp)
{
  double complex e[3],h[3];
  double x[3],dr;
  int i,j,d;
  
  dr=id->rang*2.0/(double)(id->m-1);
  
  for(i=0;i<3;i++){
    id->me[i]=0.0;
    id->mh[i]=0.0;
  }
  
  // z=0 plane  
  x[2]=0.0;
  #pragma omp parallel for schedule(dynamic) firstprivate(x) private(j,d,e,h) 
  for(i=0;i<id->m;i++){
    x[1]=id->rang-(double)i*dr;
    for(j=0;j<id->m;j++){
      x[0]=-id->rang+(double)j*dr;
      total_EH_ms(e,h,x,sp); // total field
      
      #pragma omp critical
      for(d=0;d<3;d++){
        if(cabs(e[d])>id->me[d]) id->me[d]=cabs(e[d]);
        if(cabs(h[d])>id->mh[d]) id->mh[d]=cabs(h[d]);
      }
      
      for(d=0;d<3;d++){
        id->ve[i*id->m*3+j*3+d]=e[d];
        id->vh[i*id->m*3+j*3+d]=h[d];
      }
    }
  }
}

void output_field(char *pl,IMD *id,MSPD *sp)
{
  void output_png(int nt,double complex cet,char *pl,IMD *id);
  void output_color_bar(IMD *id);
  
  FILE *fp;
  char fn[128];
  double dt;
  int n;
  
  dt=sp->bm.lambda_0/(double)id->ts;
  
  #pragma omp parallel for schedule(dynamic) 
  for(n=0;n<id->ts;n++){
    output_png(n,cexp(-I*sp->bm.omega*dt*(double)n),pl,id);
  }

  // print info
  sprintf(fn,"%s/%s_info.txt",id->dir_name,pl);
  fp=fopen(fn,"wt");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fn);
    exit(1);
  }
  fprintf(fp,"the range of color bar\n");
  fprintf(fp,"Ex is %8e to %8e\n",-id->me[0],id->me[0]);
  fprintf(fp,"Ey is %8e to %8e\n",-id->me[1],id->me[1]);
  fprintf(fp,"Ez is %8e to %8e\n",-id->me[2],id->me[2]);
  fprintf(fp,"Hx is %8e to %8e\n",-id->mh[0],id->mh[0]);
  fprintf(fp,"Hy is %8e to %8e\n",-id->mh[1],id->mh[1]);
  fprintf(fp,"Hz is %8e to %8e\n",-id->mh[2],id->mh[2]);
  fclose(fp);
  
  // output color bar image
  output_color_bar(id);
}

void output_png(int nt,double complex cet,char *pl,IMD *id)
{
  int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b); // -1 <= x <= 1
  
  FILE *fep[3],*fhp[3];
  char fname[256],sf[3]={'x','y','z'};
  int j,i,sj,si,d,m,scale;
  png_uint_32 width,height;
  png_structp png_e[3],png_h[3];
  png_infop info_e[3],info_h[3];
  png_bytepp pd_e[3],pd_h[3];
  png_byte r,g,b;

  m=id->m;
  scale=id->scale;

  width =m*(scale+1);
  height=m*(scale+1);

  for(d=0;d<3;d++){
    png_e[d] =png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_e[d]=png_create_info_struct(png_e[d]);
    sprintf(fname,"%s/%s_E%c_%03d.png",id->dir_name,pl,sf[d],nt);
    fep[d]=fopen(fname,"wb");
    if(fep[d]==NULL){
      printf("Failed to open the %s file. Exit...\n",fname);
      exit(1);
    }
    
    png_h[d] =png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_h[d]=png_create_info_struct(png_h[d]);
    sprintf(fname,"%s/%s_H%c_%03d.png",id->dir_name,pl,sf[d],nt);
    fhp[d]=fopen(fname,"wb");
    if(fhp[d]==NULL){
      printf("Failed to open the %s file. Exit...\n",fname);
      exit(1);
    }
  
    png_init_io(png_e[d],fep[d]);
    png_set_IHDR(png_e[d],info_e[d],width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
    pd_e[d]=(png_bytepp)png_malloc(png_e[d],sizeof(png_bytep)*height);
    png_set_rows(png_e[d],info_e[d],pd_e[d]);
    
    png_init_io(png_h[d],fhp[d]);
    png_set_IHDR(png_h[d],info_h[d],width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
    pd_h[d]=(png_bytepp)png_malloc(png_h[d],sizeof(png_bytep)*height);
    png_set_rows(png_h[d],info_h[d],pd_h[d]);

    for(j=0;j<height;j++){
      pd_e[d][j]=(png_bytep)png_malloc(png_e[d],sizeof(png_byte)*width*3);
      pd_h[d][j]=(png_bytep)png_malloc(png_h[d],sizeof(png_byte)*width*3);
    }
  }
  
  for(i=0;i<m;i++){
    for(j=0;j<m;j++){
      
      for(d=0;d<3;d++){
        color_rgb(creal(cet*id->ve[i*m*3+j*3+d])/id->me[d],&r,&g,&b);
        for(si=0;si<=scale;si++){
          for(sj=0;sj<=scale;sj++){
            pd_e[d][i*(scale+1)+si][(j*(scale+1)+sj)*3+0]=r;
            pd_e[d][i*(scale+1)+si][(j*(scale+1)+sj)*3+1]=g;
            pd_e[d][i*(scale+1)+si][(j*(scale+1)+sj)*3+2]=b;
          }
        }
        
        color_rgb(creal(cet*id->vh[i*m*3+j*3+d])/id->mh[d],&r,&g,&b);
        for(si=0;si<=scale;si++){
          for(sj=0;sj<=scale;sj++){
            pd_h[d][i*(scale+1)+si][(j*(scale+1)+sj)*3+0]=r;
            pd_h[d][i*(scale+1)+si][(j*(scale+1)+sj)*3+1]=g;
            pd_h[d][i*(scale+1)+si][(j*(scale+1)+sj)*3+2]=b;
          }
        }      
      }
    }
  }
  
  for(d=0;d<3;d++){
    png_write_png(png_e[d],info_e[d],PNG_TRANSFORM_IDENTITY,NULL);
    png_write_png(png_h[d],info_h[d],PNG_TRANSFORM_IDENTITY,NULL);
    
    for(j=0;j<height;j++){
      png_free(png_e[d],pd_e[d][j]);
      png_free(png_h[d],pd_h[d][j]);
    }
    png_free(png_e[d],pd_e[d]);
    png_free(png_h[d],pd_h[d]);
    
    fclose(fep[d]);
    fclose(fhp[d]);
  }
}

void output_color_bar(IMD *id)
{
  int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b); // -1 <= x <= 1
  
  FILE *fp;
  char fname[128];
  int j,i;
  
  png_uint_32 width,height;
  png_structp png;
  png_infop info;
  png_bytepp pdata;
  png_byte r,g,b;

  sprintf(fname,"%s/color_bar.png",id->dir_name);

  height=id->m*(id->scale+1);
  width=height/16;
  
  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info= png_create_info_struct(png);
  
  fp=fopen(fname,"wb");
  if(fp==NULL){
    printf("Failed to open the %s file. Exit...\n",fname);
    exit(1);
  }
  
  png_init_io(png, fp);
  png_set_IHDR(png,info,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  pdata=(png_bytepp)png_malloc(png, sizeof(png_bytep)*height);
  png_set_rows(png,info,pdata);

  for(j=0;j<height;j++){
    pdata[j]=(png_bytep)png_malloc(png,sizeof(png_byte)*width*3);
  }
  
  for(i=0;i<height;i++){
    color_rgb(1.0-(2.0/(double)height)*(double)i,&r,&g,&b);
    for(j=0;j<width;j++){
      pdata[i][j*3+0]=r;
      pdata[i][j*3+1]=g;
      pdata[i][j*3+2]=b;
    }
  }
  
  png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
  
  for(j=0;j<height;j++){
    png_free(png,pdata[j]);
  }
  png_free(png,pdata);
  fclose(fp);
}

int color_rgb(double x,png_byte *r,png_byte *g,png_byte *b) // -1 <= x <= 1
{
  double i_nc,dr,dg,db;
  unsigned int i,n,nc,nd;

  if(x<-1.0 || x>1.0){
    *r=0x00;    *g=0x00;    *b=0x00;
    return -1;
  }
  
  n=(unsigned int)floor(pow(2,23)*(x+1.0));
  nc=(unsigned int)pow(2,21);
  i_nc=1.0/(double)nc;
  
  if(n<nc*1)      i=1;
  else if(n<nc*2) i=2;
  else if(n<nc*3) i=3;
  else if(n<nc*4) i=4;
  else if(n<nc*5) i=5;
  else if(n<nc*6) i=6;
  else if(n<nc*7) i=7;
  else if(n<nc*8) i=8;
  else {
    *r=ct1[8][0];    *g=ct1[8][1];    *b=ct1[8][2];
    return 0;
  }
    
  nd=n-nc*(i-1);
  dr=(double)(ct1[i][0]-ct1[i-1][0])*i_nc;
  dg=(double)(ct1[i][1]-ct1[i-1][1])*i_nc;
  db=(double)(ct1[i][2]-ct1[i-1][2])*i_nc;
  *r=(png_byte)floor((double)ct1[i-1][0]+dr*(double)nd);
  *g=(png_byte)floor((double)ct1[i-1][1]+dg*(double)nd);
  *b=(png_byte)floor((double)ct1[i-1][2]+db*(double)nd);
  
  return 0;  
}
