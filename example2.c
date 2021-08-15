// calculation example of electric field intensity distributions
#include "emf_mie_ms.h"

int main(int argc,char *argv[]) 
{
  MSPD msp;
  FILE *fp1,*fp2;
  double complex e[3],h[3];
  double rang,dr,r[3];
  int max,i,j;

  read_dat_ms(argv[1],&msp); // read data file
  print_data_ms(&msp); // print data
  
  max=150;
  rang=4.0*msp.bm.lambda_0;
  dr=rang*2/(double)(max-1);
  
  // x=0 plane
  if((fp1=fopen("Ie_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# y z electric_field_intensity");
  if((fp2=fopen("Ih_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp2,"%s\n","# y z magnetic_field_intensity");
  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      total_EH_ms(e,h,r,&msp); // total field 
      fprintf(fp1,"%g %g %15.14e\n",r[1],r[2],creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2])));
      fprintf(fp2,"%g %g %15.14e\n",r[1],r[2],creal(h[0]*conj(h[0]))+creal(h[1]*conj(h[1]))+creal(h[2]*conj(h[2])));
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);
 
 
  // y=0 plane
  if((fp1=fopen("Ie_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# x z electric_field_intensity");
  if((fp2=fopen("Ih_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp2,"%s\n","# x z magnetic_field_intensity");
  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      total_EH_ms(e,h,r,&msp);
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[2],creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2])));
      fprintf(fp2,"%g %g %15.14e\n",r[0],r[2],creal(h[0]*conj(h[0]))+creal(h[1]*conj(h[1]))+creal(h[2]*conj(h[2])));
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);
  

  // z=0 plane  
  if((fp1=fopen("Ie_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# x y electric_field_intensity");
  if((fp2=fopen("Ih_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp2,"%s\n","# x y magnetic_field_intensity");
  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      total_EH_ms(e,h,r,&msp);
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[1],creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2])));
      fprintf(fp2,"%g %g %15.14e\n",r[0],r[1],creal(h[0]*conj(h[0]))+creal(h[1]*conj(h[1]))+creal(h[2]*conj(h[2])));
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
  }
  fclose(fp1);
  fclose(fp2);
  
  printf("Intensity plot is finished\n");
  
  free_ms(&msp);
  return 0;
}
