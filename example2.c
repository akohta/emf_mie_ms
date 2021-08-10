// calculation example of electric field intensity distributions
#include "emf_mie_ms.h"

int main(int argc,char *argv[]) 
{
  MSPD msp;
  FILE *fp;
  double complex e[3],h[3];
  double rang,dr,r[3];
  int max,i,j;

  read_dat_ms(argv[1],&msp); // read data file
  print_data_ms(&msp); // print data
  
  max=150;
  rang=4.0*msp.bm.lambda_0;
  dr=rang*2/(double)(max-1);
  
  // x=0 plane
  if((fp=fopen("Ie_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      total_EH_ms(e,h,r,&msp); // total field 
      fprintf(fp,"%g %g %15.14e\n",r[1],r[2],creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2])));
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
 
  // y=0 plane
  if((fp=fopen("Ie_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      total_EH_ms(e,h,r,&msp);
      fprintf(fp,"%g %g %15.14e\n",r[0],r[2],creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2])));
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  // z=0 plane  
  if((fp=fopen("Ie_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      total_EH_ms(e,h,r,&msp);
      fprintf(fp,"%g %g %15.14e\n",r[0],r[1],creal(e[0]*conj(e[0]))+creal(e[1]*conj(e[1]))+creal(e[2]*conj(e[2])));
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  
  printf("Intensity plot is finished\n");
  
  free_ms(&msp);
  return 0;
}
