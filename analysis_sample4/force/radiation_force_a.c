// calculation example of radiation force.
// changing the sphere radius.
#include "emf_mie_ms.h"

void set_coef(MSPD *msp);

int main(int argc,char *argv[])
{
  if(argc!=5){
    printf("Usage : %s minimum_radius radius_increase_by_step number_of_steps output_filename\n",argv[0]);
    exit(0);
  }
    
  MSPD msp;
  FILE *ff;
  double a,da,vf[3],vn[3];
  int max,i;
  
  read_data_ms(&msp);
  print_data_ms(&msp);
  setup_ms(&msp);
  
  a=atof(argv[1]);
  da=atof(argv[2]);
  max=atoi(argv[3]);
  printf("-- settings --\n");
  printf("minimum radius          = %g\n",a);
  printf("radius increase by step = %g\n",da);
  printf("number of steps         = %d\n",max);
  printf("output filename         = %s\n\n",argv[4]);
  
  if((ff=fopen(argv[4],"wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(ff,"%s\n","# a Fx[N] Fy[N] Fz[N] Nx[Nm] Ny[Nm] Nz[Nm]");

  for(i=0;i<max;i++){
    printf("\r%d/%d",i+1,max);  fflush(stdout);
    // set shpere radius
    msp.sp[0].a=a+(double)i*da;
    set_coef(&msp);
    // calculation of force and torque
    force_torque_ms(0,vf,vn,&msp);
    // print force and torqe in MKSA system of units [N]
    fprintf(ff,"%g %g %g %g %g %g %g\n",msp.sp[0].a,OSUtoMKSA_Force(vf[0]),OSUtoMKSA_Force(vf[1]),OSUtoMKSA_Force(vf[2]),OSUtoMKSA_Torque(vn[0]),OSUtoMKSA_Torque(vn[1]),OSUtoMKSA_Torque(vn[2]));
  }
  fclose(ff);
  printf("\nDone.\n"); fflush(stdout);
  
  free_ms(&msp);
  return 0;
}

void set_coef(MSPD *msp)
{
  void setup_coefficient_dispd(SPD *sp,Bobj *bm); // emf_mie_ms.c
  void initialize_eh_r(SPD *sp,Bobj *bm); // emf_mie_ms.c
  void coefficient(SPD *sp);              // emf_mie_ms.c

  setup_coefficient_dispd(&(msp->sp[0]),&(msp->bm));
  initialize_eh_r(&(msp->sp[0]),&(msp->bm));
  coefficient(&(msp->sp[0]));
}
  
