#include "emf_mie_ms.h"

int main(int argc, char *argv[])
{
  if(argc!=2){
    printf("This program needs command line argument.\n");
    printf("Usage : %s 'output datafile name'\n",argv[0]);
    exit(0);
  }
	
  MSPD msp;
  
  read_data_ms(&msp);
  print_data_ms(&msp);
  setup_ms(&msp);
  output_node_particles(argv[1],&msp);  
  
  iterative_ops_ms(&msp);
  
  write_dat_ms(argv[1],&msp);
  free_ms(&msp);

  return 0;
}
