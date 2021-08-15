#include "my_utils.h"

void vadd_d(double *r,double *v0,double *v1)
{
  int i;
  for(i=0;i<3;i++) r[i]=v0[i]+v1[i];
}

void vsub_d(double *r,double *v0,double *v1)
{
  int i;
  for(i=0;i<3;i++) r[i]=v0[i]-v1[i];
}

int vuni_d(double *r)
{
  double A=vabs_d(r);
  if(A==0.0) return -1;
  else {
    r[0]/=A;    r[1]/=A;    r[2]/=A;
    return 0;
  }
}



double vabs_d(double *v)
{
  return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double vabs_2dm(double *v0,double *v1)
{
  return sqrt((v0[0]-v1[0])*(v0[0]-v1[0])+(v0[1]-v1[1])*(v0[1]-v1[1])+(v0[2]-v1[2])*(v0[2]-v1[2]));
}

double vdot_d(double *v0,double *v1)
{
  return v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2];
}

void vcrs_d(double *r,double *v0,double *v1)
{
  r[0]=v0[1]*v1[2]-v0[2]*v1[1];
  r[1]=v0[2]*v1[0]-v0[0]*v1[2];
  r[2]=v0[0]*v1[1]-v0[1]*v1[0];
}


void prt_z_en(double complex z)
{
  printf("% 15.14e %+15.14eI\n",creal(z),cimag(z));
}

void prt_z_fen(char *txt,int n,double complex z)
{
  if(n>=0){
    printf("%s[%3d] = % 15.14e %+15.14eI\n",txt,n,creal(z),cimag(z));
  }
  else {
    printf("%s = % 15.14e %+15.14eI\n",txt,creal(z),cimag(z));
  }
}



void *m_alloc2(size_t num,size_t size, char *txt)
{
  void *tmp;
  tmp=calloc(num,size);
  if(tmp==NULL){
    printf("memory allocation error!\n");
    printf("%s. ",txt);
    printf("Exit..\n");
    exit(1);
  }
  else return tmp;
}

void continue_message()
{
  printf("continue? (y/n) : ");  if(getchar()!='y'){ printf("Exit\n");  exit(0);}
}
