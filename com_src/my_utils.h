#ifndef MY_UTILS_H_
#define MY_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void vadd_d(double *r,double *v0,double *v1);
// calc vector addition. r=v0+v1
void vsub_d(double *r,double *v0,double *v1);
// calc vector subtraction. r=v0-v1
int vuni_d(double *r);
// calc unit vector of r. r=r/|r|. returned value 0 is normal termination, returned value -1 is abnormal termination ( |r|=0 ).


double vabs_d(double *v);
// return absolute value of double precision 3d-vector
double vabs_2dm(double *v0,double *v1);
// return absolute value of double precision 3d-vector with minus operation. |v0-v1|
double vdot_d(double *v0,double *v1);
// return dot product of double precision 3d-vector. v0 dot v1
void vcrs_d(double *r,double *v0,double *v1);
// calc cross product . r=v0 times v1

void prt_z_en(double complex z);
// print function for double complex data, data format type 'e', termination '\n'
void prt_z_fen(char *txt,int n,double complex z);
// print function for double complex data, txt:name, n:matrix index(n<0 for non matrix data), z:double complex data


void *m_alloc2(size_t num,size_t size, char *txt);
// memory allocation function, num:size of unit data, size:number of data, *txt:error message when failed.


void continue_message();
// print continue message

#endif 
