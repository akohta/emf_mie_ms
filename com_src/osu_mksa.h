#if !defined OSU_MKSA_H
#define OSU_MKSA_H

#include <math.h>
#include <complex.h>
#include "physical_constant.h"

// OSU is the system of units optimized for optics 
// defined as c0 ( speed of light in vacuum ) = 1, mu0 ( permeability of vacuum ) = 1  

// definition for conversion to MKSA system of units
#define UD_UnitL 1.0e-6 // unit length [m](MKSA), UD_UnitL[m](MKSA) = 1 [L](OSU)
#define UD_UnitP 1.0e-3 // unit power  [W](MKSA), UD_UnitP[W](MKSA) = 1 [L^2 M T^(-3)](OSU)
/* ---- conversion relation ----
  a = UD_UnitL,  b = UD_UnitP

  a [m] = 1 [L]
  a*b/(c0^3) [kg] = 1 [M]
  a/c0 [s] = 1 [T]
  sqrt(b/(c0*mu0)) [A] = 1 [I]
*/

// ------OSU to MKSA conversion------
double OSUtoMKSA_length(double x);
double OSUtoMKSA_weight(double x);
double OSUtoMKSA_time(double x);
double OSUtoMKSA_current(double x);
double OSUtoMKSA_power(double x);
double complex OSUtoMKSA_ElectricField(double complex x);
double complex OSUtoMKSA_MagneticField(double complex x);
double OSUtoMKSA_Force(double x);
double OSUtoMKSA_Torque(double x);
double OSUtoMKSA_Power_per_unit_area(double x);
double complex OSUtoMKSA_dEdv(double complex x);
double complex OSUtoMKSA_dHdv(double complex x);

// ------MKSA to OSU conversion------
double MKSAtoOSU_length(double x);
double MKSAtoOSU_weight(double x);
double MKSAtoOSU_time(double x);
double MKSAtoOSU_current(double x);
double MKSAtoOSU_power(double x);
double complex MKSAtoOSU_ElectricField(double complex x);
double complex MKSAtoOSU_MagneticField(double complex x);
double MKSAtoOSU_Force(double x);
double MKSAtoOSU_Torque(double x);
double MKSAtoOSU_Power_per_unit_area(double x);
double complex MKSAtoOSU_dEdv(double complex x);
double complex MKSAtoOSU_dHdv(double complex x);

#endif
