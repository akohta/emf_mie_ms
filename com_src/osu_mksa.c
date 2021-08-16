#include "osu_mksa.h"

// OSU to MKSA conversion
double OSUtoMKSA_length(double x)
{
  return x*UD_UnitL;
}

double OSUtoMKSA_weight(double x)
{
  return x*UD_UnitL*UD_UnitP/(C_C0*C_C0*C_C0);
}

double OSUtoMKSA_time(double x)
{
  return x*UD_UnitL/C_C0;
}

double OSUtoMKSA_current(double x)
{
  return x*sqrt(UD_UnitP/(C_C0*C_Mu0));
}

double OSUtoMKSA_power(double x)
{
  return x*UD_UnitP;
}

double complex OSUtoMKSA_ElectricField(double complex x)
{
  return x*sqrt(UD_UnitP*C_C0*C_Mu0)/UD_UnitL;
}

double complex OSUtoMKSA_MagneticField(double complex x)
{
  return x*sqrt(UD_UnitP/(C_C0*C_Mu0))/UD_UnitL;
}

double OSUtoMKSA_Force(double x)
{
  return x*UD_UnitP/C_C0;
}

double OSUtoMKSA_Torque(double x)
{
  return x*UD_UnitL*UD_UnitP/C_C0;
}

double OSUtoMKSA_Power_per_unit_area(double x)
{
  return x*UD_UnitP/(UD_UnitL*UD_UnitL);
}

double complex OSUtoMKSA_dEdv(double complex x)
{
  return x*sqrt(UD_UnitP*C_C0*C_Mu0)/(UD_UnitL*UD_UnitL);
}

double complex OSUtoMKSA_dHdv(double complex x)
{
  return x*sqrt(UD_UnitP/(C_C0*C_Mu0))/(UD_UnitL*UD_UnitL);
}

// MKSA to OSU conversion
double MKSAtoOSU_length(double x)
{
  return x/UD_UnitL;
}

double MKSAtoOSU_weight(double x)
{
  return x*C_C0*C_C0*C_C0/(UD_UnitL*UD_UnitP);
}

double MKSAtoOSU_time(double x)
{
  return x*C_C0/UD_UnitL;
}

double MKSAtoOSU_current(double x)
{
  return x*sqrt(C_C0*C_Mu0/UD_UnitP);  
}

double MKSAtoOSU_power(double x)
{
  return x/UD_UnitP;
}

double complex MKSAtoOSU_ElectricField(double complex x)
{
  return x*UD_UnitL/sqrt(UD_UnitP*C_C0*C_Mu0);
}

double complex MKSAtoOSU_MagneticField(double complex x)
{
  return x*UD_UnitL*sqrt(C_C0*C_Mu0/UD_UnitP);
}

double MKSAtoOSU_Force(double x)
{
  return x*C_C0/UD_UnitP;
}

double MKSAtoOSU_Torque(double x)
{
  return x*C_C0/(UD_UnitL*UD_UnitP);
}

double MKSAtoOSU_Power_per_unit_area(double x)
{
  return x*UD_UnitL*UD_UnitL/UD_UnitP;
}

double complex MKSAtoOSU_dEdv(double complex x)
{
  return x*UD_UnitL*UD_UnitL/sqrt(UD_UnitP*C_C0*C_Mu0);
}

double complex MKSAtoOSU_dHdv(double complex x)
{
  return x*UD_UnitL*UD_UnitL*sqrt(C_C0*C_Mu0/UD_UnitP);
}
