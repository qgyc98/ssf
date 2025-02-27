/*
    File:            gardner.cpp
    Author:          Tomas Krejci, 13/04/2018
    Purpose:         retention curve accorging to Gardner model
    sources:                              
                     W. R. Gardner, Some steady-state solutions of the unsaturated moisture flow equation 
                     to evaporation from a water table, Soil Science 85 (1958), no. 4, 228-232;
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gardner.h"

gardner_reten::gardner_reten()
{
  ssat = 1.0;
  sirr = 0.33;
  alfa = 0.1;
  ksat = 0.1;
  gamaw = 10000.0;
}

gardner_reten::~gardner_reten()
{}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void gardner_reten::read(XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf", &ssat, &sirr, &alfa, &ksat);
}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void gardner_reten::print(FILE *out)
{
  fprintf (out,"\n%lf %lf %lf %lf\n", ssat, sirr, alfa, ksat);
}

/**
   function computes saturation degree

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 3/04/2018
*/
double gardner_reten::sw(double pw)
{
  double sw,kr;
  
  kr = krw(pw);
  sw = sirr + (ssat - sirr)*kr;//Warick
  //sw = sirr + (ssat - sirr)*pow(kr,1.0/3.0);//alternative expression according to Fatt and Klikof
  
  return(sw);
}


/**
   function computes partial derivative of degree of saturation with respect to pw, specific water content

   @param pw - water pressure
   @retval dsw_dpw - partial derivative of degree of saturation with respect to pw

   TKr 3/04/2018
*/
double gardner_reten::dsw_dpw(double pw)
{
  double dsw_dpw,dkr_dpw;
  
  dkr_dpw = exp(alfa*pw/gamaw)*alfa/gamaw;
  dsw_dpw = (ssat - sirr)*dkr_dpw;//Warick
  //dsw_dpw = (ssat - sirr)*alfa/3.0*pow(dkr_dpw,-2.0/3.0);//alternative expression according to Fatt and Klikof
  
  return(dsw_dpw);
}


/**
   function computes partial derivative of degree of saturation with respect to t

   @param pw - water pressure
   @retval dsw_dt - partial derivative of degree of saturation with respect to t

   TKr 3/04/2018
*/
double gardner_reten::dsw_dt(double /*pw*/)
{
  double dsw_dt;
  
  dsw_dt = 0.0;
  
  return(dsw_dt);
}


/**
   function computes water relative permeability

   @param pw - water pressure
   @retval krw - water relative permeability

   TKr 3/04/2018
*/
double gardner_reten::krw(double pw)
{
  double krw;

  krw = exp(0.1*pw/gamaw);

  return(krw);
}

/**
   function computes water conductivity of porous material 

   @retval k - conductivity

   13/04/2018, TKr
*/
double gardner_reten::get_k(double /*pw*/)
{
  double k;

  k = ksat/86400.0/gamaw;

  return(k);
}
