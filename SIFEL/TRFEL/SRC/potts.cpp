/*
    File:            potts.cpp
    Author:          Tomas Krejci, 13/04/2018
    Purpose:         retention curve accorging to Potts model
    sources:                              
                     D. M. Potts and L. Zdravkovic, Finite element analysis in geotechnical engineering theory, 
		     Thomas Telford, London, 1999;
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "potts.h"

potts_reten::potts_reten()
{ 
  sirr = 0.0;
  ssat = 0.0;
  hp_min = 0.0;
  htz = 0.14;
  r = 100.0;
  gammaw = 10000.0;
  ksat = 0.1;
}

potts_reten::~potts_reten()
{}


/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void potts_reten::read(XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf", &ssat, &sirr, &hp_min, &htz, &r, &ksat);
}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void potts_reten::print(FILE *out)
{
  fprintf (out,"\n %lf %lf %lf %lf %lf %lf \n", ssat, sirr, hp_min, htz, r, ksat);
}

/**
   function computes saturation degree

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 3/04/2018
*/
double potts_reten::sw(double pw)
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
double potts_reten::dsw_dpw(double pw)
{
  double dsw_dpw,dkr_dpw,hp;
  
  hp = pw/gammaw;
  dkr_dpw = pow(10,(hp - hp_min)*log10(r)/htz)*log(10)*log10(r)/gammaw/htz;
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
double potts_reten::dsw_dt(double /*pw*/)
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
double potts_reten::krw(double pw)
{
  double krw,hp;
  
  hp = pw/gammaw;
  krw = pow(10,(hp - hp_min)*log10(r)/htz);

  return(krw);
}


/**
   function computes water conductivity of porous material 

   @retval k - conductivity

   13/04/2018, TKr
*/
double potts_reten::get_k(double /*pw*/)
{
  double k;

  k = ksat/86400.0/gammaw;

  return(k);
}
