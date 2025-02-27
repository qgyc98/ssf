/*
    File:             sejtkrmatc.cpp
    Author:           Tomas Krejci, 12/9/2008
    Purpose:          material model for saturated one-phase flow in deforming medium
    sources:          Lewis and Schrefler pp. 75-82
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "sejtkrmatc.h"
#include "globalt.h"

sejtkrmatc::sejtkrmatc()
{
  //
  emod = 10.0e6;//Young's modulus [Pa]
  //
  nu = 0.4;//Poisson's ratio [-]
  //
  alpha = 0.9998;//Biot's constant [-] alpha = 1 - kt/ks
  //
  ks = 18.0e12;//bulk modulus of solid phase (grains) [Pa] - sejnoha to znaci jako Kz
  //
  kt = 3.6e6; //bulk modulus of porous skeleton [Pa] - sejnoha to znaci jako Ks
  //
  kw = 2.0e9;//bulk modulus of water [Pa]
  //
  //ka = 0.1e6;//bulk modulus of air [Pa]
  //
  phi0 = 0.5;//initial porosity [-]
  //
  kintr = 4.5e-13;//intrinsic permeability [m^2]
  //
  rhos = 2000.0;//solid grain density [kg/m^3]
  //
  rhow = 1000.0;//water density [kg/m^3]
  //
  muw = 1.0e-3;//water viscosity [Pa.s]
  //
  //mua = 1.8e-5;//air viscosity [Pa.s]
}

sejtkrmatc::~sejtkrmatc()
{}


/**
   function reads parameters
   
   @param in - input file

   12/9/2008, TKr
*/
void sejtkrmatc::read(XFILE *in)
{
  //xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &emod, &nu, &rhos, &alpha, &ks, &phi0, &kw, &rhow, &muw, &kintr);
}



/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kuw - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kuw(double pw)
{
  double kuw;
  
  kuw = -alpha;

  return(kuw);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kwu - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kwu(double pw)
{
  double kwu;

  kwu = 0.0;

  return(kwu);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kww - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kww(double pw)
{
  double kww;
  
  kww = kintr/muw;

  return(kww);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capuw - capacity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_capuw(double pw)
{
  double capuw;
  
  capuw = 0.0;

  return(capuw);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capwu - capacity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_capwu(double pw)
{
  double capwu;
  
  capwu = alpha;

  return(capwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capww - capacity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_capww(double pw)
{
  double n,capww;  

  n = phi;//porosity
  
  capww = (alpha-n)/ks+n/kw;

  return(capww);

}


/**
   function returns ... 
   @param pw - water pressure

   @retval fw - first part for right-hand side for continutiy equation

   12/9/2008, TKr
*/
double sejtkrmatc::get_fw1(double pw)
{
  double fw1;
  
  fw1 = kintr/muw*rhow;//doplnit!!

  return(fw1);

}



/**
   function returns ... 
   @param pw - water pressure

   @retval fu1 - first part for right-hand side for balance equation

   12/9/2008, TKr
*/
double sejtkrmatc::get_fu1(double pw)
{
  double fu1;  
  
  fu1 = rhos*(n-1) + n*rhow;//doplnit!!

  return(fu1);

}
