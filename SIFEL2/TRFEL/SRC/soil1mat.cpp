/*
  File:             soil1mat.cpp
  Author:           Tomas Krejci, 14.3.2006
  Purpose:          


  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "soil1mat.h"
#include "globalt.h"

soil1mat::soil1mat()
{
  mw = 18.01528;   //molar mass of water kg.mol-1
  ma = 28.9645;    //molar mass of dry air kg.mol-1
  gasr = 8314.41;  //universal gas constant J.mol-1.K-1
  
  t0 = 273.15;   //reference temperature
  p0 = 101325.0; //pressure
  tcr = 647.3;   //critical point of water
}
soil1mat::~soil1mat()
{}


/**
   function computes diffusion coefficient of vapour inside pores
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dg - diffusion coefficient of vapour inside pores

   16.3.2006, Tkr
*/
double soil1mat::_dg()
{
  double dg;

  dg = 2.6e-5;
  
  return(dg);
}


/**
   function returns Biot's constant
   
   @retval alpha - Biot's constant

   16.3.2006, Tkr
*/
double soil1mat::_alpha()
{
  double alpha,ks,kt;

  ks = _ks();
  kt = _kt();

  alpha = 1.0 - kt/ks;


  return(alpha);
}



/**
   function returns compresibility coefficient

   @retval ks - compresibility coefficient

   16.3.2006, Tkr
*/
double soil1mat::_ks()
{
  return(0.14e10);
}



/**
   function returns compresibility coefficient of sceleton

   @retval kt - compresibility coefficient of sceleton

   16.3.2006, Tkr
*/
double soil1mat::_kt()
{
  double kt,emod,nu;

  emod = _emod();
  nu = _nu();

  kt = emod/(3.0*(1.0 - 2.0*nu));

  return(kt);
}



/**
   function computes emod Young's modulus

   @retval emod - Young's modulus
*/
double soil1mat::_emod()
{
  return(60000000.0);
}


/**
   function computes nu Poisson's constant

   @retval nu - Poisson's constant
*/
double soil1mat::_nu()
{
  return(0.2857);
}


/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval beats - cubic thermal expansion coefficient of solid (K-1)
*/
double soil1mat::_betas()
{
  return(0.9e-6);
}

/**
   function computes thermal capacity of partially saturated concrete
   
   @retval rhocp -  thermal capacity of partially saturated concrete
*/
double soil1mat::_rhocp()
{

  //1cal = 4,1868J
  //40 kcal.K^-1.m^-3 = 167472.0 J.K^-1.m^-3
  //funguje s 40000000.0
  return(40.0);
}

/**
   function computes average thermal conductivity

   @retval lambdas - average thermal conductivity
*/
double soil1mat::_lambdaa()
{
  //0.2 cal.K^-1.m^-1.s-1 = 0.83736 W.K^-1.m^-1
  //funguje s 0.2
  return(5.55555e-05);
}



/**
   function computes degree of saturation(desorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double soil1mat::sat(double pc,double /*t*/)
{
  double sw,pb,lam,sirr;
  
  pb = 133813.0;
  sirr = 0.3216;
  lam = 2.308;
  
  sw = (1.0 - sirr)*pow((pb/pc),lam) + sirr;
  
  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double soil1mat::dsat_dpc(double pc,double /*t*/)
{
  double dsw_dpc,pb,lam,sirr;
  
  pb = 133813.0;
  sirr = 0.3216;
  lam = 2.308;
  
  dsw_dpc = pow((pb/pc),lam)*(sirr - 1.0)*lam/pc;
  
  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double soil1mat::dsat_dt(double /*pc*/,double /*t*/)
{
  return(0.0);
}


/**
   function computes porosity
   
   @retval phi - porosity
*/
double soil1mat::_phi()
{
  return(0.5);
}

/**
   function computes intrinsic permeability

   @retval kintr - intrinsic permeability
*/
double soil1mat::_kintr()
{

  //funguje s 4.0e-15
  return(4.0e-3/3600.0/1000.0);
}

/**
   function computes gas relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krg - gas relative permeability
*/
double soil1mat::_krg(double pc,double t)
{
  double krg,se,lam,sw,sirr;

  sirr = 0.3216;
  lam = 2.308;
  sw = sat(pc,t);
  se = (sw - sirr)/(1.0 - sirr);
  krg = (1.0 - se)*(1.0 - se)*(1.0 - pow(se,(2.0+lam)/lam));

  return(krg);
}

/**
   function computes water relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krw - water relative permeability
*/
double soil1mat::_krw(double pc,double t)
{
  double krw,se,lam,sw,sirr;

  sirr = 0.3216;
  lam = 2.308;
  sw = sat(pc,t);
  se = (sw - sirr)/(1.0 - sirr);
  krw = pow(se,(2.0+3.0*lam)/lam);

  return(krw);
}

/**
   function computes specific heat of solid skeleton
   
   @retval cps - specific heat capacity of solid skeleton 
*/
double soil1mat::_cps()
{
  //doplnit??!!
  return(0.0);
}


/**
   function computes  volume density of concrete skeleton

   @retval rhos - volume density of concrete skeleton
*/
double soil1mat::_rhos()
{
  return(2000.0);
}


/**
   function reads parameters

   @param in - input file
*/
void soil1mat::read(XFILE */*in*/)
{}


/**
   function prints parameters

   @param out - output file
*/
void soil1mat::print(FILE */*out*/)
{}
