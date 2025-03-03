/*
    File:             sejtkrmatc.cpp
    Author:           Tomas Krejci, 12/9/2008
    Purpose:          material model for saturated-nonsaturated one-phase flow in deforming medium
    sources:          Lewis and Schrefler pp. 86-89
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
   function computes degree of saturation(sorption curve)
   @param pw - water pressure

   @retval sw - degree of saturation

   12/9/2008, TKr
*/
double sejtkrmatc::get_sw(double pw)
{
  //predelat??!!
  double sw,pc;
  //
  pc = -pw;
  
  sw = 1.0;//fully saturated

  return(sw);
}


/**
   function computes partial derivative of degree of saturation with respect to pw
   @param pw - water pressure

   @retval dsw_dpw - partial derivative of degree of saturation with respect to pw

   12/9/2008, TKr
*/
double sejtkrmatc::get_dsw_dpw(double pw)
{
  //predelat??!!
  double dsw_dpc,pc;
  //
  pc = -pw;
  
  dsw_dpc = 0.0;

  return(dsw_dpc);
}


/**
   function computes water relative permeability
   @param pw - water pressure
   
   @retval krw - water relative permeability

   12/9/2008, TKr
*/
double sejtkrmatc::get_krw(double pw)
{
  double krw,s;
  
  s = get_sw(pw);
  krw = pow(s,10.0);
  
  return(krw);
}


/**
   function returns porosity

   @retval phi - porosity

   12/9/2008, TKr
*/
double sejtkrmatc::get_phi()
{
  double phi;

  phi = phi0;

  return(phi);
}


/**
   function computes intrinsic permeability

   @retval kintr - intrinsic permeability

   12/9/2008, TKr
*/
double sejtkrmatc::get_kintr()
{
  return(kintr);
}


/**
   function returns volume density of soil skeleton

   @retval rhos - volume density of concrete skeleton

   12/9/2008, TKr
*/
double sejtkrmatc::get_rhos()
{
  return(rhos);
}


/**
   function returns Biot's constant

   @retval alpha - Biot's constant

   12/9/2008, TKr
*/
double sejtkrmatc::get_alpha()
{
  return(alpha);
}


/**
   function returns emod Young's modulus

   @retval emod = Young's modulus

   12/9/2008, TKr
*/
double sejtkrmatc::get_emod()
{
  return(emod);
}

/**
   function returns nu = Poisson's constant

   @retval nu - Poisson's constant

   12/9/2008, TKr
*/
double sejtkrmatc::get_nu()
{
  return(nu);
}


/**
   function returns ks = 

   @retval ks - 

   12/9/2008, TKr
*/
double sejtkrmatc::get_ks()
{
  return(ks);
}


/**
   function returns kw = 

   @retval kw - 

   12/9/2008, TKr
*/
double sejtkrmatc::get_kw()
{
  return(kw);
}


/**
   function returns muw = 

   @retval muw - 

   12/9/2008, TKr
*/
double sejtkrmatc::get_muw()
{
  return(muw);
}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kuw - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kuw(double pw)
{
  double alpha,s,kuw;
  
  alpha = get_alpha();
  s = get_sw(pw);

  kuw = -alpha*s;

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
  double krw,muw,kww;
  
  krw = get_krw(pw);
  kintr = get_kintr();
  muw = get_muw();

  kww = krw*kintr/muw;

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
  double alpha,s,capwu;
  
  alpha = get_alpha();
  s = get_sw(pw);

  capwu = alpha*s;

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
  double alpha,sw,ks,dsw_dpw,kw,n,capww;  
  
  alpha = get_alpha();
  sw = get_sw(pw);
  dsw_dpw = get_dsw_dpw(pw);
  n = get_phi();
  ks = get_ks();
  kw = get_kw();
  
  capww = (alpha-n)/ks*sw*(sw+dsw_dpw/n*pw)+n*sw/kw + dsw_dpw;

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
  double krw,muw,fw1;
  
  krw = get_krw(pw);
  kintr = get_kintr();
  muw = get_muw();
  rhow = 1000.0;

  fw1 = krw*kintr/muw*rhow;//doplnit!!

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
  double sw,n,rhos,fu1;  
  
  sw = get_sw(pw);
  n = get_phi();
  rhos = get_rhos();
  rhow = 1000.0;

  
  fu1 = rhos*(n-1) + sw*n*rhow;//doplnit!!

  return(fu1);

}
