/*
    File:             concreteB.cpp
    Author:           Tomas Krejci, 1.12.2003
    Purpose:          material properties for concrete at high temperature - only for testing
    sources:           
    1. FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS BY AN ALGEBRAIC MULTIGRID METHOD
    Wang Xicheng, B.A. Schrefler
    
    2. NUMERICAL ANALYSIS OF HYGRO-THERMAL BEHAVIOUR AND DAMAGE OF CONCRETE AT HIGH TEMPERATURE
    D. Gawin, C.E. Majorana, B.A. Schrefler

    3. NONLINEAR MODELLING OF CONCRETE AS MULTIPHASE POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
    Francesco Pesavento - doctoral thesis

*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "concreteB.h"
#include "globalt.h"

concreteBmat::concreteBmat()
{
  mw = 18.01528e-3;//molar mass of water kg.mol-1
  ma = 28.9645e-3;//molar mass of dry air kg.mol-1
  gasr = 8.31441;//universal gas constant J.mol-1.K-1
  
  t0 = 273.15;   //reference temperature
  p0 = 101325.0; //pressure
  tcr = 647.3; //critical point of water
  
  // gas relative permeability
  ////from F.Pesavento's PhD-thesis,pages 
  scr = 0.9;//
  ag = 2.0; //<1;3>
  
  // water relative permeability
  //from F.Pesavento's PhD-thesis,pages 
  sir = 0.2;//
  aw = 2.0; //<1;3>
  bw = 6.0;//or 16.0
  
  /*****************************************************************/
  //input parameters

  // porosity - input parameters
  //from F.Pesavento's PhD-thesis,pages 
  //B35 - silicate concrete
  //phi0 = 0.06;
  //aphi = 0.000195;
  //limestone concrete
  //phi0 = 0.087;
  //aphi = 0.000163;
  //basalt concrete
  //phi0 = 0.0802;
  //aphi = 0.00017;
  
  //intrinsic permeability;  input parameters
  //from F.Pesavento's PhD-thesis,pages 
  //k0 = 1.5e-16;
  //ak = 0.005;
  //bk = 0.368;

  //thermal capacity of solid skeleton; input parameters
  //from F.Pesavento's PhD-thesis,pages 
  tref = 298.15;
  //ac = 0.0;
  //cps0 = 940.0;//thermal capacity of the solid skeleton at reference temperature
  
  //DEGREE OF SATURATION; input parameters
  //Saturation assuming Baroghel formulation extended for high temperature
  //from F.Pesavento's PhD-thesis,pages 198-199
  //ordinary concrete
  //ads = 18.62e6;//[Mpa]
  //bds = 2.2748;
  //nds = 1.2;
  //high performance concrete
  //ads = 46.9364e6;//[Mpa]
  //bds = 2.06;
  //nds = 1.5;
  //ordinary cement paste
  //ads = 37.5479e6;//[Mpa]
  //bds = 2.1684;
  //nds = 1.5;
  //high performance paste
  //ads = 96.2637e6;//[Mpa]
  //bds = 1.9540;
  //nds = 1.5;
  /*******************************************************************************/


  //EFFECTIVE DIFFUSION COEFFICIENT OF VAPOUR
  av = 1.0;//<1;3>
  fs = 1.0;//structure coefficient

  //Hydration energy =0.5 MJ/kg
  hydren = 0.5e+6;
  //finv= aging factor (hydration degree)
  finv = 0.65;
  //fste= Water/Cement ratio (data from Brite 1997)
  fste = 0.36;

  //initialization:
  rhos_th0 = 0.0;

}
concreteBmat::~concreteBmat()
{}

/**
   function computes degree of saturation(desorption curve)
   Saturation assuming Baroghel formulation extended for high temperature

   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double concreteBmat::concreteB_sw(double pc,double t)
  //from F.Pesavento's PhD-thesis,pages 198-199
{
  double sw,eds,gds,a,tt;
  double q0,q2,q3,e0,z;
  
  z = 0.5;
  q3= 30.0e6;

  if (t < 373.15)
    a = ads;
  else{
    tt=(t-373.15)/(647.15-373.15);
    q2 = 25.0e6;
    q0=(q3-q2)*((2.*(tt*tt*tt))-(3.0*(tt*tt))+1.0);
    a = q0 + q2;
  }

  if (t < tcr)
    eds = pow(((tcr - t0)/(tcr - t)),nds); 
  else{
    e0=pow(((tcr-293.15)/z),nds);
    eds = nds/z*e0*t + e0 - nds/z*e0*(tcr - z);
  }

  gds = pow((eds/a*pc),(bds/(bds-1.0)));
  sw = pow((gds + 1.0),(-1.0/bds));

  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   Saturation assuming Baroghel formulation extended for high temperature
   
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dpc - partial derivative of degree of saturation with respect to pc
*/
double concreteBmat::concreteB_dsw_dpc(double pc,double t)
  //from F.Pesavento's PhD-thesis,pages 198-199
{
  double dsw_dpc,eds,a,dg_dpc,tt,gds;
  double q0,q2,q3,e0,z;
  
  z = 0.5;
  q3= 30.0e6;

  if (t < 373.15)
    a = ads;
  else{
    tt=(t-373.15)/(647.15-373.15);
    q2 = 25.0e6;
    q0=(q3-q2)*((2.*(tt*tt*tt))-(3.0*(tt*tt))+1.0);
    a = q0 + q2;
  }
  
  if (t < tcr){
    eds = pow(((tcr - t0)/(tcr - t)),nds); 
  }
  else{
   e0=pow(((tcr-293.15)/z),nds);
   eds = nds/z*e0*t + e0 - nds/z*e0*(tcr - z);
  }
  gds = pow((eds/a*pc),(bds/(bds-1.0)));
  dg_dpc = (bds/(bds-1.0))*eds/a*pow((eds/a*pc),(bds/(bds-1.0)-1.0));
  dsw_dpc = (-1.0/bds)*pow((gds + 1.0),(-1.0-1.0/bds))*dg_dpc;
  
  return(dsw_dpc);
}


/**
   function computes partial derivative of degree of saturation with respect to t
   Saturation assuming Baroghel formulation extended for high temperature
   
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double concreteBmat::concreteB_dsw_dt(double pc,double t)
  //from F.Pesavento's PhD-thesis,pages 198-199
{
  double dsw_dt,eds,de_dt,gds,a,da_dt,dea_dt,dg_dt,tt,dtt_dt;
  double q0,dq0_dt,q2,q3,e0,z;
  
  z = 0.5;
  q3= 30.0e6;

  if (t < 373.15){
    a = ads;
    da_dt = 0.0;
  }
  else{
    tt=(t-373.15)/(647.15-373.15);
    q2 = 25.0e6;
    q0=(q3-q2)*((2.*(tt*tt*tt))-(3.0*(tt*tt))+1.0);
    dtt_dt = 1.0/(647.15-373.15);
    dq0_dt = (q3-q2)*(6.0*tt*tt*dtt_dt - 6.0*tt*dtt_dt + 1.0);
    a = q0 + q2;
    da_dt = dq0_dt;
  }
  
  if (t < tcr){
    eds = pow(((tcr - t0)/(tcr - t)),nds);
    de_dt = nds*pow(((tcr - t0)/(tcr - t)),(nds-1.0))*(tcr - t0)/(tcr - t)/(tcr - t);
  }
  else{
    e0=pow(((tcr-293.15)/z),nds);
    eds = nds/z*e0*t + e0 - nds/z*e0*(tcr - z);
    de_dt = nds/z*e0;
  }

  gds = pow((eds/a*pc),(bds/(bds-1.0)));
  dea_dt = (de_dt*a - da_dt*eds)/a/a;
  dg_dt = (bds/(bds-1.0))*pc*pow((eds/a*pc),(bds/(bds-1.0)-1.0))*dea_dt;
  dsw_dt = (-1.0/bds)*pow((gds + 1.0),(-1.0-1.0/bds))*dg_dt;
    
  return(dsw_dt);
}

/**
   function returns saturation solid point

   @retval ssp - saturation solid point
*/
double concreteBmat::concreteB_ssp()
{
  return(0.55);
}


/**
   function computes gas relative permeability
   @param s - degree of saturation

   @retval krg - relative permeability
*/
double concreteBmat::concreteB_krg(double s)
{
  double krg;

  krg = 1.0 - pow((s/scr),ag);

  return(krg);
}

/**
   function computes water relative permeability
   @param s - degree of saturation
   @param rh - relative humidity

   @retval krw - water relative permeability
*/
double concreteBmat::concreteB_krw(double s,double rh)
{
  double krw,help;

  if (rh < 0.75)
    krw = pow(((s-sir)/(1.0-sir)),aw);
  else{
    help = 1.0 + pow(((1.0-rh)/0.25),bw); 
    krw = pow(help,-1.0)*pow(s,aw);
  }

  return(krw);
}

/**
   function computes porosity
   @param t - temperature

   @retval phi - porosity
*/
double concreteBmat::concreteB_phi(double t)
{
  double phi;

  phi = phi0 + aphi*(t - t0);

  return(phi);
}

/**
   function computes intrinsic permeability
   @param pg - capillary gas pressure
   @param t - temperature
   @param dam - damage parameter

   @retval kintr - intrinsic permeability
*/
double concreteBmat::concreteB_kintr(double pg,double t,double dam)
{
  double kintr,ad;
      
  kintr = k0*pow(10.0,(ak*(t-t0)))*pow((pg/p0),bk);

  //damge effect - from Francesco Pesavento PhD thesis page 215
  ad = 4.0;
  kintr = kintr + k0*pow(10.0,(ad*dam));//added eq. 5.75

  //another formula eq. 5.76
  //kintr = k0*pow(10.0,(at*(t-298.15)))*pow((pg/p0),ap)*pow(10.0,(ad*dam));

  return(kintr);
}

/**
   function computes specific heat of solid skeleton
   @param t - temperature

   @retval cps - specific heat of solid skeleton
*/
double concreteBmat::concreteB_cps(double t)
{
  double cps;

  cps = cps0*(1.0 + ac*(t - tref));

  return(cps);
}

/**
   function computes thermal capacity of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval rhocp -  of thermal capacity partially saturated concrete
*/
double concreteBmat::concreteB_rhocp(double pc,double pg,double t,long /*ipp*/)
{
  double s,phi,rhocp,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos;
  state_eq tt;

  s = concreteB_sw(pc,t);
  phi = concreteB_phi(t);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  cps = concreteB_cps(t);
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();
  rhos = concreteB_rhos(t);

  rhocp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));

  return(rhocp);
}


/**
   function computes cpecific heat of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval cp - specific heat of partially saturated concrete
*/
double concreteBmat::concreteB_cp(double pc,double pg,double t,long ipp)
{
  double s,phi,cp,rho,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos;
  state_eq tt;

  s = concreteB_sw(pc,t);
  phi = concreteB_phi(t);
  rho = tt.get_rho(pc,pg,t,ipp);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  cps = concreteB_cps(t);
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();
  rhos = concreteB_rhos(t);

  cp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));
  cp = cp/rho;

  return(cp);
}


/**
   function computes structure factor
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval fs*tau - structure factor
*/ 
double concreteBmat::concreteB_fs(double pc,double t)
{
  double f,tau;
  
  tau = concreteB_tau(pc,t);
  
  f = fs*tau;
  
  return(f);
}

/**
   function computes tortuosity factor
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval tau - tortuosity factor
*/ 
double concreteBmat::concreteB_tau(double pc,double t)
{
  double phi,s,tau;
  
  s = concreteB_sw(pc,t);
  phi = concreteB_phi(t);
  
  tau = pow(phi,(1.0/3.0))*pow((1.0-s),(7.0/3.0));

  return(tau);
}

/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval deff - effective diffusion coefficient of vapour inside pores
*/
double concreteBmat::concreteB_deff(double pc,double pg,double t)
{
  //from Francesco Pesavento PhD thesis page 203
  double deff,f,phi,s,cdiff;
  state_eq tt;
  
  phi = concreteB_phi(t);
  s = concreteB_sw(pc,t);
  cdiff = tt.get_cdiff(pc,pg,t);
  f = concreteB_fs(pc,t);

  deff = phi*pow((1.0 - s),av)*f*cdiff;

  return(deff);
}

/**
   function computes effective thermal conductivity of partially saturated concrete
   @param pc - capillary pressure
   @param t - temperature

   @retval lambdaeff - effective thermal conductivity of partially saturated concrete
*/
double concreteBmat::concreteB_lambdaeff(double pc,double /*pg*/,double t)
{
  double lambdaeff,lambdad,lambdad0,alam,s,phi,rhow,rhos;
  state_eq tt;

  s = concreteB_sw(pc,t);
  phi = concreteB_phi(t);
  lambdad0 = 1.67;//W/(m.K) 
  alam = -1.017e-3;//K-1
  rhow = tt.get_rhow(t);
  rhos = concreteB_rhos(t);
  
  if (t < tcr)
    lambdad = lambdad0*(1.0 + alam*(t-298.15));
   else
     lambdad = lambdad0*(1.0 + alam*(tcr-298.15));

  lambdaeff = lambdad*(1+ 4.0*phi*rhow*s/(1.0-phi)/rhos);

  return(lambdaeff);
}

/**
   function computes bulk modulus of porous medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval kt - bulk modulus of porous medium
*/
double concreteBmat::concreteB_kt(double pc,double pg,double t)
{
  double kt,ks;
  
  ks = concreteB_ks(pc,pg,t);//not completed

  kt = ks;

  return(kt);
}

/**
   function computes bulk modulus of solid phase
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval ks - bulk modulus of solid phase
*/
double concreteBmat::concreteB_ks(double /*pc*/,double /*pg*/,double /*t*/)
  // bulk modulus of solid phase
{
  double ks;
  
  ks = emod/(3.0*(1.0 - 2.0*nu));
  
  return(ks);
}

/**
   function computes  volume density of concrete skeleton,
   changes of solid density, caused by dehydratation process

   @param t - temperature

   @retval rhos - volume density of concrete skeleton
*/
double concreteBmat::concreteB_rhos(double t)
{
  double rhos,ah,aphi,phit,phi_th0;
   
  ah = 0.0;//provisionally
  aphi = 1.0;//provisionally
  //th0 = 0.0;//provisionally
  
  phit = concreteB_phi(t);
  phi_th0 = 0.0;//provisionally
  
  rhos = rhos_th0 + ah/aphi*log((1.0-phit)/(1.0-phi_th0));
  
  return(rhos);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval beats - cubic thermal expansion coefficient of solid (K-1)
*/
double concreteBmat::concreteB_betas()
{
  return(betas);
}

/**
   function computes emod Young's modulus

   @retval emod - Young's modulus
*/
double concreteBmat::concreteB_emod()
{
  return(emod);
}

/**
   function returns Poisson's coefficient

   @retval nu - Poisson's coefficient
*/
double concreteBmat::concreteB_nu()
{
  return(nu);
}

/**
   function computes dehydrated water content
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dehydw_dt - changes of hydrated water content with temperature
*/
double concreteBmat::concreteB_dmdh_dt(double /*pc*/,double /*pg*/,double t)
{
  double dfhyt,dehydw_dt;
  
  if((t-t0)< 105.0){
    //fhy = 0.0;
    dfhyt = 0.0;
  }
  else{
    //fhy   = (1.0+sin(3.1416/2.0*(1.0-2.0*exp(-0.004*((t-t0)-105.0)))))/2.0;
    dfhyt = (3.1416*0.004/2.0)*cos(3.1416/2.0*(1.0-2.0*exp(-0.004*((t-t0)-105.0))))*exp(-0.004*((t-t0)-105.0));
  }
  
  //hydw = fste*finv*c1*fhy;
  dehydw_dt = fste*finv*c1*dfhyt;
  
  return(dehydw_dt);
}


/**
   function computes hydration energy
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dhdehydr - hydration energy
*/
double concreteBmat::concreteB_dhdehydr(double /*pc*/,double /*pg*/,double /*t*/)
{
  double dhdehydr;
  
  dhdehydr = hydren;

  return(dhdehydr);
}

/**
   function computes derivative of apparent density with respect to degree of hydration
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval drhos_dgammadh - derivative of apparent density with respect to degree of hydration
*/
double concreteBmat::concreteB_drhos_dgammadh(double /*pc*/,double /*pg*/,double /*t*/)
{
  double drhos_dgammadh;
  
  drhos_dgammadh = 0.0;//not completed

  return(drhos_dgammadh);
}

/**
   function computes changes of degre of hydration with temperature
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dgammadh_dt - derivative of degre of hydration with respect to temperature
*/
double concreteBmat::concreteB_dgammadh_dt(double /*pc*/,double /*pg*/,double /*t*/)
{
  double dgammadh_dt;
  
  dgammadh_dt = 0.0;//not completed

  return(dgammadh_dt);
}

/**
   function reads parameters
   @param in - input file
*/
void concreteBmat::read(XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	  &emod, &nu, &rhos_th0, &betas, &phi0, &aphi, &ads, &bds, &nds, &k0, &ak, &bk, &ac, &cps0);
}

/**
   function prints parameters
   @param out - output file
*/
void concreteBmat::print(FILE *out)
{
  fprintf (out,"  %e %e %e %e %e %e %e %e %e %e %e %e %e %e", 
	  emod, nu, rhos_th0, betas, phi0, aphi, ads, bds, nds, k0, ak, bk, ac, cps0);
}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   29. 11. 2013 by TKo
*/
void concreteBmat::give_reqntq(long *antq)
{
  //  damage parameter
  antq[scal_iso_damage-1] = 1;
}



