/*
  File:             C30_baroghelc.cpp
  Author:           Tomas Krejci, 1.12.2003
  Purpose:          Calculates properties of C30 concrete (Ordinary Performance Concrete)
                    Saturation assuming Baroghel formulation extended for high temperature
  sources:          USERMATsifel.f90 and SATBGMsifel.f90 from Padova
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "C30baroghelc.h"
#include "globalt.h"

C30barmatc::C30barmatc()
{
  mw = 18.01528;   //molar mass of water kg.mol-1
  ma = 28.9645;    //molar mass of dry air kg.mol-1
  gasr = 8314.41;  //universal gas constant J.mol-1.K-1
  
  t0 = 273.15;   //reference temperature
  t00 = 293.15;
  p0 = 101325.0; //pressure
  tcr = 647.3;   //critical point of water
  
  w1 = 100.0;
  c1 = 200.0;
  
  // porosity
  phi0 = 0.1368;
  aphi = 7.8e-5;
  
  //intrinsic permeability
  k0 = 3.2e-18;
  ak = 5.0e-3;
  
  // gas relative permeability
  scr = 1.0;//
  ag = 1.0; //<1;3>
  
  // water relative permeability
  sir = 0.2;//
  aw = 2.0; //<1;3>
  bw = 6.0;//or 16.0
  
  //skeleton density, experimental data by Kalifa(report BRITE 1997)
  rhos = 2625.8;
  
  fs = 1.0;//structure coefficient

  //thermal conductivity of solid skeleton
  lambdas0 = 2.0;//W/(m.K) 
  alam = -1.017e-3;//K-1
  
  //thermal capacity of solid skeleton
  ac = 0.35;
  cps0 = 940.0;//thermal capacity of the solid skeleton at reference temperature 
  
  //Hydration energy =0.5 MJ/kg
  hydren = 0.5e+6;
  //finv= aging factor (hydration degree)
  finv = 0.4;
  //fste= Water/Cement ratio (data from Brite 1997)
  fste = 0.24;

  //chracteristic length
  dld = 4.0E-2;

  //Young's mosulus
  emod0  = 3.0e+10;//in Pa at 20°C

  //Poisson's ceofficient
  vcoeff = 0.2;

  //cubic thermal expansion coeffcient
  betas  = 9.0e-6;

  //Biot's constant
  alpha = 0.5;//rigid porous materials

  //MAZAR'S COEFFICIENTS
  //tensile coefficients
  //functions not completed??!!
  at = 0.8;
  bt = 20000.0;
  //compressive coefficients
  acc = 1.5;

  ddbw0 = 1.0e-20; //diffusion of bound water at refference temperature
}

C30barmatc::~C30barmatc()
{}

/**
   function computes degree of saturation(desorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double C30barmatc::sat(double pc,double t)
{
  double sw,n,a,b,tt,q0,q2,q3,z,e,e0,g;
  
  b = 2.27; 
  q3 = 18.62e6;
  q2 = 7.0e6;
  n = 1.2;
  z = 0.5; //z is governing the transition through critical temperature
  
  if (t <= 373.15)
    a = q3;
  else{
    tt=(t-373.15)/(647.15-373.15);
    q0=(q3-q2)*((2.*(tt*tt*tt))-(3.0*(tt*tt))+1.0);
    q2 = 25.0e6;
    a = q0 + q2;
  }

  if (t <= tcr)
    e = pow(((tcr - t0)/(tcr - t)),n); 
  else{
    e0=pow(((tcr-293.15)/z),n);
    e = n/z*e0*t + e0 - n/z*e0*(tcr - z);
  }

  g = pow((e/a*pc),(b/(b-1.0)));
  sw = pow((g + 1.0),(-1.0/b));
  
  return(sw);
  
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double C30barmatc::dsat_dpc(double pc,double t)
{
  double dsw_dpc;
  double n,a,b,dg_dpc,tt,q0,q2,q3,z,e,e0,g;
  
  b = 2.27; 
  q3 = 18.62e6 ;
  q2 = 7.0e6;
  n = 1.2;
  z = 0.5; //z is governing the transition through critical temperature
  
  if (t <= 373.15)
    a = q3;
  else{
    tt=(t-373.15)/(647.15-373.15);
    q0=(q3-q2)*((2.*(tt*tt*tt))-(3.0*(tt*tt))+1.0);
    q2 = 25.0e6;
    a = q0 + q2;
  }
  
  if (t <= tcr){
    e = pow(((tcr - t0)/(tcr - t)),n); 
  }
  else{
   e0=pow(((tcr-293.15)/z),n);
   e = n/z*e0*t + e0 - n/z*e0*(tcr - z);
  }
  g = pow((e/a*pc),(b/(b-1.0)));
  dg_dpc = (b/(b-1.0))*e/a*pow((e/a*pc),(b/(b-1.0)-1.0));
  dsw_dpc = (-1.0/b)*pow((g + 1.0),(-1.0-1.0/b))*dg_dpc;

  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double C30barmatc::dsat_dt(double pc,double t)
{
  double dsw_dt,e,b,n,de_dt,g,a,da_dt,dea_dt,dg_dt,tt,dtt_dt,q0,dq0_dt,q2,q3,e0,z;
  
  b = 2.27; 
  q3 = 18.62e6 ;
  q2 = 7.0e6;
  n = 1.2;
  z = 0.5; //z is governing the transition through critical temperature
  
  if (t <= 373.15){
    a = q3;
    da_dt = 0.0;
  }
  else{
    tt=(t-373.15)/(647.15-373.15);
    q0=(q3-q2)*((2.*(tt*tt*tt))-(3.0*(tt*tt))+1.0);
    dtt_dt = 1.0/(647.15-373.15);
    dq0_dt = (q3-q2)*(6.0*tt*tt*dtt_dt - 6.0*tt*dtt_dt);
    q2 = 25.0e6;
    a = q0 + q2;
    da_dt = dq0_dt;
  }
  
  if (t <= tcr){
    e = pow(((tcr - t0)/(tcr - t)),n);
    de_dt = n*pow(((tcr - t0)/(tcr - t)),(n-1.0))*(tcr - t0)/(tcr - t)/(tcr - t);
  }
  else{
    e0=pow(((tcr-293.15)/z),n);
    e = n/z*e0*t + e0 - n/z*e0*(tcr - z);
    de_dt = n/z*e0;
  }

  g = pow((e/a*pc),(b/(b-1.0)));
  dea_dt = (de_dt*a - da_dt*e)/a/a;
  dg_dt = (b/(b-1.0))*pc*pow((e/a*pc),(b/(b-1.0)-1.0))*dea_dt;
  dsw_dt = (-1.0/b)*pow((g + 1.0),(-1.0-1.0/b))*dg_dt;

  return(dsw_dt);
}

/**
   function returns saturation solid point

   @retval ssp - saturation solid point
*/
double C30barmatc::ssp()
{
  return(0.55);
}

/**
   function computes porosity, Data by ALONSO-ANDRADE 
   (report BRITE 1997) - results valid up to 600°C
   
   @param t - temperature
   @retval phi - porosity
*/
double C30barmatc::C30bar_phi(double t)
{
  double phi;
  
  phi = phi0 + aphi*(t - t00);
  
  return(phi);
}

/**
   function computes intrinsic permeability with damage effect
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param dam - damage parameter

   @retval kintr - intrinsic permeability
*/
double C30barmatc::C30bar_kintr(double /*pc*/,double pg,double t,double /*dam*/)
{
  double kintr,bk;//ad;
  state_eq tt;

  bk = log(5.0/3.0)/log(4.0);
  
  kintr = k0*exp(log(10.0)*ak*(t-298.15))*pow((pg/p0),bk);

  //damage effect - from Francesco Pesavento PhD thesis page 215
  //ad = 4.0;
  //kintr = kintr + k0*pow(10.0,(ad*dam));//added eq. 5.75

  //another formula eq. 5.76
  //kintr = k0*pow(10.0,(at*(t-298.15)))*pow((pg/p0),ap)*pow(10.0,(ad*dam));
 
  return(kintr);
}

/**
   function computes gas relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krg - gas relative permeability
*/
double C30barmatc::C30bar_krg(double pc,double t)
{
  double krg,s;
  
  s = sat(pc,t);

  krg = 1.0 - pow((s/scr),ag);

  return(krg);
}

/**
   function computes water relative permeability
   @param pc - capillary pressure
   @param t - temperature
   @param rh - relative humidity

   @retval krw - water relative permeability
*/
double C30barmatc::C30bar_krw(double pc,double t,double rh)
{
  double krw,help,s;

  s = sat(pc,t);

  //not corrected
  if (rh < 0.75)
    krw = pow(((s-sir)/(1.0-sir)),aw);
  else{
    help = 1.0 + pow(((1.0-rh)/0.25),bw); 
    krw = pow(help,-1.0)*pow(s,aw);
  }

  return(krw);
}

/**
   function computes dd
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dd - ...
*/
double C30barmatc::C30bar_dd(double pc,double t)
{
  double tau,dd;

  tau = C30bar_tau(pc,t);
  
  dd = tau*mw*ma/gasr;//according to Frotran code
  //dd = tau;//according to Francesco's PhD thesis

  return(dd);
}

/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval deff - effective diffusion coefficient of vapour inside pores
*/
double C30barmatc::C30bar_deff(double pc,double pg,double t)
{
  //from Francesco Pesavento PhD thesis page 203
  double deff,dd,phi,s,cdiff;
  state_eq tt;
  
  dd = C30bar_dd(pc,t);
  phi = C30bar_phi(t);
  s = sat(pc,t);
  cdiff = tt.get_cdiff(pc,pg,t);

  deff = dd*phi*(1.0 - s)*fs*cdiff;

  return(deff);
}

/**
   function computes specific heat of solid skeleton
   @param t - temperature
   
   @retval cps - specific heat of solid skeleton 
*/
double C30barmatc::C30bar_cps(double t)
{
  double cps;
  
  if (t < tcr)
    cps = cps0 + ac*(t - 298.15);
  else
    cps = cps0 + ac*(tcr - 298.15);
  
  return(cps);
}

/**
   function computes thermal capacity of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   
   @retval rhocp -  of thermal capacity partially saturated concrete
*/
double C30barmatc::C30bar_rhocp(double pc,double pg,double t)
{
  double s,phi,rhos,rhocp,rhow,rhog,rhogw,cps,cpw,cpga,cpgw;
  state_eq tt;
  
  s = sat(pc,t);
  phi = C30bar_phi(t);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  rhos = C30bar_rhos();
  cps = C30bar_cps(t);
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();

  rhocp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));

  return(rhocp);
}

/**
   function computes tortuosity factor, Formulation by Baroghel 
   @param pc - capillary pressure
   @param t - temperature

   @retval tau - tortuosity factor
*/ 
double C30barmatc::C30bar_tau(double pc,double t)
{
  double phi,s,tau;
  
  s = sat(pc,t);
  phi = C30bar_phi(t);
  
  tau = pow(phi,(1.0/3.0))*pow((1.0-s),(7.0/3.0));
  
  return(tau);
}

/**
   function computes solid thermal conductivity, experimental data by Kalifa
   (report BRITE 1997) - data valid up to 600°C

   @param pc - capillary pressure
   @param t - temperature

   @retval lambdas - solid thermal conductivity
*/
double C30barmatc::C30bar_lambdas(double t)
{
  double lambdas;

  if (t < tcr)
    lambdas = lambdas0 + alam*(t-298.15);
  else
    lambdas = lambdas0 + alam*(tcr-298.15);

  return(lambdas);
}

/**
   function computes effective thermal conductivity of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   
   @retval lambdaeff - effective thermal conductivity of partially saturated concrete
*/
double C30barmatc::C30bar_lambdaeff(double pc,double /*pg*/,double t)
{
  double lambdaeff,lambdas,s,phi,rhow,rhos;
  state_eq tt;

  s = sat(pc,t);
  phi = C30bar_phi(t);
  rhow = tt.get_rhow(t);
  lambdas = C30bar_lambdas(t);
  rhos = C30bar_rhos();

  lambdaeff = lambdas*(1.0 + 4.0*phi*rhow*s/(1.0-phi)/rhos);

  return(lambdaeff);
}

/**
   function computes  volume density of concrete skeleton

   @retval rhos - volume density of concrete skeleton
*/
double C30barmatc::C30bar_rhos()
{
  return(rhos);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval beats - cubic thermal expansion coefficient of solid (K-1)
*/
double C30barmatc::C30bar_betas()
{
  return(betas);
}




/**
   function computes hydration degree
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval hydw - hydration degree
*/
double C30barmatc::C30bar_hydw(double /*pc*/,double /*pg*/,double t)
{
  double fhy,hydw;
  
  if((t-t0)< 105.0){
    fhy = 0.0;
  }
  else{
    fhy   = (1.0+sin(3.1416/2.0*(1.0-2.0*exp(-0.004*((t-t0)-105.0)))))/2.0;
  }
  
  hydw = fste*finv*c1*fhy;

  return(hydw);
}

/**
   function computes derivative of hydration degree with respect to temperature
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dehydw_dt - derivative of hydration degree with respect to temperature
*/
double C30barmatc::C30bar_dehydw_dt(double /*pc*/,double /*pg*/,double t)
{
  double dfhyt,dehydw_dt;
  
  if((t-t0)< 105.0){
    dfhyt = 0.0;
  }
  else{
    dfhyt = (3.1416*0.004/2.0)*cos(3.1416/2.0*(1.0-2.0*exp(-0.004*((t-t0)-105.0))))*exp(-0.004*((t-t0)-105.0));
  }
  
  dehydw_dt = fste*finv*c1*dfhyt;

  return(dehydw_dt);
}


/**
   function computes hydration energy
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval hydren - hydration energy
*/
double C30barmatc::C30bar_hydren(double /*pc*/,double /*pg*/,double /*t*/)
{
  return(hydren);
}

/**
   function computes Water/Cement ratio
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval fste - Water/Cement ratio
*/
double C30barmatc::C30bar_fste(double /*pc*/,double /*pg*/,double /*t*/)
{
  return(fste);
}


/**
   function computes diffusivity of bound water
   @param pc - capillary pressure
   @param t - temperature

   @retval  - diffusivity of bound water - according to Frotran code
*/
double C30barmatc::C30bar_ddbw(double /*pc*/,double /*pg*/,double t)
{
  double ddbw;
  
  ddbw = ddbw0*exp(-t/(273.15+23));
  if (t > tcr)
    ddbw = ddbw0*exp(-tcr/(273.15+23));

  return(ddbw);
}



/**
   function computes emod Young's modulus
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval emod - Young's modulus
*/
double C30barmatc::C30bar_emod(double /*pc*/,double /*pg*/,double t)
{
  double emod,ttc;
  
  ttc = t - 273.15;
  
  if (ttc <= 50.0)
    emod = emod0;
  else
    emod = emod0*(3.14e-6*pow(ttc,2.0)-3.77e-3*(ttc)+1.1806);
  if (ttc >= 600.0)
    emod = emod0*0.05;

  return(emod);
}

/**
   function computes fct tensile strenght  (Brite data - Felicetti 1999)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval fct - tensile strenght
*/
double C30barmatc::C30bar_fct(double /*pc*/,double /*pg*/,double t)
{
  double fct,ttc;
  
  ttc   = t - 273.15;
  
  if (ttc <= 600.0)
    fct = (6.0 - 8.56e-3*ttc)*1.0e6;
  else
    fct = 0.864*1.0e6;

  return(fct);
}

/**
   function computes xk0 maximum linear elastic tensile strain (elastic threshold)

   @retval xk0 - maximum linear elastic tensile strain
*/
double C30barmatc::C30bar_xk0()
{
  double xk0;
  
  xk0 = 1.0e-4;

  return(xk0);
}

/**
   function computes bcc compressive coefficient

   @retval bcc - compressive coefficient
*/
double C30barmatc::C30bar_bcc()
{
  double bcc;
  
  bcc = 2000.0;

  return(bcc);
}

/**
   function returns Biot's constant

   @retval alpha - Biot's constant
*/
double C30barmatc::C30bar_alpha()
{
  double alpha;

  alpha = 0.5;

  return(alpha);
}

/**
   function computes nu Poisson's constant

   @retval vcoeff - Poisson's constant
*/
double C30barmatc::C30bar_nu()
{
  return(vcoeff);
}

/**
   function reads parameters
   
   @param in - input file
*/
void C30barmatc::read(XFILE */*in*/)
{}
