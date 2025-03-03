/*
  File:             o30_bazant.cpp
  Author:           Tomas Krejci, 1.12.2003
  Purpose:          Calculates the properties of ordinary concrete (30MPa)
                    Saturation assuming Bazant formulation extended for high temperature
  sources:          SATUR_A.f90 and Cor1sifel.f90 from Padova
  
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "o30bazant.h"

o30bazmat::o30bazmat()
{
  mw = 18.01528;   //molar mass of water kg.mol-1
  ma = 28.9645;    //molar mass of dry air kg.mol-1
  gasr = 8314.41;  //universal gas constant J.mol-1.K-1
  
  t0 = 273.15;   //reference temperature
  t00 = 293.15;
  p0 = 101325.0; //pressure
  tcr = 647.3;   //critical point of water
  
  w1 = 98.0;
  c1 = 510.0;
  
  // porosity
  phi0 = 0.1;
  aphi = 5.0e-5;
  
  //intrinsic permeability
  k0 = 1.0e-21;
  ak = 5.0e-3;
  
  // gas relative permeability
  scr = 1.0;//
  ag = 1.0; //<1;3>
  
  // water relative permeability
  sir = 0.2;//
  aw = 2.0; //<1;3>
  bw = 6.0;//or 16.0
  
  //thermal conductivity of solid skeleton
  lambdas0 = 1.67;//W/(m.K) 
  alam = -0.25/300.0;//K-1
  
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
  dld = 3.0E-2;

  //Young's modulus
  emod0  = 3.4e+10;//in Pa at 20°C

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
o30bazmat::~o30bazmat()
{}

/**
   function computes degree of saturation(desorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double o30bazmat::sat(double pc,double t)
{
  double sw,rh,tamb,tem,tt1,mt1,fcm,mt;
  state_eq tt;

  tamb= 25.0;
  rh = tt.get_rh(pc,t);
  
  //temperature in Celsius degrees
  if (t < tcr){
    tem = t - 273.15;
    if (tem < tamb)
      tem = tamb;
  }
  else
    tem = tcr - 273.15;
  
  //function SATUR_A.f90
  tt1 = ((tem+10.0)/(25.0+10.0))*((tem+10.0)/(25.0+10.0));//modified temperature from Bazant formulation
  mt1 = 1.04 - (tt1/(24.0+tt1));//function m from Bazant formulation
  if(t < tcr){
    fcm = (tcr-t)/(tcr-298.15);
    mt  = mt1*fcm;//bazant parameter (below critical point)
    if(mt < 1.0e-4) 
      mt = 1.0e-4;
  }
  else{
    fcm = 0.0;
    mt = 1.0e-4;  //Bazant parameter (above critical point)
  }

  //conversion pow(a,b) -> exp(b*log(a))
  //sw = pow(rh,(1.0/mt));//saturation calculation
  sw = exp(log(rh)/mt);//saturation calculation
  
  if(sw < 0.0)//limit for saturation level
    {
      sw=0.0;
      fprintf (stderr,"\n\n Uprava saturace");
    } 
  
  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double o30bazmat::dsat_dpc(double pc,double t)
{
  double dsw_dpc;
  double sw,rh,tamb,tem,tt1,mt1,fcm,mt;
  double ssmh,drhdc,dsdc;
  state_eq tt;

  tamb= 25.0;
  rh = tt.get_rh(pc,t);
  drhdc = tt.get_drh_dpc(pc,t);

   //temperature in Celsius degrees
  if (t < tcr){
    tem = t - 273.15;
    if (tem < tamb)
      tem = tamb;
  }
  else
    tem = tcr - 273.15;

  tt1 = ((tem+10.0)/(25.0+10.0))*((tem+10.0)/(25.0+10.0));//modified temperature from Bazant formulation
  mt1 = 1.04 - (tt1/(24.0+tt1));//function m from Bazant formulation
  if(t < tcr){
    fcm = (tcr-t)/(tcr-298.15);
    mt  = mt1*fcm;//bazant parameter (below critical point)
    if(mt<1.0e-4) 
      mt=1.0e-4;
  }
  else{
    fcm = 0.0;
    mt = 1.0e-4;  //Bazant parameter (above critical point)
  }
  
  //conversion pow(a,b) -> exp(b*log(a))
  //sw = pow(rh,(1.0/mt));//saturation calculation
  //sw = exp(1.0/mt*log(rh));//saturation calculation
  sw = sat(pc,t);
  
  ssmh = sw/mt/rh;
  dsdc = ssmh*drhdc;

  dsw_dpc = dsdc;  

  return(dsw_dpc);
}


/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double o30bazmat::dsat_dt(double pc,double t)
{
  double dsw_dt;
  double sw,rh,tamb,tem,tt1,mt1,fcm,mt,dtt1dt,dm1tdt,dfcmdt,dmtdt;
  double lnh,xfun,drhdt,ssmh,dsdt;
  state_eq tt;

  tamb= 25.0;
  rh = tt.get_rh(pc,t);
  drhdt = tt.get_drh_dt(pc,t);

  //temperature in Celsius degrees
  if (t < tcr){
    tem = t - 273.15;
    if (tem < tamb)
      tem = tamb;
  }
  else
    tem = tcr - 273.15;

  //function SATUR_A.f90
  tt1 = ((tem+10.)/(25.0+10.))*((tem+10.)/(25.0+10.));//modified temperature from Bazant formulation
  mt1 = 1.04 - (tt1/(24.+tt1));//function m from Bazant formulation
  if(t < tcr){
    fcm = (tcr-t)/(tcr-298.15);
    mt  = mt1*fcm;//bazant parameter (below critical point)
    dtt1dt = 2.0*(tem+10.0)/((25.0+10.0)*(25.0+10.0));
    dm1tdt = -24.0*dtt1dt/((24.0+tt1)*(24.0+tt1));
    dfcmdt = -1.0/(tcr-298.15);
    if(mt<1.0e-4)  
      mt=1.0e-4;
  }
  else{
    fcm = 0.0;
    mt = 1.0e-4;  //Bazant parameter (above critical point)
    dtt1dt = 0.0;
    dm1tdt = 0.0;  
    dfcmdt = 0.0;
  }
  
  dmtdt = dm1tdt*fcm + mt1*dfcmdt;
  
  //conversion pow(a,b) -> exp(b*log(a))
  //sw = pow(rh,(1.0/mt));//saturation calculation
  //sw = exp(1.0/mt*log(rh));//saturation calculation
  sw = sat(pc,t);
  
  lnh  = log(rh);
  
  if(t < tcr){//t=temperature on gauss point ; tcr is critical temperature
    xfun = drhdt - rh/mt*lnh*dmtdt;
    ssmh = sw/mt/rh;
    dsdt = ssmh*xfun;
  }
  else{
    xfun = drhdt - rh/mt*lnh*dmtdt ;
    ssmh = sw/mt/rh;
    dsdt = 0.0;
  }
  
  dsw_dt = dsdt;

  return(dsw_dt);
}


/**
   function returns saturation solid point

   @retval ssp - saturation solid point
*/
double o30bazmat::ssp()
{
  return(0.55);
}

/**
   function computes porosity, Data by ALONSO-ANDRADE 
   (report BRITE 1997)- results valid up to 600°C
   
   @param t - temperature
   @retval phi - porosity
*/
double o30bazmat::o30baz_phi(double t)
{
  double phi;
  
  phi = phi0 + aphi*(t - 298.15);
  
  return(phi);
}

/**
   function computes intrinsic permeability
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param dam - damage parameter

   @retval kintr - intrinsic permeability
*/
double o30bazmat::o30baz_kintr(double /*pc*/,double pg,double t,double /*dam*/)
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

   @retval krg - relative permeability
*/
double o30bazmat::o30baz_krg(double pc,double t)
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
double o30bazmat::o30baz_krw(double pc,double t,double rh)
{
  double krw,help,s;

  s = sat(pc,t);

  //not corrected??!!
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

   @retval dd - ...
*/
double o30bazmat::o30baz_dd()
{
  double tau,dd;

  tau = o30baz_tau();
  
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
double o30bazmat::o30baz_deff(double pc,double pg,double t)
{
  //from Francesco Pesavento PhD thesis page 203
  double deff,dd,phi,s,cdiff,fs;
  state_eq tt;

  dd = o30baz_dd();
  phi = o30baz_phi(t);
  s = sat(pc,t);
  cdiff = tt.get_cdiff(pc,pg,t);
  fs = o30baz_fs(pc,pg,t);

  deff = dd*phi*(1.0 - s)*fs*cdiff;

  return(deff);
}

/**
   function computes specific heat of solid skeleton
   @param t - temperature
   
   @retval cps - specific heat capacity of solid skeleton 
*/
double o30bazmat::o30baz_cps(double t)
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
   @param ipp - number of integration point
   
   @retval rhocp -  of thermal capacity partially saturated concrete
*/
double o30bazmat::o30baz_rhocp(double pc,double pg,double t,long /*ipp*/)
{
  double s,phi,rhocp,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos;
  state_eq tt;
  
  s = sat(pc,t);
  phi = o30baz_phi(t);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  cps = o30baz_cps(t);
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();
  rhos = o30baz_rhos(t);

  rhocp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));

  return(rhocp);
}


/**
   function computes structure factor
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval f - structure factor
*/ 
double o30bazmat::o30baz_fs(double /*pc*/,double /*pg*/,double t)
{
  double f,fs;
  
  fs = 1.0e-3;
  f  = fs*exp(-log(10.0)*(t-298.15)/500.0);
  
  return(f);
}


/**
   function computes tortuosity factor, Formulation by Baroghel 

   @retval tau - tortuosity factor
*/ 
double o30bazmat::o30baz_tau()
{
  double tau;
    
  tau = 0.5;
  
  return(tau);
}

/**
   function computes solid thermal conductivity, experimental data by Kalifa
   (report BRITE 1997) - data valid up to 600°C

   @param pc - capillary pressure
   @param t - temperature

   @retval lambdas - solid thermal conductivity
*/
double o30bazmat::o30baz_lambdas(double t)
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
   @param t - temperature
   
   @retval lambdaeff - effective thermal conductivity of partially saturated concrete
*/
double o30bazmat::o30baz_lambdaeff(double pc,double /*pg*/,double t)
{
  double lambdaeff,lambdas,s,phi,rhow,rhos;
  state_eq tt;

  s = sat(pc,t);
  phi = o30baz_phi(t);
  rhow = tt.get_rhow(t);
  lambdas = o30baz_lambdas(t);
  rhos = o30baz_rhos(t);

  lambdaeff = lambdas*(1.0 + 4.0*phi*rhow*s/(1.0-phi)/rhos);

  return(lambdaeff);
}

/**
   function computes  volume density of concrete skeleton
   @param t - temperature

   @retval rhos - volume density of concrete skeleton
*/
double o30bazmat::o30baz_rhos(double t)
{
  double rhos;

  if (t < tcr)
    rhos = 2590.0-0.529*(t-298.15);
  else                                                              
    rhos = 2590.0-0.529*(tcr-298.15);

  return(rhos);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval beats - cubic thermal expansion coefficient of solid (K-1)
*/
double o30bazmat::o30baz_betas()
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
double o30bazmat::o30baz_hydw(double /*pc*/,double /*pg*/,double t)
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
double o30bazmat::o30baz_dehydw_dt(double /*pc*/,double /*pg*/,double t)
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
double o30bazmat::o30baz_hydren(double /*pc*/,double /*pg*/,double /*t*/)
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
double o30bazmat::o30baz_fste(double /*pc*/,double /*pg*/,double /*t*/)
{
  return(fste);
}


/**
   function computes diffusivity of bound water
   @param pc - capillary pressure
   @param t - temperature

   @retval  - diffusivity of bound water - according to Frotran code
*/
double o30bazmat::o30baz_ddbw(double /*pc*/,double /*pg*/,double t)
{
  double ddbw;
  
  ddbw = ddbw0*exp(-t/(273.15+23));
  if (t > tcr)
    ddbw = ddbw0*exp(-tcr/(273.15+23));

  return(ddbw);
}





/**
   function computes emod Young's modulus

   @retval emod - Young's modulus
*/
double o30bazmat::o30baz_emod()
{
  return(emod0);
}


/**
   function computes fct tensile strenght  (Brite data - Felicetti 1999)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval fct - tensile strenght
*/
double o30bazmat::o30baz_fct(double /*pc*/,double /*pg*/,double t)
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
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval xk0 - maximum linear elastic tensile strain
*/
double o30bazmat::o30baz_xk0(double /*pc*/,double /*pg*/,double t)
{
  double xk0;
  
  xk0 = 1.0e-4*(1.0 + (t-t00)/600.0);

  return(xk0);
}


/**
   function computes bcc compressive coefficient

   @retval bcc - compressive coefficient
*/
double o30bazmat::o30baz_bcc()
{
  double bcc;
  
  bcc = 2000.0;

  return(bcc);
}


/**
   function returns Biot's constant

   @retval alpha - Biot's constant
*/
double o30bazmat::o30baz_alpha()
{
  return(alpha);
}

/**
   function returns Poisson's coefficient

   @retval vcoeff - Poisson's coefficient
*/
double o30bazmat::o30baz_nu()
{
  return(vcoeff);
}



/**
   function reads parameters

   @param - input file
*/
void o30bazmat::read(XFILE */*in*/)
{}


/**
   function prints parameters

   @param out - output file
*/
void o30bazmat::print(FILE */*out*/)
{}



/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   29. 11. 2013, by TKo
*/
void o30bazmat::give_reqntq(long *antq)
{
  //  damage parameter
  antq[scal_iso_damage-1] = 1;
}
