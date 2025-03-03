/*
  File:             C60_baroghelc.cpp
  Author:           Tomas Krejci, 1.12.2003
  Purpose:          Calculates properties of C60 concrete (High Performance Concrete)
                    Saturation assuming Baroghel formulation extended for high temperature
  sources:          C60sifel.f90 and SATBGMsifel.f90 from Padova
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "C60baroghelc.h"
#include "globalt.h"

C60barmatc::C60barmatc()
{
  mw = 18.01528;   //molar mass of water kg.mol-1
  ma = 28.9645;    //molar mass of dry air kg.mol-1
  gasr = 8314.41;  //universal gas constant J.mol-1.K-1
  
  t0 = 273.15;   //reference temperature
  p0 = 101325.0; //pressure
  tcr = 647.3;   //critical point of water
  
  w1 = 163.0;    //free water content at 20°C for Bazant isotherms
  c1 = 450.0;    //cementconctent (for Bazant isotherms)
  
  // porosity
  phi0 = 0.0825468;
  aphi = 0.00017417;
  
  //intrinsic permeability
  k0 = 2.0e-18;
  ak = 6.29e-3;
  
  // gas relative permeability
  scr = 1.0;
  ag = 1.0; //<1;3>
  
  // water relative permeability
  sir = 0.2;
  aw = 2.0; //<1;3>
  bw = 6.0;//or 16.0
  
  //skeleton density, experimental data by Kalifa(report BRITE 1997)
  rhos = 2564.0;
  
  fs = 1.0;//structure coefficient

  //thermal conductivity of solid skeleton
  lambdas0 = 1.9215152;//W/(m.K) 
  alam = -0.0012525253;//K-1
  
  //thermal capacity of solid skeleton
  ac = -0.22626263;
  cps0 = 855.25758;//thermal capacity of the solid skeleton at reference temperature 
  
  //Hydration energy =0.5 MJ/kg
  hydren = 0.5e+6;
  //finv= aging factor (hydration degree)
  finv = 0.65;
  //fste= Water/Cement ratio (data from Brite 1997)
  fste = 0.36;

  //chracteristic length
  dld = 8.0E-2;

  //Young's mosulus
  emod0  = 3.452222e+10;//in Pa at 20°C

  //Poisson's ceofficient
  vcoeff = 0.18;

  //cubic thermal expansion coeffcient
  betas  = 3.0*10.2e-6;

  //Biot's constant
  alpha = 0.5;//rigid porous materials

  //MAZAR'S COEFFICIENTS
  //tensile coefficients
  //functions not completed??!!
  at = 1.0;
  bt = 15000.0;
  //compressive coefficients
  acc = 1.876;

  ddbw0 = 1.0e-20; //diffusion of bound water at refference temperature
}
C60barmatc::~C60barmatc()
{}

/**
   function computes degree of saturation(desorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double C60barmatc::sat(double pc,double t)
{
  double sw,n,a,b,tt,q0,q2,q3,z,e,e0,g;

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

  if (t <= tcr)
    e = pow(((tcr - t0)/(tcr - t)),n); 
  else{
    e0=pow(((tcr-293.15)/z),n);
    e = n/z*e0*t + e0 - n/z*e0*(tcr - z);
  }

  g = pow((e/a*pc),(b/(b-1.0)));
  sw = pow((g + 1.0),(-1.0/b));

  //fprintf(Outt,"sw = %e\n",sw);

  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double C60barmatc::dsat_dpc(double pc,double t)
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

  //fprintf(Outt,"dsw_dpc = %e\n",dsw_dpc);

  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double C60barmatc::dsat_dt(double pc,double t)
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

  //fprintf(Outt,"dsw_dt = %e\n",dsw_dt);

  return(dsw_dt);
}

/**
   function returns saturation solid point

   @retval ssp - saturation solid point
*/
double C60barmatc::ssp()
{
  return(0.55);
}

/**
   function computes porosity, Data by ALONSO-ANDRADE 
   (report BRITE 1997)- results valid up to 600°C
   
   @param t - temperature
   @retval phi - porosity
*/
double C60barmatc::C60bar_phi(double t)
{
  double phi;
  
  phi = phi0 + aphi*(t - t0);
  
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
double C60barmatc::C60bar_kintr(double /*pc*/,double pg,double t,double /*dam*/)
{
  double kintr,bk;//ad;
  state_eq tt;

  bk = log(5.0/3.0)/log(4.0);
  
  kintr = k0*exp(log(10.0)*ak*(t-t0))*pow((pg/p0),bk);
  
  //damage effect - from Francesco Pesavento PhD thesis page 215
  //ad = 4.0;
  //kintr = kintr + k0*pow(10.0,(ad*dam));//added eq. 5.75

  //another formula eq. 5.76
  //kintr = k0*pow(10.0,(at*(t-t0)))*pow((pg/p0),ap)*pow(10.0,(ad*dam));

  /*  printf("\n");
      printf("t = %e\n",t);
      printf("kintr = %e\n",kintr); 
  */ 
  return(kintr);
}

/**
   function computes gas relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krg - gas relative permeability
*/
double C60barmatc::C60bar_krg(double pc,double t)
{
  double krg,s;
  
  s = sat(pc,t);

  krg = 1.0 - pow((s/scr),ag);

  return(krg);
}

/**
   function computes water relative permeability
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param rh - relative humidity

   @retval krw - water relative permeability
*/
double C60barmatc::C60bar_krw(double pc,double t,double rh)
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
   @param pc - capillary pressure
   @param t - temperature

   @retval dd - ...
*/
double C60barmatc::C60bar_dd(double pc,double t)
{
  double tau,dd;

  tau = C60bar_tau(pc,t);
  
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
double C60barmatc::C60bar_deff(double pc,double pg,double t)
{
  //from Francesco Pesavento PhD thesis page 203
  double deff,dd,phi,s,cdiff;
  state_eq tt;

  dd = C60bar_dd(pc,t);
  phi = C60bar_phi(t);
  s = sat(pc,t);
  cdiff = tt.get_cdiff(pc,pg,t);

  deff = dd*phi*(1.0 - s)*fs*cdiff;

  return(deff);
}

/**
   function computes specific heat of solid skeleton
   @param t - temperature
   
   @retval cps - specific heat capacity of solid skeleton 
*/
double C60barmatc::C60bar_cps(double t)
{
  double cps;
  
  if ((t-t0) < 600.0)
    cps = cps0 + ac*(t - t0);
  else
    cps = cps0 + ac*(600.0);
  
  return(cps);
}

/**
   function computes thermal capacity of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   
   @retval rhocp -  of thermal capacity partially saturated concrete
*/
double C60barmatc::C60bar_rhocp(double pc,double pg,double t)
{
  double s,phi,rhocp,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos;
  state_eq tt;
  
  s = sat(pc,t);
  phi = C60bar_phi(t);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  cps = C60bar_cps(t);
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();
  rhos = C60bar_rhos();

  rhocp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));

  return(rhocp);
}

/**
   function computes tortuosity factor, Formulation by Baroghel 
   @param pc - capillary pressure
   @param t - temperature

   @retval tau - tortuosity factor
*/ 
double C60barmatc::C60bar_tau(double pc,double t)
{
  double phi,s,tau;
  
  s = sat(pc,t);
  phi = C60bar_phi(t);
  
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
double C60barmatc::C60bar_lambdas(double t)
{
  double lambdas;

  if ((t-t0) < 600.0)
    lambdas = lambdas0 + alam*(t-t0);
  else
    lambdas = lambdas0 + alam*(600.0);

  return(lambdas);
}

/**
   function computes effective thermal conductivity of partially saturated concrete
   @param pc - capillary pressure
   @param t - temperature
   
   @retval lambdaeff - effective thermal conductivity of partially saturated concrete
*/
double C60barmatc::C60bar_lambdaeff(double pc,double /*pg*/,double t)
{
  double lambdaeff,lambdas,s,phi,rhow,rhos;
  state_eq tt;

  s = sat(pc,t);
  phi = C60bar_phi(t);
  rhow = tt.get_rhow(t);
  lambdas = C60bar_lambdas(t);
  rhos = C60bar_rhos();

  lambdaeff = lambdas*(1.0 + 4.0*phi*rhow*s/(1.0-phi)/rhos);

  return(lambdaeff);
}

/**
   function computes  volume density of concrete skeleton

   @retval rhos - volume density of concrete skeleton
*/
double C60barmatc::C60bar_rhos()
{
  return(rhos);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval beats - cubic thermal expansion coefficient of solid (K-1)
*/
double C60barmatc::C60bar_betas()
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
double C60barmatc::C60bar_hydw(double /*pc*/,double /*pg*/,double t)
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
double C60barmatc::C60bar_dehydw_dt(double /*pc*/,double /*pg*/,double t)
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
double C60barmatc::C60bar_hydren(double /*pc*/,double /*pg*/,double /*t*/)
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
double C60barmatc::C60bar_fste(double /*pc*/,double /*pg*/,double /*t*/)
{
  return(fste);
}


/**
   function computes diffusivity of bound water
   @param pc - capillary pressure
   @param t - temperature

   @retval  - diffusivity of bound water - according to Frotran code
*/
double C60barmatc::C60bar_ddbw(double /*pc*/,double /*pg*/,double t)
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
double C60barmatc::C60bar_emod(double /*pc*/,double /*pg*/,double t)
{
  double emod,ttc;
  
  ttc   = t - 273.15;
  
  if (ttc <= 500.0)
    emod   =  3.5604e10 - 54.089*1.0e6*ttc;
  else
    emod   = 0.8560e10*exp((-54.089/8560.0)*(ttc-500.0));

  return(emod);
}


/**
   function computes fct tensile strenght  (Brite data - Felicetti 1999)
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval fct - tensile strenght
*/
double C60barmatc::C60bar_fct(double /*pc*/,double /*pg*/,double t)
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
double C60barmatc::C60bar_xk0(double pc,double pg,double t)
{
  double fct,emod,xk0;
  
  fct = C60bar_fct(pc,pg,t);
  emod = C60bar_emod(pc,pg,t);

  xk0 = fct/emod;

  return(xk0);
}


/**
   function computes bcc compressive coefficient
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval bcc - compressive coefficient
*/
double C60barmatc::C60bar_bcc(double /*pc*/,double /*pg*/,double t)
{
  double bcc,ttc;
  
  ttc   = t - 273.15;

  bcc = 1550.0*exp(-0.0015154*ttc);

  return(bcc);
}

/**
   function returns Biot's constant

   @retval alpha - Biot's constant
*/
double C60barmatc::C60bar_alpha()
{
  return(alpha);
}


/**
   function computes nu Poisson's constant

   @retval vcoeff - Poisson's constant
*/
double C60barmatc::C60bar_nu()
{
  return(vcoeff);
}

/**
   function reads parameters

   @param in - input file
*/
void C60barmatc::read(XFILE */*in*/)
{}
