/*
  File:             C30_bazant.cpp
  Author:           Tomas Krejci, 1.12.2003
  Purpose:          Calculates the properties of C30 concrete (High Performance Concrete)
                    Saturation assuming Bazant formulation extended for high temperature
  sources:          C30sifel.f90 and SATBGMsifel.f90 from Padova
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "C30bazant.h"
#include "globalt.h"
#include "globmatt.h"

C30bazmat::C30bazmat()
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
  scr = 0.8;
  ag = 2.0; //<1;3>
  
  // water relative permeability
  sir = 0.05;
  aw = 3.0; //<1;3>
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

  ddbw0 = 1.0e-20;//diffusion of bound water at refference temperature

  //chracteristic length
  dld = 8.0E-2;

  //Young's modulus
  emod0  = 3.452222e+10;//in Pa at 20°C

  //Poisson's ceofficient
  vcoeff = 0.18;

  //cubic thermal expansion coeffcient
  betas  = 3.0*10.0e-6;

  //Biot's constant
  alpha = 0.5;//rigid porous materials

  //MAZAR'S COEFFICIENTS
  //tensile coefficients
  at = 1.0;
  bt = 15000.0;
  //compressive coefficients
  acc = 1.876;
}
C30bazmat::~C30bazmat()
{}

/**
   function computes degree of saturation(desorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double C30bazmat::sat(double pc,double t)
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
  else
    mt = 1.0e-4;  //Bazant parameter (above critical point)
  
  //conversion pow(a,b) -> exp(b*log(a))
  sw = pow(rh,(1.0/mt));//saturation calculation
  //sw = exp(log(rh)/mt);//saturation calculation
  
  if(sw < 1.0e-3)//limit for saturation level
    {
      sw=1.0e-3;
      fprintf (Outt,"\n\n Uprava saturace");
    }

  /*
    fprintf (Outt,"\n\n");
    fprintf (Outt,"pc = %lf\n",pc);
    fprintf (Outt,"t  = %lf\n",t);
    fprintf (Outt,"rh = %lf\n",rh);
    fprintf (Outt,"s  = %lf\n",sw);
    fprintf (Outt,"mt = %lf\n",mt);
  */
  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double C30bazmat::dsat_dpc(double pc,double t)
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

  //function SATUR_A.f90
  tt1 = ((tem+10.0)/(25.0+10.0))*((tem+10.0)/(25.0+10.0));//modified temperature from Bazant formulation
  mt1 = 1.04 - (tt1/(24.0+tt1));//function m from Bazant formulation
  if(t < tcr){
    fcm = (tcr-t)/(tcr-298.15);
    mt  = mt1*fcm;//bazant parameter (below critical point)
    if(mt<1.0e-4)  mt=1.0e-4;
  }
  else{
    mt = 1.0e-4;  //Bazant parameter (above critical point)
  }
  
  //conversion pow(a,b) -> exp(b*log(a))
  sw = pow(rh,(1.0/mt));//saturation calculation
  //sw = exp(1.0/mt*log(rh));//saturation calculation
  //sw = sat(pc,t);
  
  if(sw < 1.0e-3)//limit for saturation level
    {
      sw=1.0e-3;
      fprintf (Outt,"\n\n Uprava saturace");
    }

  ssmh = sw/mt/rh;
  dsdc = ssmh*drhdc;
  
  //if(sw <= 1.0e-10)//limit for saturation level
  //dsdc=0.0;
  
  dsw_dpc = dsdc;  

  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double C30bazmat::dsat_dt(double pc,double t)
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
    if(mt<1.0e-4)  mt=1.0e-4;
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
  sw = pow(rh,(1.0/mt));//saturation calculation
  //sw = exp(1.0/mt*log(rh));//saturation calculation
  
  lnh  = log(rh);
  
  if(sw < 1.0e-3)//limit for saturation level
    {
      sw=1.0e-3;
      fprintf (Outt,"\n\n Uprava saturace");
    }
  
  
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
  
  //if(sw <= 1.0e-10)//limit for saturation level
  //dsdt=0.0;
  
  dsw_dt = dsdt;

  return(dsw_dt);
}


/**
   function returns saturation solid point

   @retval ssp - saturation solid point
*/
double C30bazmat::ssp()
{
  return(0.55);
}



/**
   function computes porosity
   
   @retval phi - porosity
*/
double C30bazmat::C30baz_phi()
{
  double phi;
  
  phi = phi0;
  
  return(phi);
}

/**
   function computes intrinsic permeability

   @retval kintr - intrinsic permeability
*/
double C30bazmat::C30baz_kintr()
{
  double kintr;
  
  kintr = k0;
  
  return(kintr);
}

/**
   function computes gas relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krg - gas relative permeability
*/
double C30bazmat::C30baz_krg(double pc,double t)
{
  double krg,s;
  
  s = sat(pc,t);

  //krg = 1.0 - pow((s/scr),ag);

  krg = 1.0 - s;

  return(krg);
}

/**
   function computes water relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krw - water relative permeability
*/
double C30bazmat::C30baz_krw(double pc,double t)
{
  double krw,s;

  s = sat(pc,t);

  //conversion pow(a,b) -> exp(b*log(a))
  //krw = pow(s,10.0);
  krw = exp(10.0*log(s));


  return(krw);
}

/**
   function computes dd
   @param pc - capillary pressure
   @param t - temperature

   @retval dd - ...
*/
double C30bazmat::C30baz_dd(double pc,double t)
{
  double tau,dd;

  tau = C30baz_tau(pc,t);
  
  dd = tau*mw*ma/gasr;//according to Fortran code
  //dd = tau;//according to Francesco's PhD thesis?

  return(dd);
}

/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval deff - effective diffusion coefficient of vapour inside pores
*/
double C30bazmat::C30baz_deff(double pc,double pg,double t)
{
  double deff,dd,phi,s,cdiff;
  state_eq tt;

  dd = C30baz_dd(pc,t);
  phi = C30baz_phi();
  s = sat(pc,t);
  cdiff = tt.get_cdiff(pc,pg,t);

  deff = dd*phi*(1.0 - s)*fs*cdiff;
  
  return(deff);
}

/**
   function computes specific heat of solid skeleton
   
   @retval cps - specific heat capacity of solid skeleton 
*/
double C30bazmat::C30baz_cps()
{
  double cps;
  
  cps = cps0;
  
  return(cps);
}

/**
   function computes thermal capacity of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   
   @retval rhocp -  of thermal capacity partially saturated concrete
*/
double C30bazmat::C30baz_rhocp(double pc,double pg,double t)
{
  double s,phi,rhocp,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos;
  state_eq tt;
  
  s = sat(pc,t);
  phi = C30baz_phi();
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  cps = C30baz_cps();
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();
  rhos = C30baz_rhos();

  rhocp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));

  return(rhocp);
}

/**
   function computes tortuosity factor, Formulation by Baroghel 
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval tau - tortuosity factor
*/ 
double C30bazmat::C30baz_tau(double pc,double t)
{
  double phi,s,tau;
  
  s = sat(pc,t);
  phi = C30baz_phi();
  
  //conversion pow(a,b) -> exp(b*log(a))
  tau = exp(1.0/3.0*log(phi))*exp(7.0/3.0*log(1.0-s));
  //tau = pow(phi,(1.0/3.0))*pow((1.0-s),(7.0/3.0));
  
  return(tau);
}

/**
   function computes solid thermal conductivity

   @retval lambdas - solid thermal conductivity
*/
double C30bazmat::C30baz_lambdas()
{
  double lambdas;

  lambdas = lambdas0;

  return(lambdas);
}

/**
   function computes effective thermal conductivity of partially saturated concrete
   @param t - temperature
   
   @retval lambdaeff - effective thermal conductivity of partially saturated concrete
*/
double C30bazmat::C30baz_lambdaeff(double pc,double t)
{
  double lambdaeff,lambdas,s,phi,rhow,rhos;
  state_eq tt;

  s = sat(pc,t);
  phi = C30baz_phi();
  rhow = tt.get_rhow(t);
  lambdas = C30baz_lambdas();
  rhos = C30baz_rhos();

  lambdaeff = lambdas*(1.0 + 4.0*phi*rhow*s/(1.0-phi)/rhos);
  
  return(lambdaeff);
}

/**
   function computes  volume density of concrete skeleton

   @retval rhos - volume density of concrete skeleton
*/
double C30bazmat::C30baz_rhos()
{
  return(rhos);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval beats - cubic thermal expansion coefficient of solid (K-1)
*/
double C30bazmat::C30baz_betas()
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
double C30bazmat::C30baz_hydw(double /*pc*/,double /*pg*/,double t)
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
double C30bazmat::C30baz_dehydw_dt(double /*pc*/,double /*pg*/,double t)
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
double C30bazmat::C30baz_hydren(double /*pc*/,double /*pg*/,double /*t*/)
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
double C30bazmat::C30baz_fste(double /*pc*/,double /*pg*/,double /*t*/)
{
  return(fste);
}


/**
   function computes diffusivity of bound water
   @param pc - capillary pressure
   @param t - temperature

   @retval  - diffusivity of bound water
*/
double C30bazmat::C30baz_ddbw(double /*pc*/,double /*pg*/,double t)
{
  double ddbw;
  
  ddbw = ddbw0*exp(-1.0*t/(295.0));
  if (t > tcr)
    ddbw = ddbw0*exp(-1.0*tcr/(295.0));

  return(ddbw);
}


/**
   function computes emod Young's modulus
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval emod - Young's modulus
*/
double C30bazmat::C30baz_emod(double /*pc*/,double /*pg*/,double t)
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
double C30bazmat::C30baz_fct(double /*pc*/,double /*pg*/,double t)
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
double C30bazmat::C30baz_xk0(double pc,double pg,double t)
{
  double fct,emod,xk0;
  
  fct = C30baz_fct(pc,pg,t);
  emod = C30baz_emod(pc,pg,t);

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
double C30bazmat::C30baz_bcc(double /*pc*/,double /*pg*/,double t)
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
double C30bazmat::C30baz_alpha()
{
  return(alpha);
}


/**
   function computes nu Poisson's constant

   @retval vcoeff - Poisson's constant
*/
double C30bazmat::C30baz_nu()
{
  return(vcoeff);
}

/**
   function reads parameters

   @param in - input file
*/
void C30bazmat::read(XFILE */*in*/)
{}


/**
   function prints parameters

   @param out - output file
*/
void C30bazmat::print(FILE */*out*/)
{}


/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Krejci according to Tomas Koudelka, 11/01/2023
*/
void C30bazmat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_capillary_press;
  dofname[1] = trf_gas_press;
  dofname[2] = trf_temperature;
}
