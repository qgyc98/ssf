/*
    File:             baroghelB.cpp
    Author:           Tomas Krejci, 1.12.2003
    Purpose:          material properties for BO and BH concrete by Baroghel at normal temperatures - for testing
    sources:          Cor1sifel.f90 and SATBGMsifel.f90 from Padova
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "constrel.h"
#include "baroghelB.h"
#include "globalt.h"
#include "globmatt.h"

baroghelmat::baroghelmat()
{
  mw = 18.01528;   //molar mass of water kg.mol-1
  ma = 28.9645;    //molar mass of dry air kg.mol-1
  gasr = 8314.41;  //universal gas constant J.mol-1.K-1
  
  t0 = 273.15;   //reference temperature
  p0 = 101325.0; //pressure
  tcr = 647.3;   //critical point of water
  
  // gas relative permeability
  scr = 1.0;//provisionally
  ag = 1.0; //<1;3>
  
  //EFFECTIVE DIFFUSION COEFFICIENT OF VAPOUR
  av = 1.0;//<1;3>
  fs = 1.0;//structure coefficient
  
  //cubic thermal expansion coeffcient
  betas  = 3.0*10.2e-6;

  //Biot's constant
  alpha = 0.5;//rigid porous materials

  /*   parameters from input file
       //DEGREE OF SATURATION
       //parametres:
       abo = 18.62*1000000.0;
       abh = 46.93*1000000.0;
       bbo = 2.27;
       bbh = 2.06;
       
       //INTRINSIC PERMEABILITY
       bo = 3.0e-21;
       bh = 5.22e-22;
       co = 1.0e-20;
       ch = 4.0e-21;  
       
       //THERMAL CONDUCTIVITY
       lambdabo = 1.67;
       lambdabh = 2.0;
       
       //SPECIFIC HEAT
       cpbo = 940.0;
       cpbh = 860.0;
       
       //DRY DENSITY (APPARENT)
       rhosbo = 2286.0;
       rhosbh = 2385.0;
       
       //POROSITY
       phibo = 0.122;
       phibh = 0.082;
  */

  //free water content at 20°C
  w1 = 163.0;
  //cement content
  c1 = 450.0;
  //Hydration energy =0.5 MJ/kg
  hydren = 0.0;//hydren = 0.5e+6;??!!
  //finv= aging factor (hydration degree)
  finv = 0.65;
  //fste= Water/Cement ratio (data from Brite 1997)
  fste = 0.36;
  
  ddbw0 = 1.0e-20; //diffusion of bound water at refference temperature
}

baroghelmat::~baroghelmat()
{}

/**
   function computes degree of saturation(sorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double baroghelmat::baroghel_sw(double pc,double /*t*/)
{
  double sw,k,n,m,help,a,b;
  
  a = ab;
  b = bb;
  
  m = 1.0/(1.0-1.0/b);
  n = -1.0/b;
  k = 1.0/a;

  //conversion pow(a,b) -> exp(b*log(a))
  //help = pow((k*pc),m);
  help = exp(m*log(k*pc));

  //conversion pow(a,b) -> exp(b*log(a))
  //sw = pow((1.0 + help),n);
  sw = exp(n*log(1.0+help));

  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dpc - partial derivative of degree of saturation with respect to pc
*/
double baroghelmat::baroghel_dsw_dpc(double pc,double /*t*/)
{
  double dsw_dpc;
  double sw1,k,n,m,help,help2,a,b;
  
  a = ab;
  b = bb;
  
  m = 1.0/(1.0-1.0/b);
  n = -1.0/b;
  k = 1.0/a;
  help = pow((k*pc),m);
  sw1 = n*(pow((1.0 + help),(n-1.0)));
  help2 = m*(pow((k*pc),(m-1.0)))*k;
  dsw_dpc = sw1*help2;

  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double baroghelmat::baroghel_dsw_dt(double /*pc*/, double /*t*/)
{
  double dsw_dt;
    
  dsw_dt = 0.0;

  return(dsw_dt);
}

/**
   function returns saturation solid point

   @retval ssp - saturation solid point
*/
double baroghelmat::baroghel_ssp()
{
  return(0.55);
}

/**
   function computes gas relative permeability
   @param s - degree of saturation

   @retval krg - gas relative permeability
*/
double baroghelmat::baroghel_krg(double s)
{
  double krg;

  if (s < scr)
    krg = 1.0 - pow((s/scr),ag);
  else
    krg = 0.0;
  
  return(krg);
}

/**
   function computes water relative permeability
   @param pc - capillary pressure
   @param t - temperature
   
   @retval krw - water relative permeability
*/
double baroghelmat::baroghel_krw(double pc, double t)
{
  double krw,s;
  
  s = baroghel_sw(pc,t);
  krw = pow(s,10.0);
  
  return(krw);
}

/**
   function computes porosity

   @retval phi - porosity
*/
double baroghelmat::baroghel_phi()
{
  double phi;

  phi = phib;

  return(phi);
}

/**
   function computes intrinsic permeability

   @retval kintr - intrinsic permeability
*/
double baroghelmat::baroghel_kintr()
{
  double kintr;

  kintr = c;

  return(kintr);
}

/**
   function computes specific heat of solid skeleton 

   @retval cps - specific heat of solid skeleton
*/
double baroghelmat::baroghel_cps()
{
  double cps;

  cps = cpb;

  return(cps);
}

/**
   function computes thermal capacity of partially saturated concrete
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhocp - thermal capacity of partially saturated concrete
*/
double baroghelmat::baroghel_rhocp(double pc,double pg,double t)
{
  double s,pgw,phi,rhocp,rhow,cps,cpw,rhocpg,cpgw,rhos;
  state_eq tt;

  s = baroghel_sw(pc,t);
  phi = baroghel_phi();
  pgw = tt.get_pgw(pc,t);
  rhow = tt.get_rhow(t);
  cps = baroghel_cps();
  cpw = tt.get_cpw();
  rhocpg = tt.get_rhocpg(pc,pg,t);
  cpgw = tt.get_cpgw();
  rhos = baroghel_rhos();

  rhocp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*rhocpg);
  
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
double baroghelmat::baroghel_cp(double pc,double pg,double t,long ipp)
{
  double s,pgw,phi,cp,rho,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos;
  state_eq tt;

  s = baroghel_sw(pc,t);
  phi = baroghel_phi();
  pgw = tt.get_pgw(pc,t);
  rho = tt.get_rho(pc,pg,t,ipp);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);
  rhogw = tt.get_rhogw(pc,t);
  cps = baroghel_cps();
  cpw = tt.get_cpw();
  cpga = tt.get_cpga();
  cpgw = tt.get_cpgw();
  rhos = baroghel_rhos();

  cp = (1.0-phi)*rhos*cps + phi*(s*rhow*cpw + (1.0-s)*(rhog*cpga + rhogw*(cpgw-cpga)));
  cp = cp/rho;

  return(cp);
}

/**
   function computes tortuosity factor
   @param pc - capillary pressure
   @param t - temperature

   @retval tau - tortuosity factor
*/ 
double baroghelmat::baroghel_tau(double pc,double t)
{
  double phi,s,tau;
  
  s = baroghel_sw(pc,t);
  phi = baroghel_phi();
  
  tau = pow(phi,(1.0/3.0))*pow((1.0-s),(7.0/3.0));

  return(tau);
}

/**
   function computes dd
   @param pc - capillary pressure
   @param t - temperature

   @retval dd - ...
*/
double baroghelmat::baroghel_dd(double pc,double t)
{
  double tau,dd;

  tau = baroghel_tau(pc,t);
  
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
double baroghelmat::baroghel_deff(double pc,double pg,double t)
{
  double dd,deff,phi,s,cdiff;
  state_eq tt;

  dd = baroghel_dd(pc,t);
  phi = baroghel_phi();
  s = baroghel_sw(pc,t);
  cdiff = tt.get_cdiff(pc,pg,t);

  deff = dd*phi*(1.0 - s)*fs*cdiff;

  return(deff);
}

/**
   function computes effective thermal conductivity of partially saturated concrete
   @param pc - capillary pressure
   @param t - temperature

   @retval lambdaeff - effective thermal conductivity of partially saturated concrete
*/
double baroghelmat::baroghel_lambdaeff(double pc,double /*pg*/,double t)
{
  double lambdaeff,lambdad,lambdad0,alam,s,phi,rhow,rhos;
  state_eq tt;

  s = baroghel_sw(pc,t);
  phi = baroghel_phi();
  lambdad0 = lambdab;
  alam = 0.000;
  rhow = tt.get_rhow(t);
  rhos = baroghel_rhos();

  lambdad = lambdad0*(1.0 + alam*(t-t0));

  lambdaeff = lambdad*(1.0 + 4.0*phi*rhow*s/(1.0-phi)/rhos);

  return(lambdaeff);
}

/**
   function computes  volume density of concrete skeleton

   @retval rhos - volume density of concrete skeleton
*/
double baroghelmat::baroghel_rhos()
{
  double rhos;

  rhos = rhosb;

  return(rhos);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval betas - cubic thermal expansion coefficient of solid (K-1)
*/
double baroghelmat::baroghel_betas()
{
  return(betas);
}

/**
   function returns Biot's constant

   @retval alpha - Biot's constant
*/
double baroghelmat::baroghel_alpha()
{
  return(alpha);
}



/**
   function computes hydration degree
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval hydw - hydration degree
*/
double baroghelmat::baroghel_hydw(double /*pc*/,double /*pg*/,double t)
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
double baroghelmat::baroghel_dehydw_dt(double /*pc*/,double /*pg*/,double t)
{
  double dfhyt,dehydw_dt;
  
  if((t-t0)< 105.0){
    dfhyt = 0.0;
  }
  else{
    dfhyt = (3.1416*0.004/2.0)*cos(3.1416/2.0*(1.0-2.0*exp(-0.004*((t-t0)-105.0))))*exp(-0.004*((t-t0)-105.0));
  }
  
  dehydw_dt = fste*finv*c1*dfhyt;

  dehydw_dt = 0.0;//!!??

  return(dehydw_dt);
}


/**
   function computes hydration energy
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval hydren - hydration energy
*/
double baroghelmat::baroghel_hydren(double /*pc*/,double /*pg*/,double /*t*/)
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
double baroghelmat::baroghel_fste(double /*pc*/,double /*pg*/,double /*t*/)
{
  return(fste);
}


/**
   function computes diffusivity of bound water
   @param pc - capillary pressure
   @param t - temperature

   @retval  - diffusivity of bound water - according to Frotran code
*/
double baroghelmat::baroghel_ddbw(double /*pc*/,double /*pg*/,double t)
{
  double ddbw;
  
  ddbw = ddbw0*exp(-t/(273.15+23));
  if (t > tcr)
    ddbw = ddbw0*exp(-tcr/(273.15+23));

  ddbw = 0.0;//!!??

  return(ddbw);
}

/**
   function computes emod Young's modulus

   @retval emod - Young's modulus
*/
double baroghelmat::baroghel_emod()
{
  return(emod);
}

/**
   function computes nu Poisson's constant

   @retval nu - Poisson's constant
*/
double baroghelmat::baroghel_nu()
{
  return(nu);
}


/**
   function reads parameters
   
   @param in - input file
*/
void baroghelmat::read(XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &ab, &bb, &c, &lambdab, &cpb, &rhosb, &phib, &emod, &nu);
}


/**
   function prints parameters
   
   @param out - output file
*/
void baroghelmat::print(FILE *out)
{
  fprintf (out,"  %e %e %e %e %e %e %e %e %e", ab, bb, c, lambdab, cpb, rhosb, phib, emod, nu);
}




/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Krejci according to Tomas Koudelka, 09/12/2022
*/
void baroghelmat::give_dof_names(namevart *dofname, long ntm)
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
