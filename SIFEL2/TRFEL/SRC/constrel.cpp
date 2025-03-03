/*
    File:             constrel.cpp
    Author:           Tomas Krejci,  1.12.2002
    Purpose:          constitutive relations for thermo-hydro-mechanical problem

    FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
    --------------------------------------------------
    sources: 
    1. FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS BY AN ALGEBRAIC MULTIGRID METHOD
    Wang Xicheng, B.A. Schrefler
    
    2. NUMERICAL ANALYSIS OF HYGRO-THERMAL BEHAVIOUR AND DAMAGE OF CONCRETE AT HIGH TEMPERATURE
    D. Gawin, C.E. Majorana, B.A. Schrefler

    3. NONLINEAR MODELLING OF CONCRETE AS MULTIPHASE POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
    Francesco Pesavento - doctoral thesis
    
    LIST OF USED VARIABLES:
    -----------------------
    
    cp    ... effective specific heat of porous medium (J.kg-1.K-1)    
    cpg   ... effective specific heat of gaz mixture (J.kg-1.K-1)    
    cpw   ... effective specific heat of liquid (water) (J.kg-1.K-1)    
    cps   ... effective specific heat of solid matrix (J.kg-1.K-1)    
    ct    ... tangent matrix (Pa)
    deff  ... effective diffusivity of gas mixture (m2.s-1)
    g     ... acceleration of gravity (m.s-2)
    kh    ... hydraulic conductivity
    kap   ... absolute permeability (m2)
    krg   ... relative permeability of gas phase
    krw   ... relative permeability of liquid phase
    ks    ... bulk modulus of solid phase
    kt    ... bulk modulus of porous medium
    m     ... molar mass of gas mixture (moist air) (kg.kmol-1)
    ma    ... molar mass of dry air (kg.kmol-1)
    mw    ... molar mass of water vapour (kg.kmol-1)
    gasr  ... universal gas constant 8.3144e3 J.mol-1.K-1
    pav   ... average pressure of mixture (Pa)
    patm  ... atmospheric pressure (Pa)
    pc    ... capillary pressure (Pa)
    pg    ... presure of gas phase (Pa)
    pgw   ... water vapour partial presure (Pa)
    pgws  ... water vapour saturation presure (Pa)
    pw    ... liquid water presure (Pa)
    s     ... liquid phase volumic saturation (=Sw = liquid volume/pore volume)
    t     ... temperature (K)

    alpha ... Biot's constant
    alphac... convective heat transfer coefficient on the boundary
    betac ... mass transfer coefficient on the boundary
    betas ... cubic thermal expansion coefficient of solid (K-1)
    betaswg ... volume thermal dilatation of solid - water vapor for incompressible solid grains
    betasw ... volume thermal dilatation of solid - water for incompressible solid grains
    betaw ... volume thermal dilatation of water
    dhvap ... enthalpy of vaporization per unit mass (J.kg-1)
    eps   ... strain tensor
    epsv  ... volumetric strain
    epst  ... thermoelastic strain
    lambdaeff ... effective thermal conductivity (W.m-1.K-1)
    rho   ... effective density of porous medium (kg.m-3)
    rhog  ... gas phase density (kg.m-3)
    rhoga ... mass concentration of dry air in gas phase (kg.m-3)
    rhogw ... mass concentration of water vapour in gas phase (kg.m-3)
    rhow  ... liquid phase density (kg.m-3)
    phi   ... porosity (pore volume/total volume)
    mug   ... dynamic viscosity of gas phase (kg.m-1.s-1)
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "globalt.h"
#include "constrel.h"

state_eq::state_eq()
{    
  mw = 18.01528;//molar mass of water kg.mol-1
  ma = 28.9645;//molar mass of dry air kg.mol-1
  gasr = 8314.41;//universal gas constant J.kmol-1.K-1

  t0 = 273.15;
  p0 = 101325.0;

  // PHYSICAL PROPERTIES OF DRY AIR
  muga0 = 17.17e-6;//[Pa.s]
  alphaa = 4.73e-8;//[Pa.s.K-1]
  betaa = 2.222e-11;//[Pa.s.K-2]

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  // from D.Gawin, F.Pesavento (PRVAP.f90)
  //Hyland and Wexler equation (1983) for the saturation water vapor pressure
  dv0 = 2.58e-5; //effective diffusion coefficient of vapour m^2/s at reference temperature 273.15 
  bv = 1.667;
  c8 = -5.8002206e+03;
  c9 = 1.3914993;
  c10 =-4.8640239e-02;
  c11 = 4.1764768e-05;
  c12 = -1.4452093e-08;
  c13 = 6.5459673;
  mugw0 = 8.85e-6;//[Pa.s]
  alphaw = 3.633e-8;//[Pa.s.K-1]
  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // alphaw = 3.53e-8    (eq. 50)
  
  // PHYSICAL PROPERTIES OF WATER
  // from Dariusz Gawin (WATPROP.f90)
  rhow0 = 999.84;
  tcr = 647.3;
  cwat  =  0.2e9;
  betawat = -0.414e-3;
  hvap0 =  2.7e+5;
  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // hvap0 = 2.672e+5    (eq. 49)

  a0 =  4.8863e-7; a1 = -1.6528e-9; a2 =  1.8621e-12;
  a3 =  2.4266e-13; a4 = -1.5996e-15; a5 =  3.3703e-18;
  b0 =  1.0213e3; b1 = -7.7377e-1; b2 =  8.7696e-3;
  // M. Starnoni 10-11-2010
  // Gaw-Pes-Sch "Modelling of hygro-thermal behaviour and damage of concrete at T above the cr point of water"
  // b0 = 1.02e-3    (eq. 35)
  
  b3 = -9.2118e-5; b4 =  3.3534e-7; b5 = -4.4034e-10;
  pr1 = 1.0e7; prif = 2.0e7;
  muw0 =  0.6612;
  conb  =  229.0;
  conc  = -1.562;
  cpw = 4181.0;
  lambdaw = 0.6;  
  kw0 = 0.43e10;//compresibility coefficient of water

  //minimal pc value
  pcmin = 100.0;
}
state_eq::~state_eq()
{}

//CAPILLARY PRESSURE
/**
   function computes capillary pressure
   @param pg - capillary gas pressure
   @param pw - capillary water pressure

   @retval pc - capillary pressure
*/
double state_eq::get_pc(double pg,double pw)
{
  double pc;

  pc = pg - pw;

  return(pc);
}

//DEGREE OF SATURATION AND DERIVATIVES
/**
   function computes degree of saturation
   @param pc - capillary pressure
   @param t - temperature

   @retval s - degree of saturation
*/
double state_eq::get_s(double pc,double t,long ipp)
{
  long i;
  double s;
  
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    s = Tm->concrete[i].concreteB_sw(pc,t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    s = Tm->baroghel[i].baroghel_sw(pc,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    s = Tm->C60baroghel[i].sat(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    s = Tm->C30baroghel[i].sat(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    s = Tm->o30bazant[i].sat(pc,t);//function from o30bazant.cpp
    break;
  }
  case C60bazantB:{
    s = Tm->C60bazant[i].sat(pc,t);//function from C60bazant.cpp
    break;
  }
  case C30bazantB:{
    s = Tm->C30bazant[i].sat(pc,t);//function from C30bazant.cpp
    break;
  }

  case soilmat1:{
    s = Tm->soil1[i].sat(pc,t);//function from soil1mat.cpp
    break;
  }  
    

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(s);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double state_eq::get_ds_dpc(double pc,double t,long ipp)
{
  long i;
  double ds_dpc;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    ds_dpc = Tm->concrete[i].concreteB_dsw_dpc(pc,t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    ds_dpc = Tm->baroghel[i].baroghel_dsw_dpc(pc,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    ds_dpc = Tm->C60baroghel[i].dsat_dpc(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    ds_dpc = Tm->C30baroghel[i].dsat_dpc(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    ds_dpc = Tm->o30bazant[i].dsat_dpc(pc,t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantB:{
    ds_dpc = Tm->C60bazant[i].dsat_dpc(pc,t);//function from C60bazant.cpp
    break;
  }    
  case C30bazantB:{
    ds_dpc = Tm->C30bazant[i].dsat_dpc(pc,t);//function from C30bazant.cpp
    break;
  }    

  case soilmat1:{
    ds_dpc = Tm->soil1[i].dsat_dpc(pc,t);//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(ds_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dt - partial derivative of degree of saturation with respect to t
*/
double state_eq::get_ds_dt(double pc,double t,long ipp)
{
  long i;
  double ds_dt;
  
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    ds_dt = Tm->concrete[i].concreteB_dsw_dt(pc,t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    ds_dt = Tm->baroghel[i].baroghel_dsw_dt(pc,t);//function from baroghelB.cpp
    break;
  }    
  case C60baroghelB:{
    ds_dt = Tm->C60baroghel[i].dsat_dt(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    ds_dt = Tm->C30baroghel[i].dsat_dt(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    ds_dt = Tm->o30bazant[i].dsat_dt(pc,t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantB:{
    ds_dt = Tm->C60bazant[i].dsat_dt(pc,t);//function from C60bazant.cpp
    break;
  }    
  case C30bazantB:{
    ds_dt = Tm->C30bazant[i].dsat_dt(pc,t);//function from C30bazant.cpp
    break;
  }    

  case soilmat1:{
    ds_dt = Tm->soil1[i].dsat_dt(pc,t);//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(ds_dt);
}

/**
   function computes saturation solid point
   @param ipp - number of integration point

   @retval ssp - saturation solid point
*/
double state_eq::get_ssp(long ipp)
{
  long i;
  double ssp;
  
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    ssp = Tm->concrete[i].concreteB_ssp();//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    ssp = Tm->baroghel[i].baroghel_ssp();//function from baroghelB.cpp
    break;
  }    
  case C60baroghelB:{
    ssp = Tm->C60baroghel[i].ssp();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    ssp = Tm->C30baroghel[i].ssp();//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    ssp = Tm->o30bazant[i].ssp();//function from o30bazant.cpp
    break;
  }    
  case C60bazantB:{
    ssp = Tm->C60bazant[i].ssp();//function from C60bazant.cpp
    break;
  }    
  case C30bazantB:{
    ssp = Tm->C30bazant[i].ssp();//function from C30bazant.cpp
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(ssp);
}


/**
   function computes diffusivity of bound water
   @param pc - capillary pressure
   @param t - temperature

   @retval  - diffusivity of bound water - according to Frotran code
*/
double state_eq::get_ddbw(double pc,double pg,double t,long ipp)
{
  long i;
  double ddbw;

  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//type of material
    /*  case concreteB:{
	ddbw = Tm->concrete[i].concreteB_sw(pc,t);//function from concreteB.cpp
	break;
	}
    */
  case o30bazantB:{
    ddbw = Tm->o30bazant[i].sat(pc,t);//function from o30bazant.cpp
    break;
  }
  case C30baroghelB:{
    ddbw = Tm->C30baroghel[i].sat(pc,t);//function from C30baroghel.cpp
    break;
  }
  case C60baroghelB:{
    ddbw = Tm->C60baroghel[i].sat(pc,t);//function from C60baroghel.cpp
    break;
  }
  case baroghelB:{
    ddbw = Tm->baroghel[i].baroghel_ddbw(pc,pg,t);//function from baroghelB.cpp
    break;
  }
  case C60bazantB:{
    ddbw = Tm->C60bazant[i].C60baz_ddbw(pc,pg,t);//function from C60bazant.cpp
    break;
  }
  case C30bazantB:{
    ddbw = Tm->C30bazant[i].C30baz_ddbw(pc,pg,t);//function from C30bazant.cpp
    break;
  }
    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
    
  return(ddbw);
}



//GAS - DRY AIR+WATER VAPOR; DENSITY AND DERIVATIVES
/**
   function computes relative humidity of moist air inside the pores = Kelvin-Laplace law
   @param pc - capillary pressure
   @param t - temperature

   @retval rh - relative humidity of moist air inside the pores = Kelvin-Laplace law
*/
double state_eq::get_rh(double pc,double t)
{
  double rhow,rh;

  //critical point of water check
  if(t >= tcr)
    t = tcr;

  rhow = get_rhow(t);
  rh = exp(-1.0*pc/rhow*mw/gasr/t);
    
  //if (rh > 1.0){
  //rh = 1.0;
  //fprintf (Outt,"\n\n Uprava rel. hum > 1.0");
  //}
  if (rh <= 0.001){
    fprintf (Outt,"\n\n Uprava rel. hum < 0.0");
    rh = 0.001;
  }
  return(rh);
}

/**
   function computes capillary pressure from relative humidity = inverse Kelvin-Laplace law
   @param rh - relative humidity
   @param t - temperature

   @retval pc - capillary pressure from relative humidity = inverse Kelvin-Laplace law
*/
double state_eq::get_pcrh(double rh,double t)
{
  double rhow,pc;

  //critical point of water check
  if(t >= tcr)
    t = tcr;
  
  rhow = get_rhow(t);
  pc = -rhow*gasr*t/mw*log(rh);   
  
  
  //minimum capillary pressure check
  //  JK, 25.11.2010
  if (pc < 100.0)
    pc = 1.0e-10;
  
  
  return(pc);
}


/**
   function computes derivative of relative humidity with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval drh_dpc - derivative of relative humidity with respect to pc
*/
double state_eq::get_drh_dpc(double pc,double t)
{
  double drh_dpc,rhow;

  //critical point of water check
  if(t >= tcr)
    t = tcr;

  rhow = get_rhow(t);
  drh_dpc = -1.0*mw/rhow/gasr/t*exp(-1.0*pc*mw/rhow/gasr/t);

  return(drh_dpc);
}

/**
   function computes derivative of relative humidity with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval drh_dt - derivative of relative humidity with respect to t
*/
double state_eq::get_drh_dt(double pc,double t)
{
  double drh_dt,rhow,drhow_dt;

  //critical point of water check
  if(t >= tcr)
    t = tcr;
  
  rhow = get_rhow(t);
  drhow_dt = get_drhow_dt(pc,t);
  drh_dt = (pc/rhow*mw/gasr/t/t+pc/rhow*mw/gasr/t*drhow_dt/rhow)*exp(-1.0*pc/rhow*mw/gasr/t);

  if(t >= tcr)
    drh_dt =0.0;

  return(drh_dt);
}


/**
   function computes water content w [kg/kg]
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval w - computes water content w [kg/kg]
*/
double state_eq::get_w(double pc,double pg,double t,long ipp)
{
  double w,por,s,rhow,rhog,rhos;

  por = get_phi(t,ipp);
  s = get_s(pc,t,ipp);
  rhow = get_rhow(t);
  rhog = get_rhog(pc,pg,t);
  rhos = get_rhos(t,ipp);

  w = (por*s*rhow + por*(1.0 - s)*rhog)/(1.0 - por)/rhos;

  return(w);
}



/**
   function computes gas pressure
   @param pga - capillary dry air pressure
   @param pgw - capillary water vapour pressure
   @param t - temperature

   @retval pg - gas pressure
*/
double state_eq::get_pg(double pga,double pgw,double /*t*/)
{
  double pg;

  pg = pga + pgw;

  return(pg);
}

/**
   function computes gas phase density
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhog - gas phase density
*/
double state_eq::get_rhog(double pc,double pg,double t)
{
  double rhog,pgw;

  pgw = get_pgw(pc,t);

  rhog = (pg*ma + (mw - ma)*pgw)/gasr/t;
  
  return(rhog);
}


/**
   function computes partial derivative of rhog (gas phase density) with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval drhog_dpc - partial derivative of rhog with respect to pc
*/
double state_eq::get_drhog_dpc(double pc,double t)
{
  double drhog_dpc,dpgw_dpc;

  dpgw_dpc = get_dpgw_dpc(pc,t);
  
  drhog_dpc = (mw - ma)*dpgw_dpc/(gasr*t);
  
  return(drhog_dpc);
}

/**
   function computes partial derivative of rhog (gas phase density) with respect to pg
   @param t - temperature

   @retval drhog_dpg - partial derivative of rhog with respect to pg
*/
double state_eq::get_drhog_dpg(double t)
{
  double drhog_dpg;

  drhog_dpg = ma/(gasr*t);
  
  return(drhog_dpg);
}


/**
   function computes partial derivative of rhog (gas phase density) with respect to t
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval drhog_dt - partial derivative of rhog with respect to t
*/
double state_eq::get_drhog_dt(double pc,double pg,double t)
{
  double drhog_dt,dpgw_dt,rhog;

  rhog = get_rhog(pc,pg,t);
  dpgw_dt = get_dpgw_dt(pc,t);

  drhog_dt = (mw - ma)*dpgw_dt/(gasr*t) - rhog/t;
  
  return(drhog_dt);
}

/**
   function computes dynamic viscosity of moist air
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval mug - dynamic viscosity of moist air
*/
double state_eq::get_mug(double pc,double pg,double t)
{
  double mug,mugw,muga,pga,pgw;

  pgw = get_pgw(pc,t);

  pga = pg - pgw;
  muga = get_muga(t);
  mugw = get_mugw(t);

  //gas pressure check
  if(pgw <= pg)
    mug = mugw + (muga - mugw)*pow((1.0 - pgw/pg),0.6083);
  else
    mug = mugw;

  return(mug);
}

/**
   function computes molar mass of moist air
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval mg - molar mass of moist air
*/
double state_eq::get_mg(double pc,double pg,double t)
{
  /*  double rhogw,rhog,rhoga,mg;   
      rhogw = get_rhogw(pc,pg,t);
      rhog = get_rhog(pc,pg,t);
      rhoga = get_rhoga(pc,pg,t);
      
      mg = 1.0/(rhogw/rhog/mw + rhoga/rhog/ma);
  */

  double mg,pgw;

  pgw = get_pgw(pc,t);

  mg = ma + (mw-ma)*pgw/pg;
  
  return(mg);
}


/**
   function computes thermal capacity of air mixture
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhocpg - air mixture thermal capcity
*/
double state_eq::get_rhocpg(double pc,double pg,double t)
{
  double rhocpg,rhog,rhogw,cpga,cpgw;

  rhog = get_rhog(pc,pg,t);
  cpga = get_cpga();
  rhogw = get_rhogw(pc,t);
  cpgw = get_cpgw();

  rhocpg = rhog*cpga + rhogw*(cpgw-cpga);

  return(rhocpg);
}


/**
   function computes air specific heat
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval cpg - air specific heat
*/
double state_eq::get_cpg(double pc,double pg,double t)
{
  double cpg,rhog,rhogw,cpga,cpgw;

  rhog = get_rhog(pc,pg,t);
  cpga = get_cpga();
  rhogw = get_rhogw(pc,t);
  cpgw = get_cpgw();
  
  cpg = (rhog*cpga + rhogw*(cpgw-cpga))/rhog;//zkontrolovat??!!
  
  return(cpg);
}


//DRY AIR; DENSITY AND DERIVATIVES
/**
   function computes air pressure
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval pga - air pressure
*/
double state_eq::get_pga(double pc,double pg,double t)
{
  double rhoga,pga;
  
  rhoga = get_rhoga(pc,pg,t);
  
  pga = rhoga*t*gasr/ma;

  return(pga);
}

/**
   function computes mass concentration of dry air in gas phase
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhoga - mass concentration of dry air in gas phase
*/
double state_eq::get_rhoga(double pc,double pg,double t)
{
  double rhoga,pgw;

  //gas pressure check
  pgw = get_pgw(pc,t);

  if (pgw <= pg)
    rhoga = (pg - pgw)*ma/gasr/t;
  else
    rhoga = 0.0;
    //rhoga = (pg - pgw)*ma/gasr/t;

  return(rhoga);
}

/**
   function computes partial derivative of rhoga with respect to pg
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval drhoga_dpg - partial derivative of rhoga with respect to pg
*/
double state_eq::get_drhoga_dpg(double pc,double pg,double t)
{
  double drhoga_dpg,pgw;

  //gas pressure check
  pgw = get_pgw(pc,t);

  if(pgw <= pg)
    drhoga_dpg = ma/gasr/t;
  else
    drhoga_dpg = 0.0;

  return(drhoga_dpg);
}

/**
   function computes partial derivative of rhoga with respect to pc
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval drhoga_dpc - partial derivative of rhoga with respect to pc
*/
double state_eq::get_drhoga_dpc(double pc,double pg,double t)
{
  double drhoga_dpc,dpgw_dpc,pgw;
  
  //gas pressure check
  pgw = get_pgw(pc,t);

  if(pgw <= pg){
    dpgw_dpc = get_dpgw_dpc(pc,t);
    drhoga_dpc = -ma*dpgw_dpc/gasr/t;
  }
  else
    drhoga_dpc = 0.0;
  
  return(drhoga_dpc);
}

/**
   function computes partial derivative of rhoga with respect to t
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval drhoga_dt - partial derivative of rhoga with respect to t
*/
double state_eq::get_drhoga_dt(double pc,double pg,double t)
{
  double drhoga_dt,dpgw_dt,rhoga,pgw;

  //gas pressure check
  pgw =get_pgw(pc,t); 

  if(pgw <= pg){
    rhoga = get_rhoga(pc,pg,t);
    dpgw_dt = get_dpgw_dt(pc,t);
    drhoga_dt = -ma*dpgw_dt/gasr/t - rhoga/t;
  }
  else
    drhoga_dt = 0.0;
  
  return(drhoga_dt);
}

/**
   function computes dynamic viscosity of dry air
   @param t - temperature

   @retval muga - dynamic viscosity of dry air
*/
double state_eq::get_muga(double t)
{
  double muga;

  muga = muga0 + alphaa*(t-t0) - betaa*pow((t-t0),2.0);

  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // muga = muga0 + alphaa*(t-t0) + betaa*pow((t-t0),2.0)    (eq. 50)

  return(muga);
}

/**
   function computes dry air specific heat

   @retval cpga - dry air specific heat
*/
double state_eq::get_cpga()
{
  double cpga;

  cpga = 1005.7;

  return(cpga);
}


//WATER VAPOR; DENSITY AND DERIVATIVES
//evaluate the Water Vapor Pressure in                
//equilibrium with curved surface liquid phase and
//its derivatives, relative humidity and the saturation value
//of vapour pressure
//from D.Gawin, F.Pesavento (PRVAP.f90)

/**
   function computes water vapour diffusivity in air
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval cdiff - water vapour diffusivity in air
*/
double state_eq::get_cdiff(double /*pc*/,double pg,double t)
{
  double cdiff;

  //conversion pow(a,b) -> exp(b*log(a))
  //cdiff = dv0*p0/pg*pow((t/t0),bv);
  if (t < tcr)
    cdiff = dv0*p0/pg*exp(bv*log(t/t0));
  else
    cdiff = dv0*p0/pg*exp(bv*log(tcr/t0));

  return(cdiff);
}

/**
   function computes water vapour pressure = Kelvin equation
   @param pc - capillary pressure
   @param t - temperature

   @retval pgw - water vapour pressure = Kelvin equation
*/
double state_eq::get_pgw(double pc,double t)
{
  double pgw,rhow,pgws,tt;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  
  pgw = pgws*exp(-1.0*pc*mw/rhow/gasr/tt);

  return(pgw);
}

/**
   function computes capillary pressure from water vapour pressure = inverse Kelvin equation
   @param pgw - capillary pressure of water vapour
   @param t - temperature

   @retval pc - capillary pressure from water vapour pressure = inverse Kelvin equation
*/
double state_eq::get_pcpgw(double pgw,double t)
{
  double pc,rhow,pgws,tt;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  
  pc = -rhow*gasr*tt/mw*log(pgw/pgws);

  //capillary pressure minimum check
  if (pc < 100.0)
    pc = 1.0e-10;

  return(pc);
}

/**
   function computes partial derivative of pgw with respect to pc (Kelvin equation)
   @param pc - capillary pressure
   @param t - temperature

   @retval dpgw_dpc - partial derivative of pgw with respect to pc (Kelvin equation)
*/
double state_eq::get_dpgw_dpc(double pc,double t)
{
  double dpgw_dpc,pgws,rhow,tt;
  
  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
    
  dpgw_dpc = -pgws*exp(-pc*mw/rhow/gasr/tt)*mw/rhow/gasr/tt;
  
  return(dpgw_dpc);
}

/**
   function computes partial derivative of pgw with respect to t (Kelvin equation)
   @param pc - capillary pressure
   @param t - temperature

   @retval dpgw_dt - partial derivative of pgw with respect to t (Kelvin equation)
*/
double state_eq::get_dpgw_dt(double pc,double t)
{
  double dpgw_dt,dpgws_dt,tt,pgws,rhow;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  dpgws_dt = get_dpgws_dt(tt);
  
  if(t < tcr)
    dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt) + pgws*exp(-pc*mw/rhow/gasr/tt)*(pc*mw/rhow/gasr/tt/tt);
  else
    dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt);

  return(dpgw_dt);
}

/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure
*/
double state_eq::get_pgws(double t)
{
  /* //Clausius-Clapeyron equation
     double pgws,pgws0,t0,dhvap;

     switch (Tm->ip[ipp].tm){//type of material
     case concreteB:{
     
     Tm->concrete[i].concreteB_pgws0(pgws0,t0);//funkce z concreteB.cpp
     
     break;
     }
     
     default:{
     fprintf (stderr,"\n\n unknown material type is required in function in function (%s, line %d).\n",__FILE__,__LINE__);
     }
     }
     
     dhvap = get_dhvap(pc,pg,t);
     
     pgws = pgws0*exp(-mw*dhvap/gasr*(1.0/t - 1.0/t0));
  */

  double t1,t2,t3,pgws,psl;
  
  t1 = 1.0/t;
  t2 = t*t;
  t3 = t*t*t;

  //critical point of water check
  if (t < tcr){
    psl = c8*t1 + c9 + c10*t + c11*t2 + c12*t3 + c13*log(t);
    pgws = exp(psl);
  }
  else
    pgws = 21780137.37214;

  return(pgws);
}

/**
   function computes partial derivative of water vapour saturation pressure with respect to t
   @param t - temperature

   @retval dpgws_dt - partial derivative of water vapour saturation pressure with respect to t
*/
double state_eq::get_dpgws_dt(double t)
{
  double t1,t12,t2,dpgws_dt,dpsl_dt,pgws;
  
  t1 = 1.0/t;
  t2 = t*t;
  t12 = 1.0/t2;

  pgws = get_pgws(t);

  //critical point of water check
  if (t < tcr){
    dpsl_dt = -c8*t12 + c10 + c11*2.0*t + c12*3.0*t2 + c13*t1;
    dpgws_dt = pgws*dpsl_dt;
  }
  else
    dpgws_dt = 0.0;
  
  return(dpgws_dt);
}

/**
   function computes mass concentration of water vapour air in gas phase
   @param pc - capillary pressure
   @param t - temperature

   @retval rhogw - mass concentration of water vapour air in gas phase
*/
double state_eq::get_rhogw(double pc,double t)
{
  double rhogw,pgw;

  pgw = get_pgw(pc,t);

  rhogw = pgw/t/gasr*mw;

  return(rhogw);
}

/**
   function computes capillary pressure from mass concentration of water vapour air in gas phase
   @param pc - capillary pressure
   @param t - temperature

   @retval pc - capillary pressure from mass concentration of water vapour air in gas phase
*/
double state_eq::get_pcrhogw(double rhogw,double t)
{
  double pc,pgw;

  pgw = rhogw*t*gasr/mw;
  pc = get_pcpgw(pgw,t);

  //minimum capillary pressure check
  if (pc < 100.0)
    pc = 1.0e-10;

  return(pc);
}

/**
   function computes partial derivative of rhogw with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval drhogw_dpc - partial derivative of rhogw with respect to pc
*/
double state_eq::get_drhogw_dpc(double pc,double t)
{
  double drhogw_dpc,dpgw_dpc;
  
  dpgw_dpc = get_dpgw_dpc(pc,t);

  drhogw_dpc = dpgw_dpc/t/gasr*mw;
  
  return(drhogw_dpc);
}

/**
   function computes partial derivative of rhogw with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval drhogw_dt - partial derivative of rhogw with respect to t
*/
double state_eq::get_drhogw_dt(double pc,double t)
{
  double drhogw_dt,rhogw,dpgw_dt;

  rhogw = get_rhogw(pc,t);
  dpgw_dt = get_dpgw_dt(pc,t);

  drhogw_dt = dpgw_dt/t/gasr*mw - rhogw/t;

  return(drhogw_dt);
}

/**
   function computes dynamic viscosity of water vapour
   @param t - temperature

   @retval mugw - dynamic viscosity of water vapour
*/
double state_eq::get_mugw(double t)
{
  double mugw;
 
  mugw = mugw0 + alphaw*(t-t0);

  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // alphaw = 3.633e-8

  return(mugw);
}

/**
   function computes water vapour specific heat
   @param t - temperature

   @retval cpgw - water vapour specific heat
*/
double state_eq::get_cpgw()
{
  double cpgw;

  cpgw = 1805.0;

  return(cpgw);
}

//WATER; DENSITY AND DERIVATIVES - DEFINES THE PHYSICAL PROPERTIES OF WATER
// from Dariusz Gawin (WATPROP.f90)
//WATER PRESSURE
/**
   function computes capillary pressure
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval pw - water pressure
*/
double state_eq::get_pw(double pc,double pg,double /*t*/)
{
  double pw;

  pw = pg - pc;

  return(pw);
}

/**
   function computes water density
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhow - water density
*/
double state_eq::get_rhow(double t)
{
  double rhow,tem;
  //double pw;

  /*   
  //nuw = 1.0/rhow; //water specific volume
  //betaw = 1.0/nuw*dnuw_dt;//volume thermal expansion coefficient
  //alphaw = 1.0/nuw*dnuw_dp;//isothermal compression modulus of water
  
  betaw = 0.68e-4;//[K-1] at t=273.15
  //betaw = 10.1e-4;//[K-1] at t=420.0
  alphaw = 4.3e-9;//[Pa-1]
  pw = pc - pg;
  rhow = rhow0*(1.0 - betaw*(t-t0) + alphaw*(pw - p0));
  */
  
  if (t < tcr){
    //rhow =  rhow0 *(1.+(pg - p0)/cwat + betawat*(t-t0));
    //drhow_dpg =  rhow0/cwat;
    //drhow_dpc =  0.0;
    //drhow_dt =  rhow0*betawat;
    //am_rhow =  rhow0 *(1.+ betawat*(t-t0));
    tem = t - t0;
    rhow =  (b0+(b1+(b2+(b3+(b4+b5*tem)*tem)*tem)*tem)*tem) + (pr1-prif)*
      (a0+(a1+(a2+(a3+(a4+a5*tem)*tem)*tem)*tem)*tem);
   
    //drhow_dpg =  0.0;
    //drhow_dpc =  0.0;
    //am_drhow_dt =  rhow0*betawat;
    //drhow_dt =  (b1+3(2*b2+(3*b3+(4*b4+5*b5*tem)*tem)*tem)*tem) + (pr1-prif)*
    //(a1+(2*a2+(3*a3+(4*a4+5*a5*tem)*tem)*tem)*tem);
    //am_d2rhow_d2t = 0.0;
    //d2rhow_d2t = (2*b2+(6*b3+(12*b4+20*b5*tem)*tem)*tem) + (pr1-prif)*
    //(2*a2+(6*a3+(12*a4+20*a5*tem)*tem)*tem);
    //d2rhow_d2pc = 0.0;
    //d2rhow_dtdpc = 0.0;
  }
  else{
    //rhow = rhow0 *(1.+(pg-p0)/cwat+betawat*(tcr-t0));
    //drhow_dpg =  rhow0/cwat;
    //am_rhow =  rhow0 *(1.+ betawat*(t-t0));
    tem = tcr - t0;
    rhow =  (b0+(b1+(b2+(b3+(b4+b5*tem)*tem)*tem)*tem)*tem) + (pr1-prif)*
      (a0+(a1+(a2+(a3+(a4+a5*tem)*tem)*tem)*tem)*tem);

    //drhow_dpg = 0.0;
    //drhow_dpc = 0.0;
    //drhow_dt = 0.0;
    //d2rhow_d2t = 0.0;
    //d2rhow_d2pc = 0.0;
    //d2rhow_dtdpc = 0.0;
  }
  
  return(rhow);
}

/**
   function computes derivative of water density with respect to temperature
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval drhow_dt - derivative of rhow with respect to t
*/
double state_eq::get_drhow_dt(double /*pc*/,double t)
{
  double drhow_dt,tem;
  
  if (t < tcr){
    tem = t - t0;
    //drhow_dt = (b1+(2.0*b2+(3.0*b3+(4.0*b4+5.0*b5*tem)*tem)*tem)*tem) + (pr1-prif)*(a1+(2.0*a2+(3.0*a3+(4.0*a4+5.0*a5*tem)*tem)*tem)*tem);
    drhow_dt = (b1+(2*b2+(3*b3+(4*b4+5*b5*tem)*tem)*tem)*tem) + (pr1-prif)*(a1+(2*a2+(3*a3+(4*a4+5*a5*tem)*tem)*tem)*tem);
  }
  else{
    drhow_dt = 0.0;
  }
    
  return(drhow_dt);
}

/**
   function computes enthalpy of evaporation (latent heat of vaporization)
   @param t - temperature

   @retval - enthalpy of evaporation (latent heat of vaporization)
*/
double state_eq::get_dhvap(double t)
{
  double dhvap,tem;

  tem = tcr - t;
  if (t < tcr)
    dhvap = hvap0*pow(tem,0.38);
  else
    dhvap = 0.0;

  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // dhvap = 2.672e+5 * pow((t-tcr),0.38)    (eq. 49)

  return(dhvap);
}

/**
   function computes dynamic viscosity of water = 1000e-6 Pa*s at 20 C
   @param t - temperature

   @retval muw dynamic viscosity of water = 1000e-6 Pa*s at 20 deg. C
*/
double state_eq::get_muw(double t)
{
  double muw;
  
  muw = muw0*pow((t - conb),conc);

  return(muw);
}

/**
   function computes water specific heat

   @retval cpw - water specific heat
*/
double state_eq::get_cpw()
{
  double c;

  c = cpw;

  return(c);
}

/**
   function computes water heat conductivity
   @param t - temperature
   
   @retval lambdaw - water heat conductivity
*/
double state_eq::get_lambdaw(double /*t*/)
{
  double lam;

  lam = lambdaw;

  return(lam);
}

/**
   function computes volume thermal expansion coefficient of water
   @param t - temperature
   
   @retval betaw - volume thermal expansion coefficient of water
*/
double state_eq::get_betaw(double /*t*/)
{
  double betaw;
  
  betaw = 0.68e-4;//[K-1] at t=273.15
  //betaw = 10.1e-4;//[K-1] at t=420.0

  //16.3.2006 docasne??!!
  betaw = 0.63e-5;
  //betaw = 0.68 + (10.1 - 0.68)/(420.0 - 273.15)*(t - 273.15);//linear expression (rough)

  return(betaw);
}


/**
   function computes compresibility coefficient of water
   
   @retval kw - compresibility coefficient of water

   16.3.2006, Tkr
*/
double state_eq::get_kw()
{
  double kw;
  
  kw = kw0;

  return(kw);
}




//EFFECTIVE PROPERTIES AND POROUS MATRIX PROPETIES
//AVERAGED DENSITY

/**
   function computes averaged density of multi-phase medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rho - averaged density of multi-phase medium
*/
double state_eq::get_rho(double pc,double pg,double t,long ipp)
{
  long i;
  double phi,rho,rhos,s,rhog;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    phi = Tm->concrete[i].concreteB_phi(t);//function from concreteB.cpp
    rhos = Tm->concrete[i].concreteB_rhos(t);//function from concreteB.cpp
    s = Tm->concrete[i].concreteB_sw(pc,t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    phi = Tm->baroghel[i].baroghel_phi();//function from baroghelB.cpp
    rhos = Tm->baroghel[i].baroghel_rhos();//function from baroghelB.cpp
    s = Tm->baroghel[i].baroghel_sw(pc,t);//function from baroghelB.cpp
    break;
  }    
  case C60baroghelB:{
    phi = Tm->C60baroghel[i].C60bar_phi(t);//function from C60baroghel.cpp
    rhos = Tm->C60baroghel[i].C60bar_rhos();//function from C60baroghel.cpp
    s = Tm->C60baroghel[i].sat(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    phi = Tm->C60baroghel[i].C60bar_phi(t);//function from C30baroghel.cpp
    rhos = Tm->C30baroghel[i].C30bar_rhos();//function from C30baroghel.cpp
    s = Tm->C30baroghel[i].sat(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    phi = Tm->o30bazant[i].o30baz_phi(t);//function from o30bazant.cpp
    rhos = Tm->o30bazant[i].o30baz_rhos(t);//function from o30bazant.cpp
    s = Tm->o30bazant[i].sat(pc,t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantB:{
    phi = Tm->C60bazant[i].C60baz_phi();//function from C60bazant.cpp
    rhos = Tm->C60bazant[i].C60baz_rhos();//function from C60bazant.cpp
    s = Tm->C60bazant[i].sat(pc,t);//function from C60bazant.cpp
    break;
  }    
  case C30bazantB:{
    phi = Tm->C30bazant[i].C30baz_phi();//function from C30bazant.cpp
    rhos = Tm->C30bazant[i].C30baz_rhos();//function from C30bazant.cpp
    s = Tm->C30bazant[i].sat(pc,t);//function from C30bazant.cpp
    break;
  }    

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  rhog = get_rhog(pc,pg,t);
  
  rho = (1.0 - phi)*rhos + phi*s + phi*(1.0 - s)*rhog;
  
  return(rho);
}

// BIOT'S CONSTANT
/**
   function computes Biot's constant
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval alpha - Biot's constant
*/
double state_eq::get_alpha(double pc, double pg, double t,long ipp)
{
  long i;
  double alpha;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//type of material
    /* case concreteB:{
       alpha = Tm->concrete[i].concreteB_alpha();//function from concreteB.cpp
       break;
       }
    */
  case baroghelB:{
    alpha = Tm->baroghel[i].baroghel_alpha();//function from concreteB.cpp
    break;
  }    
  case C60baroghelB:{
    alpha = Tm->C60baroghel[i].C60bar_alpha();//value from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    alpha = Tm->C30baroghel[i].C30bar_alpha();//value from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    alpha = Tm->o30bazant[i].o30baz_alpha();//value from o30bazant.cpp
    break;
  }
  case C60bazantB:{
    alpha = Tm->C60bazant[i].C60baz_alpha();//value from C60bazant.cpp
    break;
  }
  case C30bazantB:{
    alpha = Tm->C30bazant[i].C30baz_alpha();//value from C30bazant.cpp
    break;
  }
    
  case soilmat1:{
    alpha = Tm->soil1[i]._alpha();//function from soil1mat.cpp
    break;
  }  

  default:{
    
    double alpha,kt,ks;

    kt = get_kt(pc,pg,t,ipp);
    ks = get_ks(pc,pg,t,ipp);
    
    alpha = 1.0 - kt/ks;
    
    if(alpha > 1.0)
      alpha = 1.0;
    
    alpha = 1.0;//provisionally
    //fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(alpha);
}
    
//MATERIAL PROPERTIES

/**
   function computes density of solid phase
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhos - density of solid phase
*/
double state_eq::get_rhos(double t,long ipp)
{
  long i;
  double rhos;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    rhos = Tm->concrete[i].concreteB_rhos(t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    rhos = Tm->baroghel[i].baroghel_rhos();//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    rhos = Tm->C60baroghel[i].C60bar_rhos();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    rhos = Tm->C30baroghel[i].C30bar_rhos();//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    rhos = Tm->o30bazant[i].o30baz_rhos(t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantB:{
    rhos = Tm->C60bazant[i].C60baz_rhos();//function from C60bazant.cpp
    break;
  }    
  case C30bazantB:{
    rhos = Tm->C30bazant[i].C30baz_rhos();//function from C30bazant.cpp
    break;
  }    

  case soilmat1:{
    rhos = Tm->soil1[i]._rhos();//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return(rhos);
}

/**
   function computes bulk modulus of porous medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval kt - bulk modulus of porous medium
*/
double state_eq::get_kt(double pc,double pg,double t,long ipp)
{
  long i;
  double kt;
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    kt = Tm->concrete[i].concreteB_kt(pc,pg,t);//function from concreteB.cpp
    break;
  }
   default:{
   fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 

  return(kt);
}

/**
   function computes bulk modulus of solid phase
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval ks - bulk modulus of solid phase
*/
double state_eq::get_ks(double pc,double pg,double t,long ipp)
{
  long i;
  double ks;
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    ks = Tm->concrete[i].concreteB_ks(pc,pg,t);//function from concreteB.cpp
    break;
  }

  case soilmat1:{
    ks = Tm->soil1[i]._ks();//function from soil1mat.cpp
    break;
  }  

   default:{
   fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
  return(ks);
}

/**
   function computes gas relative permeability
   @param pc - capillary pressure
   @param t - temperature

   @retval krg - gas relative permeability
*/
double state_eq::get_krg(double pc,double t,long ipp)
{
  long i;
  double krg,s;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    s = Tm->concrete[i].concreteB_sw(pc,t);//function from concreteB.cpp
    krg = Tm->concrete[i].concreteB_krg(s);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    s = Tm->baroghel[i].baroghel_sw(pc,t);//function from baroghelB.cpp
    krg = Tm->baroghel[i].baroghel_krg(s);//function from baroghelB.cpp
    break;
  }    
  case C60baroghelB:{
    krg = Tm->C60baroghel[i].C60bar_krg(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    krg = Tm->C30baroghel[i].C30bar_krg(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    krg = Tm->o30bazant[i].o30baz_krg(pc,t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantB:{
    krg = Tm->C60bazant[i].C60baz_krg(pc,t);//function from C60bazant.cpp
    break;
  }    
  case C30bazantB:{
    krg = Tm->C30bazant[i].C30baz_krg(pc,t);//function from C30bazant.cpp
    break;
  }    
    
  case soilmat1:{
    krg = Tm->soil1[i]._krg(pc,t);//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(krg);
}

/**
   function computes water relative permeability
   @param pg - capillary gas pressure
   @param t - temperature

   @retval krw - water relative permeability
*/
double state_eq::get_krw(double pc,double t,long ipp)
{
  long i;
  double krw,s,rh;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    s = Tm->concrete[i].concreteB_sw(pc,t);//function from concreteB.cpp
    rh = get_rh(pc,t);
    krw = Tm->concrete[i].concreteB_krw(s,rh);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    krw = Tm->baroghel[i].baroghel_krw(pc,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    rh = get_rh(pc,t);
    krw = Tm->C60baroghel[i].C60bar_krw(pc,t,rh);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    rh = get_rh(pc,t);
    krw = Tm->C30baroghel[i].C30bar_krw(pc,t,rh);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    rh = get_rh(pc,t);
    krw = Tm->o30bazant[i].o30baz_krw(pc,t,rh);//function from o30bazant.cpp
    break;
  }        
  case C60bazantB:{
    krw = Tm->C60bazant[i].C60baz_krw(pc,t);//function from C60bazant.cpp
    break;
  }        
  case C30bazantB:{
    krw = Tm->C30bazant[i].C30baz_krw(pc,t);//function from C30bazant.cpp
    break;
  }        

  case soilmat1:{
    krw = Tm->soil1[i]._krw(pc,t);//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(krw);
}

/**
   function computes porosity
   @param t - temperature

   @retval phi - porosity
*/
double state_eq::get_phi(double t,long ipp)
{
  long i;
  double phi;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    phi = Tm->concrete[i].concreteB_phi(t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    phi = Tm->baroghel[i].baroghel_phi();//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    phi = Tm->C60baroghel[i].C60bar_phi(t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    phi = Tm->C30baroghel[i].C30bar_phi(t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    phi = Tm->o30bazant[i].o30baz_phi(t);//function from o30bazant.cpp
    break;
  }            
  case C60bazantB:{
    phi = Tm->C60bazant[i].C60baz_phi();//function from o30bazant.cpp
    break;
  }            
  case C30bazantB:{
    phi = Tm->C30bazant[i].C30baz_phi();//function from o30bazant.cpp
    break;
  }            
  case soilmat1:{
    phi = Tm->soil1[i]._phi();//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(phi);
}

/**
   function computes porosity
   @param t - temperature

   @retval dphi_dt - derivative of porosity with respect to temperature
*/
double state_eq::get_dphi_dt(double pc, double pg,double t,long ipp)
{
  long i;
  double dphi_dt,dehydw_dt,rhos;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
    /*  case concreteB:{
	phi = Tm->concrete[i].concreteB_phi(t);//function from concreteB.cpp
	break;
	}
    */ 
  case o30bazantB:{
    dehydw_dt = Tm->o30bazant[i].o30baz_dehydw_dt(pc,pg,t);//function from C30baroghel.cpp
    rhos = Tm->o30bazant[i].o30baz_rhos(t);//function from C30baroghel.cpp
    break;
  }
  case C30baroghelB:{
    dehydw_dt = Tm->C30baroghel[i].C30bar_dehydw_dt(pc,pg,t);//function from C30baroghel.cpp
    rhos = Tm->C30baroghel[i].C30bar_rhos();//function from C30baroghel.cpp
    break;
  }
  case C60baroghelB:{
    dehydw_dt = Tm->C60baroghel[i].C60bar_dehydw_dt(pc,pg,t);//function from C60baroghel.cpp
    rhos = Tm->C60baroghel[i].C60bar_rhos();//function from C60baroghel.cpp
    break;
  }
  case baroghelB:{
    dehydw_dt = Tm->baroghel[i].baroghel_dehydw_dt(pc,pg,t);//function from baroghelB.cpp
    rhos = Tm->baroghel[i].baroghel_rhos();//function from baroghelB.cpp
    break;
  }            
  case C60bazantB:{
    dehydw_dt = Tm->C60bazant[i].C60baz_dehydw_dt(pc,pg,t);//function from C60bazant.cpp
    rhos = Tm->C60bazant[i].C60baz_rhos();//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    dehydw_dt = Tm->C30bazant[i].C30baz_dehydw_dt(pc,pg,t);//function from C30bazant.cpp
    rhos = Tm->C30bazant[i].C30baz_rhos();//function from C30bazant.cpp
    break;
  }            
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  dphi_dt = dehydw_dt/rhos;

  return(dphi_dt);
}

/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval deff - effective diffusion coefficient of vapour inside pores
*/
double state_eq::get_deff(double pc,double pg,double t,long ipp)
{
  long i;
  double deff;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    deff = Tm->concrete[i].concreteB_deff(pc,pg,t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    deff = Tm->baroghel[i].baroghel_deff(pc,pg,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    deff = Tm->C60baroghel[i].C60bar_deff(pc,pg,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    deff = Tm->C30baroghel[i].C30bar_deff(pc,pg,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    deff = Tm->o30bazant[i].o30baz_deff(pc,pg,t);//function from o30bazant.cpp
    break;
  }            
  case C60bazantB:{
    deff = Tm->C60bazant[i].C60baz_deff(pc,pg,t);//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    deff = Tm->C30bazant[i].C30baz_deff(pc,pg,t);//function from C30bazant.cpp
    break;
  }            
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(deff);
}


/**
   function computes diffusion coefficient of vapour inside pores
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dg - diffusion coefficient of vapour inside pores
*/
double state_eq::get_dg(double /*pc*/,double /*pg*/,double /*t*/,long ipp)
{
  long i;
  double dg;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case soilmat1:{
    dg = Tm->soil1[i]._dg();//function from soil1mat.cpp
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(dg);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
*/
double state_eq::get_betas(long ipp)
{
  long i;
  double betas;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    betas = Tm->concrete[i].concreteB_betas();//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    betas = Tm->baroghel[i].baroghel_betas();//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    betas = Tm->C60baroghel[i].C60bar_betas();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    betas = Tm->C30baroghel[i].C30bar_betas();//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    betas = Tm->o30bazant[i].o30baz_betas();//function from o30bazant.cpp
    break;
  }                
  case C60bazantB:{
    betas = Tm->C60bazant[i].C60baz_betas();//function from C60bazant.cpp
    break;
  }                
  case C30bazantB:{
    betas = Tm->C30bazant[i].C30baz_betas();//function from C30bazant.cpp
    break;
  }                

  case soilmat1:{
    betas = Tm->soil1[i]._betas();//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(betas);
}

/**
   function computes intrinsic permeability
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval kintr - intrinsic permeability
*/
double state_eq::get_kintr(double pc,double pg,double t,long ipp)
{
  long i;
  double kintr,dam;

  i = Tm->ip[ipp].idm;

  // if damage was required by the material model then restore it
  if (Tm->givestatusntq(scal_iso_damage))
    dam = Tm->givenontransq(scal_iso_damage, ipp);
  else
    dam = 0.0;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    kintr = Tm->concrete[i].concreteB_kintr(pg,t,dam);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    kintr = Tm->baroghel[i].baroghel_kintr();//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    kintr = Tm->C60baroghel[i].C60bar_kintr(pc,pg,t,dam);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    kintr = Tm->C30baroghel[i].C30bar_kintr(pc,pg,t,dam);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    kintr = Tm->o30bazant[i].o30baz_kintr(pc,pg,t,dam);//function from o30bazant.cpp
    break;
  }            
  case C60bazantB:{
    kintr = Tm->C60bazant[i].C60baz_kintr();//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    kintr = Tm->C30bazant[i].C30baz_kintr();//function from C30bazant.cpp
    break;
  }            

  case soilmat1:{
    kintr = Tm->soil1[i]._kintr();//function from soil1mat.cpp
    break;
  }  
    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(kintr);
}

/**
   function computes specific heat of solid skeleton
   @param t - temperature

   @retval cps - specific heat of solid skeleton
*/
double state_eq::get_cps(double t,long ipp)
{
  long i;
  double cps;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    cps = Tm->concrete[i].concreteB_cps(t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    cps = Tm->baroghel[i].baroghel_cps();//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    cps = Tm->C60baroghel[i].C60bar_cps(t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    cps = Tm->C30baroghel[i].C30bar_cps(t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    cps = Tm->o30bazant[i].o30baz_cps(t);//function from o30bazant.cpp
    break;
  }            
  case C60bazantB:{
    cps = Tm->C60bazant[i].C60baz_cps();//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    cps = Tm->C30bazant[i].C30baz_cps();//function from C30bazant.cpp
    break;
  }            
    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(cps);
}

/**
   function computes effective thermal capacity of partially saturated medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval rhocp - effective thermal capacity of partially saturated medium
*/
double state_eq::get_rhocp(double pc,double pg,double t,long ipp)
{
  long i;
  double rhocp;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    rhocp = Tm->concrete[i].concreteB_rhocp(pc,pg,t,ipp);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    rhocp = Tm->baroghel[i].baroghel_rhocp(pc,pg,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    rhocp = Tm->C60baroghel[i].C60bar_rhocp(pc,pg,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    rhocp = Tm->C30baroghel[i].C30bar_rhocp(pc,pg,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    rhocp = Tm->o30bazant[i].o30baz_rhocp(pc,pg,t,ipp);//function from o30bazant.cpp
    break;
  }                
  case C60bazantB:{
    rhocp = Tm->C60bazant[i].C60baz_rhocp(pc,pg,t);//function from C60bazant.cpp
    break;
  }                
  case C30bazantB:{
    rhocp = Tm->C30bazant[i].C30baz_rhocp(pc,pg,t);//function from C30bazant.cpp
    break;
  }                

  case soilmat1:{
    rhocp = Tm->soil1[i]._rhocp();//function from soil1mat.cpp
    break;
  }  

  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(rhocp);
}


/**
   function computes specific heat of partially saturated medium
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval cp - specific heat of partially saturated medium
*/
double state_eq::get_cp(double pc,double pg,double t,long ipp)
{
  long i;
  double cp;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    cp = Tm->concrete[i].concreteB_cp(pc,pg,t,ipp);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    cp = Tm->baroghel[i].baroghel_cp(pc,pg,t,ipp);//function from baroghelB.cpp
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(cp);
}

/**
   function computes effective thermal conductivity
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval lambdaeff - effective thermal conductivity
*/
double state_eq::get_lambdaeff(double pc,double pg,double t,long ipp)
{
  long i;
  double lambdaeff;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteB:{
    lambdaeff = Tm->concrete[i].concreteB_lambdaeff(pc,pg,t);//function from concreteB.cpp
    break;
  }
  case baroghelB:{
    lambdaeff = Tm->baroghel[i].baroghel_lambdaeff(pc,pg,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelB:{
    lambdaeff = Tm->C60baroghel[i].C60bar_lambdaeff(pc,pg,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelB:{
    lambdaeff = Tm->C30baroghel[i].C30bar_lambdaeff(pc,pg,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantB:{
    lambdaeff = Tm->o30bazant[i].o30baz_lambdaeff(pc,pg,t);//function from o30bazant.cpp
    break;
  }            
  case C60bazantB:{
    lambdaeff = Tm->C60bazant[i].C60baz_lambdaeff(pc,t);//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    lambdaeff = Tm->C30bazant[i].C30baz_lambdaeff(pc,t);//function from C30bazant.cpp
    break;
  }            
  case soilmat1:{
    lambdaeff = Tm->soil1[i]._lambdaa();//function from C30bazant.cpp
    break;
  }            
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(lambdaeff);
}

/**
   function volume computes thermal expansion coefficient of solid - water vapor for incompressible solid grains
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval betaswg - volume thermal expansion coefficient of solid - water vapor for incompressible solid grains
*/
double state_eq::get_betaswg(double pc,double /*pg*/,double t,long ipp)
{
  double betaswg,betas,phi,s,rhogw,rhow,betaw;

  betas = get_betas(ipp);
  phi = get_phi(t,ipp);
  s = get_s(pc,t,ipp);
  rhogw = get_rhogw(pc,t);
  rhow = get_rhow(t);
  betaw = get_betaw(t);
  
  betaswg = betas*(1.0 - phi)*((1.0 - s)*rhogw + s*rhow) + phi*betaw*rhow*s;

  return(betaswg);
}


/**
   function volume computes thermal expansion coefficient of solid - water vapor for compressible solid grains
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval betaswg - volume thermal expansion coefficient of solid - water vapor for compressible solid grains
*/
double state_eq::get_betaswg_c(double pc,double pg,double t,long ipp)
{
  double alpha,betaswg,betas,phi,sw,sg,rhogw,rhow,betaw;

  alpha = get_alpha(pc,pg,t,ipp);
  betas = get_betas(ipp);
  phi = get_phi(t,ipp);
  sw = get_s(pc,t,ipp);
  sg = 1.0 -sw;
  rhogw = get_rhogw(pc,t);
  rhow = get_rhow(t);
  betaw = get_betaw(t);
  
  betaswg = betas*(alpha - phi)*(sg*rhogw + sw*rhow) + phi*betaw*rhow*sw;

  return(betaswg);
}


/**
   function volume computes thermal expansion coefficient of solid - water for incompressible solid grains
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval betasw - volume thermal expansion coefficient of solid - water for incompressible solid grains
*/
double state_eq::get_betasw(double pc,double /*pg*/,double t,long ipp)
{
  double betasw,betas,phi,s,betaw;

  betas = get_betas(ipp);
  phi = get_phi(t,ipp);
  s = get_s(pc,t,ipp);
  betaw = get_betaw(t);
  
  betasw = s*((1.0 - phi)*betas + phi*betaw);

  return(betasw);
}


/**
   function volume computes thermal expansion coefficient of solid - water for compressible solid grains
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval betasw - volume thermal expansion coefficient of solid - water for compressible solid grains
*/
double state_eq::get_betasw_c(double pc,double pg,double t,long ipp)
{
  double betasw,betas,phi,sw,betaw,alpha,rhow;

  alpha = get_alpha(pc,pg,t,ipp);
  betas = get_betas(ipp);
  phi = get_phi(t,ipp);
  sw = get_s(pc,t,ipp);
  betaw = get_betaw(t);
  rhow = get_rhow(t);
  
  betasw = rhow*sw*((alpha - phi)*betas + phi*betaw);

  return(betasw);
}


/**
   function volume computes thermal expansion coefficient of solid - gas for compressible solid grains
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval betasg - volume thermal expansion coefficient of solid - gas for compressible solid grains
*/
double state_eq::get_betasg_c(double pc,double pg,double t,long ipp)
{
  double betasw,betas,phi,sg,alpha;

  alpha = get_alpha(pc,pg,t,ipp);
  betas = get_betas(ipp);
  phi = get_phi(t,ipp);
  sg = 1.0 - get_s(pc,t,ipp);
  
  betasw = sg*(alpha - phi)*betas;

  return(betasw);
}




/* this part is only for concrte at high temperatures */
/**
   function computes derivative of hydration degree with respect to temperature
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval dehydw_dt - derivative of hydration degree with respect to temperature
*/
double state_eq::get_dehydw_dt(double pc,double pg,double t,long ipp)
{
  double dehydw_dt;
  long i;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
    /*  case concreteB:{
	dhdehydr = Tm->concrete[i].concreteB_dhdehydr(pc,pg,t);//function from concreteB.cpp
	break;
	}
    */ 
  case o30bazantB:{
    dehydw_dt = Tm->o30bazant[i].o30baz_dehydw_dt(pc,pg,t);//function from o30bazant.cpp
    break;
  }           
  case C30baroghelB:{
    dehydw_dt = Tm->C30baroghel[i].C30bar_dehydw_dt(pc,pg,t);//function from C30baroghel.cpp
    break;
  }
  case C60baroghelB:{
    dehydw_dt = Tm->C60baroghel[i].C60bar_dehydw_dt(pc,pg,t);//function from C60baroghel.cpp
    break;
  }
  case baroghelB:{
    dehydw_dt = Tm->baroghel[i].baroghel_dehydw_dt(pc,pg,t);//function from baroghelB.cpp
    break;
  }            
  case C60bazantB:{
    dehydw_dt = Tm->C60bazant[i].C60baz_dehydw_dt(pc,pg,t);//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    dehydw_dt = Tm->C30bazant[i].C30baz_dehydw_dt(pc,pg,t);//function from C30bazant.cpp
    break;
  }            
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(dehydw_dt);
}


/**
   function computes hydration energy
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval hydren - hydration energy
*/
double state_eq::get_hydren(double pc,double pg,double t,long ipp)
{
  double hydren;
  long i;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
    /*  case concreteB:{
	dhdehydr = Tm->concrete[i].concreteB_dhdehydr(pc,pg,t);//function from concreteB.cpp
	break;
	}
    */ 
  case o30bazantB:{
    hydren = Tm->o30bazant[i].o30baz_hydren(pc,pg,t);//function from o30bazant.cpp
    break;
  }           
  case C30baroghelB:{
    hydren = Tm->C30baroghel[i].C30bar_hydren(pc,pg,t);//function from C30baroghel.cpp
    break;
  }
  case C60baroghelB:{
    hydren = Tm->C60baroghel[i].C60bar_hydren(pc,pg,t);//function from C60baroghel.cpp
    break;
  }
  case baroghelB:{
    hydren = Tm->baroghel[i].baroghel_hydren(pc,pg,t);//function from baroghelB.cpp
    break;
  }            
  case C60bazantB:{
    hydren = Tm->C60bazant[i].C60baz_hydren(pc,pg,t);//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    hydren = Tm->C30bazant[i].C30baz_hydren(pc,pg,t);//function from C30bazant.cpp
    break;
  }            
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(hydren);
}



/**
   function computes Water/Cement ratio
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval fste - Water/Cement ratio
*/
double state_eq::get_fste(double pc,double pg,double t,long ipp)
{
  double fste;
  long i;

  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
    /*  case concreteB:{
	dhdehydr = Tm->concrete[i].concreteB_dhdehydr(pc,pg,t);//function from concreteB.cpp
	break;
	}
    */ 
  case o30bazantB:{
    fste = Tm->o30bazant[i].o30baz_fste(pc,pg,t);//function from o30bazant.cpp
    break;
  }           
  case C30baroghelB:{
    fste = Tm->C30baroghel[i].C30bar_fste(pc,pg,t);//function from C30baroghel.cpp
    break;
  }
  case C60baroghelB:{
    fste = Tm->C60baroghel[i].C60bar_fste(pc,pg,t);//function from C60baroghel.cpp
    break;
  }
  case baroghelB:{
    fste = Tm->baroghel[i].baroghel_fste(pc,pg,t);//function from baroghelB.cpp
    break;
  }            
  case C60bazantB:{
    fste = Tm->C60bazant[i].C60baz_fste(pc,pg,t);//function from C60bazant.cpp
    break;
  }            
  case C30bazantB:{
    fste = Tm->C30bazant[i].C30baz_fste(pc,pg,t);//function from C30bazant.cpp
    break;
  }            
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(fste);
}

