/*
  File:             consol_hawf3.cpp
  Author:           Tomas Krejci, 30/05/2016 revised 03/10/2023
  Purpose:          computes conductivity and capacity matrices in a material point for consol_hawf3 porous media;
                    material model for saturated-nonsaturated air and water flow and heat transfer 
                    in a deforming porous medium (soils), liquid water and gas (moist air) advection, diffusion of vapour, heat conduction plus convective terms of water and gas
  unknowns:         number of unknowns=3, pw = pore water pressure, pg = pore gas(air) pressure, t = temperature
  sources:          
  
  FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
  --------------------------------------------------
  
  sources: 
  1. THE FINITE ELEMENT METHOD IN THE STATIC AND DYNAMIC DEFORMATION AND CONSOLIDATION OF POROUS MEDIA
  R. W. Lewis, B.A. Schrefler, pp. 354-396 - some mistakes were found
  
  2. NONLINEAR MODELLING OF CONCRETE AS CONSOL_HAWF3 POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
  Francesco Pesavento - doctoral thesis                      
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "constrel.h"
#include "consol_hawf3.h"
#include "globalt.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "globmatt.h"

#include <errno.h>

con_hawf3mat::con_hawf3mat()
{
  compress = 0;          //compressible grains: 0=no; 1=yes
  por_type = 0;          //porosity calculation type
  kintr_type = 0;        //intrinsic permability calculation type
  krw_type = 0;          //relative permeability calculation type
  krg_type = 0;          //relative permability calculation type
  deff_type = 0;         // diffusion calculation type
  sr_type = 1;           //retention curve calculation type
  xi_type = 1;           //effective stress parameter type
  lambda_type =0;        //heat conduction calculation type
  cps_type = 0;          //specific heat calculation type
  thermal_capacity_type = 0; //thermal capacity type
  betas_type = 0;        //thermal expansion calculation type

  vol_strain_effect = 0; //volumetric strain rate influence: 0=no; 1=yes 
  mefel_units = 1.0; //basic units for pressures = Pa (Pascals)
  wrc_vol_strain_effect = 0; //volumetric strain rate influence on water retention curve: 0=no; 1=yes

  //STATE VARIABLES

  mw = 18.01528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

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
  bv0 = 1.667;
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
  //rhow0 = 999.84;//testing??!!
  rhow0 = 1000.0;//testing??!!
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
  muw0const = 1.0e-3;
  conb  =  229.0;
  conc  = -1.562;
  cpw0 = 4181.0;
  lambdaw = 0.6;  
  //kw0 = 0.43e10;//compresibility coefficient of water - this is nonsense //09/06/2017
  kw0 = 2.0e9;//bulk modulus of water
 
  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha0 = 0.0;
  ks0 = 0.0;
  kt0 = 0.0;
  phi0 = 0.0;
  kintr0 = 0.0;
  kintrw0 = 0.0;
  kintrg0 = 0.0;
  betas0 = 0.0;
  betas_dry = 0.0;
  betas_wet = 0.0;
  deff0 = 0.0;
  cdiff0 = 25.4e-6;//initial water vapour diffusivity in air
  rhos0 = 0.0;
  cps0 = 0.0;
  lambda_eff0 = 0.0;
  lambda_dry = 0.0;
  lambda_wet = 0.0;
  sr_dry = 0.0;
  sr_wet = 0.0;
  rhocp_dry = 0.0;
  rhocp_wet = 0.0;
  tau0 = 0.0;

  gamma = 0.0;
  lambda0 = 0.0;
  s_entry = 0.0;

  sirr = 0.0;
  ssat = 1.0;

  krw0 = 1.0;
  lambda_krw = 1.9;  //parameter for exponential function of relative permeability
  beta_krw = 0.0;    //parameter for power function of relative permeability

  bb1 = 0.0;
  phi01 = 0.0;
  ng = 0.0;

  //gas relative permeability parameters:
  krg0 = 1.0;
  s_crit = 0.8; 
  ag = 2.0; 

  scale_pw = Tp->scale[0];
  scale_pg = Tp->scale[1];
  scale_t = Tp->scale[2];

  pw_bc = 0.0; //free boundary pressure 
  rel_gas_press = 0;

  //parameters for the artificial material:
  kww0 = kwg0 = kwt0 = kgw0 = kgg0 = kgt0 = ktw0 = ktg0 = ktt0 = 0.0;
  capww0 = capwg0 = capwt0 = capgw0 = capgg0 = capgt0 = captw0 = captg0 = captt0 =0.0;

  kggmin = 1.0e-20;
  kg0 = kgn = 0.0;

  cps_dry = cps_wet = cps0 = cps_lin = 0.0;
}
con_hawf3mat::~con_hawf3mat()
{}


/**
   function reads parameters
   
   @param in - input file

   20/05/2017, TKr, revised 03/10/2023
*/
void con_hawf3mat::read(XFILE *in)
{
  xfscanf (in,"%k%m","heatairwaterflowtype",&heatairwaterflowtype_kwdset, &model_type);
  xfscanf (in,"%d", &compress);

  // common material parameters
  xfscanf (in,"%le %le %le %le %le %le %d %d %d %d %d %d %d %d %d %d %d", &alpha0, &ks0, &rhos0, &pw_bc, &tau0, &kintr0, &por_type, &kintr_type, &krw_type, &krg_type, &deff_type, &sr_type, &xi_type, &lambda_type, &cps_type, &thermal_capacity_type, &betas_type);

  switch (model_type){
  case artificial3:{//artificial isotropic material for non-isotherma air-water flow
    xfscanf (in,"%le %le %le %le %le %le %le %le %le", &kww0, &kwg0, &kwt0, &kgw0, &kgg0, &kgt0, &ktw0, &ktg0, &ktt0);
    xfscanf (in,"%le %le %le %le %le %le %le %le %le", &capww0, &capwg0, &capwt0, &capgw0, &capgg0, &capgt0, &captw0, &captg0, &captt0);
  }
  case lewis_and_schrefler3:{//Lewis and Schrefler's model approach
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's model approach coupled with mechanics
    xfscanf (in,"%le %d %d", &mefel_units, &vol_strain_effect, &wrc_vol_strain_effect);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  if(model_type != artificial3){
    //porosity calculation:
    switch (por_type){
    case 0:{//constant
      xfscanf (in,"%le",&phi0);
      break;
    }
    case 1:{//dependent on mefel calculation
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //intrinsic permebility calculation:
    switch (kintr_type){
    case 0:{//constant - for permeability of both phases
      break;
    }
    case 1:{//dependent on porosity - for permeability of both phases
      xfscanf (in,"%le %le", &bb1, &phi01);
      break;
    }
    case 2:{//dependent on porosity - cubic and quadratic - for permeability of both phases
      xfscanf (in,"%le", &phi01);
      break;
    }
    case 3:{//constant separate permeabilities of liquid water and gas
      xfscanf (in,"%le", &kintrw0);
      xfscanf (in,"%le", &kintrg0);
      break;
    }
    case 4:{//dependent on porosity - separate permeabilities of liquid water and gas
      xfscanf (in,"%le %le %le", &kintrw0, &bb1, &phi01);
      xfscanf (in,"%le %le", &kintrg0, &ng);
      break;
    }      
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //relative permebility calculation:
    switch (krw_type){
    case 0:{//constant
      xfscanf (in,"%le", &krw0);
      break;
    }
    case 1:
    case 2:{//dependent on saturation degree:
      xfscanf (in,"%le %le", &sirr, &ssat);
      break;
    }
    case 3:{//exponential
      xfscanf (in,"%le", &lambda_krw);
      break;
    }
    case 4:{//liakopoulos
      break;
    }
    case 5:{//double exponential
      xfscanf (in,"%le %le", &sirr, &ssat);
      xfscanf (in,"%le", &beta_krw);
      break;
    }
    case 6:{//FEBEX granit
      //xfscanf (in,"%le %le %le %le", &m1, &m2, &m3, &m4);
      break;
    }
    case 7:{//van_genuchten
      break;
    }
    case 9:{//bazant
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //gas relative permebility calculation:
    switch (krg_type){
    case 0:{//constant
      xfscanf (in,"%le", &krg0);
      break;
    }
    case 1:{//dependent on saturation degree:
      xfscanf (in,"%le %le", &s_crit, &ag);
      break;
    }
    case 2:{//dependent on saturation degree:
      xfscanf (in,"%le %le", &s_crit, &ag);
      break;
    }
    case 3:{//exponential dependence on void ratio and saturation for Kg
      xfscanf (in,"%le %le", &kg0, &kgn); //not finished
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //retention curve type:
    switch (sr_type){
    case baroghel_sr:{//Baroghle-Bouny approach
      baroghel_ret.read(in);
      break;
    }
    case bazant_sr:{//Bazant approach
      break;
    }
    case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
      lewis_ret.read(in);
      break;
    }
    case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
      van_genuchten_ret.read(in);
      break;
    }
    case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
      //van_genuchten_ret.read2(in);
      break;
    }
    case mefel_sr:{//from MEFEL
      break;
    }
    case table_sr:{//from table
      //reading of retention curve:
      data.read (in);
      break;
    }
    case masin_sr:{//extended formulation from Brooks and Correy according to Masin
      masin_ret.read(in);
      break;
    }
    case febex_granit_sr:{//FEBEX granit
      febex_granit_ret.read(in);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    switch (xi_type){
    case biot_xi:
      break;
    case biot_reduced_xi:{//coefficient for effective stress factor reduction
      xfscanf (in,"%le  %le", &gamma,&lambda0);
      break;
    }
    case biot_masin_xi:{
      xfscanf (in,"%le  %le", &gamma,&s_entry);
      break;
    }
    case masin_xi:
      break;
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //conductivity calculation:
    switch (lambda_type){
    case 0:{//constant
      xfscanf (in,"%le", &lambda_eff0);
      break;
    }
    case 1:{//dependent on moisture
      xfscanf (in,"%le %le", &lambda_dry, &lambda_wet);
      break;
    }
    case 2:{//dependent on moisture
      xfscanf (in,"%le %le %le %le", &lambda_dry, &lambda_wet, &sr_dry, &sr_wet);
      break;
    }
    case 3:{//dependent on moisture
      xfscanf (in,"%le %le", &lambda_dry, &lambda_wet);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //specific heat calculation:
    switch (cps_type){
    case 0:{//constant
      xfscanf (in,"%le", &cps0);
      break;
    }
    case 1:{//dependent on moisture
      xfscanf (in,"%le %le %le %le", &cps_dry, &cps_wet, &sr_dry, &sr_wet);
      break;
    }
    case 2:{//dependent on moisture
      xfscanf (in,"%le %le", &cps_dry, &cps_wet);
      break;
    }
    case 3:{//dependent on temperature
      xfscanf (in,"%le %le", &cps0, &cps_lin);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    //thermal capacity calculation:    
    switch (thermal_capacity_type){
    case 0:
    case 1:
    case 2:{
      //nothing to read
      break;
    }
    case 3:{//effective heat capacity is directly measured with respec to saturation degree
      //linear approximation
      xfscanf (in,"%le %le %le %le", &rhocp_dry, &rhocp_wet, &sr_dry, &sr_wet);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    }
    
    //thermal expansion calculation:
    switch (betas_type){
    case 0:{//constant
      xfscanf (in,"%le", &betas0);
      break;
    }
    case 1:{//dependent on moisture
      xfscanf (in,"%le %le", &betas_dry, &betas_wet);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
  }
}


/**
   function prints parameters
   
   @param out - output file

   20/05/2017, TKr
*/
void con_hawf3mat::print(FILE */*out*/)
{
  /*   fprintf (out,"\n %d ", int(model_type));
       fprintf (out,"\n %d ", compress);
       
       switch (model_type){
       case lewis_and_schrefler3:{//Lewis and Schrefler's book
       fprintf (out,"\n %le %le %le %le %le %le %le %le ", alpha0, ks0, phi0, kintr0, betas0, rhos0, cps0, lambda0);
       lewis_ret.print(out);
       break;
       }
       case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book
       fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le %le %le %le ", mefel_units, alpha0, ks0, phi0, kintr0, betas0, rhos0, cps0, lambda_dry, lambda_wet, sirr0, tau0, pw_bc);
       break;
       }
       case van_genuchten3:{//partially saturated medium =Van Genuchten model
       fprintf (out,"\n %le %le %le %le %le %le %le %le %le ", mefel_units, alpha0, ks0, phi0, kintr0, betas0, rhos0, cps0, lambda0);
       van_genuchten_ret.print(out);
       break;
       }
       case lewis_and_schrefler3_2:{//Lewis and Schrefler's book p. 381
       fprintf (out,"\n %le %le %le %le %le %le %le %le ", alpha0, ks0, phi0, kintr0, betas0, rhos0, cps0, lambda0);
       break;
       }
       default:{
       print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
       abort();
       }
       } 
  */
}


/****************************************************/
/** The basic functions and constitutive relations **/
/****************************************************/


/**
   function computes cubic thermal expansion coefficient of solid (K-1)
   @param ipp - number of integration point


   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
   
   03/10/2023 TKr
*/
double con_hawf3mat::get_betas(long /*ipp*/)
{
  double betas;

  betas = betas0;

  switch (betas_type){
  case 0:{ //constant
    betas = betas0;
    break;
  }
  case 1:{ //moisture dependent - not finished
    break;
  }
 default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(betas);
}


/**
   function computes specific heat of solid skeleton
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval cps - specific heat of soil skeleton

   03/10/2023 TKr
*/
double con_hawf3mat::get_cps(double pw, double pg, double t, long ipp)
{
  double cps;
  double sw;

  cps = cps0;

  switch (cps_type){
  case 0:{ //constant
    cps = cps0;
    break;
  }
  case 1:{ //moisture dependent
    sw = get_sw(pw,pg,t,ipp);
    cps = cps_dry + (cps_wet - cps_dry)*(sw - sr_dry)/(sr_wet - sr_dry);
    break;
  }
  case 2:{ //moisture dependent
    sw = get_sw(pw,pg,t,ipp);
    cps = cps_dry + (cps_wet - cps_dry)*sw; //this is maybe linear approximation for bentonite:
    break;
  }
  case 3:{ //temperature dependent
    cps = cps_lin*t + cps_0; //this is maybe linear approximation for bentonite:
    break;
  }
 default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(cps);
}


/**
   function computes  volume density of concrete skeleton,
   changes of solid density, caused by dehydratation process

   @param t - temperature

   @retval rhos - volume density of soil skeleton

   03/10/2023 TKr
*/
double con_hawf3mat::get_rhos(double /*t*/)
{
  double rhos;
   
  rhos = rhos0;//temp.
  
  return(rhos);
}

/**
   function computes water vapour diffusivity in air
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval cdiff - water vapour diffusivity in air

   03/10/2023 TKr
*/
double con_hawf3mat::get_cdiff(double pw,double pg,double t)
{
  double cdiff,pgw;
  
  if(rel_gas_press == 1)
    pg = pg + p0;

  //actualized - check if pgw>pg
  pgw = 0.0;
  get_pgw(pw,pg,t);
  if (pgw > pg)
    pg = pgw;

  //conversion pow(a,b) -> exp(b*log(a))
  //cdiff = dv0*p0/pg*pow((t/t0),bv0);
  if (t < tcr)
    cdiff = dv0*p0/pg*exp(bv0*log(t/t0));
  else
    cdiff = dv0*p0/pg*exp(bv0*log(tcr/t0));

  //this equation below is from Patek's thesis:
  //dv0 = 5.9e-6;
  //cdiff = dv0*pow(t,2.3)/pg;

  return(cdiff);
}

/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval deff - effective diffusion coefficient of vapour inside pores

   03/10/2023 TKr
*/
double con_hawf3mat::get_dg(double pw,double pg,double t,long ipp)
{
  double deff;
  double phi,cdiff,sw,tau;

  deff = deff0;

  switch (deff_type){
  case 0:{//constant
    deff = cdiff0;
    break;
  }
  case 1:{//moisture dependent
    phi = get_porosity(ipp);
    cdiff = get_cdiff(pw,pg,t);
    sw = get_sw(pw,pg,t,ipp);
    tau = tau0;
    deff = phi*(1.0 - sw)*tau*cdiff;
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(deff);
}



/**
   function computes effective stress factor xi
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval xi - factor xi

   04/04/2023, TKr
*/
double con_hawf3mat::get_xi(double pw, double pg, double t, long ipp)
{
  double suc=0.0,sr=0.0,xi=0.0;
  
  switch (xi_type){
  case biot_xi:
    xi = get_sw(pw,pg,t,ipp);
    break;
  case biot_reduced_xi:
    //xi = gamma*get_sw(pw,ipp);
    sr = get_sw(pw,pg,t,ipp);
    //xi = pow(sr,(gamma/lambda0));
    xi = pow(sr,gamma);
    xi = (1-gamma)*xi;
    break;
  case biot_masin_xi:{//according to masin for testing
    sr = get_sw(pw,pg,t,ipp);
    suc = pg - p0 - pw;//pore gas pressure without atmospheric pressure
    if (suc>=s_entry)
      xi = pow((s_entry/suc),gamma);
    else
      xi = 1.0;
    if (suc>=s_entry)
      xi = (1-gamma)*xi;
    else
      xi = 1.0;
    break;
  }
  case masin_xi:{
    double e=0.0,dpw=0.0,dpg=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }
    
    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpg = Tm->ip[ipp].av[1]-Tm->ip[ipp].pv[1];
    dpc = dpg - dpw;

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;

    xi = masin_ret.psi(pc,dpc,e,t);//positive value of suction
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(xi);
}



/**
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval sw - degree of saturation

   20/05/2017, TKr, revised 03/10/2023 TKr
*/
double con_hawf3mat::get_sw(double pw, double pg, double t, long ipp)
{
  double sw,pc;
  sw = 0.0;  
  pc = 0.0;
  
  switch (sr_type){
  case bazant_sr:{//Bazant
    pc = pg - pw;
    sw = bazant_ret.sat(pc,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;

    sw = lewis_ret.sw(-pc);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    sw = Tm->givenontransq(saturation_deg, ipp); //actual saturation degree
    break;
  }
  case table_sr:{//saturation degree and its derivative are obtained from table;
    //actual saturation degree    
    if (data.tfunc == tab)
      {
	if ((pw < data.tabf->x[0]) || (pw > data.tabf->x[data.tabf->asize-1]))
	  {
	    print_err("required value %le is out of table range <%le;%le> on ip=%ld\n", 
		      __FILE__, __LINE__, __func__, pw, data.tabf->x[0], data.tabf->x[data.tabf->asize-1], ipp);
	    abort();
	  }
      }
    sw = data.getval (pw);
    
    break;
  }
    
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    sw = van_genuchten_ret.sw(-pc,t);//suction = capillary pressure
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model type 2
    //sw = van_genuchten_ret.sw2(pw);
    break;
  }
  case masin_sr:{//extended formulation from Brooks and Correy according to Masin
    double e=0.0,dpw=0.0,dpg=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }

    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpg = Tm->ip[ipp].av[1]-Tm->ip[ipp].pv[1];
    dpc = dpg - dpw;
    
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;
    sw = masin_ret.sw(pc,dpc,e,t);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    //pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    sw = febex_granit_ret.sw(pc);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  
  return(sw);
}


/**
   function computes specific moisture content (partial derivative of degree of saturation with respect to pc)
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ds_dpc - specific moisture content = partial derivative of degree of saturation with respect to pc


   20/05/2017, TKr, revised 03/10/2023 TKr
*/
double con_hawf3mat::get_ds_dpc(double pw, double pg, double t, long ipp)
{
  double dsw_dpc,pc;
  dsw_dpc = pc = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    pc = pg - pw;
    dsw_dpc = bazant_ret.dsat_dpc(pc,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    pc = pg - pw;

    dsw_dpc = -lewis_ret.dsw_dpw(-pc);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    dsw_dpc = -Tm->givenontransq(der_saturation_deg, ipp); //actual derivative of saturation degree
    dsw_dpc = dsw_dpc/mefel_units; //basic units = Pa
    break;
  }

  case table_sr:{//saturation degree and its derivative are obtained from table;
    pc = pg - pw;

    dsw_dpc = -data.getderiv (-pc); //actual derivative of saturation degree
    break;
  }

  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    dsw_dpc = -van_genuchten_ret.dsw_dpw(-pc,t);//suction = capillary pressure
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
    //dsw_dpc = -van_genuchten_ret.dsw_dpw2(pw);
    break;
  }
  case masin_sr:{//extended formulation from Brooks and Correy according to Masin
    double e=0.0,dpw=0.0,dpg=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }

    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpg = Tm->ip[ipp].av[1]-Tm->ip[ipp].pv[1];
    dpc = dpg - dpw;
    
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;

    dsw_dpc = masin_ret.dsw_dpw(pc,dpc,e,t);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    //pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;
    
    dsw_dpc = febex_granit_ret.dsw_ds(pc);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ds_dt - partial derivative of degree of saturation with respect to temperature


   20/05/2017, TKr, revised 03/10/2023
*/
double con_hawf3mat::get_ds_dt(double pw, double pg, double t, long ipp)
{
  double dsw_dt,pc;
  dsw_dt = pc = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    pc = pg - pw;
    dsw_dt = bazant_ret.dsat_dt(pc,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    pc = pg - pw;

    dsw_dt = lewis_ret.dsw_dt(-pc);//suction = capillary pressure
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    dsw_dt = Tm->givenontransq(der_saturation_deg_dtemp, ipp); //actual derivative of saturation degree
    dsw_dt = dsw_dt/mefel_units; //basic units = Pa
    break;
  }
    
  case table_sr:{//saturation degree and its derivative are obtained from table;
    pc = pg - pw;

    //dsw_dt = -data.getderiv_dt (pc); //actual derivative of saturation degree
    break;
  }
    
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    dsw_dt = van_genuchten_ret.dsw_dt(-pc,t);
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
    //dsw_dt = van_genuchten_ret.dsw_dt2(pw);
    break;
  }
  case masin_sr:{//extended formulation from Brooks and Correy according to Masin
    double e=0.0,dpw=0.0,dpg=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }
    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpg = Tm->ip[ipp].av[1]-Tm->ip[ipp].pv[1];
    dpc = dpg - dpw;

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;

    dsw_dt = masin_ret.dsw_dt(pc,dpc,e,t);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    //pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    dsw_dt = febex_granit_ret.dsw_dt(pc);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(dsw_dt);
}


/**
   function returns porosity
   @param ipp - number of integration point

   @retval phi - porosity

   04/04/2023, TKr
*/
double con_hawf3mat::get_porosity(long ipp)
{
  double por;
  
  //porosity calculation:
  switch (por_type){
  case 0:{//constant
    por = phi0;
    break;
  }
  case 1:{//dependent on mefel calculation
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
    }
    else
      por = phi0;//for testing
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(por);
}


/**
   function computes intrinsic permeability assumed for water and gas (some models have two separate values)
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval kintr - intrinsic permeability

   27/07/2017, TKr, revised 03/10/2023 TKr
*/
double con_hawf3mat::get_kintr(double pw, double pg, double t, long ipp)
{
  double kintr;
  double phi;

  kintr = kintr0;

  switch (kintr_type){
  case 0:{//constant
    kintr = kintr0;
    break;
  }
  case 1:{//dependent on porosity
    phi = get_porosity(ipp);
    kintr = kintr0*exp(bb1*(phi - phi01));
    break;
  }
  case 2:{//dependent on porosity - cubic and quadratic
    //Kozeny's approach for bentonite:
    phi = get_porosity(ipp);
    kintr = kintr0*phi*phi*phi/(1 - phi01)*(1 - phi01)*(1 - phi)*(1 - phi)/(phi01*phi01*phi01);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(kintr);
}



/**
   function computes intrinsic permeability assumed only for liquid (some models have two separate values)
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval kintr - intrinsic permeability of water

   02/02/2024, TKr
*/
double con_hawf3mat::get_kintrw(double pw, double pg, double t, long ipp)
{
  double kintr;
  double phi;

  kintr = kintrw0;
  
  switch (kintr_type){
  case 0:
  case 3:{//constant
    kintr = kintrw0;
    break;
  }
  case 1:
  case 4:{//dependent on porosity
    phi = get_porosity(ipp);
    kintr = kintrw0*exp(bb1*(phi - phi01));
    break;
  }
  case 2:{//dependent on porosity - cubic and quadratic
    //Kozeny's approach for bentonite:
    phi = get_porosity(ipp);
    kintr = kintrw0*phi*phi*phi/(1 - phi01)*(1 - phi01)*(1 - phi)*(1 - phi)/(phi01*phi01*phi01);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(kintr);
}



/**
   function computes intrinsic permeability assumed only for gas (some models have two separate values)
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval kintr - intrinsic permeability of gas

   02/02/2024, TKr
*/
double con_hawf3mat::get_kintrg(double pw, double pg, double t, long ipp)
{
  double kintr;
  double phi,e;

  kintr = kintrg0;

  switch (kintr_type){
  case 0:
  case 3:{//constant
    kintr = kintrg0;
    break;
  }
  case 1:
  case 4:{//dependent on porosity
    phi = get_porosity(ipp);
    e = phi/(1 - phi);
    kintr = kintrg0*pow(e,ng);
    break;
  }
  case 2:{//dependent on porosity - cubic and quadratic
    //Kozeny's approach for bentonite:
    //phi = get_porosity(ipp);
    //kintr = kintrw0*phi*phi*phi/(1 - phi01)*(1 - phi01)*(1 - phi)*(1 - phi)/(phi01*phi01*phi01);
    kintr = kintrg0;//temp.
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(kintr);
}

/**
   function computes water relative permeability
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval krw - water relative permeability

   20/05/2017, TKr, revised 03/10/2023 TKr
*/
double con_hawf3mat::get_krw(double pw, double pg, double t, long ipp)
{
  double pc,sw,sef,krw;
  krw = 1.0;

  switch (krw_type){
  case 0:{//constant
    krw = krw0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,pg,t,ipp);
    krw = (sw-sirr)/(ssat-sirr);
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,pg,t,ipp);
    krw = (sw-sirr)*(sw-sirr)*(sw-sirr)/((ssat-sirr)*(ssat-sirr)*(ssat-sirr));
    break;
  }
  case 3:{//exponential
    sw = get_sw(pw,pg,t,ipp);
    krw = pow(sw,lambda_krw);
    break;
  }
  case 4:{//liakopoulos
    sw = get_sw(pw,pg,t,ipp);
    krw = 1.0-2.207*pow((1.0-sw),1.0121);
    break;
  }
  case 5:{//double exponential
    sw = get_sw(pw,pg,t,ipp);
    sef = (sw-sirr)/(ssat-sirr); //effective saturation degree
    krw = pow(sef,(1.0/beta_krw));
    krw = pow((1.0-krw),beta_krw);
    krw = pow((1.0-krw),2.0);
    krw = pow(sef,0.5)*krw;
    break;
  }
  case 6:{//FEBEX granit
    sw = get_sw(pw,pg,t,ipp);
    krw = febex_granit_ret.get_krw(sw);
    break;
  }
  case 7:{
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;
    
    krw = van_genuchten_ret.get_krw(-pc,t);
  }
  case 9:{//bazant
    sw = get_sw(pw,pg,t,ipp);
    krw = exp(10.0*log(sw));
    break;
  }

  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(krw);
}



/**
   function computes air relative permeability
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval krg - gas relative permeability

   20/05/2017, TKr, revised 03/10/2023 TKr
*/
double con_hawf3mat::get_krg(double pw, double pg, double t, long ipp)
{
  double krg;
  double sw,n,kg,mug,g,rhog,kintr,e;
  
  krg=1.0;
  
  switch (krg_type){
  case 0:{//constant
    krg = krg0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,pg,t,ipp);
    krg = 1.0 - pow((sw/s_crit),ag);
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,pg,t,ipp);
    krg = pow((1.0-sw),ag);
    break;
  }
  case 3:{//exponential dependence of Kg on void ration and saturation degree - FEBEX bentonite
    sw = get_sw(pw,pg,t,ipp);
    if (sw < 0.9999){ //because of numerical issues
      n = get_porosity(ipp);
      e = n/(1.0-n);
      mug = get_mug(pw,pg,t);
      g = 9.81;
      rhog = get_rhog(pw,pg,t);
      kintr = get_kintr(pw,pg,t,ipp);
      kg = kg0*pow((e*(1.0-sw)),kgn);
      krg = kg*mug/(g*rhog*kintr);
    }
    else
      krg = 1.0e-10;
    break;
  }

  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(krg);
}


/**
   function computes Biot's constant
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval alpha - Biot's constant

   03/10/2023 TKr
*/
double con_hawf3mat::get_alpha(double /*pw*/, double /*pg*/, double /*t*/,long /*ipp*/)
{
  double alpha;
  //double kt,ks;
  
  //kt = get_kt(pw,pg,t,ipp);
  //ks = get_ks(pw,pg,t,ipp);
  
  //alpha = 1.0 - kt/ks;  

  alpha = alpha0;//temp.
  
  return(alpha);
}

/**
   function computes bulk modulus of solid phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ks - bulk modulus of solid phase

   03/10/2023 TKr
*/
double con_hawf3mat::get_ks(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double ks;

  ks = ks0;//temp.

  return(ks);
}


/**
   function computes bulk modulus of porous medium
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval kt - bulk modulus of porous medium

   03/10/2023 TKr
*/
double con_hawf3mat::get_kt(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double kt;
  
  kt = kt0;//temp.
  
  return(kt);
}


/**
   function computes mass concentration of water vapour air in gas phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhogw - mass concentration of water vapour air in gas phase
*/
double con_hawf3mat::get_rhogw(double pw,double pg,double t)
{
  double rhogw,pgw;

  pgw = get_pgw(pw,pg,t);

  rhogw = pgw/t/gasr*mw;

  return(rhogw);
}

/**
   function computes water vapour specific heat
   @param t - temperature

   @retval cpgw - water vapour specific heat
*/
double con_hawf3mat::get_cpgw()
{
  double cpgw;

  cpgw = 1805.0;

  return(cpgw);
}

/**
   function computes water vapour pressure = Kelvin equation
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval pgw - water vapour pressure = Kelvin equation
*/
double con_hawf3mat::get_pgw(double pw,double pg,double t)
{
  double pgw,rhow,pgws,tt,arg,pc;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);

  if(rel_gas_press == 1)
    pg = pg + p0;
  
  pc = pg - pw;
  
  arg = -1.0*pc*mw/rhow/gasr/tt;
  
  if(arg < -7){
    pgw = pgws*0.001;
    //printf("relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
    //fprintf(Outt,"relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
    //print_err("relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
  }
  else{
    //pgw = pgws*exp(-1.0*pc*mw/rhow/gasr/tt);
    pgw = pgws*exp(arg);
  }
  
  check_math_err();
  
  return(pgw);
}



/**
   function computes partial derivative of pgw with respect to pc (Kelvin equation)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   
   @retval dpgw_dpc - partial derivative of pgw with respect to pc (Kelvin equation)
*/
double con_hawf3mat::get_dpgw_dpc(double pw,double pg,double t)
{
  double dpgw_dpc,pgws,rhow,tt,pc,arg;
  
  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;
  
  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  
  if(rel_gas_press == 1)
    pg = pg + p0;
  
  pc = pg - pw;
  
  arg = -1.0*pc*mw/rhow/gasr/tt;
  
  if(arg < -7){
    dpgw_dpc = 0.0;
    //printf("relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
    //fprintf(Outt,"relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
    //print_err("relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
  }
  else{
    //dpgw_dpc = -pgws*exp(-pc*mw/rhow/gasr/tt)*mw/rhow/gasr/tt;
    dpgw_dpc = -pgws*exp(arg)*mw/rhow/gasr/tt;
  }
  
  check_math_err();
  
  return(dpgw_dpc);
}


/**
   function computes partial derivative of pgw with respect to t (Kelvin equation)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval dpgw_dt - partial derivative of pgw with respect to t (Kelvin equation)
*/
double con_hawf3mat::get_dpgw_dt(double pw,double pg,double t)
{
  double dpgw_dt,dpgws_dt,tt,pgws,rhow,pc,drhow_dt,arg;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  pc = pg - pw;
  dpgws_dt = get_dpgws_dt(tt);
  
  drhow_dt = get_drhow_dt(tt);
  
  arg = -1.0*pc*mw/rhow/gasr/tt;
  
  if(arg < -7){
    dpgw_dt = 0.0;
    if(t < tcr){
      dpgw_dt = dpgws_dt*0.0001 + pgws*0.0001*(7/tt);
    }
    else{
      dpgw_dt = dpgws_dt*0.0001;
    }
    //printf("relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
    //fprintf(Outt,"relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
    //print_err("relative humidity correction to limit value 0.001",__FILE__,__LINE__,__func__);
  }
  else{
    if(t < tcr){
      dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt) + pgws*exp(-pc*mw/rhow/gasr/tt)*(pc*mw/rhow/gasr/tt/tt);
      //actualized form:
      //dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt) + pgws*exp(-pc*mw/rhow/gasr/tt)*(pc*mw/rhow/gasr/tt)*(drhow_dt/rhow + 1/tt);
    }   
    else
      dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt);
  }

  return(dpgw_dt);
}



/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure
*/
double con_hawf3mat::get_pgws(double t)
{
  /* //Clausius-Clapeyron equation
     double pgws,pgws0,t0,dhvap;

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


  check_math_err();

  return(pgws);
}



/**
   function computes partial derivative of water vapour saturation pressure with respect to t
   @param t - temperature

   @retval dpgws_dt - partial derivative of water vapour saturation pressure with respect to t
*/
double con_hawf3mat::get_dpgws_dt(double t)
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
   function computes water density
   @param t - temperature

   @retval rhow - water density

   03/10/2023 TKr
*/
double con_hawf3mat::get_rhow(double t)
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

  }
  
  return(rhow);
}



/**
   function computes compresibility coefficient of water
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval kw - compresibility coefficient of water

   03/10/2023 TKr
*/
double con_hawf3mat::get_kw(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double kw;
  
  kw = kw0;//temp.

  return(kw);
}



/**
   function volume computes thermal expansion coefficient of solid - gas for compressible solid grains
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval betasg - volume thermal expansion coefficient of solid - gas for compressible solid grains

   03/10/2023 TKr
*/
double con_hawf3mat::get_betasg(double pw,double pg,double t,long ipp)
{
  double alpha,betasg,betas,n,sw,sg;

  alpha = get_alpha(pw,pg,t,ipp);
  betas = get_betas(ipp);
  n = get_porosity(ipp);
  sw = get_sw(pw,pg,t,ipp);
  sg = 1.0 -sw;

  betasg = betas*(alpha - n)*sg;

  return(betasg);
}



/**
   function volume computes thermal expansion coefficient of solid - water for compressible solid grains
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval betasw - volume thermal expansion coefficient of solid - water for compressible solid grains

   03/10/2023 TKr
*/
double con_hawf3mat::get_betasw(double pw,double pg,double t,long ipp)
{
  double alpha,betasw,betas,n,sw,betaw;

  alpha = get_alpha(pw,pg,t,ipp);
  betas = get_betas(ipp);
  n = get_porosity(ipp);
  sw = get_sw(pw,pg,t,ipp);
  betaw = get_betaw(t);

  betasw = sw*(betas*(alpha - n) + n*betaw);

  return(betasw);
}



/**
   function volume computes thermal expansion coefficient of solid - water vapor for compressible solid grains
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval betasgw - volume thermal expansion coefficient of solid - water vapor for compressible solid grains

   03/10/2023 TKr
*/
double con_hawf3mat::get_betasgw(double pw,double pg,double t,long ipp)
{
  double alpha,betasgw,betas,n,sw,sg,rhogw,rhow,betaw;

  alpha = get_alpha(pw,pg,t,ipp);
  betas = get_betas(ipp);
  n = get_porosity(ipp);
  sw = get_sw(pw,pg,t,ipp);
  sg = 1.0 -sw;
  rhogw = get_rhogw(pw,pg,t);
  rhow = get_rhow(t);
  betaw = get_betaw(t);
  
  betasgw = betas*(alpha - n)*(sg*rhogw + sw*rhow) + n*betaw*rhow*sw;

  return(betasgw);
}


/**
   function computes derivative of water density with respect to temperature
   @param t - temperature

   @retval drhow_dt - derivative of rhow with respect to t
*/
double con_hawf3mat::get_drhow_dt(double t)
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
    
  drhow_dt = 0.0;//temp.

  return(drhow_dt);
}

/**
   function computes enthalpy of evaporation (latent heat of vaporization)
   @param t - temperature

   @retval - enthalpy of evaporation (latent heat of vaporization)
*/
double con_hawf3mat::get_dhvap(double t)
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
   
   03/10/2023 TKr
*/
double con_hawf3mat::get_muw(double t)
{
  double muw;
  
  muw = muw0*pow((t - conb),conc);

  //check_math_errel(0);

  //muw = muw0const;//constant value

  return(muw);
}

/**
   function computes water specific heat

   @retval cpw - water specific heat

   03/10/2023 TKr
*/
double con_hawf3mat::get_cpw()
{
  double c;

  c = cpw0;//temp.

  return(c);
}

/**
   function computes water heat conductivity
   @param t - temperature
   
   @retval lambdaw - water heat conductivity

   03/10/2023 TKr
*/
double con_hawf3mat::get_lambdaw(double /*t*/)
{
  double lam;

  lam = lambdaw;//temp.

  return(lam);
}

/**
   function computes volume thermal expansion coefficient of water
   @param t - temperature
   
   @retval betaw - volume thermal expansion coefficient of water


   03/10/2023 TKr
*/
double con_hawf3mat::get_betaw(double t)
{
  double betaw;
  
  //betaw = 0.68e-4;//[K-1] at t=273.15
  //betaw = 10.1e-4;//[K-1] at t=420.0
  
  //betaw = 0.63e-5;//temp.
  betaw = 0.68 + (10.1 - 0.68)/(420.0 - 273.15)*(t - 273.15);//linear expression (rough)
  betaw = betaw*1.0e-4;

  return(betaw);
}



/**
   function computes molar mass of moist air
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval mg - molar mass of moist air
*/
double con_hawf3mat::get_mg(double pw,double pg,double t)
{
  double mg,pgw;

  pgw = get_pgw(pw,pg,t);

  if(rel_gas_press == 1)
    pg = pg + p0;

  mg = ma + (mw-ma)*pgw/pg;
  
  return(mg);
}


/**
   function computes gas phase density
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhog - gas phase density
*/
double con_hawf3mat::get_rhog(double pw,double pg,double t)
{
  double rhog,pgw;

  pgw = get_pgw(pw,pg,t);

  if(rel_gas_press == 1)
    pg = pg + p0;

  rhog = (pg*ma + (mw - ma)*pgw)/gasr/t;
  
  return(rhog);
}



/**
   function computes dynamic viscosity of moist air
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval mug - dynamic viscosity of moist air
*/
double con_hawf3mat::get_mug(double pw,double pg,double t)
{
  double mug,mugw,muga,pga,pgw;

  pgw = get_pgw(pw,pg,t);

  if(rel_gas_press == 1)
    pg = pg + p0;

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
   function computes dynamic viscosity of dry air
   @param t - temperature

   @retval muga - dynamic viscosity of dry air
*/
double con_hawf3mat::get_muga(double t)
{
  double muga;

  muga = muga0 + alphaa*(t-t0) - betaa*pow((t-t0),2.0);

  return(muga);
}

/**
   function computes dynamic viscosity of water vapour
   @param t - temperature

   @retval mugw - dynamic viscosity of water vapour
*/
double con_hawf3mat::get_mugw(double t)
{
  double mugw;
 
  mugw = mugw0 + alphaw*(t-t0);

  return(mugw);
}



/**
   function computes mass concentration of dry air in gas phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhoga - mass concentration of dry air in gas phase
*/
double con_hawf3mat::get_rhoga(double pw,double pg,double t)
{
  double rhoga,pgw;

  //gas pressure check
  pgw = get_pgw(pw,pg,t);

  if(rel_gas_press == 1)
    pg = pg + p0;

  if (pgw <= pg)
    rhoga = (pg - pgw)*ma/gasr/t;
  else
    rhoga = 0.0;

  return(rhoga);
}


/**
   function computes dry air specific heat

   @retval cpga - dry air specific heat
*/
double con_hawf3mat::get_cpga()
{
  double cpga;

  cpga = 1005.7;

  return(cpga);
}


/**
   function computes effective thermal capacity of partially saturated medium
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval rhocp - effective thermal capacity of partially saturated medium

   03/10/2023 TKr
*/
double con_hawf3mat::get_rhocp(double pw,double pg,double t,long ipp)
{
  double sw,n,rhocp,rhow,rhog,rhogw,cps,cpw,cpga,cpgw,rhos,cp;

  sw = get_sw(pw,pg,t,ipp);
  n = get_porosity(ipp);
  
  rhow = get_rhow(t);
  rhog = get_rhog(pw,pg,t);
  rhogw = get_rhogw(pw,pg,t);
  
  cps = get_cps(pw,pg,t,ipp);
  
  cpw = get_cpw();
  cpga = get_cpga();
  cpgw = get_cpgw();
  
  rhos = get_rhos(t);

  switch (thermal_capacity_type){
  case 0:{//constant value without water and gas phases
    rhocp = (1.0-n)*rhos*cps;// constant effective heat capacity only for solid phase
    break;
  }
  case 1:{//effective heat capacity from all components
    rhocp = (1.0-n)*rhos*cps;// solid phase = rhod*cps
    rhocp = rhocp + n*(sw*rhow*cpw + (1.0-sw)*(rhog*cpga + rhogw*(cpgw-cpga)));//liquid and gas phases
    break;
  }
  case 2:{
    rhocp = rhos*cps;// special case - constant heat capacity only for solid phase
    break;
  }
  case 3:{//effective heat capacity is directly measured with respec to saturation degree
    sw = get_sw(pw,pg,t,ipp);
    //linear approximation
    rhocp = rhocp_dry + (rhocp_wet - rhocp_dry)*(sw - sr_dry)/(sr_wet - sr_dry);
    if(rhocp < rhocp_dry)
      rhocp = rhocp_dry;
    if(rhocp > rhocp_wet)
      rhocp = rhocp_wet;  
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  //this below is temporarilly for mock_up experiment testing
  switch (model_type){
  case lewis_and_schrefler3:{//Lewis and Schrefler's book
    //rhocp = rhos0*cps0;//temporarily for testing??!!
    break;
  }
  case lewis_and_schrefler3_mefel:{//from mefel;
    //rhocp = rhos0*cps0;//temporarily for testing??!!
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(rhocp);
}


/**
   function computes air specific heat
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval cpg - air specific heat
*/
double con_hawf3mat::get_cpg(double pw,double pg,double t)
{
  double cpg,rhog,rhogw,cpga,cpgw;

  rhog = get_rhog(pw,pg,t);
  cpga = get_cpga();
  rhogw = get_rhogw(pw,pg,t);
  cpgw = get_cpgw();
  
  cpg = (rhog*cpga + rhogw*(cpgw-cpga))/rhog;
  
  return(cpg);
}



/**
   function computes effective thermal conductivity of partially saturated concrete
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval lambdaeff - effective thermal conductivity of partially saturated soil

   03/10/2023 TKr
*/
double con_hawf3mat::get_lambdaeff(double pw,double pg,double t,long ipp)
{
  double lambdaeff;
  double sw;

  lambdaeff = lambda_eff0;
 
  switch (lambda_type){
  case 0:{//constant
    lambdaeff = lambda_eff0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,pg,t,ipp);
    lambdaeff = lambda_dry + (lambda_wet - lambda_dry)*sw; //this is linear approximation for bentonite:
    break;
  }
  case 2:{// linear approximation
    sw = get_sw(pw,pg,t,ipp);
    lambdaeff = lambda_dry + (lambda_wet - lambda_dry)*(sw - sr_dry)/(sr_wet - sr_dry);
    if(lambdaeff < lambda_dry)
      lambdaeff = lambda_dry;
    if(lambdaeff > lambda_wet)
      lambdaeff = lambda_wet; 
    //fprintf(Outt,"lambda = %e, sr = %e\n",lambdaeff,sw);
    //fflush(Outt);
    break;
  }
  case 3:{//
    sw = get_sw(pw,pg,t,ipp);
    lambdaeff = pow(lambda_wet,sw)*pow(lambda_dry,(1-sw)); //this is from Patek's thesis p. 18, but it seems to be strange
    break;
  }
  case 4:{//complex
    //lambdaeff = lambdas + lambdaw + lambdag; //not finished
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(lambdaeff);
}


/***********************************************************/
/** The end of basic functions and constitutive relations **/
/***********************************************************/




/**
   function computes conductivity matrix of the material  in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
*/
void con_hawf3mat::matcond (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;
  
  switch (n){
  case 1:{
    matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err ("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}


/**
   function creates conductivity matrix of the material for 1D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pw,pg,t,ipp);// *scale_t;//scaling
  
  d[0][0] = k;
  check_math_errel(ipp);
}

/**
   function creates conductivity matrix of the material for 2D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pw,pg,t,ipp);// *scale_t;//scaling

  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
 
  check_math_errel(ipp);
}

/**
   function creates conductivity matrix of the material for 3D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void con_hawf3mat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    k = get_kwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    k = get_kgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    k = get_ktw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    k = get_ktg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    k = get_ktt1(pw,pg,t,ipp);// *scale_t;//scaling
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;

  check_math_errel(ipp);
}


/**
   function creates capacity matrix of the material
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::matcap (double &c,long ri,long ci,long ipp)
{
  double pw,pg,t;
  c = 0.0;

  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capww(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = get_capwg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 0) && (ci == 2))
    c = get_capwt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    c = get_capgw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    c = get_capgg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 1) && (ci == 2))
    c = get_capgt(pw,pg,t,ipp);// *scale_t;//scaling
  
  if((ri == 2) && (ci == 0))
    c = get_captw(pw,pg,t,ipp);// *scale_pw;//scaling
  if((ri == 2) && (ci == 1))
    c = get_captg(pw,pg,t,ipp);// *scale_g;//scaling
  if((ri == 2) && (ci == 2))
    c = get_captt(pw,pg,t,ipp);// *scale_t;//scaling
  
  check_math_errel(ipp);
}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void con_hawf3mat::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    matcond1d_2 (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d_2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d_2 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err ("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}




/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::matcond1d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c,dd;
  double pw,pg,t;
      
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 2) && (ci == 0)){
  }

  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pw,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,pg,t,ipp);// *scale_t;//scaling
    dd = get_ktt2d(pw,pg,t,ipp);// *scale_t;//scaling
    
    d[0][0] = -1.0*a*Tm->ip[ipp].grad[0][0] - b*Tm->ip[ipp].grad[1][0] + (a*c+b*dd)*Tp->gr[0];//zkontrolovat??!!
  }  
}


/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::matcond2d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c,dd;
  double pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 2) && (ci == 0)){
  }  
  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pw,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,pg,t,ipp);// *scale_t;//scaling
    dd = get_ktt2d(pw,pg,t,ipp);// *scale_t;//scaling
    
    d[0][0] = -1.0*a*Tm->ip[ipp].grad[0][0] - b*Tm->ip[ipp].grad[1][0] + (a*c+b*dd)*Tp->gr[0];//zkontrolovat??!!
    d[0][1] = -1.0*a*Tm->ip[ipp].grad[0][1] - b*Tm->ip[ipp].grad[1][1] + (a*c+b*dd)*Tp->gr[1];//zkontrolovat??!!
  }

  check_math_errel(ipp);
}



/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::matcond3d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c,dd;
  double pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  fillm(0.0,d);
    
  if((ri == 2) && (ci == 0)){
  }  
  if((ri == 2) && (ci == 1)){
  }
  
  if((ri == 2) && (ci == 2)){
    a = get_ktt2a(pw,pg,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,pg,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,pg,t,ipp);// *scale_t;//scaling
    dd = get_ktt2d(pw,pg,t,ipp);// *scale_t;//scaling

    d[0][0] = -1.0*a*Tm->ip[ipp].grad[0][0] - b*Tm->ip[ipp].grad[1][0] + (a*c+b*dd)*Tp->gr[0];//zkontrolovat??!!
    d[0][1] = -1.0*a*Tm->ip[ipp].grad[0][1] - b*Tm->ip[ipp].grad[1][1] + (a*c+b*dd)*Tp->gr[1];//zkontrolovat??!!
    d[0][2] = -1.0*a*Tm->ip[ipp].grad[0][2] - b*Tm->ip[ipp].grad[1][2] + (a*c+b*dd)*Tp->gr[2];//zkontrolovat??!!
  }
}



/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   11.2.2003
*/
void con_hawf3mat::rhs_volume (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;
  
  switch (m){
  case 1:{
    rhs1d1 (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    rhs2d1 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    rhs3d1 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err ("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}



/**
   function computes volume part 2 of right-hand side matrix
   in the required integration point
   
   @param cc - right-hand side coefficient of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::rhs_volume2 (double &cc,long ri,long /*ci*/,long ipp)
{
  double pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t  = Tm->ip[ipp].av[2];

  if(ri == 0){
    cc = get_fwu(pw,pg,t,ipp);
  }
  if(ri == 1){
    cc = get_fgu(pw,pg,t,ipp);
  }
  if(ri == 1){
    cc = get_ftu(pw,pg,t,ipp);
  }  
}




/**
   function creates volume right-hand side matrix of the material for 1D problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::rhs1d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fw1(pw,pg,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
  if(ri == 1){
    f = get_fg(pw,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
  if(ri == 2){
    f = get_ft1(pw,pg,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
}

/**
   function creates volume right-hand side matrix of the material for 2D problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::rhs2d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fw1(pw,pg,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
  if(ri == 1){
    f = get_fg(pw,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
  if(ri == 2){
    f = get_ft1(pw,pg,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
}

/**
   function creates volume right-hand side matrix of the material for 3D problems

   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hawf3mat::rhs3d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg,t;

  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  pg = Tm->ip[ipp].av[1];// *scale_pg;//scaling
  t  = Tm->ip[ipp].av[2];// *scale_t;//scaling

  if(ri == 0){
    f = get_fw1(pw,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
  if(ri == 1){
    f = get_fg(pw,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
  if(ri == 2){
    f = get_ft1(pw,pg,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
}




/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval fwu - volumetric strain rate effect on pore water pressure; only for partially coupled METR version

   13/07/2018, TKr revised 03/10/2023 TKr
*/
double con_hawf3mat::get_fwu(double pw, double pg, double t, long ipp)
{
  double depsv_r,sw,fwu,dsr_depsv,n,sg,rhow,rhogw,alpha;
    
  fwu = 0.0;
  
  switch (model_type){
  case artificial3:
  case lewis_and_schrefler3:{
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book
    if(vol_strain_effect == 1){
      n = get_porosity(ipp);
      depsv_r = Tm->givenontransq(strain_vol_rate, ipp);      //actual rate of volumetric strain from MEFEL
      Tm->ip[ipp].eqother[0] = depsv_r; //this is not necessary
      
      
      /*
	if (sr_type == mefel_sr){
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	Tm->ip[ipp].eqother[1] = dsr_depsv; //this is not necessary
	}
      */
      
      switch (sr_type){
      case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
	dsr_depsv = 0.0;
	break;
      }
      case gardner_exponential_sr:{//dependent on saturation degree:
	dsr_depsv = 0.0;
	break;
      }
      case potts_log_linear_sr:{//exponential
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case mefel_sr:{//from MEFEL
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	break;
      }
      case table_sr:{//from table
	dsr_depsv = 0.0;
	break;
      }
      case masin_sr:{//extended formulation from Brooks and Correy according to Masin
	double e=0.0,dpw=0.0,por=0.0;
	if(Tm->nontransq != NULL){
	  por = Tm->givenontransq(porosity, ipp);// from mefel
	  e = por/(1-por);
	}
	else{
	  e = phi0/(1-phi0);
	}
	dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
	dsr_depsv = masin_ret.dsw_depsv(-pw,-dpw,e,t);//positive value of suction
	break;
      }
      default:{
	print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
	abort();
      }
      } 

      
      sw = get_sw(pw,pg,t,ipp);
      sg = 1.0 - sw;
      rhow = get_rhow(t);
      rhogw = get_rhogw(pw,pg,t);
      alpha = get_alpha(pw,pg,t,ipp);
      
      fwu = -1.0*(sg*rhogw + sw*rhow)*alpha*depsv_r;//volumetric strain effect        

      if(wrc_vol_strain_effect == 1)
	fwu = fwu - n*(rhow - rhogw)*dsr_depsv*depsv_r;//volumetric strain effect on retention curve
    }
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }    
  return(fwu);
}




/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval fgu - volumetric strain rate effect on pore gas pressure; only for partially coupled METR version

   13/07/2018, TKr revised 03/10/2023 TKr
*/
double con_hawf3mat::get_fgu(double pw, double pg, double t, long ipp)
{
  double depsv_r,sw,fgu,dsr_depsv,n,sg,rhoga,alpha;
  
  fgu = 0.0;
  
  switch (model_type){
  case artificial3:
  case lewis_and_schrefler3:{
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book
    if(vol_strain_effect == 1){
      n = get_porosity(ipp);
      depsv_r = Tm->givenontransq(strain_vol_rate, ipp);      //actual rate of volumetric strain from MEFEL
      Tm->ip[ipp].eqother[0] = depsv_r; //this is not necessary
      
      
      /*
	if (sr_type == mefel_sr){
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	Tm->ip[ipp].eqother[1] = dsr_depsv; //this is not necessary
	}
      */
      
      switch (sr_type){
      case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
	dsr_depsv = 0.0;
	break;
      }
      case gardner_exponential_sr:{//dependent on saturation degree:
	dsr_depsv = 0.0;
	break;
      }
      case potts_log_linear_sr:{//exponential
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case mefel_sr:{//from MEFEL
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	break;
      }
      case table_sr:{//from table
	dsr_depsv = 0.0;
	break;
      }
      case masin_sr:{//extended formulation from Brooks and Correy according to Masin
	double e=0.0,dpw=0.0,por=0.0;
	if(Tm->nontransq != NULL){
	  por = Tm->givenontransq(porosity, ipp);// from mefel
	  e = por/(1-por);
	}
	else{
	  e = phi0/(1-phi0);
	}
	dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
	dsr_depsv = masin_ret.dsw_depsv(-pw,-dpw,e,t);//positive value of suction
	break;
      }
      default:{
	print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
	abort();
      }
      } 

      sw = get_sw(pw,pg,t,ipp);
      sg = 1.0 - sw;
      rhoga = get_rhoga(pw,pg,t);
      alpha = get_alpha(pw,pg,t,ipp);
      fgu = -alpha*sg*rhoga*depsv_r;

      if(wrc_vol_strain_effect == 1)
	fgu = fgu + n*rhoga*dsr_depsv*depsv_r;//volumetric strain rate effect
    }
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }    
  return(fgu);
}




/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ftu - volumetric strain rate effect on temperature; only for partially coupled METR version

   13/07/2018, TKr revised 03/10/2023 TKr
*/
double con_hawf3mat::get_ftu(double pw, double pg, double t, long ipp)
{
  double depsv_r,sw,ftu,dsr_depsv,n,rhow,dhvap,alpha;
  
  ftu = 0.0;
  
  
  switch (model_type){
  case artificial3:
  case lewis_and_schrefler3:{
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book
    if(vol_strain_effect == 1){
      n = get_porosity(ipp);
      depsv_r = Tm->givenontransq(strain_vol_rate, ipp);      //actual rate of volumetric strain from MEFEL
      Tm->ip[ipp].eqother[0] = depsv_r; //this is not necessary
      
      
      /*
	if (sr_type == mefel_sr){
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	Tm->ip[ipp].eqother[1] = dsr_depsv; //this is not necessary
	}
      */
      
      switch (sr_type){
      case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
	dsr_depsv = 0.0;
	break;
      }
      case gardner_exponential_sr:{//dependent on saturation degree:
	dsr_depsv = 0.0;
	break;
      }
      case potts_log_linear_sr:{//exponential
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case mefel_sr:{//from MEFEL
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	break;
      }
      case table_sr:{//from table
	dsr_depsv = 0.0;
	break;
      }
      case masin_sr:{//extended formulation from Brooks and Correy according to Masin
	double e=0.0,dpw=0.0,por=0.0;
	if(Tm->nontransq != NULL){
	  por = Tm->givenontransq(porosity, ipp);// from mefel
	  e = por/(1-por);
	}
	else{
	  e = phi0/(1-phi0);
	}
	dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
	dsr_depsv = masin_ret.dsw_depsv(-pw,-dpw,e,t);//positive value of suction
	break;
      }
      default:{
	print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
	abort();
      }
      } 

      sw = get_sw(pw,pg,t,ipp);
      rhow = get_rhow(t);
      dhvap = get_dhvap(t);
      alpha = get_alpha(pw,pg,t,ipp);
      
      ftu = dhvap*sw*rhow*alpha*depsv_r;
      if(wrc_vol_strain_effect == 1)
	ftu = ftu + dhvap*rhow*n*dsr_depsv*depsv_r;//volumetric strain rate effect
    }
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }    
  return(ftu);
}



/**
   function computes new transmission coefficient (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double con_hawf3mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,pg,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_transcoeff_tt(pw,pg,t,bc,ipp);// *scale_t;//scaling

  c = c*trc;

  return (c);
}


/**
   function computes new transmission coefficient (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double con_hawf3mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp, int flag)
{
  double c,pw,pg,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,pg,t,bc,ipp,flag);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_transcoeff_tt(pw,pg,t,bc,ipp);// *scale_t;//scaling

  c = c*trc;

  return (c);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double con_hawf3mat::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);
  
  if((ri == 0) && (ci == 0))
    c = get_transmission_nodval_ww(nodval,pw,pg,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_nodval_tt(nodval,trc2,pw,pg,t,bc,ipp);// *scale_t;//scaling

  return (c);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double con_hawf3mat::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg,t;
  c = 0.0;

  pw = nodalval (nn,0);
  pg = nodalval (nn,1);
  t = nodalval (nn,2);
  
  if((ri == 0) && (ci == 0))
    c = get_transmission_flux_ww(nodval,pw,pg,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  if((ri == 0) && (ci == 2))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = 0.0;
  if((ri == 1) && (ci == 2))
    c = 0.0;
  
  if((ri == 2) && (ci == 0))
    c = 0.0;
  if((ri == 2) && (ci == 1))
    c = 0.0;
  if((ri == 2) && (ci == 2))
    c = get_transmission_flux_tt(nodval,trc2,pw,pg,t,bc,ipp);// *scale_t;//scaling

  return (c);
}


/**
   function checks if gas pressure is greater than vapour pressure

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   TKr, 18.5.2005
   TKr 06/12/2022 - actualized
*/
void con_hawf3mat::gaspress_check(double pw,double &pg,double t,long ipp)
{
  double pgw;
  
  //general
  //gas pressure check
  pgw = get_pgw(pw,pg,t);

  if (pgw > pg){
    pg = pgw;
    //fprintf(Outt,"gas pressure was modified in integration point No. %ld\n",ipp);
    //print_err("gas pressure was modified in integration point No. %ld",__FILE__,__LINE__,__func__,ipp);
  }

  //general
  //fully saturated state
  if(pgw < 100.0){
    pg = 101325.0;
  }

  // M. Starnoni 24-11-2010
  if (pg < 0.0){
    print_err("gas pressure was modified in integration point No. %ld",__FILE__,__LINE__,__func__,ipp);
    pg = 0.0;
  }

  check_math_err();
}


/**
   function checks pore water pressure

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
   TKr 06/12/2022 - actualized
*/
void con_hawf3mat::waterpress_check(double &pw,double pg,double t,long ipp)
{
  if(pw < -5.0e8){
    //printf("water pressure was modified from value pw = %e to -5.0e8 in integration point No. %ld\n",pw,ipp);
    //fprintf(Outt,"water pressure was modified from value pw = %e to -5.0e8 in integration point No. %ld\n",pw,ipp);
    //print_err("water pressure was modified from value pw = %e to -5.0e8 in integration point No. %ld",__FILE__,__LINE__,__func__,pw,ipp);
    pw = -5.0e8;
  }

  check_math_err();
}

/**
   function corrects values of variables
   
   @param nv - array with variables
   
*/
void con_hawf3mat::values_correction (vector &nv, long ipp)
{
  //  pore water pressure
  waterpress_check(nv[0],nv[1],nv[2],ipp);
  
  //  gas pressure
  gaspress_check(nv[0],nv[1],nv[2],ipp);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capww - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_capww(double pw,double pg,double t,long ipp)
{
  double capww;
  double alpha,n,ks,sw,rhogw,sg,rhow,kw,dpgw_dpc,ds_dpc;
  
  if(model_type == artificial3)
    capww = capww0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);

    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    rhogw = get_rhogw(pw,pg,t);
    sg = 1.0 - sw;
    rhow = get_rhow(t);
    kw = get_kw(pw,pg,t,ipp);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    ds_dpc = get_ds_dpc(pw,pg,t,ipp);//this is equal to minus derivative with respect to suction
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      capww = (alpha - n)/ks*sw*(rhogw*sg + rhow*sw) + rhow*sw*n/kw;
      capww = capww - sg*n*mw/t/gasr*dpgw_dpc;
      capww = capww - ((alpha - n)/ks*(rhogw*sg*(pw-pg) + rhow*sw*(pw-pg)) + n*(rhow - rhogw))*ds_dpc;
    }
    else{
      //incompressible grains:
      capww = -sg*n*mw/t/gasr*dpgw_dpc;
      capww = capww - n*(rhow - rhogw)*ds_dpc;    
    }
  }

  return(capww);
}



/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capwg - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_capwg(double pw,double pg,double t,long ipp)
{
  double capwg;
  double alpha,n,ks,sw,rhogw,sg,rhow,kw,dpgw_dpc,ds_dpc;

  if(model_type == artificial3)
    capwg = capwg0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    rhogw = get_rhogw(pw,pg,t);
    sg = 1.0 - sw;
    rhow = get_rhow(t);
    kw = get_kw(pw,pg,t,ipp);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    ds_dpc = get_ds_dpc(pw,pg,t,ipp);//this is equal to minus derivative with respect to suction
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      capwg = (alpha - n)/ks*sg*(rhogw*sg + rhow*sw);
      capwg = capwg + sg*n*mw/t/gasr*dpgw_dpc;
      capwg = capwg + ((alpha - n)/ks*(rhogw*sg*(pw-pg) + rhow*sw*(pw-pg)) + n*(rhow - rhogw))*ds_dpc;
    }    
    else{
      //incompressible grains:
      capwg = sg*n*mw/t/gasr*dpgw_dpc;
      capwg = capwg + n*(rhow - rhogw)*ds_dpc;    
    }
  }  

  return(capwg);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capwt - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_capwt(double pw,double pg,double t,long ipp)
{
  double capwt;
  double betasgw,sw,sg,dpgw_dt,pgw,alpha,n,ks,rhogw,rhow,ds_dt;

  if(model_type == artificial3)
    capwt = capwt0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    betasgw = get_betasgw(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    sg = 1.0 - sw;
    dpgw_dt = get_dpgw_dt(pw,pg,t);
    pgw = get_pgw(pw,pg,t);
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    rhogw = get_rhogw(pw,pg,t);
    rhow = get_rhow(t);
    ds_dt = get_ds_dt(pw,pg,t,ipp);
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      capwt = -betasgw;
      capwt = capwt + sg*n*mw/t/gasr*(dpgw_dt-pgw/t);
      capwt = capwt + ((alpha - n)/ks*(rhogw*sg*(pw-pg) + rhow*sw*(pw-pg)) + n*(rhow - rhogw))*ds_dt;
    }
    else{
      //incompressible grains:
      capwt = -betasgw;
      capwt = capwt + sg*n*mw/t/gasr*(dpgw_dt-pgw/t);
      capwt = capwt + n*(rhow - rhogw)*ds_dt;    
    }
  }
  return(capwt);
}

/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kww - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_kww(double pw,double pg,double t,long ipp)
{
  double kww;
  double rhow,krw,muw,rhog,mg,dg,dpgw_dpc,kintr;

  if(model_type == artificial3)
    kww = kww0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhow = get_rhow(t);
    krw = get_krw(pw,pg,t,ipp);
    muw = get_muw(t);
    rhog =get_rhog(pw,pg,t);
    mg = get_mg(pw,pg,t);
    dg = get_dg(pw,pg,t,ipp);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrw(pw,pg,t,ipp);
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    kww = rhow*kintr*krw/muw;//water
    kww = kww - rhog*ma*mw/mg/mg*dg/pg*dpgw_dpc;//water vapour
  }
  return(kww);
}


/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kwg - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_kwg(double pw,double pg,double t,long ipp)
{
  double kwg;
  double rhogw,krg,mug,dg,pgw,rhow,krw,muw,rhog,mg,dpgw_dpc,kintr;
  
  if(model_type == artificial3)
    kwg = kwg0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhogw = get_rhogw(pw,pg,t);
    krg = get_krg(pw,pg,t,ipp);
    mug = get_mug(pw,pg,t);
    dg = get_dg(pw,pg,t,ipp);
    pgw = get_pgw(pw,pg,t);
    rhow = get_rhow(t);
    krw = get_krw(pw,pg,t,ipp);
    muw = get_muw(t);
    rhog = get_rhog(pw,pg,t);
    mg = get_mg(pw,pg,t);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrg(pw,pg,t,ipp);

    if(rel_gas_press == 1)
      pg = pg + p0;
    
    kwg = rhogw*kintr*krg/mug;//air
    kwg = kwg + rhog*ma*mw/mg/mg*dg*(dpgw_dpc/pg - pgw/pg/pg);//water vapour
  }
  
  check_math_err();

  return(kwg);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgw - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_capgw(double pw,double pg,double t,long ipp)
{
  double capgw,alpha,n,ks,sw,sg,pc,rhoga,ds_dpc,dpgw_dpc;
  
  if(model_type == artificial3)
    capgw = capgw0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    sg = 1.0 - sw;
    pc = pg - pw;
    rhoga = get_rhoga(pw,pg,t);
    ds_dpc = get_ds_dpc(pw,pg,t,ipp);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);//this is equal to minus derivative with respect to suction
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      capgw = (alpha - n)/ks*sw*sg*rhoga;
      capgw = capgw + ((alpha - n)/ks*sg*pc + n)*rhoga*ds_dpc;
      capgw = capgw + sg*n*mw/t/gasr*dpgw_dpc;
    }
    else{
      //incompressible grains:
      capgw =  n*rhoga*ds_dpc;
      capgw = capgw + sg*n*mw/t/gasr*dpgw_dpc;
    }
  }
  return(capgw);
}



/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgg - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_capgg(double pw,double pg,double t,long ipp)
{
  double capgg,alpha,n,ks,sw,sg,pc,rhoga,ds_dpc,dpgw_dpc;
  
  if(model_type == artificial3)
    capgg = capgg0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    sg = 1.0 - sw;
    rhoga = get_rhoga(pw,pg,t);
    pc = pg - pw;
    ds_dpc = get_ds_dpc(pw,pg,t,ipp);//this is equal to minus derivative with respect to suction
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      capgg = (alpha - n)/ks*sg*sg*rhoga;
      capgg = capgg - ((alpha - n)/ks*sg*pc + n)*rhoga*ds_dpc;
      capgg = capgg + sg*n*ma/t/gasr;
      capgg = capgg - sg*n*mw/t/gasr*dpgw_dpc;
    }
    else{
      //incompressible grains:
      capgg = -n*rhoga*ds_dpc;
      capgg = capgg + sg*n*ma/t/gasr;
      capgg = capgg - sg*n*mw/t/gasr*dpgw_dpc;    
    }
  }
  return(capgg);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capgt - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_capgt(double pw,double pg,double t,long ipp)
{
  double capgt,rhoga,betasg,alpha,n,ks,sw,sg,pc,ds_dt,dpgw_dt,pgw;
  
  if(model_type == artificial3)
    capgt = capgt0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhoga = get_rhoga(pw,pg,t);
    betasg = get_betasg(pw,pg,t,ipp);
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    sg = 1.0 - sw;
    pc = pg - pw;
    ds_dt = get_ds_dt(pw,pg,t,ipp);
    dpgw_dt = get_dpgw_dt(pw,pg,t);
    pgw = get_pgw(pw,pg,t);
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      capgt = -1.0*rhoga*betasg;
      capgt = capgt - ((alpha - n)/ks*sg*pc + n)*rhoga*ds_dt;
      capgt = capgt - (sg*n*ma/t/t/gasr*pg + sg*n*ma/t/gasr*(dpgw_dt-pgw/t));
    }
    else{
      //incompressible grains:
      capgt = -1.0*rhoga*betasg;
      capgt = capgt - n*rhoga*ds_dt;
      capgt = capgt - (sg*n*ma/t/t/gasr*pg + sg*n*ma/t/gasr*(dpgw_dt-pgw/t));    
    }
  }
  return(capgt);
}



/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgw - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_kgw(double pw,double pg,double t,long ipp)
{
  double kgw,rhog,dg,mg,dpgw_dpc;
  
  if(model_type == artificial3)
    kgw = kgw0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhog = get_rhog(pw,pg,t);
    dg = get_dg(pw,pg,t,ipp);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    mg = get_mg(pw,pg,t);
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    kgw = rhog*ma*mw/mg/mg*dg/pg*dpgw_dpc;//water vapour
  }
  return(kgw);
}



/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgg - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_kgg(double pw,double pg,double t,long ipp)
{
  double kgg,rhoga,krg,mug,rhog,dg,mg,dpgw_dpc,pgw,kintr;
  
  if(model_type == artificial3)
    kgg = kgg0;
  else{

    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhoga = get_rhoga(pw,pg,t);
    krg = get_krg(pw,pg,t,ipp);
    mug = get_mug(pw,pg,t);
    rhog = get_rhog(pw,pg,t);
    dg = get_dg(pw,pg,t,ipp);
    mg = get_mg(pw,pg,t);
    dpgw_dpc = get_dpgw_dpc(pw,pg,t);
    pgw = get_pgw(pw,pg,t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrg(pw,pg,t,ipp);

    if(rel_gas_press == 1)
      pg = pg + p0;
    
    kgg = rhoga*kintr*krg/mug;//air
    kgg = kgg - rhog*ma*mw/mg/mg*dg*(dpgw_dpc/pg-pgw/pg/pg);//water vapour

    // minimal value check in case of fully saturation state:
    if(kgg <= kggmin)
      kgg = kggmin;
  }

  return(kgg);
}




/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kwt - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_kwt(double pw,double pg,double t,long ipp)
{
  double kwt,rhog,mg,dg,dpgw_dt;
  if(model_type == artificial3)
    kwt = kwt0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhog =get_rhog(pw,pg,t);
    mg = get_mg(pw,pg,t);
    dg = get_dg(pw,pg,t,ipp);
    dpgw_dt = get_dpgw_dt(pw,pg,t);
    
    kwt = rhog*ma*mw/mg/mg*dg/pg*dpgw_dt;
  }
  return(kwt);
}



/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captw - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_captw(double pw,double pg,double t,long ipp)
{
  double captw,dhvap,rhow,alpha,n,ks,sw,kw,dsw_dpc;
  
  if(model_type == artificial3)
    captw = captw0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    kw = get_kw(pw,pg,t,ipp);
    dsw_dpc = get_ds_dpc(pw,pg,t,ipp);//this is equal to minus derivative with respect to suction
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      captw = -1.0*dhvap*rhow*((alpha - n)/ks*sw*sw + sw*n/kw);
      captw = captw + dhvap*rhow*((alpha - n)/ks*sw*pw - (alpha - n)/ks*sw*pg + n)*dsw_dpc;
    }
    else{
      //incompressible grains:
      captw = dhvap*rhow*n*dsw_dpc;
    }
  }
  return(captw);
}


/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kgt - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_kgt(double pw,double pg,double t,long ipp)
{
  double kgt,rhog,mg,dg,dpgw_dt;
  if(model_type == artificial3)
    kgt = kgt0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhog =get_rhog(pw,pg,t);
    mg = get_mg(pw,pg,t);
    dg = get_dg(pw,pg,t,ipp);
    dpgw_dt = get_dpgw_dt(pw,pg,t);
    
    
    dpgw_dt = get_dpgw_dt(pw,pg,t);
    
    kgt = -rhog*ma*mw/mg/mg*dg/pg*dpgw_dt;
  }
  return(kgt);
}


/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktg - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktg(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double ktg;

  ktg = 0.0;

  return(ktg);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captg - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_captg(double pw,double pg,double t,long ipp)
{
  double captg,dhvap,rhow,alpha,n,ks,sw,sg,kw,dsw_dpc;
 
  if(model_type == artificial3)
    captg = captg0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    sg = 1.0 - sw;
    kw = get_kw(pw,pg,t,ipp);
    dsw_dpc = get_ds_dpc(pw,pg,t,ipp);//this is equal to minus derivative with respect to suction
    
    if(rel_gas_press == 1)
      pg = pg + p0;

    if(compress == 1){
      //compressible grains:
      captg = -dhvap*rhow*(alpha - n)/ks*sg*sw;//water
      captg = captg - dhvap*rhow*((alpha - n)/ks*sw*pw - (alpha - n)/ks*sw*pg + n)*dsw_dpc;
    }
    else{
      //incompressible grains:
      captg = -dhvap*rhow*n*dsw_dpc;
    }
  }
  return(captg);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captt - capacity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_captt(double pw,double pg,double t,long ipp)
{
  double captt,rhocp,dhvap,betasw,rhow,alpha,n,ks,sw,dsw_dt;
   
  if(model_type == artificial3)
    captt = captt0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhocp = get_rhocp(pw,pg,t,ipp);
    dhvap = get_dhvap(t);
    betasw = get_betasw(pw,pg,t,ipp);
    rhow = get_rhow(t);
    alpha = get_alpha(pw,pg,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,pg,t,ipp);
    sw = get_sw(pw,pg,t,ipp);
    dsw_dt = get_ds_dt(pw,pg,t,ipp);
    
    if(rel_gas_press == 1)
      pg = pg + p0;
    
    if(compress == 1){
      //compressible grains:
      captt = rhocp;
      captt = captt + dhvap*betasw*rhow;
      captt = captt - dhvap*(rhow*((alpha - n)/ks*sw*pw - (alpha - n)/ks*sw*pg + n)*dsw_dt);
    }
    else{
      //incompressible grains:
      captt = rhocp;
      captt = captt + dhvap*betasw*rhow;
      captt = captt - dhvap*n*dsw_dt;
    }
  }
  return(captt);
}



/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval ktw - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktw(double pw,double pg,double t,long ipp)
{
  double dhvap,rhow,krw,muw,ktw,kintr;
  
  if(model_type == artificial3)
    ktw = ktw0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    krw = get_krw(pw,pg,t,ipp);
    muw = get_muw(t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrw(pw,pg,t,ipp);


    //correction 02/12/2022
    ktw = -dhvap*(rhow*kintr*krw/muw);//water
  }

  return(ktw);
}




/**
   function creates conductivity coefficient of a general material - first part
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt1 - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktt1(double pw,double pg,double t,long ipp)
{
  double lambdaeff,ktt1;
   
  if(model_type == artificial3)
    ktt1 = ktt0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    lambdaeff = get_lambdaeff(pw,pg,t,ipp);
    
    ktt1 = lambdaeff;
  }

  return(ktt1);
}


/**
   function creates conductivity coefficient of a general material - second (A) part (convective term)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2a - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktt2a(double pw,double pg,double t,long ipp)
{
  double n,sw,cpw,rhow,krw,muw,ktt2a,kintr;
    
  if(model_type == artificial3)
    ktt2a = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    n = get_porosity(ipp);
    sw = get_sw(pw,pg,t,ipp);
    rhow = get_rhow(t);
    cpw = get_cpw();
    krw = get_krw(pw,pg,t,ipp);
    muw = get_muw(t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrw(pw,pg,t,ipp);

    ktt2a = n*sw*rhow*cpw*kintr*krw/muw;
  }
  return(ktt2a);
}

/**
   function creates conductivity coefficient of a general material - second (C) part (convective term)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2b - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktt2b(double pw,double pg,double t,long ipp)
{
  double n,sw,sg,rhog,cpg,krg,mug,kintr;
  double ktt2b;
  
  if(model_type == artificial3)
    ktt2b = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    n = get_porosity(ipp);
    sw = get_sw(pw,pg,t,ipp);
    sg = 1.0 - sw;
    rhog = get_rhog(pw,pg,t);
    cpg = get_cpg(pw,pg,t);
    krg = get_krg(pw,pg,t,ipp);
    mug = get_mug(pw,pg,t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrg(pw,pg,t,ipp);

    ktt2b = n*sg*rhog*cpg*kintr*krg/mug;
  }

  return(ktt2b);
}


/**
   function creates conductivity coefficient of a general material - second (B) part (convective term)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2c - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktt2c(double pw,double pg,double t,long ipp)
{
  double rhow,ktt2c;
  
  if(model_type == artificial3)
    ktt2c = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhow = get_rhow(t);
    
    ktt2c = rhow;
  }

  return(ktt2c);
}


/**
   function creates conductivity coefficient of a general material - second (D) part (convective term)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2d - conductivity coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_ktt2d(double pw,double pg,double t,long ipp)
{
  double rhog,ktt2d;
    
  if(model_type == artificial3)
    ktt2d = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhog = get_rhog(pw,pg,t);
    
    ktt2d = rhog;
  }

  return(ktt2d);
}



/**
   function creates right-hand side coefficient of a general material for c medium
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   03/10/2023 TKr
*/
double con_hawf3mat::get_fw1(double pw,double pg,double t,long ipp)
{
  double fc1;
  double rhow,rhogw,rhog,krw,muw,krg,mug,kintr,kintrg,kintrw;
  
  if(model_type == artificial3)
    fc1 = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhow = get_rhow(t);
    rhogw = get_rhogw(pw,pg,t);
    rhog = get_rhog(pw,pg,t);  
    krw = get_krw(pw,pg,t,ipp);
    muw = get_muw(t);
    krg = get_krg(pw,pg,t,ipp);
    mug = get_mug(pw,pg,t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3){
      kintr = get_kintr(pw,pg,t,ipp);
      fc1 = rhogw*kintr*krg*rhog/mug + rhow*kintr*krw*rhow/muw;
    }
    else{
      kintrw = get_kintrw(pw,pg,t,ipp);
      kintrg = get_kintrg(pw,pg,t,ipp);
      fc1 = rhogw*kintrg*krg*rhog/mug + rhow*kintrw*krw*rhow/muw;
    }
  }
  return(fc1);
}


/**
   function creates right-hand side coefficient of a general material for g medium
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   03/10/2023 TKr
*/
double con_hawf3mat::get_fg(double pw,double pg,double t,long ipp)
{
  double fg;
  double rhoga,rhog,krg,mug,kintr;
  
  if(model_type == artificial3)
    fg = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    rhoga = get_rhoga(pw,pg,t);
    rhog = get_rhog(pw,pg,t);
    krg = get_krg(pw,pg,t,ipp);
    mug = get_mug(pw,pg,t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrg(pw,pg,t,ipp);

    fg = rhoga*kintr*krg*rhog/mug;
    fg = 0.0;//no gravity force is included for the moist air 
  }

  return(fg);
}


/**
   function creates right-hand side coefficient of a general material for t medium
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   03/10/2023 TKr
*/
double con_hawf3mat::get_ft1(double pw,double pg,double t,long ipp)
{
  double ft1;
  double dhvap,rhow,krw,muw,kintr;

  if(model_type == artificial3)
    ft1 = 0.0;
  else{
    //gas pressure check
    gaspress_check(pw,pg,t,ipp);
    //water pressure check
    waterpress_check(pw,pg,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    krw = get_krw(pw,pg,t,ipp);
    muw = get_muw(t);
    //kintr = get_kintr(pw,pg,t,ipp);
    if(kintr_type < 3)
      kintr = get_kintr(pw,pg,t,ipp);
    else
      kintr = get_kintrw(pw,pg,t,ipp);
    
    ft1 = -dhvap*rhow*kintr*krw*rhow/muw;
  }

  return(ft1);
}





/**function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval new_nodval - nodal value of transmission coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_transcoeff_ww(double pw,double pg,double t,long bc,long ipp)
{
  double trc,pgws,rhow;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{

    //relative humidity - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = get_pgws(t);
    rhow = get_rhow(t);
    
    trc=1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 31:{
    //mass concentration - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = get_pgws(t);
    rhow = get_rhow(t);
    
    trc=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 32:{
    //water vapour pressure - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = get_pgws(t);
    rhow = get_rhow(t);
    
    trc=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 33:{//relative humidity - flux
    trc=0.0;
    break;
  }
  case 34:{//mass concentration
    trc=0.0;
    break;
  } 
  case 35:{//water vapour pressure
    trc=0.0;
    break;
  } 
  case 36:{//relative humidity - pokus
    trc=0.0;
    break;
  }
  case 90:{//fire = heat transfer + radiation
    trc = 1.0;
    break;
  }
  case 40:{//simulation of free soil surface
    if(pw < pw_bc)
      trc = 0.0;
    else
      trc = 1.0e+03;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(trc);
}



/**function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval new_nodval - nodal value of transmission coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_transcoeff_ww(double pw,double pg,double t,long bc,long ipp,int flag)
{
  double trc,pgws,rhow;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{

    //relative humidity - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    //pgws = get_pgws(t);
    //rhow = get_rhow(t);

    trc = -1.0;//get_pwrh(bv,t);//tady??!!
    //trc=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie

    break;
  }
  case 31:{
    //mass concentration - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = get_pgws(t);
    rhow = get_rhow(t);
    
    trc=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 32:{
    //water vapour pressure - for testing = Kelvin eq. approximated by Taylor's serie in point 0
    pgws = get_pgws(t);
    rhow = get_rhow(t);
    
    trc=-1.0*pgws*mw*mw/gasr/gasr/t/t/rhow;//1st member of T. serie
    
    break;
  }
  case 33:{//relative humidity - flux
    if(flag == 1)
      trc=0.0;//into right hand side (matrix)
    else{
      pgws = get_pgws(t);
      trc=1.0;//into left hand side (flux)
      //trc=gasr*t/mw*pgws;//pokus
    }
    break;
  }
  case 34:{//mass concentration
    trc=0.0;
    break;
  } 
  case 35:{//water vapour pressure
    trc=0.0;
    break;
  } 
  case 36:{//relative humidity - pokus
    trc=0.0;
    break;
  }
  case 90:{//fire = heat transfer + radiation
    trc = 1.0;
    break;
  }
  case 40:{//simulation of free soil surface
    if(pw < pw_bc)
      trc = 0.0;
    else
      trc = 1.0e+03;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(trc);
}

/**function creates correct new nodal value on the boundary (transmission) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval new_nodval - nodal value of transmission coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_nodval_ww(double bv,double pw,double pg,double t,long bc,long ipp)
{
  double new_nodval;
  state_eq tt;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{
    new_nodval = tt.get_pcrh(bv,t);
    new_nodval = pg - new_nodval;
    break;
  }
  case 40:{//simulation of free soil surface
    new_nodval = bv;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(new_nodval);
}



/**function creates flux on the boundary (transmission - convective mass transfer) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval flux - bounday flux

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_flux_ww(double bv,double pw,double pg,double t,long bc,long ipp)
{
  double flux,trc,pc;
  state_eq tt;
  
  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//for testing - boundary flux
    pc = pg - pw;

    flux = tt.get_pcrh(bv,t);
    flux = -1.0*(flux - pc);

    break;
  }
  case 40:{//simulation of free soil surface
    if(pw < pw_bc)
      trc = 0.0;
    else
      trc = 1.0e+03;
    flux = trc*(bv - pw);//minus sign
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function creates correct transfer coefficient on the boundary (transmission) for t medium

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval new_nodval - nodal value of transmission coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_transcoeff_tt(double pw,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    ft3=1.0;
    break;
  }
  case 100:{//combined condition for flux and dirichlet b.c.
    ft3=1.0;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  
  return(ft3);
}




/**
   function creates correct new nodal value on the boundary (transmission) for t medium

   @param bv - value of prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval new_nodval - nodal value of transmission coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_nodval_tt(double bv,double /*trr*/,double pw,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    ft3 = bv;
    break;
  }
  case 100:{//combined condition for flux and dirichlet b.c.
    ft3 = bv;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }

  return(ft3);
}


/**
   function creates flux on the boundary (transmission - convective mass transfer) for c medium

   @param bv - prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   @retval new_nodval - nodal value of transmission coefficient

   03/10/2023 TKr
*/
double con_hawf3mat::get_transmission_flux_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp)
{
  double ft3;

  //gas pressure check
  gaspress_check(pw,pg,t,ipp);
  //water pressure check
  waterpress_check(pw,pg,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission - boundary flux
    ft3 = (bv - t);
    break;
  }
  case 100:{//combined condition for flux and dirichlet b.c.
    if(trr > 1.0){
      ft3 = (bv - t);//heat transfer b.c. for dirichlet b.c.
    }
    else
      ft3 = bv; //heat transfer b.c. for prescribed flux
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }

  return(ft3);
}



/**
   function computes all variables in nodes
   
   @param compother - number of other components
   @param ipp - first integration point on element
   @param pw - pore water pressure on actual node
   @param pg - pore gas pressure on actual node
   @param t - temperature on actual node

   @retval other - other variables

   03/10/2023 TKr
*/

double con_hawf3mat::get_othervalue(long compother,long ipp,double *r)
{
  double other;

  switch (compother){
  case 0:{//capillary pressure
    other = r[1] - r[0];
      break;
  }
  case 1:{//gas pressure
    other = r[1];
    break;
  }
  case 2:{//temperature in deg. C
    other = r[2] - 273.15;
    break;
  }
  case 3:{//saturation
    other = get_sw(r[0],r[1],r[2],ipp);
    break;
  }
  case 4:{//vapour pressure
    other = get_pgw(r[0],r[1],r[2]);
    break;
  }
  case 5:{//liquid water pressure
    other = r[0];
    break;
  }
    //case 6:{//moisture content
    //other = get_w(r[0],r[1],r[2],ipp);
    //break;
    //}    
  case 6:{//relative humidity
    other = get_pgw(r[0],r[1],r[2])/get_pgws(r[2]);
    break;
  }    
  case 7:{//dsw_dpw
    other = -get_ds_dpc(r[0],r[1],r[2],ipp);
    break;
  }    
  case 8:{//dsw_dt
    other = get_ds_dt(r[0],r[1],r[2],ipp);
    break;
  }
  default:{
    print_err ("unknown type of component is required in function", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  return (other);

}

/**
     function prints names of all variables in nodes

     @param out - output file
     @param compother - number of other components
*/
void con_hawf3mat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//capillary pressure
    fprintf (out,"Capilary pressure (Pa)        ");
    break;
  }
  case 1:{//gas pressure
    fprintf (out,"Gas pressure (Pa)             ");
    break;
  }
  case 2:{//temperature
    fprintf (out,"Temperature (C)               ");
    break;
  }
  case 3:{//saturation
    fprintf (out,"Degree of saturation ()       ");
    break;
  }
  case 4:{//vapour pressure
    fprintf (out,"Pore water vapor pressure (Pa)");
    break;
  }
  case 5:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
    //case 6:{//moisture content
    //fprintf (out,"Moisture content (kg/kg)      ");
    //break;
    //}    
  case 6:{//relative humidity
    fprintf (out,"Relative humidity ()          ");
    break;
  }    
  case 7:{//dS_dpw
    fprintf (out," dS_dpw Derivative of Saturation degree with respect to pore water pressure      ");
    break;
  }    
  case 8:{//dS_dt
    fprintf (out," dS_dt Derivative of Saturation degree with respect to temperature      ");
    break;
  }    
  default:{
    print_err ("unknown type of component is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.

   23/05/2016, TKr
*/
void con_hawf3mat::updateval (long /*ipp*/)
{
}



/**
   This function initializes material model data

   @param ipp - integration point number

   23/05/2016, TKr
*/
void con_hawf3mat::initval(long /*ipp*/)
{
}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Krejci according to Tomas Koudelka, 23/05/2016
*/
void con_hawf3mat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_water_press;
  dofname[1] = trf_gas_press;
  dofname[2] = trf_temperature;
}


/**
   function returns temperature

   @param ipp - integration point number

   @retval t - temperature

   16/07/2018, TKr
*/
double con_hawf3mat::give_temperature(long ipp)
{
  double t;

  t = Tm->ip[ipp].av[2];

  return(t);
}


/**
   function returns effective pore pressure (pressure has negative sign - mechanical convention)

   @param ipp - integration point number

   @retval pw - pore water pressure

   03/10/2023 TKr
*/
double con_hawf3mat::give_effective_pore_pressure(long ipp)
{
  double ps,pw,pg,t,xi=1.0;

  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];


  //needs to be corrected??!!
  switch (model_type){
  case lewis_and_schrefler3:{//Lewis and Schrefler's book
    xi =  get_xi(pw,pg,t,ipp);//effective stress factor

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    ps = xi*pw + (1.0 - xi)*pg;
    ps = -ps;//returns effective pore pressure
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;

    xi = get_xi(pw,pg,t,ipp); //effective stress factor

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    ps = xi*pw + (1.0 - xi)*pg;
    ps = -ps/mefel_units;//returns effective pore pressure //corrected units for mefel //basic units = Pa  
    
    //debug???!!!
    //pw = -pw/mefel_units;//returns effective pore pressure //corrected units for mefel //basic units = Pa  
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }   
  return(ps);
}


/**
   function returns water pressure

   @param ipp - integration point number

   @retval pw - pore water pressure

   27/05/2016, TKr
*/
double con_hawf3mat::give_water_pressure(long ipp)
{
  double pw;

  pw = Tm->ip[ipp].av[0];

  return(pw);
}

/**
   function returns pore pressure

   @param ipp - integration point number

   @retval pp - pore pressure

   27/05/2016, TKr
*/
double con_hawf3mat::give_pore_pressure(long ipp)
{
  double pp,pg,pw;

  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];

  pp = pg - pw;

  return(pp);
}

/**
   function returns gas pressure

   @param ipp - integration point number

   @retval pg - pore gas pressure

   27/05/2016, TKr
*/
double con_hawf3mat::give_gas_pressure(long ipp)
{
  double pg;

  pg = Tm->ip[ipp].av[1];

  return(pg);
}


/**
   function computes suction stress s = -(pg - pw) = -pc;
   @param ipp - integration point number

   @retval suction - suction stress [Pa]

   27/05/2016, TKr
*/
double con_hawf3mat::give_suction(long ipp)
{
  double pw,pg,suction;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if(rel_gas_press == 1)
    pg = pg + p0;

  switch (model_type){
  case lewis_and_schrefler3:{//Lewis and Schrefler's book
    suction = -1.0*(pg - pw);//this is correct, because pore water pressure is negative
    break;
  }
  case lewis_and_schrefler3_2:{//Lewis and Schrefler's book p. 381
    suction = -1.0*(pg - pw);//this is correct, because pore water pressure is negative
    break;
  }
  case van_genuchten3:{//partially saturated medium = Van Genuchten model
    suction = -1.0*(pg - pw);//this is correct, because pore water pressure is negative
    suction = suction/mefel_units;//corrected units for mefel //basic units = Pa
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel
    pg = pg - p0;//pore gas pressure without atmospheric pressure

    suction = -1.0*(pg - pw);//this is correct, because capillary water pressure is negative
    //suction = pw; //debug??!!
    suction = suction/mefel_units;//corrected units for mefel //basic units = Pa
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(suction);
}


/**
   function gives saturation degree
   @param ipp - integration point number

   @retval saturation degree []

   27/05/2016, TKr
*/
double con_hawf3mat::give_saturation_degree(long ipp)
{
  double pw,pg,t,s;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];
  

  s = get_sw(pw,pg,t,ipp);
  
  return(s);
}

/**
   function computes water content w [kg/kg]
   @param pw - pore water pressure
   @param pg - gas pore pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval w - computes water content w [kg/kg]
*/
double con_hawf3mat::get_w(double pw,double pg,double t,long ipp)
{
  double w,n,s,rhow;

  n = get_porosity(ipp);
  s = get_sw(pw,pg,t,ipp);
  rhow = get_rhow(t);

  w = (n*s*rhow)/(1.0 - n)/rhos0;

  return(w);
}




/**
  The funtion marks required non-transport quantities in the array antq.

  @param antq - array with flags for used material types
                antq[i] = 1 => quantity type nontransquant(i+1) is required
                antq[i] = 0 => quantity type nontransquant(i+1) is not required

  @return The function does not return anything, but it may change content of antq array.
  
  27/05/2016, TKr
*/
void con_hawf3mat::give_reqntq(long *antq)
{
  switch (model_type){
  case artificial3:
  case lewis_and_schrefler3:{//Lewis and Schrefler's book
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    antq[porosity-1] = 1;
    if (sr_type == mefel_sr){
      antq[saturation_deg-1] = 1;
      antq[der_saturation_deg-1] = 1;
      antq[der_saturation_deg_depsv-1] = 1;
      antq[der_saturation_deg_dtemp-1] = 1;
    }
    if(vol_strain_effect == 1)
      antq[strain_vol_rate-1] = 1;
    break;
  }  
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
  }
  } 
}
