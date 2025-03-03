/*
  File:             consol_hwf2.cpp
  Author:           Tomas Krejci, 06/04/2018
  Purpose:          computes conductivity and capacity matrices in a material point for consol_hwf2 porous media;
                    material model for saturated-nonsaturated water flow and heat transfer 
                    in a deforming porous medium (soils). It combines water flow (Richards equation) and water vapour diffusion and heat conduction

  unknowns:         number of unknowns=2, pw = pore water pressure, t = temperature
  sources:          

  1. THE FINITE ELEMENT METHOD IN THE STATIC AND DYNAMIC DEFORMATION AND CONSOLIDATION OF POROUS MEDIA
  R. W. Lewis, B.A. Schrefler
  
  2. NONLINEAR MODELLING OF CONCRETE AS POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
  Francesco Pesavento - doctoral thesis                      
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "consol_hwf2.h"
#include "globalt.h"
#include "C60bazant.h" //testing
#include "C60baroghel.h" //testing
#include "globmatt.h"

con_hwf2mat::con_hwf2mat()
{
  compress = 0; //compressible grains: 0=no; 1=yes
  kintr_type = 0;        //intrinsic permability calculation type
  krw_type = 0;          //relative permeability calculation type
  deff_type = 0;         // diffusion calculation type
  sr_type = 1;           //retention curve calculation type
  xi_type = 1;           //effective stress parameter type
  lambda_type =0;        //heat conduction calculation type
  cps_type = 0;          //specific heat calculation type
  betas_type = 0;        //thermal expansion calculation type

  vol_strain_effect = 0; //volumetric strain rate influence: 0=no; 1=yes 
  mefel_units = 1.0;//basic units for pressures = Pa (Pascals)

  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha0 = 0.0;
  ks0 = 0.0;
  kt0 = 0.0;
  phi0 = 0.0;
  kintr0 = 0.0;
  betas0 = 0.0;
  rhos0 = 0.0;
  cps0 = 0.0;
  lambda_eff0 = 0.0;
  lambda_dry = 0.0;
  lambda_wet = 0.0;
  sr_dry = 0.0;
  sr_wet = 0.0;

  sirr = 0.0;
  ssat = 1.0;
  lambda_krw = 1.9;  //parameter for exponential function of relative permeability
  beta_krw = 0.0;    //parameter for power function of relative permeability

  bb1 = 0.0;
  phi01 = 0.0;
  

  // PHYSICAL PROPERTIES OF WATER
  kw0 = 2.0e9;//bulk modulus of water
  rhow0 = 1000.0;
  muw0 = 1.0e-3;
  cpw0 = 4181.0;
  hvap0 =  2.7e+5;

  pw_bc = 0.0; //free boundary pressure

  tcr = 647.3;
  mw = 18.01528;//molar mass of water kg.mol-1
  gasr = 8314.41;//universal gas constant J.kmol-1.K-1

  c8 = -5.8002206e+03;
  c9 = 1.3914993;
  c10 =-4.8640239e-02;
  c11 = 4.1764768e-05;
  c12 = -1.4452093e-08;
  c13 = 6.5459673;

  dv0 = 2.16e-5;//initial water vapour diffusivity in air
  tau0 = 0.8;

  cpgw0 = 1805.0;

  //parameters for the artificial material:
  kww0 = kwt0 = ktw0 = ktt0 = 0.0;
  capww0 = capwt0 = captw0 = captt0 =0.0;

  gamma = 0.0;
  lambda0 = 0.0;
  s_entry = 0.0;

  cps_dry = cps_wet = cps0 = cps_lin = 0.0;
}

con_hwf2mat::~con_hwf2mat()
{}


/**
   function reads parameters
   
   @param in - input file

   20/05/2017, TKr
*/
void con_hwf2mat::read(XFILE *in)
{
  xfscanf (in,"%k%m","heatwaterflowtype",&heatwaterflowtype_kwdset, &model_type);
  xfscanf (in,"%d", &compress);

  // common material parameters
  xfscanf (in,"%le %le %le %le %le %le %d %d %d %d %d %d %d %d %d", &alpha0, &ks0, &rhos0, &pw_bc, &tau0, &kintr0, &por_type, &kintr_type, &krw_type, &deff_type, &sr_type, &xi_type, &lambda_type, &cps_type, &betas_type);

  switch (model_type){
  case artificial2:{//artificial isotropic material for non-isotherma air-water flow
    xfscanf (in,"%le %le %le %le", &kww0, &kwt0, &ktw0, &ktt0);
    xfscanf (in,"%le %le %le %le", &capww0, &capwt0, &captw0, &captt0);
  }
  case lewis_and_schrefler2hw:{//Lewis and Schrefler's model
    break;
  }
  case lewis_and_schrefler2hw_mefel:{//Lewis and Schrefler's model approach coupled with MEFEL
    xfscanf (in,"%le %d", &mefel_units, &vol_strain_effect);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 


  if(model_type != artificial2){
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
    case 0:{//constant
      break;
    }
    case 1:{//dependent on porosity
      xfscanf (in,"%le %le", &bb1, &phi01);
      break;
    }
    case 2:{//dependent on porosity - cubic and quadratic
      xfscanf (in,"%le", &phi01);
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
      xfscanf (in,"%le %le", &cps_dry, &cps_wet);
      break;
    }
    case 2:{//dependent on temperature
      xfscanf (in,"%le %le", &cps0, &cps_lin);
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
void con_hwf2mat::print(FILE *out)
{
  /*  fprintf (out,"\n %d ", int(model_type));
      fprintf (out,"\n %d ", compress);
      
      switch (model_type){
      case lewis_and_schrefler2hw:{//Lewis and Schrefler's model
      fprintf (out,"\n %le %le %le %le %le %le %le %le ",alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda0);
      lewis_ret.print(out);
      break;
      }
      case lewis_and_schrefler2hw_mefel:{//Lewis and Schrefler's book
      fprintf (out,"\n %le %le %le %le %le %le %le %le %le  %le %le %le %le  %d \n",mefel_units,alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda1,lambda2,sr1,sr2,pw_bc, vol_strain_effect);
      break;
      }
      case van_genuchten2hw:{//partially saturated medium =Van Genuchten model
      fprintf (out,"\n %le %le %le %le %le %le %le %le %le ",mefel_units,alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda0);
      van_genuchten_ret.print(out);
      break;
      }
      default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
      }
      } 
  */
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
*/
double con_hwf2mat::get_betas(long /*ipp*/)
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
   @param t - temperature

   @retval cps - specific heat of soil skeleton
*/
double con_hwf2mat::get_cps(double pw, double t, long ipp)
{
  double cps;
  double sw;

  cps = cps0;

  switch (cps_type){
  case 0:{ //constant
    cps = cps0;
    break;
  }
  case 1:
  case 3:{ //moisture dependent
    sw = get_sw(pw,t,ipp);
    cps = cps_dry + (cps_wet - cps_dry)*sw; //this is maybe linear approximation for bentonite:
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
*/
double con_hwf2mat::get_rhos(double /*t*/)
{
  double rhos;
   
  rhos = rhos0;
  
  return(rhos);
}



/**
   function computes effective stress factor xi
   @param pw - water pressure
   @param t - temperature

   @retval xi - factor xi

   03/04/2023, TKr
*/
double con_hwf2mat::get_xi(double pw, double t, long ipp)
{
  double suc=0.0,sr=0.0,xi=0.0;
  
  switch (xi_type){
  case biot_xi:
    xi = get_sw(pw,t,ipp);
    break;
  case biot_reduced_xi:
    //xi = gamma*get_sw(pw,ipp);
    sr = get_sw(pw,t,ipp);
    //xi = pow(sr,(gamma/lambda0));
    xi = pow(sr,gamma);
    xi = (1-gamma)*xi;
    break;
  case biot_masin_xi:{//according to masin for testing
    sr = get_sw(pw,t,ipp);
    suc = -pw;//pore gas pressure without atmospheric pressure
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
    double e=0.0,dpw=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }
    
    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpc = -dpw;
    
    pc = -pw;
    
    xi = masin_ret.psi(pc,dpc,e,t);//positive value of suction
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  //old:
  /* switch (xi_type){
     case biot_xi:
     xi = get_sw(pw,t,ipp);
     break;
     case masin_xi:{
     double e=0.0,dpw=0.0,por=0.0;
     if(Tm->nontransq != NULL){
     por = Tm->givenontransq(porosity, ipp);// from mefel
     e = por/(1-por);
     }
     else{
     e = phi0/(1-phi0);
     }
     
     dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
     xi = masin_ret.psi(-pw,-dpw,e,t);//positive value of suction
     break;
     }
     default:{
     print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
     abort();
     }
     } 
  */

  return(xi);
}



/**
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param t - temperature

   @retval sw - degree of saturation

   20/05/2017, TKr
*/
double con_hwf2mat::get_sw(double pw, double t, long ipp)
{
  double sw;
  sw = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    sw = bazant_ret.sat(-pw,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    sw = lewis_ret.sw(pw);
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
    
  case gardner_exponential_sr:{//partially saturated medium = Exponential model, Gardner 
    sw = gardner_ret.sw(pw);
    break;
  }
  case potts_log_linear_sr:{//partially saturated medium = Log linear model, Potts
    sw = potts_ret.sw(pw);
    break;
  }
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    sw = van_genuchten_ret.sw(pw,t);
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
    //sw = van_genuchten_ret.sw2(pw);
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
    sw = masin_ret.sw(-pw,-dpw,e,t);//positive value of suction
    break;
  }
    case febex_granit_sr:{//FEBEX granit
    sw = febex_granit_ret.sw(-pw);
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
   function computes specific moisture content (partial derivative of degree of saturation with respect to pw)
   @param pw - water pressure
   @param t - temperature

   @retval ds_dpw - specific moisture content = partial derivative of degree of saturation with respect to pw


   20/05/2017, TKr
*/
double con_hwf2mat::get_dsw_dpw(double pw, double t, long ipp)
{
  double dsw_dpw;
  dsw_dpw = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    dsw_dpw = -bazant_ret.dsat_dpc(-pw,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    dsw_dpw = lewis_ret.dsw_dpw(pw);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    dsw_dpw = Tm->givenontransq(der_saturation_deg, ipp); //actual derivative of saturation degree
    dsw_dpw = dsw_dpw/mefel_units; //basic units = Pa
    break;
  }

  case table_sr:{//saturation degree and its derivative are obtained from table;
    dsw_dpw = data.getderiv (pw); //actual derivative of saturation degree
    break;
  }

  case gardner_exponential_sr:{//partially saturated medium = Exponential model, Gardner 
    dsw_dpw = gardner_ret.dsw_dpw(pw);
    break;
  }
  case potts_log_linear_sr:{//partially saturated medium = Log linear model, Potts
    dsw_dpw = potts_ret.dsw_dpw(pw);
    break;
  }
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    dsw_dpw = van_genuchten_ret.dsw_dpw(pw,t);
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
    //dsw_dpw = van_genuchten_ret.dsw_dpw2(pw);
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
    dsw_dpw = -masin_ret.dsw_dpw(-pw,-dpw,e,t);//positive value of suction
    //dsw_dpw = dsw_dpw/mefel_units; //basic units = Pa
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    dsw_dpw = -febex_granit_ret.dsw_ds(-pw);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(dsw_dpw);
}


/**
   function computes partial derivative of degree of saturation with respect to t
   @param pw - water pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to temperature


   20/05/2017, TKr
*/
double con_hwf2mat::get_dsw_dt(double pw, double t, long ipp)
{
  double dsw_dt;
  dsw_dt = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    dsw_dt = bazant_ret.dsat_dt(-pw,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    dsw_dt = lewis_ret.dsw_dt(pw);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    dsw_dt = Tm->givenontransq(der_saturation_deg_dtemp, ipp); //actual derivative of saturation degree
    dsw_dt = dsw_dt/mefel_units; //basic units = Pa
    break;
  }

  case table_sr:{//saturation degree and its derivative are obtained from table;
    //dsw_dt = data.getderiv (pw); //actual derivative of saturation degree
    break;
  }

  case gardner_exponential_sr:{//partially saturated medium = Exponential model, Gardner 
    dsw_dt = gardner_ret.dsw_dt(pw);
    break;
  }
  case potts_log_linear_sr:{//partially saturated medium = Log linear model, Potts
    dsw_dt = potts_ret.dsw_dt(pw);
    break;
  }
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    dsw_dt = van_genuchten_ret.dsw_dt(pw,t);
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
    //dsw_dt = van_genuchten_ret.dsw_dt2(pw);
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
    dsw_dt = masin_ret.dsw_dt(-pw,-dpw,e,t);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    dsw_dt = febex_granit_ret.dsw_dt(-pw);
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
   function computes intrinsic permeability
   @param pw - water pressure
   @param t - temperature
   
   @retval kintr - intrinsic permeability

   27/07/2017, TKr
*/
double con_hwf2mat::get_kintr(double pw, double t, long ipp)
{
  double kintr,phi;

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
   function computes water relative permeability
   @param pw - water pressure
   @param t - temperature
   
   @retval krw - water relative permeability

   20/05/2017, TKr
*/
double con_hwf2mat::get_krw(double pw, double t, long ipp)
{
  double sw,krw,sef=0.0;
  
  switch (krw_type){
  case 0:{//constant
    krw = 1.0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,t,ipp);
    krw = (sw-sirr)/(ssat-sirr);
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,t,ipp);
    krw = (sw-sirr)*(sw-sirr)*(sw-sirr)/((ssat-sirr)*(ssat-sirr)*(ssat-sirr));
    break;
  }
  case 3:{//exponential
    sw = get_sw(pw,t,ipp);
    krw = pow(sw,lambda_krw);
    break;
  }
  case 4:{//liakopoulos
    sw = get_sw(pw,t,ipp);
    krw = 1.0-2.207*pow((1.0-sw),1.0121);
    break;
  }
  case 5:{//double exponential
    sw = get_sw(pw,t,ipp);
    sef = (sw-sirr)/(ssat-sirr); //effective saturation degree
    krw = pow(sef,(1.0/beta_krw));
    krw = pow((1.0-krw),beta_krw);
    krw = pow((1.0-krw),2.0);
    krw = pow(sef,0.5)*krw;
    break;
  }
  case 6:{//FEBEX granit
    sw = get_sw(pw,t,ipp);
    krw = febex_granit_ret.get_krw(sw);
    break;
  }
  case 7:{
    krw = van_genuchten_ret.get_krw(pw,t);
  }
  case 9:{//bazant
    sw = get_sw(pw,t,ipp);
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
   function returns porosity

   @retval phi - porosity

   03/04/2023, TKr
*/
double con_hwf2mat::get_porosity(long ipp)
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
   function computes Biot's constant
   @param pw - pore water pressure
   @param t - temperature
   
   @retval alpha - Biot's constant
*/
double con_hwf2mat::get_alpha(double /*pw*/, double /*t*/,long /*ipp*/)
{
  double alpha;
  //double kt,ks;
  
  //kt = get_kt(pw,t,ipp);
  //ks = get_ks(pw,t,ipp);
  
  //alpha = 1.0 - kt/ks;  

  alpha = alpha0;//provisionally
  
  return(alpha);
}

/**
   function computes bulk modulus of solid phase
   @param pw - pore water pressure
   @param t - temperature

   @retval ks - bulk modulus of solid phase
*/
double con_hwf2mat::get_ks(double /*pw*/,double /*t*/,long /*ipp*/)
{
  double ks;

  ks = ks0;

  return(ks);
}


/**
   function computes water density
   @param t - temperature

   @retval rhow - water density
*/
double con_hwf2mat::get_rhow(double /*t*/)
{
  double rhow;

  rhow = rhow0;//provisionally
  
  return(rhow);
}


/**
   function computes compresibility coefficient of water
   
   @retval kw - compresibility coefficient of water

*/
double con_hwf2mat::get_kw(double /*pw*/,double /*t*/,long /*ipp*/)
{
  double kw;
  
  kw = kw0;

  return(kw);
}


/**
   function volume computes thermal expansion coefficient of solid - water for compressible solid grains
   @param pw - pore water pressure
   @param t - temperature

   @retval betasw - volume thermal expansion coefficient of solid - water for compressible solid grains
*/
double con_hwf2mat::get_betasw(double pw,double t,long ipp)
{
  double alpha,betasw,betas,n,sw,betaw;

  alpha = get_alpha(pw,t,ipp);
  betas = get_betas(ipp);
  n = get_porosity(ipp);
  sw = get_sw(pw,t,ipp);
  betaw = get_betaw(t);

  betasw = sw*(betas*(alpha - n) + n*betaw);

  return(betasw);
}


/**
   function computes dynamic viscosity of water = 1000e-6 Pa*s at 20 C
   @param t - temperature

   @retval muw dynamic viscosity of water = 1000e-6 Pa*s at 20 deg. C
*/
double con_hwf2mat::get_muw(double /*t*/)
{
  double muw;

  muw = muw0;

  return(muw);
}

/**
   function computes water specific heat

   @retval cpw - water specific heat
*/
double con_hwf2mat::get_cpw()
{
  double c;

  c = cpw0;

  return(c);
}

/**
   function computes volume thermal expansion coefficient of water
   @param t - temperature
   
   @retval betaw - volume thermal expansion coefficient of water
*/
double con_hwf2mat::get_betaw(double /*t*/)
{
  double betaw;
  
  betaw = 0.68e-4;//[K-1] at t=273.15
  //betaw = 10.1e-4;//[K-1] at t=420.0
  
  //temporarily??!!
  betaw = 0.63e-5;
  //betaw = 0.68 + (10.1 - 0.68)/(420.0 - 273.15)*(t - 273.15);//linear expression (rough)
  
  return(betaw);
}


/**
   function computes effective thermal capacity of partially saturated medium
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval rhocp - effective thermal capacity of partially saturated medium
*/
double con_hwf2mat::get_rhocp(double pw,double t,long ipp)
{
  double rhocp;
  double sw,n;

  sw = get_sw(pw,t,ipp);
  n = get_porosity(ipp);
  
  rhocp = (1-n)*rhos0*cps0;// solid phase
  rhocp = rhocp + n*sw*rhow0*cpw0;//liquid phases (constant values for rhow and cpw)
  
  switch (model_type){
  case lewis_and_schrefler2hw:{//Lewis and Schrefler's model
    break;
  }
  case lewis_and_schrefler2hw_mefel:{//from mefel;
    rhocp = rhos0*cps0;//temporarily??!!
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
   function computes effective thermal conductivity of partially saturated concrete
   @param pw - pore water pressure
   @param t - temperature

   @retval lambdaeff - effective thermal conductivity of partially saturated soil
*/
double con_hwf2mat::get_lambdaeff(double pw,double t,long ipp)
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
    sw = get_sw(pw,t,ipp);
    lambdaeff = lambda_dry + (lambda_wet - lambda_dry)*sw; //this is linear approximation for bentonite:
    break;
  }
  case 2:{// linear approximation
    sw = get_sw(pw,t,ipp);
    lambdaeff = lambda_dry + (lambda_wet - lambda_dry)*(sw - sr_dry)/(sr_wet - sr_dry);
    break;
  }
  case 3:{//
    sw = get_sw(pw,t,ipp);
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


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
*/
void con_hwf2mat::matcond (matrix &d,long ri,long ci,long ipp)
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
void con_hwf2mat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwt(pw,t,ipp);// *scale_g;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_ktw(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_ktt1(pw,t,ipp);// *scale_g;//scaling
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hwf2mat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwt(pw,t,ipp);// *scale_g;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_ktw(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_ktt1(pw,t,ipp);// *scale_g;//scaling
  
  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;

}

/**
   function creates conductivity matrix of the material for 3D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void con_hwf2mat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    k = get_kwt(pw,t,ipp);// *scale_g;//scaling
  
  if((ri == 1) && (ci == 0))
    k = get_ktw(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    k = get_ktt1(pw,t,ipp);// *scale_g;//scaling
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}


/**
   function creates capacity matrix of the material
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hwf2mat::matcap (double &c,long ri,long ci,long ipp)
{
  double pw,t;
  c = 0.0;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capww(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = get_capwt(pw,t,ipp);// *scale_t;//scaling
  
  if((ri == 1) && (ci == 0))
    c = get_captw(pw,t,ipp);// *scale_pw;//scaling
  if((ri == 1) && (ci == 1))
    c = get_captt(pw,t,ipp);// *scale_g;//scaling

  check_math_errel(0);

}



/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void con_hwf2mat::matcond2 (matrix &d,long ri,long ci,long ipp)
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
void con_hwf2mat::matcond1d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c;
  double pw,t;
      
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 1) && (ci == 0)){
  }

  if((ri == 1) && (ci == 1)){
    a = get_ktt2a(pw,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,t,ipp);// *scale_t;//scaling
    
    d[0][0] = -1.0*a*Tm->ip[ipp].grad[0][0] - b*Tm->ip[ipp].grad[1][0] + (a*c)*Tp->gr[0];
  }  
}


/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hwf2mat::matcond2d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c;
  double pw,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  fillm(0.0,d);
  
  if((ri == 1) && (ci == 0)){
  }  
  if((ri == 1) && (ci == 1)){
    a = get_ktt2a(pw,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,t,ipp);// *scale_t;//scaling
    
    d[0][0] = -1.0*a*Tm->ip[ipp].grad[0][0] - b*Tm->ip[ipp].grad[1][0] + (a*c)*Tp->gr[0];
    d[0][1] = -1.0*a*Tm->ip[ipp].grad[0][1] - b*Tm->ip[ipp].grad[1][1] + (a*c)*Tp->gr[1];
  }
}



/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_hwf2mat::matcond3d_2 (matrix &d,long ri,long ci,long ipp)
{
  double a,b,c;
  double pw,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  fillm(0.0,d);
    
  if((ri == 1) && (ci == 0)){
  }  
  if((ri == 1) && (ci == 1)){
    a = get_ktt2a(pw,t,ipp);// *scale_t;//scaling
    b = get_ktt2b(pw,t,ipp);// *scale_t;//scaling
    c = get_ktt2c(pw,t,ipp);// *scale_t;//scaling
    
    d[0][0] = -1.0*a*Tm->ip[ipp].grad[0][0] - b*Tm->ip[ipp].grad[1][0] + (a*c)*Tp->gr[0];
    d[0][1] = -1.0*a*Tm->ip[ipp].grad[0][1] - b*Tm->ip[ipp].grad[1][1] + (a*c)*Tp->gr[1];
    d[0][2] = -1.0*a*Tm->ip[ipp].grad[0][2] - b*Tm->ip[ipp].grad[1][2] + (a*c)*Tp->gr[2];
  }
}



/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   ...
   @param ipp - number of integration point
   
   11.2.2003
*/
void con_hwf2mat::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
void con_hwf2mat::rhs_volume2 (double &cc,long ri,long /*ci*/,long ipp)
{
  double pw,t;
  
  pw = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];

  if(ri == 0){
    cc = get_fwu(pw,t,ipp);
  }
  if(ri == 1){
    cc = 0.0;
  }
}



/**
   function creates volume right-hand side matrix of the material for 1D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void con_hwf2mat::rhs1d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fw1(pw,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
  if(ri == 1){
    f = get_ft1(pw,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void con_hwf2mat::rhs2d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,t;
  
  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling
  
  if(ri == 0){
    f = get_fw1(pw,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
  if(ri == 1){
    f = get_ft1(pw,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   ...
   @param ipp - number of integration point
*/
void con_hwf2mat::rhs3d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,t;

  pw = Tm->ip[ipp].av[0];// *scale_pw;//scaling
  t  = Tm->ip[ipp].av[1];// *scale_t;//scaling

  if(ri == 0){
    f = get_fw1(pw,t,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
  if(ri == 1){
    f = get_ft1(pw,t,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
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
double con_hwf2mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  t = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_transcoeff_tt(pw,t,bc,ipp);// *scale_t;//scaling

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
double con_hwf2mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp, int flag)
{
  double c,pw,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  t = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,t,bc,ipp,flag);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_transcoeff_tt(pw,t,bc,ipp);// *scale_t;//scaling
  
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
double con_hwf2mat::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  t = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_nodval_ww(nodval,pw,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_nodval_tt(nodval,trc2,pw,t,bc,ipp);// *scale_t;//scaling

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
double con_hwf2mat::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,t;
  c = 0.0;
  
  pw = nodalval (nn,0);
  t = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_flux_ww(nodval,pw,t,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_flux_tt(nodval,trc2,pw,t,bc,ipp);// *scale_t;//scaling

  return (c);
}



/**
   function checks pore water pressure

   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void con_hwf2mat::waterpress_check(double &/*pw*/,double /*t*/,long /*ipp*/)
{
}

/**
   function corrects values of variables
   
   @param nv - array with variables
   
*/
void con_hwf2mat::values_correction (vector &nv,long /*ipp*/)
{
  //  pore water pressure
  waterpress_check(nv[0],nv[1],0);
  
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capww - capacity coefficient
*/
double con_hwf2mat::get_capww(double pw,double t,long ipp)
{
  double capww;
  double alpha,n,ks,sg,sw,rhow,kw,ds_dpc,rhogw,drhogw_dpw,dpgw_dpc;

  if(model_type == artificial2)
    capww = capww0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    alpha = get_alpha(pw,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,t,ipp);
    sw = get_sw(pw,t,ipp);
    rhow = get_rhow(t);
    kw = get_kw(pw,t,ipp);
    ds_dpc = -get_dsw_dpw(pw,t,ipp);
    drhogw_dpw = get_drhogw_dpw(pw,t);
    rhogw = get_rhogw(pw,t);
    dpgw_dpc = get_dpgw_dpc(pw,t);
    sg = 1.0 - sw;
    
    if(compress == 1){
      //compressible grains:
      capww = (alpha - n)/ks*sw*(rhogw*sg + rhow*sw) + rhow*sw*n/kw;
      capww = capww - sg*n*mw/t/gasr*dpgw_dpc;
      capww = capww - ((alpha - n)/ks*(rhogw*sg*(pw) + rhow*sw*(pw)) + n*(rhow - rhogw))*ds_dpc;
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
   @param t - temperature
   @param ipp - number of integration point
   
   @retval capwt - capacity coefficient
*/
double con_hwf2mat::get_capwt(double pw,double t,long ipp)
{
  double capwt;
  double sg,sw,alpha,n,ks,rhow,dsw_dt,pgw,dpgw_dt,rhogw;
  double betasw;

  if(model_type == artificial2)
    capwt = capwt0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    betasw = get_betasw(pw,t,ipp);
    sw = get_sw(pw,t,ipp);
    sg = 1.0 - sw;
    dpgw_dt = get_dpgw_dt(pw,t);
    pgw = get_pgw(pw,t);
    alpha = get_alpha(pw,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,t,ipp);
    rhogw = get_rhogw(pw,t);
    rhow = get_rhow(t);
    dsw_dt = get_dsw_dt(pw,t,ipp);

    if(compress == 1){
      //compressible grains:
      capwt = -betasw*rhow;
      capwt = capwt + (1.0-sw)*n*mw/t/gasr*(dpgw_dt-pgw/t);
      capwt = capwt + ((alpha - n)/ks*(rhogw*(1.0-sw)*(pw) + rhow*sw*(pw)) + n*(rhow - rhogw))*dsw_dt;
    }
    else{
      //incompressible grains:
      capwt = -betasw*rhow;
      capwt = capwt + (1.0-sw)*n*mw/t/gasr*(dpgw_dt-pgw/t);
      capwt = capwt + n*(rhow - rhogw)*dsw_dt;    
    }
  }
  return(capwt);
}



/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kww - conductivity coefficient
*/
double con_hwf2mat::get_kww(double pw,double t,long ipp)
{
  double kww;
  double rhow,krw,muw,kintr;
  double drhogw_dpw,dv;

  if(model_type == artificial2)
    kww = kww0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    rhow = get_rhow(t);
    krw = get_krw(pw,t,ipp);
    muw = get_muw(t);
    kintr = get_kintr(pw,t,ipp);
    drhogw_dpw = get_drhogw_dpw(pw,t);
    dv = get_dv(pw,t,ipp);
    
    kww = rhow*kintr*krw/muw;//water
    kww = kww + dv*drhogw_dpw;//water vapour
  }
  return(kww);
}



/**
   function creates conductivity coefficient of a general material
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kwt - conductivity coefficient
*/
double con_hwf2mat::get_kwt(double pw,double t,long ipp)
{
  double kwt;
  double drhogw_dt,dv;

  if(model_type == artificial2)
    kwt = kwt0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    drhogw_dt = get_drhogw_dt(pw,t);
    dv = get_dv(pw,t,ipp);
    
    kwt = dv*drhogw_dt;//water vapour
    //kwt = 0.0;
  }
  return(kwt);
}


/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captw - capacity coefficient

   03/10/2023 TKr
*/
double con_hwf2mat::get_captw(double pw,double t,long ipp)
{
  double captw,dhvap,rhow,alpha,n,ks,sw,kw,dsw_dpc;
  
  if(model_type == artificial2)
    captw = captw0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    alpha = get_alpha(pw,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,t,ipp);
    sw = get_sw(pw,t,ipp);
    kw = get_kw(pw,t,ipp);
    dsw_dpc = -get_dsw_dpw(pw,t,ipp);//this is equal to minus derivative with respect to suction    
    if(compress == 1){
      //compressible grains:
      captw = -1.0*dhvap*rhow*((alpha - n)/ks*sw*sw + sw*n/kw);
      captw = captw + dhvap*rhow*((alpha - n)/ks*sw*pw + n)*dsw_dpc;
    }
    else{
      //incompressible grains:
      captw = dhvap*rhow*n*dsw_dpc;
    }
  }

  return(captw);
}



/**
   function creates capacity coefficient of a general material
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval captt - capacity coefficient

   03/10/2023 TKr
*/
double con_hwf2mat::get_captt(double pw,double t,long ipp)
{
  double captt,rhocp,dhvap,betasw,rhow,alpha,n,ks,sw,dsw_dt;
   
  if(model_type == artificial2)
    captt = captt0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    rhocp = get_rhocp(pw,t,ipp);
    dhvap = get_dhvap(t);
    betasw = get_betasw(pw,t,ipp);
    rhow = get_rhow(t);
    alpha = get_alpha(pw,t,ipp);
    n = get_porosity(ipp);
    ks = get_ks(pw,t,ipp);
    sw = get_sw(pw,t,ipp);
    dsw_dt = get_dsw_dt(pw,t,ipp);
    
    if(compress == 1){
      //compressible grains:
      captt = rhocp;
      captt = captt + dhvap*betasw*rhow;
      captt = captt - dhvap*(rhow*((alpha - n)/ks*sw*pw + n)*dsw_dt);
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
   @param t - temperature
   @param ipp - number of integration point
   
   @retval ktw - conductivity coefficient
*/
double con_hwf2mat::get_ktw(double pw,double t,long ipp)
{
  double dhvap,rhow,krw,muw,ktw,kintr;
  
  if(model_type == artificial2)
    ktw = ktw0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    krw = get_krw(pw,t,ipp);
    muw = get_muw(t);
    kintr = get_kintr(pw,t,ipp);
    
    ktw = -dhvap*(rhow*kintr*krw/muw);//water
  }

  return(ktw);
}




/**
   function creates conductivity coefficient of a general material - first part
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt1 - conductivity coefficient
*/
double con_hwf2mat::get_ktt1(double pw,double t,long ipp)
{
  double lambdaeff,ktt1;
   
  if(model_type == artificial2)
    ktt1 = ktt0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    lambdaeff = get_lambdaeff(pw,t,ipp);
    
    ktt1 = lambdaeff;
  }
  return(ktt1);
}


/**
   function creates conductivity coefficient of a general material - second (A) part (convective term)
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2a - conductivity coefficient
*/
double con_hwf2mat::get_ktt2a(double pw,double t,long ipp)
{
  double n,sw,sg,cpw,rhow,krw,muw,ktt2a,kintr,rhogw,drhogw_dpw,dv;
    
  if(model_type == artificial2)
    ktt2a = 0.0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    n = get_porosity(ipp);
    sw = get_sw(pw,t,ipp);
    sg = 1.0 - sw;
    rhow = get_rhow(t);
    cpw = get_cpw();
    krw = get_krw(pw,t,ipp);
    muw = get_muw(t);
    kintr = get_kintr(pw,t,ipp);
    rhogw = get_rhogw(pw,t);
    drhogw_dpw = get_drhogw_dpw(pw,t);
    dv = get_dv(pw,t,ipp);
    
    ktt2a = n*sw*rhow*cpw*kintr*krw/muw;
    ktt2a = ktt2a + n*sg*rhogw*cpgw0*dv*drhogw_dpw;
  }
  return(ktt2a);
}


/**
   function creates conductivity coefficient of a general material - second (C) part (diffusion term)
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2b - conductivity coefficient

   17/10/2023 TKr
*/
double con_hwf2mat::get_ktt2b(double pw,double t,long ipp)
{
  double n,sw,sg,rhogw,drhogw_dt,dv;
  double ktt2b;
  
  if(model_type == artificial2)
    ktt2b = 0.0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    n = get_porosity(ipp);
    sw = get_sw(pw,t,ipp);
    sg = 1.0 - sw;
    rhogw = get_rhogw(pw,t);
    drhogw_dt = get_drhogw_dt(pw,t);
    dv = get_dv(pw,t,ipp);
    
    ktt2b = n*sg*rhogw*cpgw0*dv*drhogw_dt;
  }

  return(ktt2b);
}



/**
   function creates conductivity coefficient of a general material - second (B) part (convective term)
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval ktt2c - conductivity coefficient
*/
double con_hwf2mat::get_ktt2c(double pw,double t,long ipp)
{
  double rhow,ktt2c;
  
  if(model_type == artificial2)
    ktt2c = 0.0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    rhow = get_rhow(t);
    
    ktt2c = rhow;
  }
  return(ktt2c);
}


/**
   function creates right-hand side coefficient of a general material for c medium
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double con_hwf2mat::get_fw1(double pw,double t,long ipp)
{
  double fw1;
  double rhow,krw,muw,kintr;
  
  if(model_type == artificial2)
    fw1 = 0.0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    rhow = get_rhow(t);
    krw = get_krw(pw,t,ipp);
    muw = get_muw(t);
    kintr = get_kintr(pw,t,ipp);
    
    fw1 = rhow*kintr*krw*rhow/muw;
  }
  return(fw1);
}


/**
   function creates right-hand side coefficient of a general material for t medium
   @param pw - pore water pressure
   @param t - temperature
   @param ipp - number of integration point
*/

double con_hwf2mat::get_ft1(double pw,double t,long ipp)
{
  double ft1;
  double dhvap,rhow,krw,muw,kintr;

  if(model_type == artificial2)
    ft1 = 0.0;
  else{
    //water pressure check
    waterpress_check(pw,t,ipp);
    
    dhvap = get_dhvap(t);
    rhow = get_rhow(t);
    krw = get_krw(pw,t,ipp);
    muw = get_muw(t);
    kintr = get_kintr(pw,t,ipp);
    
    ft1 = -dhvap*rhow*kintr*krw*rhow/muw;
  }
  return(ft1);
}



/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval fwu - volumetric strain rate effect on pore water pressure; only for partially coupled METR version

   22/05/2018, TKr
*/
double con_hwf2mat::get_fwu(double pw, double t, long ipp)
{
  double depsv_r,sw,sg,fwu,alpha,dsr_depsv,n,rhow,rhogw;
  
  fwu = 0.0;
  
  switch (model_type){
  case artificial2:
  case lewis_and_schrefler:{
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's book
    if(vol_strain_effect == 1) {
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

      Tm->ip[ipp].eqother[1] = dsr_depsv; //this is not necessary

      sw = get_sw(pw,t,ipp);
      sg = 1.0 - sw;
      rhow = get_rhow(t);
      rhogw = get_rhogw(pw,t);
      alpha = get_alpha(pw,t,ipp);
      
      fwu = -1.0*(sg*rhogw + sw*rhow)*alpha*depsv_r;//volumetric strain effect        
      
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
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_transcoeff_ww(double pw,double t,long bc,long ipp)
{
  double trc;
  
  //water pressure check
  waterpress_check(pw,t,ipp);

  //other conditions types will be added according to the neeed
  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure

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



/**
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_transcoeff_ww(double pw,double t,long bc,long ipp,int flag)
{
  double trc;
  
  //water pressure check
  waterpress_check(pw,t,ipp);
  
  //other conditions types will be added according to the neeed
  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure
    
    trc = 1.0;
    
    break;
  }
  case 33:{//pore water pressure - flux
    if(flag == 1)
      trc=0.0;//into right hand side (matrix)
    else{
      trc=1.0;//into left hand side (flux)
    }
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

/**
   function creates correct new nodal value on the boundary (transmission) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_nodval_ww(double bv,double pw,double t,long bc,long ipp)
{
  double new_nodval;
  
  //water pressure check
  waterpress_check(pw,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure
    new_nodval = bv;
    break;
  }
  case 33:{//pore water pressure - flux
    new_nodval = bv - pw;//minus sign
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



/**
   function creates flux on the boundary (transmission - convective mass transfer) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_flux_ww(double bv,double pw,double t,long bc,long ipp)
{
  double flux,trc;
  
  //water pressure check
  waterpress_check(pw,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure
    flux = bv - pw;//minus sign
    break;
  }
  case 33:{//pore water pressure - flux
    flux = bv - pw;//minus sign
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
   @param t - temperature
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_transcoeff_tt(double pw,double t,long bc,long ipp)
{
  double trc;

  //water pressure check
  waterpress_check(pw,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    trc = 1.0;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  
  return(trc);
}




/**
   function creates correct new nodal value on the boundary (transmission) for t medium

   @param bv - value of prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual pore water pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_nodval_tt(double bv,double /*trr*/,double pw,double t,long bc,long ipp)
{
  double new_nodval;

  //water pressure check
  waterpress_check(pw,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
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


/**
   function creates flux on the boundary (transmission - convective mass transfer) for c medium

   @param bv - prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual pore water pressure on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_hwf2mat::get_transmission_flux_tt(double bv,double /*trr*/,double pw,double t,long bc,long ipp)
{
  double flux;

  //water pressure check
  waterpress_check(pw,t,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission - boundary flux
    flux = (bv - t);//minus sign
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
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param pw - pore water pressure on actual node
   @param t - temperature on actual node
*/

double con_hwf2mat::get_othervalue(long compother,long ipp,double *r)
{
  double other;

  switch (compother){
  case 0:{//capillary pressure
    other = r[0];
    break;
  }
  case 1:{//temperature in deg. C
    other = r[1] - 273.15;
    break;
  }
  case 2:{//relative humidity from Kelvin-Laplace law
    other = get_rh(r[0],r[1]);
    break;
  }
  case 3:{//saturation
    other = get_sw(r[0],r[1],ipp);
    break;
  }
  case 4:{//liquid water pressure
    other = r[0];
    break;
  }
  case 5:{//moisture content
    other = get_w(r[0],r[1],ipp);
    break;
  }    
  case 6:{//dsw_dpw
    other = get_dsw_dpw(r[0],r[1],ipp);
    break;
  }    
  case 7:{//dsw_dt
    other = get_dsw_dt(r[0],r[1],ipp);
    break;
  }    
  case 8:{//depsv
    other = Tm->ip[ipp].eqother[1];
    break;
  }   
  default:{
    print_err ("unknown type of component is required in function", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  check_math_errel(Tm->elip[ipp]);

  return (other);

}

/**
     function prints names of all variables in nodes
     @param out - output file
     @param compother - number of other components
*/
void con_hwf2mat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//capillary pressure
    fprintf (out,"Capilary pressure (Pa)        ");
    break;
  }
  case 1:{//temperature
    fprintf (out,"Temperature (K)               ");
    break;
  }
  case 2:{//relative humidity
    fprintf (out,"Relative humidity ()          ");
    break;
  }
  case 3:{//saturation
    fprintf (out,"Degree of saturation ()       ");
    break;
  }
  case 4:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
  case 5:{//moisture content
    fprintf (out,"Moisture (water) content (kg/kg)      ");
    break;
  }    
  case 6:{//dS_dpw
    fprintf (out,"dS_dpw Derivative of Saturation degree with respect to pore water pressure      ");
    break;
  }    
  case 7:{//dS_dt
    fprintf (out,"dS_dt Derivative of Saturation degree with respect to temperature      ");
    break;
  }    
  case 8:{//depsv
    fprintf (out,"Volumetric strain rate (-)       ");
    break;
  }    
  default:{
    print_err ("unknown type of component is required in function", __FILE__, __LINE__, __func__);
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
void con_hwf2mat::updateval (long /*ipp*/)
{
}



/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number

   23/05/2016, TKr
*/
void con_hwf2mat::initval(long /*ipp*/)
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
void con_hwf2mat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_water_press;
  dofname[1] = trf_temperature;
}


/**
   function returns temperature

   @param ipp - integration point number

   @retval t - temperature

   16/07/2018, TKr
*/
double con_hwf2mat::give_temperature(long ipp)
{
  double t;

  t = Tm->ip[ipp].av[1];

  return(t);
}


/**
   function returns effective pore pressure (pressure has negative sign - mechanical convention)

   @param ipp - integration point number

   @retval pw - pore water pressure
*/
double con_hwf2mat::give_effective_pore_pressure(long ipp)
{
  double pw,t,xi=1.0;
  
  pw = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];  
  
  switch (model_type){
  case lewis_and_schrefler:{//Lewis and Schrefler's book
    xi =  get_xi(pw,t,ipp);
    pw = -xi*pw;//returns effective pore pressure
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;

    xi = get_xi(pw,t,ipp); //effective stress factor
    pw = -xi*pw/mefel_units;//returns effective pore pressure //corrected units for mefel //basic units = Pa 
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }   
  return(pw);
}


/**
   function returns water pressure

   @param ipp - integration point number

   @retval pw - pore water pressure

   27/05/2016, TKr
*/
double con_hwf2mat::give_water_pressure(long ipp)
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
double con_hwf2mat::give_pore_pressure(long ipp)
{
  double pp,pw;

  pw = Tm->ip[ipp].av[0];

  pp = -pw;

  return(pp);
}


/**
   function computes suction stress s = -(pg - pw) = -pc;
   @param ipp - integration point number

   @retval suction - suction stress [Pa]

   27/05/2016, TKr
*/
double con_hwf2mat::give_suction(long ipp)
{
  double pw,suction;
  
  pw = Tm->ip[ipp].av[0];
  
  switch (model_type){
  case lewis_and_schrefler2hw:{
    suction = pw; //this is correct for Lewis and Schrefler's notation
    break;
  }
  case lewis_and_schrefler2hw_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel
    suction = pw;
    suction = suction/mefel_units;//corrected units for mefel //basic units = Pa
    break;
  }
  case van_genuchten2hw:{//partially saturated medium = Van Genuchten model
    pw = pw/mefel_units;
    suction = pw;//this is correct, because pore water pressure is negative
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
double con_hwf2mat::give_saturation_degree(long ipp)
{
  double pw,t,s;
  
  pw = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  

  s = get_sw(pw,t,ipp);
  
  return(s);
}

/**
   function computes relative humidity from Kelvin-Laplace law [-]
   @param pw - water capillary pressure
   @param t - temperature
   

   @retval rh - computes relative humidity rh [-]
*/
double con_hwf2mat::get_rh(double pw,double t)
{
  double rhow,rh;

  rhow = get_rhow(t);
  //numerical issue
  if(pw < -2.0e-09)
    rh = 0.001;
  else
    rh = exp(pw/rhow*mw/gasr/t);

  return(rh);
}


/**
   function computes water content w [kg/kg]
   @param pw - water capillary pressure
   @param t - temperature
   

   @retval w - computes water content w [kg/kg]
*/
double con_hwf2mat::get_w(double pw,double t,long ipp)
{
  double w,n,s,rhow;

  n = get_porosity(ipp);
  s = get_sw(pw,t,ipp);
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
void con_hwf2mat::give_reqntq(long *antq)
{
  switch (model_type){
  case lewis_and_schrefler2hw:{//Lewis and Schrefler's book
    break;
  }
  case lewis_and_schrefler2hw_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
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
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
}

/**
   function computes water vapour pressure = Kelvin equation
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval pgw - water vapour pressure = Kelvin equation
*/
double con_hwf2mat::get_rhogw(double pw,double t)
{
  double rhogw,pgw;
  
  pgw = get_pgw(pw,t);
  rhogw = pgw*mw/gasr/t;
  
  check_math_err();

  return(rhogw);
}


/**
   function computes partial derivative of rhogw with respect to pw (Kelvin equation)
   @param pw - pore water pressure
   @param t - temperature

   @retval drhogw_dpw - partial derivative of pgw with respect to pw (Kelvin equation)
*/
double con_hwf2mat::get_drhogw_dpw(double pw,double t)
{
  double drhogw_dpw,dpgw_dpw;

  dpgw_dpw = -get_dpgw_dpc(pw,t);
  
  drhogw_dpw = dpgw_dpw*mw/gasr/t;

  return(drhogw_dpw);
}


/**
   function computes partial derivative of pgw with respect to t (Kelvin equation)
   @param pw - pore water pressure
   @param t - temperature

   @retval get_drhogw_dt - partial derivative of rhogw with respect to t (Kelvin equation)
*/
double con_hwf2mat::get_drhogw_dt(double pw,double t)
{
  double drhogw_dt,dpgw_dt;

  dpgw_dt = get_dpgw_dt(pw,t);

  drhogw_dt = dpgw_dt*mw/gasr/t;

  return(drhogw_dt);
}

/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure
*/
double con_hwf2mat::get_pgws(double t)
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

  return(pgws);
}



/**
   function computes partial derivative of water vapour saturation pressure with respect to t
   @param t - temperature

   @retval dpgws_dt - partial derivative of water vapour saturation pressure with respect to t
*/
double con_hwf2mat::get_dpgws_dt(double t)
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
   function computes water vapour pressure = Kelvin equation
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval pgw - water vapour pressure = Kelvin equation
*/
double con_hwf2mat::get_pgw(double pw,double t)
{
  double pgw,rhow,pgws,tt,pc;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);

  pc = -pw;
  
  pgw = pgws*exp(-1.0*pc*mw/rhow/gasr/tt);

  check_math_err();

  return(pgw);
}


/**
   function computes partial derivative of pgw with respect to pc (Kelvin equation)
   @param pw - pore water pressure
   @param t - temperature

   @retval dpgw_dpc - partial derivative of pgw with respect to pc (Kelvin equation)
*/
double con_hwf2mat::get_dpgw_dpc(double pw, double t)
{
  double dpgw_dpc,pgws,rhow,tt,pc;
  
  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  
  pc = -pw;
    
  dpgw_dpc = -pgws*exp(-pc*mw/rhow/gasr/tt)*mw/rhow/gasr/tt;
  
  return(dpgw_dpc);
}



/**
   function computes partial derivative of pgw with respect to t (Kelvin equation)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval dpgw_dt - partial derivative of pgw with respect to t (Kelvin equation)
*/
double con_hwf2mat::get_dpgw_dt(double pw,double t)
{
  double dpgw_dt,dpgws_dt,tt,pgws,rhow,pc;
  //double drhow_dt;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  pc = -pw;
  dpgws_dt = get_dpgws_dt(tt);
  
  //drhow_dt = get_drhow_dt(tt);

  if(t < tcr)
    dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt) + pgws*exp(-pc*mw/rhow/gasr/tt)*(pc*mw/rhow/gasr/tt/tt);
    //actualized form:
    //dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt) + pgws*exp(-pc*mw/rhow/gasr/tt)*(pc*mw/rhow/gasr/tt)*(drhow_dt/rhow + 1/tt);
  else
    dpgw_dt = dpgws_dt*exp(-pc*mw/rhow/gasr/tt);

  return(dpgw_dt);
}



/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval deff - effective diffusion coefficient of vapour inside pores
*/
double con_hwf2mat::get_dv(double pw,double t,long ipp)
{
  double dv;
  double phi,sw,tau;

  dv = dv0;

  switch (deff_type){
  case 0:{//constant
    dv = dv0;
    break;
  }
  case 1:{//moisture dependent
    phi = get_porosity(ipp);
    sw = get_sw(pw,t,ipp);
    tau = tau0;
    dv = phi*(1.0 - sw)*tau*dv0*pow((t/273.15),1.8);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(dv);
}


/**
   function computes enthalpy of evaporation (latent heat of vaporization)
   @param t - temperature

   @retval - enthalpy of evaporation (latent heat of vaporization)
*/
double con_hwf2mat::get_dhvap(double t)
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
