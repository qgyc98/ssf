/*
    File:            masin_reten.cpp
    Author:          Tomas Krejci, 28/02/2023
    Purpose:         retention curve - extended formulation from Brooks and Correy according to Masin
    sources:         
    sr_type: 1       Masin, D. (2010), Predicting the dependency of a degree of saturation on void ratio and suction using effective stress principle for unsaturated soils. Int. J. Numer. Anal. Meth. Geomech., 34: 73-90. https://doi.org/10.1002/nag.808
    sr_type: 2       Masin, D, (2013), Double structure hydromechanical coupling formalism and a model for unsaturated expansive clays, Engineering Geology, 165, 73-88, ISSN 0013-7952, https://doi.org/10.1016/j.enggeo.2013.05.026.
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalt.h"
#include "masin_reten.h"

masin_reten::masin_reten()
{
  sr_type = 1;

  se0 = 2700.0;
  e0 = 0.5;
  tref = 294.0;
  lambdap0 = 0.15;
  at = 0.118;
  bt = -0.000154;
  gamma0 = 0.55;
  ae = 1.0;
}

masin_reten::~masin_reten()
{}

/**
   function reads parameters
   
   @param in - input file

   TKr 28/02/2023
*/
void masin_reten::read(XFILE *in)
{
  xfscanf (in,"%ld %lf %lf %lf %lf %lf %lf %lf %lf", &sr_type, &se0,&e0,&ae,&lambdap0,&at,&bt,&tref,&gamma0);
}

/**
   function reads parameters
   
   @param in - input file

   TKr 28/02/2023
*/
void masin_reten::print(FILE *out)
{
  fprintf (out,"\n%ld %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", sr_type,se0,e0,ae,lambdap0,at,bt,tref,gamma0);
}

/**
   function computes effective stress factor psi

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 28/02/2023
*/
double masin_reten::psi(double suction, double dsuc, double e, double temp)
{
  double psi,se,sr,lambdap,chi;
  
  if(e == 0.0)//intial porosity
     e = e0;

  se = get_se(suction,dsuc,e,temp);
  sr = sw(suction,dsuc,e,temp);

  switch (sr_type){
  case 0:
  case 1:{
    if(suction < se)
      psi = 1.0-gamma0;
    else{
      lambdap = get_lambdap(suction,dsuc,e,temp);
      chi = pow(sr,(gamma0/lambdap));
      psi = (1-gamma0)*chi;
    }
    break;
  }
  case 2:{
    if(suction < se)
      psi = 1.0;
    else{
      psi = sr;
    }
    break;
  }
  default:{
    print_err("Unknown WRC type is required",__FILE__,__LINE__,__func__);
  }
  }

  check_math_errel(0);

  return(psi);
}


/**
   function computes parameter se - suction air-entry value

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 28/02/2023
*/
double masin_reten::get_se(double suction, double dsuc, double e, double temp)
{
  double se;

  switch (sr_type){
  case 0:{
    se = se0;
    break;
  }
  case 1:{
    se = se0*e0/e;
    break;
  }
  case 2:{
    se = se0*e0/e;//drying branch
    if(dsuc >= 0.0)
      se = ae*se;//wetting branch
    break;
  }
  default:{
    print_err("Unknown WRC type is required",__FILE__,__LINE__,__func__);
  }
  }

  return(se);
}


/**
   function computes derivative of parameter se  with repsect to e prosity

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 28/02/2023
*/
double masin_reten::get_dse_de(double suction, double dsuc, double e, double temp)
{
  double dse_de;

  switch (sr_type){
  case 0:{
    dse_de = 0.0;
    break;
  }
  case 1:{
    dse_de = -se0*e0/e/e;
    break;
  }
  case 2:{
    dse_de = -se0*e0/e/e;//drying branch
    if(dsuc >= 0.0)
      dse_de = ae;//wetting branch
    break;
  }
  default:{
    print_err("Unknown WRC type is required",__FILE__,__LINE__,__func__);
  }
  }

  return(dse_de);
}


/**
   function computes parameter lambdap

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 28/02/2023
*/
double masin_reten::get_lambdap(double suction, double dsuc, double e, double temp)
{
  double lambdap,chi0;
  
  chi0 = pow((se0/suction),gamma0);

  lambdap = gamma0/log(chi0)*log((pow(chi0,(lambdap0/gamma0)) - chi0)*pow((e/e0),(gamma0-1)) + chi0);

  return(lambdap);
}




/**
   function computes saturation degree

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 28/02/2023
*/
double masin_reten::sw(double suction, double dsuc, double e, double temp)
{
  double sr,se,gamma,lambdap;
  
  if(e == 0.0)//intial porosity
    e = e0;
  
  se = get_se(suction,dsuc,e,temp);
  
  if(suction < se)
    sr = 1.0;
  else{
    switch (sr_type){
    case 0:
    case 1:{
      lambdap = get_lambdap(suction,dsuc,e,temp);
      gamma = lambdap;
      break;
    }
    case 2:{
      gamma = gamma0;
      break;
    }
    default:{
      print_err("Unknown WRC type is required",__FILE__,__LINE__,__func__);
    }
    }
    
    sr = pow((se/suction),gamma);
  }
  
  check_math_errel(0);
  
  return(sr);
}


/**
   function computes partial derivative of degree of saturation with respect to pw, specific water content

   @param pw - water pressure
   @retval dsw_dpw - partial derivative of degree of saturation with respect to pw

   TKr 28/02/2023
*/
double masin_reten::dsw_dpw(double suction, double dsuc, double e, double temp)
{
  double dsr_ds,se,lambdap,gamma;
  
  if(e == 0.0)//intial porosity
    e = e0;

  se = get_se(suction,dsuc,e,temp);

  if(suction < se)
    dsr_ds = 0.0;
  else{
    switch (sr_type){
    case 0:
    case 1:{
      lambdap = get_lambdap(suction,dsuc,e,temp);
      gamma = lambdap;
      break;
    }
    case 2:{
      gamma = gamma0;
      break;
    }
    default:{
      print_err("Unknown WRC type is required",__FILE__,__LINE__,__func__);
    }
    }
    
    dsr_ds = pow((se/suction),gamma-1)*(-se)/(suction*suction);//not completed??!!, derivative of lambdap must be added
  }
  
  check_math_errel(0);

  return(dsr_ds);
}


/**
   function computes partial derivative of degree of saturation with respect to t

   @param pw - water pressure
   @retval dsw_dt - partial derivative of degree of saturation with respect to t

   TKr 28/02/2023
*/
double masin_reten::dsw_dt(double suction, double dsuc, double e, double temp)
{
  double dsr_dt;
  
  dsr_dt = 0.0;//p. 167
  
  return(dsr_dt);
}


/**
   function computes partial derivative of degree of saturation with respect to epsv, volumetric strain

   @param pw - water pressure
   @retval dsw_depsv - partial derivative of degree of saturation with respect to pw

   TKr 28/02/2023
*/
double masin_reten::dsw_depsv(double suction, double dsuc, double e, double temp)
{
  double dsr_depsv,dsr_de,dsr_dse,se,dse_de,lambdap,gamma;
  
  se = get_se(suction,dsuc,e,temp);

  if(suction < se)
    dsr_depsv = 0.0;
  else{
    if(e == 0.0)//intial porosity
      e = e0;
    
    dse_de = get_dse_de(suction,dsuc,e,temp);
    
    switch (sr_type){
    case 0:
    case 1:{
      lambdap = get_lambdap(suction,dsuc,e,temp);
      gamma = lambdap;
      break;
    }
    case 2:{
      gamma = gamma0;
      break;
    }
    default:{
      print_err("Unknown WRC type is required",__FILE__,__LINE__,__func__);
    }
    }
    
    dsr_dse = gamma*pow((se/suction),(gamma-1))/suction;
    dsr_de = dsr_dse*dse_de;
    dsr_depsv = dsr_de*(1+e0);//not completed??!!, derivative of lambdap must be added
  }

  check_math_errel(0);

  return(dsr_depsv);
}
