/*
    File:            van_genuchten.cpp
    Author:          Tomas Krejci, 13/04/2018
    Purpose:         retention curve accorging to Van_Genuchten model
    sources:                              
    classical expr:  M. Th. van Genuchten, A closed equation for predicting the hydraulic conductivity of
		     unsaturated soils, Journal Soil Science Society of America 44 (1980), 892-898;
		     pore water pressure (suction) is assumed in kPa (kilo Pascals)

    extension 1:     A full-scale in situ heating test for high-level nuclear waste disposal: observations, analysis and interpretation
                     A. Gens, M. Sanchez, L. Do N. Guimaraes, E. E. Alonso, A. Lloret, S. Olivella, M. V. Villar, and F. Huertas
                     Geotechnique 2009 59:4, 377-399
  

    extension 2:     ??

    extension 3:     Abel Carlos Jacinto, Maria Victoria Villar, Roberto Gomez-Espina, Alberto Ledesma,
                     Adaptation of the van Genuchten expression to the effects of temperature and density for compacted bentonites, Applied Clay Science,
                     Volume 42, Issues 3–4, 2009,Pages 575-582,ISSN 0169-1317,https://doi.org/10.1016/j.clay.2008.04.001.     
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalt.h"
#include "van_genuchten.h"

van_genuchten_reten::van_genuchten_reten()
{ 
  vg_ret_type = 0; //retention curve type

  ssat = 1.0;      //saturation degree [-]
  sirr = 0.2;      //irreversible saturation degree [-]
  ksat = 0.1;      //saturation permeability [??]

  //first fomulation
  gamaw = 10.0;    //water weigth [kN/m^3]
  delta = 20.0;    //parameter in m^(-1); delta = 20.0 [1/m]; delta = 0.2 [1/cm]
  expn = 2.0;      //parameter n - exponent
  expm = 0.0;      //parameter m - exponent

  //second formulation
  p0 = 0.0;        //suction - pore-air entry value
  pd = 0.0;        //suction value at fully dry conditions
  lambda0 = 1.0;   //parameter - exponent
  lambdap = 1.0;   //parameter - exponent
  sig0 = 0.072;    //surface tension water-gas at reference temperature t0
  t0 = 293.15;     //reference temperature t0 (20 deg.C)
}

van_genuchten_reten::~van_genuchten_reten()
{}


/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void van_genuchten_reten::read(XFILE *in)
{
  //retention curve type:
  xfscanf (in,"%ld",&vg_ret_type);
  

  switch (vg_ret_type){
  case 0:{//classical expression
    xfscanf (in,"%lf %lf %lf %lf %lf %lf", &ssat, &sirr, &expn, &expm, &ksat, &delta);
    break;
  }
  case 1:{//other expression
    xfscanf (in,"%lf %lf %lf %lf", &ssat, &sirr, &p0, &lambda0);
    break;
  }
  case 2:{//extended for hight suction
    xfscanf (in,"%lf %lf %lf %lf %lf %lf", &ssat, &sirr, &p0, &lambda0, &pd, &lambdap);
    break;
  }
  case 3:{//extended for hight suction and temperature effect
    xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf", &ssat, &sirr, &p0, &lambda0, &pd, &lambdap, &sig0, &t0);
    break;
  }
  default:{
    print_err ("unknown van genuchten retention curve type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void van_genuchten_reten::print(FILE *out)
{
  switch (vg_ret_type){
  case 0:{//classical expression
    fprintf (out,"\n%lf %lf %lf %lf %lf %lf\n", ssat, sirr, expn, expm, ksat, delta);
    break;
  }
  case 1:{//other expression
    fprintf (out,"\n%lf %lf %lf %lf\n", ssat, sirr, p0, lambda0);
    break;
  }
  case 2:{//extended for hight suction
    fprintf (out,"\n%lf %lf %lf %lf %lf %lf\n", ssat, sirr, p0, lambda0, pd, lambdap);
    break;
  }
  case 3:{//extended for hight suction and temperature effect
    fprintf (out,"\n%lf %lf %lf %lf %lf %lf %lf %lf\n", ssat, sirr, p0, lambda0, pd, lambdap, sig0, t0);
    break;
  }
  default:{
    print_err ("unknown van genuchten retention curve type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
}

/**
   function computes saturation degree

   @param pw - pore water pressure
   @param t - temperature
   @retval sw - saturation degree
   
   TKr 3/04/2018
*/
double van_genuchten_reten::sw(double pw, double t)
{
  double sw,hp,theta,fd=1.0,ft=1.0,sigt,a;

  switch (vg_ret_type){
  case 0:{//classical expression
    pw = pw/1000.0;//in kPa
    
    if(pw >= 0.0){
      sw = 1.0;//for positive pressure
    }
    else{
      hp = pw/gamaw;
      if(expm == 0.0)
	expm = 1.0 - 1.0/expn;
      theta = pow(1+(pow(delta*fabs(hp),expn)),-expm);
      sw = sirr + (ssat - sirr)*theta;
    }
    break;
  }
  case 1:{//other expression (in Pa)
    if(pw >= 0.0){
      sw = 1.0;//for positive pressure (suction)
    }
    else{
      expm = 1.0/(1.0-lambda0);
      theta = pow(1+(pow((-pw/p0),expm)),-lambda0);
      sw = sirr + (ssat - sirr)*theta;
      //sw = theta;
    }
    break;
  }
  case 2:{//extended for hight suction (in Pa)
    if(pw >= 0.0){
      sw = 1.0;//for positive pressure (suction)
    }
    else{
      expm = 1.0/(1.0-lambda0);
      fd = pow((1.0-(-pw)/pd),lambdap);
      theta = pow(1+(pow((-pw/p0),expm)),-lambda0)*fd;
      sw = sirr + (ssat - sirr)*theta;
    }
    break;
  }
  case 3:{//extended for hight suction and temperature effect (in Pa)
    if(pw >= 0.0){
      sw = 1.0;//for positive pressure (suction)
    }
    else{
      if (t < 633.15){
	a = (374.15-(t-273.15))/647.3;// temperature in deg.C
	sigt = (1.0 - 0.625*a)*(0.2358*pow(a,1.256));
      }	
      else{
	sigt = 0.0019106*exp(0.05*(633.15-t));
      }
      p0 = p0*sigt/sig0;
      expm = 1.0/(1.0-lambda0);
      fd = pow((1.0-(-pw)/pd),lambdap);
      ft = 1.0;//temp
      theta = pow(1+(pow((-pw/p0),expm)),-lambda0)*fd*ft;
      sw = sirr + (ssat - sirr)*theta;
    }
    break;
  }
  default:{
    print_err ("unknown van genuchten retention curve type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  if(sw > 1.0)
    sw = 1.0;

  return(sw);
}


/**
   function computes partial derivative of degree of saturation with respect to pw, specific water content
   
   @param pw - water pressure
   @param t - temperature
   @retval dsw_dpw - partial derivative of degree of saturation with respect to pw
   
   TKr 3/04/2018
*/
double van_genuchten_reten::dsw_dpw(double pw, double t)
{
  double sr,dsw_dpw,hp,fd,theta,dtheta_dpw,dfd_dpw,sigt,a;

  sr = sw(pw,t);
  
  switch (vg_ret_type){
  case 0:{//classical expression
    pw = pw/1000.0;//in kPa
    
    if(pw >= 0.0 || sr == 1.0){
      dsw_dpw = 0.0;
    }
    else{
      hp = pw/gamaw;
      if(expm == 0.0)
	expm = 1.0 - 1.0/expn;
      dsw_dpw = (expm*expn*delta*(ssat - sirr)*pow(delta*fabs(hp),(expn-1)))/(pow(1+(pow(delta*fabs(hp),expn)),(expm+1)))/gamaw;

      dsw_dpw =  dsw_dpw/1000.0;//basic units are Pa
    }
    break;
  }
  case 1:{//other expression (in Pa)
    if(pw >= 0.0 || sr == 1.0){
      dsw_dpw = 0.0;
    }
    else{
      expm = 1.0/(1.0-lambda0);
      theta = pow(1+(pow((-pw/p0),expm)),-lambda0);
      dtheta_dpw = -lambda0*pow(1+(pow((-pw/p0),expm)),-lambda0-1.0)*expm*pow((-pw/p0),expm-1)*(-1.0/p0);
      dsw_dpw = (ssat - sirr)*dtheta_dpw;
    }
    break;
  }
  case 2:{//extended for hight suction (in Pa)
    if(pw >= 0.0 || sr == 1.0){
      dsw_dpw = 0.0;
    }
    else{
      expm = 1.0/(1.0-lambda0);
      fd = pow((1.0-pw/pd),lambdap);
      theta = pow(1+(pow((-pw/p0),expm)),-lambda0);
      dtheta_dpw = -lambda0*pow(1+(pow((-pw/p0),expm)),-lambda0-1.0)*expm*pow((-pw/p0),expm-1)*(-1.0/p0);
      dfd_dpw = lambdap*pow((1.0-pw/pd),lambdap-1.0)*(1.0/pd);
      dsw_dpw = dtheta_dpw*fd + theta*dfd_dpw;
      dsw_dpw = (ssat - sirr)*dsw_dpw;
    }
    break;
  }
  case 3:{//extended for hight suction and temperature effect (in Pa)
    if(pw >= 0.0 || sr == 1.0){
      dsw_dpw = 0.0;
    }
    else{
      expm = 1.0/(1.0-lambda0);
      fd = pow((1.0-pw/pd),lambdap);

      if (t < 633.15){
	a = (374.15-(t-273.15))/647.3;// temperature in deg.C
	sigt = (1.0 - 0.625*a)*(0.2358*pow(a,1.256));
      }	
      else{
	sigt = 0.0019106*exp(0.05*(633.15-t));
      }
      p0 = p0*sigt/sig0;
      theta = pow(1+(pow((-pw/p0),expm)),-lambda0);
      dtheta_dpw = -lambda0*pow(1+(pow((-pw/p0),expm)),-lambda0-1.0)*expm*pow((-pw/p0),expm-1)*(-1.0/p0);
      dfd_dpw = lambdap*pow((1.0-pw/pd),lambdap-1.0)*(1.0/pd);
      dsw_dpw = dtheta_dpw*fd + theta*dfd_dpw;
      dsw_dpw = (ssat - sirr)*dsw_dpw;
    }
    break;
  }
  default:{
    print_err ("unknown van genuchten retention curve type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(dsw_dpw);
}


/**
   function computes partial derivative of degree of saturation with respect to t

   @param pw - water pressure
   @param t - temperature
   @retval dsw_dt - partial derivative of degree of saturation with respect to t

   TKr 3/04/2018
*/
double van_genuchten_reten::dsw_dt(double pw, double t)
{
  double sr,dsw_dt,fd,a,da_dt,dsigt_dt,dp0_dt,dtheta_dt;
  
  sr = sw(pw,t);

  switch (vg_ret_type){
  case 0://classical expression (in kPa)
  case 1://other expression (in Pa)
  case 2:{//extended for hight suction (in Pa)
    dsw_dt = 0.0;
    break;
  }
  case 3:{//extended for hight suction and temperature effect (in Pa) - not implemented yet
    dsw_dt = 0.0;
    
    if(pw >= 0.0 || sr == 1.0){
      dsw_dt = 0.0;
    }
    else{
      expm = 1.0/(1.0-lambda0);
      fd = pow((1.0-pw/pd),lambdap);
      
      if (t < 633.15){
	a = (374.15-(t-273.15))/647.3;// temperature in deg.C
	da_dt = -1.0/647.3;
	//sigt = (1.0 - 0.625*a)*(0.2358*pow(a,1.256));
	dsigt_dt = -0.332478*pow(a,0.256)*(a-0.89078)*da_dt;
      }	
      else{
	//sigt = 0.0019106*exp(0.05*(633.15-t));
	dsigt_dt = -5.335571e9*pow(0.951229,t);
      }
      dp0_dt = p0*dsigt_dt/sig0;
      
      dtheta_dt = -lambda0*pow(1+(pow((-pw/p0),expm)),-lambda0-1.0)*expm*pow((-pw/p0),expm-1)*(pw/p0/p0)*dp0_dt;
      dsw_dt = dtheta_dt*fd;
      dsw_dt = (ssat - sirr)*dsw_dt;
    }
    break;
  }
  default:{
    print_err ("unknown van genuchten retention curve type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(dsw_dt);
}


/**
   function computes water relative permeability

   @param pw - water pressure
   @param t - temperature
   @retval krw - water relative permeability

   TKr 3/04/2018
*/
double van_genuchten_reten::get_krw(double pw, double t)
{
  double krw,hp;

  switch (vg_ret_type){
  case 0:{//classical expression
    pw = pw/1000.0;
    
    if(pw >= 0.0)
      pw = 0.0;//for positive pressure
    if(expm == 0.0)
      expm = 1.0 - 1.0/expn;
    hp = pw/gamaw;
    krw = pow(1.0-pow(delta*fabs(hp),expn)*pow(1.0+pow(delta*fabs(hp),expn),-expm),2.0)/
      pow(1.0+pow(delta*fabs(hp),expn),(expm/2.0));
    //krw = 1.0/(pow(1.0+(pow(delta*fabs(hp),expn)),(3.0/expm)));//alternative expression according to Fatt and klikof
    //krw = 1.0/(pow(1.0+(pow(delta*fabs(hp),expn)),(1.0/expm)));//alternative expression according Warrick
    break;
  }
  case 1:{//other expression 
    double sr = 0.0;

    sr = sw(pw,t);    
    krw = pow(sr,(1.0/lambda0));
    krw = pow((1.0-krw),lambda0);
    krw = pow(krw,2.0);
    krw = sqrt(sr)*krw;
    break;
  }
  default:{
    print_err ("unknown van genuchten retention curve type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(krw);
}


/**
   function computes water conductivity of porous material 

   @retval k - conductivity

   13/04/2018, TKr
*/
double van_genuchten_reten::get_k(double /*pw*/)
{
  double k;

  k = ksat/86400.0/gamaw;

  return(k);
}
