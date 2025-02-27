/*
  File:             bazant_reten.cpp
  Author:           Tomas Krejci, 09/12/2022
  Purpose:          Calculates the retention curve assuming Bazant formulation extended for high temperature
  sources:          C60sifel.f90 and SATBGMsifel.f90 from Padova
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "constrel.h"
#include "bazant_reten.h"
#include "globalt.h"

bazant_reten::bazant_reten()
{
  t0 = 273.15;   //reference temperature
  p0 = 101325.0; //pressure
  tcr = 647.3;   //critical point of water
}
bazant_reten::~bazant_reten()
{}

/**
   function computes degree of saturation (desorption curve)

   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double bazant_reten::sat(double pc,double t)
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
      fprintf (Outt,"\n\n Uprava saturace, function sat (C60bazant.cpp)");
    }

  /*
    fprintf (Outt,"\n\n");
    fprintf (Outt,"pc = %lf\n",pc);
    fprintf (Outt,"t  = %lf\n",t);
    fprintf (Outt,"rh = %lf\n",rh);
    fprintf (Outt,"s  = %lf\n",sw);
    fprintf (Outt,"mt = %lf\n",mt);
  */
  
  //fprintf(Outt,"sw = %e, pc = %e, t = %e, ipp = %ld\n",sw,pc,t,0);
  //fflush(Outt);

  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double bazant_reten::dsat_dpc(double pc,double t)
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
      fprintf (Outt,"\n\n Uprava saturace, function dsat_dpc (C60bazant.cpp)");
    }

  ssmh = sw/mt/rh;
  dsdc = ssmh*drhdc;
  
  //if(sw <= 1.0e-10)//limit for saturation level
  //dsdc=0.0;
  
  dsw_dpc = dsdc;  

  //fprintf(Outt,"dsw_dpc = %e, pc = %e, t = %e, ipp = %ld\n",dsw_dpc,pc,t,0);
  //fflush(Outt);

  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double bazant_reten::dsat_dt(double pc,double t)
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
      fprintf (Outt,"\n\n Uprava saturace, function dsat_dt (C60bazant.cpp)");
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
  //dsdt=0.0;//debug??!!
  
  dsw_dt = dsdt;

  //fprintf(Outt,"dsw_dt = %e, pc = %e, t = %e, ipp = %ld\n",dsw_dt,pc,t,0);
  //fflush(Outt);

  return(dsw_dt);
}
