/*
  File:             baroghel_reten.cpp
  Author:           Tomas Krejci, 09/12/2022
  Purpose:          Calculates the retention curve assuming Baroghel formulation extended for high temperature
  sources:          C60sifel.f90 and SATBGMsifel.f90 from Padova
  
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "baroghel_reten.h"
#include "globalt.h"

baroghel_reten::baroghel_reten()
{
  t0 = 273.15;   //reference temperature
  tcr = 647.3;   //critical point of water

  //for concrete BO:
  b = 2.2748;     //retention curve coefficient b - exponent
  q3 = 18.6237e6; //retention curve coefficient q3[Pa]
  q2 = 7.0e6;     //retention curve coefficient q2[Pa]
  n = 1.2;        //retention curve coefficient N - exponent
  z = 0.5;        //z is governing the transition through critical temperature

  //for concrete BH:
  //b = 2.0601;     //retention curve coefficient b - exponent
  //q3 = 46.9364e6; //retention curve coefficient q3[Pa]
  //q2 = 7.0e6;     //retention curve coefficient q2[Pa]
  //n = 1.2;        //retention curve coefficient N - exponent
  //z = 0.5;        //z is governing the transition through critical temperature

  //for concrete CO:
  //b = 2.1684;     //retention curve coefficient b - exponent
  //q3 = 37.5479e6; //retention curve coefficient q3[Pa]
  //q2 = 7.0e6;     //retention curve coefficient q2[Pa]
  //n = 1.2;        //retention curve coefficient N - exponent
  //z = 0.5;        //z is governing the transition through critical temperature

  //for concrete CH:
  //b = 1.9570;     //retention curve coefficient b - exponent
  //q3 = 96.2837e6; //retention curve coefficient q3[Pa]
  //q2 = 7.0e6;     //retention curve coefficient q2[Pa]
  //n = 1.2;        //retention curve coefficient N - exponent
  //z = 0.5;        //z is governing the transition through critical temperature
}

baroghel_reten::~baroghel_reten()
{
}


/**
   function reads parameters
   
   @param in - input file
*/
void baroghel_reten::read(XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf", &b, &q2, &q3, &n, &z);
}


/**
   function prints parameters
   
   @param out - output file
*/
void baroghel_reten::print(FILE *out)
{
  fprintf (out,"  %e %e %e %e %e\n", b, q2, q3, n, z);
}


/**
   function computes degree of saturation (sorption curve)
   @param pc - capillary pressure
   @param t - temperature

   @retval sw - degree of saturation
*/
double baroghel_reten::sw(double pc,double t)
{
  double sw,n,a,b,tt,q0,q2,q3,z,e,e0,g;
  
  b = 2.27; 
  q3 = 18.62e6;
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
  
  return(sw);
}

/**
   function computes partial derivative of degree of saturation with respect to pc
   @param pc - capillary pressure
   @param t - temperature

   @retval ds_dpc - partial derivative of degree of saturation with respect to pc
*/
double baroghel_reten::dsw_dpc(double pc,double t)
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

  return(dsw_dpc);
}

/**
   function computes partial derivative of degree of saturation with respect to t
   @param pc - capillary pressure
   @param t - temperature

   @retval dsw_dt - partial derivative of degree of saturation with respect to t
*/
double baroghel_reten::dsw_dt(double pc,double t)
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

  return(dsw_dt);
}
