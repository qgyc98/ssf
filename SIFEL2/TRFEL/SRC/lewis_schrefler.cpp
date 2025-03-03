/*
    File:            lewis_schrefler.cpp
    Author:          Tomas Krejci, 13/04/2018
    Purpose:         retention curve accorging to Lewis and Schrefler's book p. 167
    sources:         Lewis and Schrefler's book p. 167
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalt.h"
#include "lewis_schrefler.h"

lewis_reten::lewis_reten()
{
  sr = 0.2; //0.06689; //residual saturation 
}

lewis_reten::~lewis_reten()
{}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void lewis_reten::read(XFILE *in)
{
  xfscanf (in,"%lf", &sr);
}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void lewis_reten::print(FILE *out)
{
  fprintf (out,"\n%lf\n", sr);
}

/**
   function computes saturation degree

   @param pw - pore water pressure
   @retval sw - saturation degree
   
   TKr 3/04/2018
*/
double lewis_reten::sw(double pw)
{
  double sw;

  if(pw > 0.0)
    {
      sw = 1.0;
    }
  else
    {
      sw = 1.0 - (1.9722e-11)*pow(-pw,2.4279);//p. 167 - Liakopoulos test
    }
  return(sw);
}


/**
   function computes partial derivative of degree of saturation with respect to pw, specific water content

   @param pw - water pressure
   @retval dsw_dpw - partial derivative of degree of saturation with respect to pw

   TKr 3/04/2018
*/
double lewis_reten::dsw_dpw(double pw)
{
  double dsw_dpw;
  
  if(pw > 0.0)
    {
      dsw_dpw = 0.0;
    }
  else
    {
      dsw_dpw = 2.4279*(1.9722e-11)*pow(-pw,(2.4279-1.0));//p. 167
    }
  return(dsw_dpw);
}


/**
   function computes partial derivative of degree of saturation with respect to t

   @param pw - water pressure
   @retval dsw_dt - partial derivative of degree of saturation with respect to t

   TKr 3/04/2018
*/
double lewis_reten::dsw_dt(double /*pw*/)
{
  double dsw_dt;
  
  dsw_dt = 0.0;//p. 167
  
  return(dsw_dt);
}


/**
   function computes water relative permeability

   @param sw - saturation degree   
   @retval krw - water relative permeability

   TKr 3/04/2018
*/
double lewis_reten::krw(double sw)
{
  double krw;

  krw = 1.0-2.207*pow((1.0-sw),1.0121);//liakopoulos

  return(krw);
}

/**
   function computes degree of saturation se (internal material variable)
   @param sw - saturation degree

   @retval se - degree of saturation se

   TKr 3/04/2018
*/
double lewis_reten::get_se(double sw)
{
  double se;
  
  se = (sw - sr)/(1.0 - sr);//p. 169

  return(se);
}


/**
   function computes gas relative permeability

   @param pw - pore water pressure
   @retval krg - water relative permeability

   TKr 3/04/2018
*/
double lewis_reten::get_krg(double pw)
{
  double krg,s,se;

  s = sw(pw);
  se = get_se(s);

  krg = pow((1.0-se),2.0)*(1.0 - pow(se,(5.0/3.0)));//p. 169//Liakopoulos
  
  if((krg < 0.0001) || (s >= 0.998))
    krg = 0.0001;

  return(krg);
}
