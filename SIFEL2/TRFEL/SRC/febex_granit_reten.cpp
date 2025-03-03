/*
    File:            febex_granit_retn.cpp
    Author:          Tomas Krejci, 19/06/23
    Purpose:         retention curve accorging to FEBEX experiment - granit Grimsel
    sources:         DECOVALEX project report
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globalt.h"
#include "febex_granit_reten.h"

febex_granit_reten::febex_granit_reten()
{ 
  consta = 1.74e6; //reference suction Pa
  n1 = 1.68;       //exponent 1
  n2 = 0.405;      //exponent 2

  //FEBEX granit krw:
  m1 = 0.5;
  m2 = 1.68;
  m3 = 0.595;
  m4 = 2.0;
}

febex_granit_reten::~febex_granit_reten()
{}


/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void febex_granit_reten::read(XFILE *in)
{
  xfscanf (in,"%lf %lf %lf", &consta, &n1, &n2);
  xfscanf (in,"%lf %lf %lf %lf", &m1, &m2, &m3, &m4);
}

/**
   function reads parameters
   
   @param in - input file

   TKr 3/04/2018
*/
void febex_granit_reten::print(FILE *out)
{
  fprintf (out,"  %lf  %lf  %lf", consta, n1, n2);
  fprintf (out,"  %lf  %lf  %lf  %lf", m1, m2, m3, m4);
}

/**
   function computes saturation degree

   @param s - suction
   @retval sw - saturation degree
   
   TKr 3/04/2018
*/
double febex_granit_reten::sw(double s)
{
  double sw=0.0;
  
  if(s <= 0.0)
    sw = 1.0;
  else{
    sw = pow((s/consta),(1.0/n2));
    sw = pow((sw + 1),(-1.0/n1));
  }
  return(sw);
}


/**
   function computes partial derivative of degree of saturation with respect to pw, specific water content
   
   @param s - suction
   @retval dsw_ds - partial derivative of degree of saturation with respect to suction
   
   TKr 3/04/2018
*/
double febex_granit_reten::dsw_ds(double s)
{
  double dsw_ds=0.0;
  
  if(s <= 0.0)
    dsw_ds = 0.0;
  else{
    dsw_ds = pow((s/consta),(1.0/n2));
    dsw_ds = pow((dsw_ds + 1),(-1.0/n1 - 1.0));
    dsw_ds = -1.0*pow((s/consta),(1.0/n2))*dsw_ds;
    dsw_ds = dsw_ds/(n2*n1*s);
  }

  return(dsw_ds);
}


/**
   function computes partial derivative of degree of saturation with respect to t

   @param s - suction
   @retval dsw_dt - partial derivative of degree of saturation with respect to temperature

   TKr 3/04/2018
*/
double febex_granit_reten::dsw_dt(double pw)
{
  double dsw_dt;

  dsw_dt = 0.0;
  
  return(dsw_dt);
}




/**
   function computes relative permeability

   @param sw - saturation degree
   @retval krw - relative permeability for water

   TKr 3/04/2018
*/
double febex_granit_reten::get_krw(double sw)
{
  double krw;
  
  krw = pow(sw,m2);
  krw = pow((1.0 - krw),m3);
  krw = pow((1.0 - krw),m4);
  krw = krw*pow(sw,m1);
  
  return (krw);
}
