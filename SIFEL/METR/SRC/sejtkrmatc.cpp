/*
    File:             sejtkrmatc.cpp
    Author:           Tomas Krejci, 12/9/2008
    Purpose:          material model for saturated one-phase flow in deforming medium
    sources:          Numericke metody mechaniky 2, Z. Bittnar - J. Sejnoha, 47 - 56
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "sejtkrmatc.h"
#include "globalt.h"

sejtkrmatc::sejtkrmatc()
{
  //parametry pro soudrznou zeminu (hlinito-piscita?)
  emod = 10.0e6;//Younguv modul pruznosti zeminy [Pa]
  //
  nu = 0.4;//Poissonova kostanta [-]
  //
  alpha = 0.5;//0.9998;//Biot's constant [-] alpha = 1 - ks/kz
  //
  kz = 18.0e12;//modul objemove pruznosti zakladni latky [Pa]
  //
  ks = 3.6e6; //modul objemove pruznosti porezniho skeletu [Pa]
  //
  kk = 2.0e9;//modul objemove pruznosti kapaliny-vody [Pa]
  //
  phi0  = 0.5;//pocatecni porozita [-] (znaceno v literature n)
  //
  k     = 1.13e-10;// soucinitel filtrace [m/s]
  //
  rhok = 998.0;// hustota kapaliny [kg/m^3]
  //
  g = 9.806;// tihove zrychleni [m.s^-2]
  //
  lambda = 0.0;// soucinitel
  //
  c = 4.1667e-8;// soucinitel konsolidace [m^2/s]
}

sejtkrmatc::~sejtkrmatc()
{}


/**
   function reads parameters
   
   @param in - input file

   12/9/2008, TKr
*/
void sejtkrmatc::read(XFILE */*in*/)
{
  //xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &emod, &nu, &kz, &ks, &kk, &phi0, &k, &rhok, &g, &c);
  alpha = 1.0 - ks/kz;
}


/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmatc::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates volume right-hand side matrix of the material for 1D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmatc::rhs1d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  
  f = get_fw1(pw);
  
  fillm(0.0,d);
  d[0][0] = f*Tp->gr[0];
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmatc::rhs2d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  
  f = get_fw1(pw);
  
  fillm(0.0,d);
  d[0][0] = f*Tp->gr[0];
  d[1][0] = f*Tp->gr[1];
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmatc::rhs3d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;

  pw = Tm->ip[ipp].av[0];// *scale_pc;//scaling
  
  f = get_fw1(pw);
  
  fillm(0.0,d);
  d[0][0] = f*Tp->gr[0];
  d[1][0] = f*Tp->gr[1];
  d[2][0] = f*Tp->gr[2];
}



/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kuw - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kuw(double /*pw*/)
{
  double kuw;
  
  kuw = -alpha;

  return(kuw);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kwu - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kwu(double /*pw*/)
{
  double kwu;

  kwu = 0.0;

  return(kwu);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kww - conductivity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_kww(double /*pw*/)
{
  double kww;
  
  kww = c/ks;

  return(kww);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capuw - capacity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_capuw(double /*pw*/)
{
  double capuw;
  
  capuw = 0.0;

  return(capuw);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capwu - capacity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_capwu(double /*pw*/)
{
  double capwu;
  
  capwu = alpha;

  return(capwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capww - capacity coefficient

   12/9/2008, TKr
*/
double sejtkrmatc::get_capww(double /*pw*/)
{
  double capww;  
  
  lambda = phi0/kk + (1-phi0)/kz - (1 - alpha)/kz;

  capww = lambda;

  return(capww);

}


/**
   function returns ... 
   @param pw - water pressure

   @retval fw - first part for right-hand side for continutiy equation

   12/9/2008, TKr
*/
double sejtkrmatc::get_fw1(double /*pw*/)
{
  double fw1;
  
  fw1 = c/ks;

  return(fw1);

}



/**
   function returns ... 
   @param pw - water pressure

   @retval fu1 - first part for right-hand side for balance equation

   12/9/2008, TKr
*/
double sejtkrmatc::get_fu1(double /*pw*/)
{
  double fu1;  
  
  fu1 = 0.0;

  return(fu1);

}
