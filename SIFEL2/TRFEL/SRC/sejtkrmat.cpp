/*
    File:             sejtkrmat.cpp
    Author:           Tomas Krejci, 12/9/2008
    Purpose:          material model for saturated one-phase flow in deforming medium
    sources:          Numericke metody mechaniky 2, Z. Bittnar - J. Sejnoha, 47 - 56
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "sejtkrmat.h"
#include "globalt.h"

sejtkrmat::sejtkrmat()
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

sejtkrmat::~sejtkrmat()
{}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void sejtkrmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of components of conductivity tensor is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function creates conductivity %matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  kk = get_kww(pw);
  
  fillm(0.0,d);

  d[0][0] = kk;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  kk = get_kww(pw);
  
  fillm(0.0,d);
  
  d[0][0] = kk;   d[0][1] = 0.0;
  d[1][0] = 0.0;  d[1][1] = kk;
}

/**
   function creates conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void sejtkrmat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  kk = get_kww(pw);
  
  fillm(0.0,d);
  
  d[0][0]=kk;   d[0][1]=0.0;  d[0][2]=0.0;
  d[1][0]=0.0;  d[1][1]=kk;   d[1][2]=0.0;
  d[2][0]=0.0;  d[2][1]=0.0;  d[2][2]=kk;
}


/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void sejtkrmat::matcap (double &cc,long /*ri*/,long /*ci*/,long ipp)
{
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  cc = get_capww(pw);
}

/**
   function reads parameters
   
   @param in - input file

   12/9/2008, TKr
*/
void sejtkrmat::read(XFILE */*in*/)
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
void sejtkrmat::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
void sejtkrmat::rhs1d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
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
void sejtkrmat::rhs2d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
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
void sejtkrmat::rhs3d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
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
double sejtkrmat::get_kuw(double /*pw*/)
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
double sejtkrmat::get_kwu(double /*pw*/)
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
double sejtkrmat::get_kww(double /*pw*/)
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
double sejtkrmat::get_capuw(double /*pw*/)
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
double sejtkrmat::get_capwu(double /*pw*/)
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
double sejtkrmat::get_capww(double /*pw*/)
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
double sejtkrmat::get_fw1(double /*pw*/)
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
double sejtkrmat::get_fu1(double /*pw*/)
{
  double fu1;  
  
  fu1 = 0.0;

  return(fu1);

}

/**
   function computes all other variables at nodes
   @param compother - number of other components
   @param pw - water capillary pressure on actual node
   @param pg - gas(air) pressure on actual node
   @param ipp - first integration point on element

   @retval other - other variable

   03/03/2011, TKr
*/

double sejtkrmat::get_othervalue(long compother,double pw, long /*ipp*/)
{
  double other;
  state_eq tt;

  switch (compother){
  case 0:{//capillary pressure
    other = -pw;
      break;
  }
  case 1:{//saturation
    other = 1.0;
    break;
  }
  case 2:{//liquid water pressure
    other = pw;
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);

}



/**
     function prints names of all other variables at nodes
     @param out - output rhle
     @param compother - number of other components

     03/03/2011, TKr
*/
void sejtkrmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//capillary pressure
    fprintf (out,"Capillary pressure (Pa)");
    break;
  }
  case 1:{//saturation
    fprintf (out,"Degree of saturation ()       ");
    break;
  }
  case 2:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}



/**
   function checks if computed unknowns are physically reasonable
   @param nv  - vector of unknowns
   @param ipp - number of integration point


   28/02/2011, TKr
*/
void sejtkrmat::values_correction (vector &nv, long ipp)
{
  //  pore water pressure control
  water_pressure_check(nv[0],ipp);
}


/**
   function checks if water pressure is non-positive
   @param pw - pore water pressure
   @param ipp - number of integration point


   28/02/2011, TKr
*/
void sejtkrmat::water_pressure_check(double &pw,long ipp)
{
  if (pw >= 0.0) 
    pw = -1.0;
  
  //storing into gauss point
  Tm->ip[ipp].av[0] = pw;
}
