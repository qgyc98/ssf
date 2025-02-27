/*
    File:             consol_awf2c.cpp
    Author:           Tomas Krejci, 29/10/2009, revised 03/03/2010
    Purpose:          material model for saturated-nonsaturated air and water flow in a deforming porous medium
    sources:          Lewis and Schrefler pp. 93-97, material parameters are set for benchmark on page n. 168
    unknowns:         number of unknowns=2, pw = liqud water pressure, pg = gas(air) pressure
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "twomediac.h"
#include "coupmatu.h"
#include "coupmatl.h"
#include "globalc.h"
#include "global.h"
#include "globalt.h"
#include "consol_awf2c.h"
#include "intpoints.h"


con_awf2matc::con_awf2matc()
{
  compress = 0; //compressible grains: 0=no; 1=yes
  //
  alpha = 1.0;//Biot's constant [-] alpha = 1 - kt/ks
  //
  ks = 2.167e6;//bulk modulus of solid phase (grains) [Pa]
  //
  kw = 2.0e9;//bulk modulus of water [Pa]
  //
  phi0 = 0.297;//initial porosity [-]
  //
  kintr = 4.5e-13;//intrinsic permeability [m^2]
  //
  rhow = 1000.0;//water density [kg/m^3]
  //
  muw0 = 1.0e-3;//water viscosity [Pa.s]
  //
  mug0 = 1.8e-5;//air viscosity [Pa.s]
  //
  p_atm = 101325.0; //atmospheric pressure [Pa]
  //
  rhos0 = 2500.0; //solid density
  //
  mefel_units = 1.0;//basic units for pressures = Pa (Pascals)
}

con_awf2matc::~con_awf2matc()
{}


/**
   function reads parameters
   
   @param in - input file

   29/10/2009, TKr
*/
void con_awf2matc::read(XFILE *in)
{
  xfscanf (in,"%k%m","airwaterflowmechtype",&airwaterflowmechtype_kwdset, &model_type);
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    xfscanf (in,"%le %le %le %le %le %le %le %le %le", &alpha, &ks, &phi0, &kw, &rhow, &muw0, &mug0, &kintr, &rhos0);
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    xfscanf (in,"%le %le %le %le %le %le %le %le %le %le", &mefel_units, &alpha, &ks, &phi0, &kw, &rhow, &muw0, &mug0, &kintr, &rhos0);
    break;
  }
  case van_genuchten2_coup:{//partially saturated medium =Van Genuchten model
    xfscanf (in,"%le %le %le %le %le %le %le %le %le %le", &mefel_units, &alpha, &ks, &phi0, &kw, &rhow, &muw0, &mug0, &kintr, &rhos0);
    van_genuchten_ret.read(in);
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }  
  
  // incompressible grains:
  if(compress == 0)
    alpha = 1.0;
  
}


/**
   function reads parameters
   
   @param out - output file

   29/10/2009, TKr
*/
void con_awf2matc::print(FILE *out)
{
  fprintf (out,"\n %d ", int(model_type));
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    fprintf (out,"\n %le %le %le %le %le %le %le %le %le \n",alpha,ks,phi0,kw,rhow,muw0,mug0,kintr,rhos0);
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le \n",mefel_units,alpha,ks,phi0,kw,rhow,muw0,mug0,kintr,rhos0);
    break;
  }
  case van_genuchten2_coup:{//partially saturated medium =Van Genuchten model
    fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le \n",mefel_units,alpha,ks,phi0,kw,rhow,muw0,mug0,kintr,rhos0);
    van_genuchten_ret.print(out);
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function computes conductivity matrix of the material
   in the required integration point for upper block
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond_u (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;

  switch (m){
  case 1:{
    matcond1d_u (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcond2d_u (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcond2d_ax_u (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcond3d_u (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates conductivity matrix of the material for 1D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond1d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond2d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,ipp);
  
  fillm(0.0,d);

  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = 0.0;  
}

/**
   function creates conductivity matrix of the material for 2D axisymmetric problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   11/06/2018, TKr
*/
void con_awf2matc::matcond2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,ipp);
  
  fillm(0.0,d);

  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = k;  
  d[3][0] = 0.0;  
}


/**
   function creates conductivity matrix of the material for 3D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond3d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[1][0]=k;   
  d[2][0]=k;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}


/**
   function computes capacity matrix of the material
   in the required integration point for upper block
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   09/03/2011, TKr
*/
void con_awf2matc::matcap_u (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.m;

  switch (m){
  case 1:{
    matcap1d_u (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcap2d_u (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcap2d_ax_u (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcap3d_u (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates capacity matrix of the material for 1D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcap1d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,ipp);
  
  d[0][0] = c;
}


/**
   function creates capacity matrix of the material for 2D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcap2d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,ipp);
  
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = 0.0;
}



/**
   function creates capacity matrix of the material for 2D axisymmetric problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   11/06/2018, TKr
*/
void con_awf2matc::matcap2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,ipp);
  
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = c;
  d[3][0] = 0.0;
}


/**
   function creates capacity matrix of the material for 3D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcap3d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,ipp);
  
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = c;
  d[3][0] = 0.0;
  d[4][0] = 0.0;
  d[5][0] = 0.0;
}


/**
   function computes conductivity matrix of the material
   in the required integration point for lower block
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond_l (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.n;

  switch (m){
  case 1:{
    matcond1d_l (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcond2d_l (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcond2d_ax_l (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcond3d_l (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates conductivity matrix of the material for 1D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond1d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kgu(pw,pg,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond2d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kgu(pw,pg,ipp);
  
  fillm(0.0,d);
  
  d[0][0] = k;
  d[0][1] = k;
  d[0][2] = 0.0;  
}

/**
   function creates conductivity matrix of the material for 2D axisymmetric problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   11/06/2018, TKr
*/
void con_awf2matc::matcond2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kgu(pw,pg,ipp);
  
  fillm(0.0,d);
  
  d[0][0] = k;
  d[0][1] = k;
  d[0][3] = k;  
  d[0][2] = 0.0;  
}

/**
   function creates conductivity matrix of the material for 3D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcond3d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kgu(pw,pg,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[0][1]=k;   
  d[0][2]=k;
  d[0][3]=0.0;   
  d[0][4]=0.0;  
  d[0][5]=0.0;
}


/**
   function computes capacity matrix of the material
   in the required integration point for lower block
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   09/03/2011, TKr
*/
void con_awf2matc::matcap_l (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.n;

  switch (m){
  case 1:{
    matcap1d_l (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    matcap2d_l (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    matcap2d_ax_l (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }   
  case 6:{
    matcap3d_l (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates capacity matrix of the material for 1D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcap1d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capgu(pw,pg,ipp);
  
  d[0][0] = c;
}


/**
   function creates capacity matrix of the material for 2D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcap2d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capgu(pw,pg,ipp);
  
  
  fillm(0.0,d);

  
  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = 0.0;  
}


/**
   function creates capacity matrix of the material for 2D axisymmetric problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   11/06/2018, TKr
*/
void con_awf2matc::matcap2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capgu(pw,pg,ipp);
  
  
  fillm(0.0,d);
  
  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = c;  
  d[0][3] = 0.0;  
}



/**
   function creates capacity matrix of the material for 3D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   09/03/2011, TKr
*/
void con_awf2matc::matcap3d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capgu(pw,pg,ipp);
  
  
  fillm(0.0,d);
  
  d[0][0]=c;  
  d[0][1]=c;   
  d[0][2]=c;
  d[0][3]=0.0;   
  d[0][4]=0.0;  
  d[0][5]=0.0;
}



/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs_u1 (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;
  
  switch (m){
  case 1:{
    rhs1d1 (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    rhs2d1 (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    rhs2d_ax_1 (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }
  case 6:{
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
void con_awf2matc::rhs1d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,pg;
  
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fu1(pw,pg,ipp);
  
  fillm(0.0,d);
  d[0][0] = f*Cp->gr[0];
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs2d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,pg;
    
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fu1(pw,pg,ipp);
  
  fillm(0.0,d);
  d[0][0] = f*Cp->gr[0];
  d[1][0] = f*Cp->gr[1];
}

/**
   function creates volume right-hand side matrix of the material for 2D axisymmetric problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs2d_ax_1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,pg;
    
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fu1(pw,pg,ipp);
  
  fillm(0.0,d);
  d[0][0] = f*Cp->gr[0];
  d[1][0] = f*Cp->gr[1];
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs3d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,pg;
  
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];// *scale_pc;//scaling

  f = get_fu1(pw,pg,ipp);
  
  fillm(0.0,d);
  d[0][0] = f*Cp->gr[0];
  d[1][0] = f*Cp->gr[1];
  d[2][0] = f*Cp->gr[2];
}


/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs_u2 (matrix &d,long ri,long ci,long ipp)
{
  long m;
  m = d.m;
  
  switch (m){
  case 1:{
    rhs1d2 (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    rhs2d2 (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    rhs2d_ax_2 (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }
  case 6:{
    rhs3d2 (d,ri,ci,ipp);//3D
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
void con_awf2matc::rhs1d2 (matrix &/*d*/,long /*ri*/,long /*ci*/,long /*ipp*/)
{
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs2d2 (matrix &/*d*/,long /*ri*/,long /*ci*/,long /*ipp*/)
{
}

/**
   function creates volume right-hand side matrix of the material for 2D axisymmetric problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs2d_ax_2 (matrix &/*d*/,long /*ri*/,long /*ci*/,long /*ipp*/)
{
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2matc::rhs3d2 (matrix &/*d*/,long /*ri*/,long /*ci*/,long /*ipp*/)
{
}


/**
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param pg - air pressure

   @retval sw - degree of saturation

   29/10/2009, TKr
*/
double con_awf2matc::get_sw(double pw, double /*pg*/,long ipp)
{
  double sw;
    
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    sw = lewis_ret.sw(pw);
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    sw = Tm->givenontransq(saturation_deg, ipp); //actual saturation degree from models included in TRFEL
    break;
  }
  case van_genuchten2_coup:{//partially saturated medium = Van Genuchten model
    pw = pw/mefel_units;
    //sw = van_genuchten_ret.sw(pw,t0);
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }  
  
  return(sw);
}

/**
   function returns porosity

   @retval phi - porosity

   13/04/2018, TKr
*/
double con_awf2matc::get_porosity(long ipp)
{
  double por;

  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    por = phi0; //constant value
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book and moisture storage functions is included in MEFEL
    por = Tm->givenontransq(porosity, ipp);//actual porosity from models included in TRFEL
    break;
  }
  case van_genuchten2_coup:{//partially saturated medium = Van Genuchten model
    por = phi0;
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 

  return por;
}


/**
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kuw - conductivity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_kuw(double pw, double pg, long ipp)
{
  double sw,kuw;
  
  switch (model_type){
  case van_genuchten2_coup://partially saturated medium = Van Genuchten model
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,pg,ipp);
    kuw = -sw*alpha;//p. 95
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    sw = get_sw(pw,pg,ipp);
    kuw = -sw*alpha;//p. 95
    kuw = 0.0;// this part is in MEFEL already included
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   
  
  return(kuw);
}



/**
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kug - conductivity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_kug(double pw, double pg, long ipp)
{
  double sw,sg,kug;
  
  switch (model_type){
  case van_genuchten2_coup://partially saturated medium = Van Genuchten model
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,pg,ipp);
    sg = 1.0 - sw;
    kug = -sg*alpha;//p. 95
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    sw = get_sw(pw,pg,ipp);
    sg = 1.0 - sw;
    kug = -sg*alpha;//p. 95
    kug = 0.0;// this part is in MEFEL already included
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   
  
  //kug = 0.0;//debug??!!
  return(kug);
  
}



/**
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kgu - conductivity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_kgu(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  double kgu;

  kgu = 0.0;

  return(kgu);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kwu - conductivity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_kwu(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  double kwu;

  kwu = 0.0;

  return(kwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capuw - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_capuw(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  double capuw;
  
  capuw = 0.0;//p. 97
  
  return(capuw);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capwu - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_capwu(double pw, double pg, long ipp)
{
  double sw,capwu,n,dsr_depsv;
  
  switch (model_type){
  case van_genuchten2_coup://partially saturated medium = Van Genuchten model
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,pg,ipp);
    capwu = sw*alpha;//p. 95
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    sw = get_sw(pw,pg,ipp);
    capwu = sw*alpha;//p. 95

    n = get_porosity(ipp);
    dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
    capwu = capwu + n*dsr_depsv;//influence of volumetric strain on saturation degree
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   
  
  return(capwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capug - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_capug(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  double capug;
  
  capug = 0.0;//p. 97
  
  return(capug);

}



/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capgu - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2matc::get_capgu(double pw, double pg, long ipp)
{
  double sg,sw,capgu,n,dsr_depsv;
 
  switch (model_type){
  case van_genuchten2_coup://partially saturated medium = Van Genuchten model
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,pg,ipp);
    sg = 1.0-sw;
    capgu = sg*alpha;//p. 95
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    sw = get_sw(pw,pg,ipp);
    sg = 1.0-sw;
    capgu = sg*alpha;//p. 95

    n = get_porosity(ipp);
    dsr_depsv = 0.0;//Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
    capgu = capgu - n*dsr_depsv;//influence of volumetric strain on saturation degree
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   

  return(capgu);

}



/**
   function returns coefficient for righ-hand side of the general material - influence of gravity accel.
   @param pw - water pressure
   @param pg - air pressure

   @retval fu1 - first part for right-hand side for balance equation

   29/10/2009, TKr
*/
double con_awf2matc::get_fu1(double pw, double pg, long ipp)
{
  double rhog,sw,sg,n,fu1;  

  sw = get_sw(pw,pg,ipp);
  sg = 1.0 - sw;
  n = get_porosity(ipp);

  rhog = 1.25;//temporarilly

  switch (model_type){
  case lewis_and_schrefler_coup:{//Lewis and Schrefler's book
    fu1 = rhos0*(n-1) + sw*n*rhow + sg*n*rhog;//p. 95
    //fu1 = sw*n*rhow + sg*n*rhog;//debug
    break;
  }
  case lewis_and_schrefler_mefel_coup:{//Lewis and Schrefler's model
    fu1 = rhos0*(n-1) + sw*n*rhow + sg*n*rhog;//p. 95
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
  return(fu1);
}


/**
   function returns coefficient for righ-hand side of the general material - influence of initial values
   @param pw - water pressure
   @param pg - air pressure

   @retval fu2 - first part for right-hand side for balance equation

   12/06/2018, TKr
*/
double con_awf2matc::get_fu2(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  return(0.0);
}
