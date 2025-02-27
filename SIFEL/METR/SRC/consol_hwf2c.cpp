/*
    File:             consol_hwf2c.cpp
    Author:           Tomas Krejci, 29/10/2009, revised 03/03/2010
    Purpose:          material model for saturated-nonsaturated heat and water flow in a deforming porous medium
    sources:          Lewis and Schrefler pp. ??-??, 
    unknowns:         number of unknowns=2, pw = liqud water pressure, t = temperature
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
#include "consol_hwf2c.h"
#include "intpoints.h"


con_hwf2matc::con_hwf2matc()
{
  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha0 = 0.0;
  ks0 = 0.0;
  phi0 = 0.0;
  kintr0 = 0.0;
  betas0 = 0.0;
  rhos0 = 0.0;
  cps0 = 0.0;
  lambda0 = 0.0;
  rhow0 = 1000.0;
  emod = 1.3e6;
  nu = 0.25;
}

con_hwf2matc::~con_hwf2matc()
{}


/**
   function reads parameters
   
   @param in - input file

   TKr 13/04/2018
*/
void con_hwf2matc::read(XFILE *in)
{
  xfscanf (in,"%k%m","heatwaterflowmechtype",&heatwaterflowmechtype_kwdset, &model_type);
  switch (model_type){
  case lewis_and_schrefler2hw_coup:
  case lewis_and_schrefler2hw_mefel_coup:{//Lewis and Schrefler's model
    xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &alpha0, &ks0, &phi0, &kintr0, &betas0, &rhos0, &cps0, &lambda0, &emod, &nu);
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
}


/**
   function prints parameters
   
   @param out - output file

   TKr 13/04/2018
*/
void con_hwf2matc::print(FILE *out)
{
  fprintf (out,"\n %d ", int(model_type));
  switch (model_type){
  case lewis_and_schrefler2hw_coup:
  case lewis_and_schrefler2hw_mefel_coup:{//Lewis and Schrefler's model
    fprintf (out,"\n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda0,emod,nu);
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
void con_hwf2matc::matcond_u (matrix &d,long ri,long ci,long ipp)
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
  
  if (ci == 1){
    matrix s(d.m,d.m);
    
    Cmu->matstiff (s,ipp);//only for temperature
    mxm(s,d,d);
    
    destrm (s);
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
void con_hwf2matc::matcond1d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kut(pw,t,ipp);
  
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
void con_hwf2matc::matcond2d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kut(pw,t,ipp);
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::matcond2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kut(pw,t,ipp);
  
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
void con_hwf2matc::matcond3d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kut(pw,t,ipp);
  
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
void con_hwf2matc::matcap_u (matrix &d,long ri,long ci,long ipp)
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
void con_hwf2matc::matcap1d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_caput(pw,t,ipp);
  
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
void con_hwf2matc::matcap2d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_caput(pw,t,ipp);
  
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::matcap2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_caput(pw,t,ipp);
  
  
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
void con_hwf2matc::matcap3d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  
  pw = Cmu->ip[ipp].av[0];
  t = Cmu->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_caput(pw,t,ipp);
  
  
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
void con_hwf2matc::matcond_l (matrix &d,long ri,long ci,long ipp)
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
void con_hwf2matc::matcond1d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_ktu(pw,t,ipp);
  
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
void con_hwf2matc::matcond2d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_ktu(pw,t,ipp);
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::matcond2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_ktu(pw,t,ipp);
  
  fillm(0.0,d);
  
  d[0][0] = k;
  d[0][1] = k;
  d[0][2] = k;  
  d[0][3] = 0.0;  
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
void con_hwf2matc::matcond3d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,t;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_ktu(pw,t,ipp);
  
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
void con_hwf2matc::matcap_l (matrix &d,long ri,long ci,long ipp)
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
void con_hwf2matc::matcap1d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_captu(pw,t,ipp);
  
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
void con_hwf2matc::matcap2d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_captu(pw,t,ipp);
  
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::matcap2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_captu(pw,t,ipp);
  
  
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
void con_hwf2matc::matcap3d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,t,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  t = Cml->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_captu(pw,t,ipp);
  
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::rhs_u1 (matrix &d,long ri,long ci,long ipp)
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

   12/06/2018, TKr
*/
void con_hwf2matc::rhs1d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
  
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t = Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fuw1(pw,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = f*Cp->gr[0];
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hwf2matc::rhs2d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
    
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t= Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fuw1(pw,t,ipp);
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::rhs2d_ax_1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
    
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t= Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fuw1(pw,t,ipp);
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::rhs3d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
  
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t = Cmu->ip[ipp].av[1];// *scale_pc;//scaling

  f = get_fuw1(pw,t,ipp);
  
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

   12/06/2018, TKr
*/
void con_hwf2matc::rhs_u2 (matrix &d,long ri,long ci,long ipp)
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
    rhs2d_ax_2 (d,ri,ci,ipp);//2D axisymmetric
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
  
  if (ci == 1){
    matrix s(d.m,d.m);
    
    Cmu->matstiff (s,ipp);//only for temperature
    mxm(s,d,d);
    
    destrm (s);
  }
}




/**
   function creates volume right-hand side matrix of the material for 1D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hwf2matc::rhs1d2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
  
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t = Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fut2(pw,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = f;
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hwf2matc::rhs2d2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
    
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t= Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fut2(pw,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = f;
  d[1][0] = f;
  d[2][0] = 0.0;
}


/**
   function creates volume right-hand side matrix of the material for 2D axisymmetric problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hwf2matc::rhs2d_ax_2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
    
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t= Cmu->ip[ipp].av[1];// *scale_pc;//scaling
  
  f = get_fut2(pw,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = f;
  d[1][0] = f;
  d[2][0] = f;
  d[3][0] = 0.0;
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hwf2matc::rhs3d2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw,t;
  
  pw = Cmu->ip[ipp].av[0];// *scale_pc;//scaling
  t = Cmu->ip[ipp].av[1];// *scale_pc;//scaling

  f = get_fut2(pw,t,ipp);
  
  fillm(0.0,d);
  d[0][0]=f;  
  d[1][0]=f;   
  d[2][0]=f;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}





/**
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param t - temperature

   @retval sw - degree of saturation

   TKr 13/04/2018
*/
double con_hwf2matc::get_sw(double pw, double /*t*/,long ipp)
{
  double sw;
    
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    sw = lewis_ret.sw(pw);
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//saturation degree ind its derivative are obtained from mefel;
    sw = Tm->givenontransq(saturation_deg, ipp); //actual saturation degree from models included in MEFEL
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
double con_hwf2matc::get_porosity(long ipp)
{
  double por;

  switch (model_type){
  case lewis_and_schrefler:{//Lewis and Schrefler's book
    por = phi0; //constant value
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's book
    por = Tm->givenontransq(porosity, ipp);// from mefel
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
   function computes Biot's constant
   @param pw - pore water pressure
   @param t - temperature
   
   @retval alpha - Biot's constant
*/
double con_hwf2matc::get_alpha(double /*pw*/, double /*t*/,long /*ipp*/)
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
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param t - temperature

   @retval kuw - conductivity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_kuw(double pw, double t, long ipp)
{
  double sw,kuw,alpha;
  
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,t,ipp);
    alpha = get_alpha(pw,t,ipp);
    kuw = -sw*alpha;//p. 362, but non-incremental eq. for displ.
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    sw = get_sw(pw,t,ipp);
    alpha = get_alpha(pw,t,ipp);
    kuw = -sw*alpha;//p. 362, but non-incremental eq. for displ.
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
   @param t - temperature

   @retval kut - conductivity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_kut(double /*pw*/, double /*t*/, long ipp)
{
  double kut,betas;
  
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    betas = get_betas(ipp);  
    kut = -1.0*betas/3.0;//p. 348
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    betas = get_betas(ipp);  
    kut = -1.0*betas/3.0;//p. 348
    kut = 0.0;// this part is in MEFEL already included
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   
  return(kut);
  
}



/**
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param t - temperature

   @retval ktu - conductivity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_ktu(double /*pw*/, double /*t*/, long /*ipp*/)
{
  double ktu;

  ktu = 0.0;//p. 362

  return(ktu);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param t - temperature

   @retval kwu - conductivity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_kwu(double /*pw*/, double /*t*/, long /*ipp*/)
{
  double kwu;

  kwu = 0.0;;//p. 362

  return(kwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param t - temperature

   @retval capuw - capacity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_capuw(double /*pw*/, double /*t*/, long /*ipp*/)
{
  double capuw;
  
  capuw = 0.0;//p. 97 and p. 362, but non-incremental eq. for displ.
    
  return(capuw);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param t - temperature

   @retval capwu - capacity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_capwu(double pw, double t, long ipp)
{
  double sw,capwu,alpha,n,dsr_depsv;
  
  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,t,ipp);
    alpha = get_alpha(pw,t,ipp);
    capwu = sw*alpha;//p. 95
    break;
  }    
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's book
    sw = get_sw(pw,t,ipp);
    alpha = get_alpha(pw,t,ipp);
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
   @param t - temperature

   @retval capug - capacity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_caput(double /*pw*/, double /*t*/, long /*ipp*/)
{
  double caput;
  
  caput = 0.0;//p. 97
  
  return(caput);

}



/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param t - temperature

   @retval captu - capacity coefficient

   TKr 13/04/2018
*/
double con_hwf2matc::get_captu(double /*pw*/, double /*t*/, long /*ipp*/)
{
  double captu;
 
  captu = 0.0;//no evaporation effect

  return(captu);

}



/**
   function returns coefficient for righ-hand side of the general material 
   @param pw - water pressure
   @param t - temperature

   @retval fuw1 - first part for right-hand side for balance equation

   TKr 13/04/2018
*/
double con_hwf2matc::get_fuw1(double pw, double t, long ipp)
{
  double sw,n,fuw1;  

  sw = get_sw(pw,t,ipp);
  n = get_porosity(ipp);

  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    fuw1 = rhos0*(n-1) + sw*n*rhow0;//p. 95
    //fuw1 = sw*n*rhow;//debug
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's model
    fuw1 = rhos0*(n-1) + sw*n*rhow0;//p. 95
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
  
  return(fuw1);
}

/**
   function returns coefficient for righ-hand side of the general material 
   @param pw - water pressure
   @param t - temperature

   @retval fut1 - first part for right-hand side for balance equation

   TKr 12/06/2018
*/
double con_hwf2matc::get_fut2(double /*pw*/, double /*t*/, long ipp)
{
  double betas,fut2;  

  betas = get_betas(ipp);

  switch (model_type){
  case lewis_and_schrefler2_coup:{//Lewis and Schrefler's book
    fut2 = -1.0*betas/3.0;
    break;
  }
  case lewis_and_schrefler2_mefel_coup:{//Lewis and Schrefler's model
    fut2 = -1.0*betas/3.0;
    fut2 = 0.0;// this part is in MEFEL already included
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
  
  return(fut2);
}



/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
*/
double con_hwf2matc::get_betas(long /*ipp*/)
{
  double betas;

  betas = betas0;

  return(betas);
}



/**
   function returns Young's modulus of solid [Pa]

   @retval e -  Young's modulus of solid [Pa]
*/
double con_hwf2matc::get_e(long /*ipp*/)
{
  return(emod); //provisionally
}



/**
   function returns Poisson's ratio [-]

   @retval nu -  Poisson's ratio [-]
*/
double con_hwf2matc::get_nu(long /*ipp*/)
{
  return(nu); //provisionally
}

