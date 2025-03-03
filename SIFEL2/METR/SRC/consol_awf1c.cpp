/*
    File:            consol_awf1c.cpp
    Author:          Tomas Krejci, 12/09/2008, revised 23/08/2012, and revised 13/11/2023
    Purpose:         material model for saturated-nonsaturated one-phase flow (water) in deforming medium plus vapour diffusion
    sources:         Lewis and Schrefler's book
    unknowns:        TRFEL: number of unknowns=1, pw = liqud water pressure; MEFEL: according to problem_type
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constrel.h"
#include "onemediumc.h"
#include "coupmatu.h"
#include "coupmatl.h"
#include "globalc.h"
#include "global.h"
#include "globalt.h"
#include "consol_awf1c.h"
#include "intpoints.h"


con_awf1matc::con_awf1matc()
{
  alpha = 1.0;//Biot's constant [-] alpha = 1 - kt/ks
  //
  vol_strain_effectc = 0; //volumetric strain rate influence: 0=no; 1=yes
  //
  wrc_vol_strain_effectc = 0; //volumetric strain rate influence on water retention curve: 0=no; 1=yes
  // 
  pore_press_effectc = 0; //pore pressure effect acording to effective stress concept:  0=no; 1=yes
  //
  rhos0 = 2000.0;//initial solid grain density [kg/m^3]
  //
  phi0 = 0.297;//initial porosity [-]
  //
  rhow0 = 1000.0;//reference water density [kg/m^3]
  //
  t0 = 293.15; //reference temperature

  //STATE VARIABLES
  mw = 18.01528; //molar mass of water kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  c8 = -5.8002206e+03;
  c9 = 1.3914993;
  c10 =-4.8640239e-02;
  c11 = 4.1764768e-05;
  c12 = -1.4452093e-08;
  c13 = 6.5459673;
}

con_awf1matc::~con_awf1matc()
{}


/**
   function reads parameters
   
   @param in - input file

   13/11/2023, TKr
*/
void con_awf1matc::read(XFILE *in)
{ 
  xfscanf (in,"%k%m","waterflowmechtype",&waterflowmechtype_kwdset, &model_type); //water flow model type
  xfscanf (in,"%lf %d  %d %d %lf %lf %lf %d %d", &alpha, &vol_strain_effectc, &pore_press_effectc, &wrc_vol_strain_effectc, &phi0, &rhos0, &t0, &sr_type, &xi_type);
  
  switch (model_type){
  case lewis_and_schrefler_coup:
  case lewis_and_schrefler_mefel_coup: {//Lewis and Schrefler's model  
    
    //retention curve type:
    switch (sr_type){
    case baroghel_sr:{//Baroghle-Bouny approach
      baroghel_ret.read(in);
      break;
    }
    case bazant_sr:{//Bazant approach
      break;
    }
    case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
      lewis_ret.read(in);
      break;
    }
    case gardner_exponential_sr:{//dependent on saturation degree:
      gardner_ret.read(in);
      break;
    }
    case potts_log_linear_sr:{//exponential
      potts_ret.read(in);
      break;
    }
    case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
      van_genuchten_ret.read(in);
      break;
    }
    case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
      //van_genuchten_ret.read2(in);
      break;
    }
    case mefel_sr:{//from MEFEL
      break;
    }
    case table_sr:{//from table
      //reading of retention curve:
      data.read (in);
      break;
    }
    case masin_sr:{//extended formulation from Brooks and Correy according to Masin
      masin_ret.read(in);
      break;
    }
    case febex_granit_sr:{//FEBEX granit
      febex_granit_ret.read(in);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    switch (xi_type){
    case biot_xi:
      break;
    case biot_reduced_xi:{//coefficient for effective stress factor reduction
      xfscanf (in,"%le  %le", &gamma,&lambda0);
      break;
    }
    case biot_masin_xi:{
      xfscanf (in,"%le  %le", &gamma,&s_entry);
      break;
    }
    case masin_xi:
      break;
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 

    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}


/**
   function prints parameters
   
   @param out - output file

   13/11/2023, TKr
*/
void con_awf1matc::print(FILE *out)
{
}


/**
   function computes conductivity matrix of the material
   in the required integration point for upper block
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   01/05/2018, TKr
*/
void con_awf1matc::matcond_u (matrix &d,long ri,long ci,long ipp)
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

   01/05/2018, TKr
*/
void con_awf1matc::matcond1d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   01/05/2018, TKr
*/
void con_awf1matc::matcond2d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,ipp);
  
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

   10/05/2018, TKr
*/
void con_awf1matc::matcond2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,ipp);
  
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

   01/05/2018, TKr
*/
void con_awf1matc::matcond3d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,ipp);
  
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
   
   01/05/2018, TKr
*/
void con_awf1matc::matcap_u (matrix &d,long ri,long ci,long ipp)
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

   01/05/2018, TKr
*/
void con_awf1matc::matcap1d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,ipp);
  
  d[0][0] = c;
}


/**
   function creates capacity matrix of the material for 2D problems
   for upper block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   01/05/2018, TKr
*/
void con_awf1matc::matcap2d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,ipp);
  
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

   10/05/2018, TKr
*/
void con_awf1matc::matcap2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  c = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,ipp);

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

   01/05/2018, TKr
*/
void con_awf1matc::matcap3d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  
  pw = Cmu->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,ipp);
  
  
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

   01/05/2018, TKr
*/
void con_awf1matc::matcond_l (matrix &d,long ri,long ci,long ipp)
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

   01/05/2018, TKr
*/
void con_awf1matc::matcond1d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   01/05/2018, TKr
*/
void con_awf1matc::matcond2d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];

  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,ipp);
  
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

   10/05/2018, TKr
*/
void con_awf1matc::matcond2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];

  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,ipp);
  
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

   01/05/2018, TKr
*/
void con_awf1matc::matcond3d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,ipp);
  
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
   
   01/05/2018, TKr
*/
void con_awf1matc::matcap_l (matrix &d,long ri,long ci,long ipp)
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

   01/05/2018, TKr
*/
void con_awf1matc::matcap1d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,ipp);
  
  d[0][0] = c;
}


/**
   function creates capacity matrix of the material for 2D problems
   for lower block

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   01/05/2018, TKr
*/
void con_awf1matc::matcap2d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,ipp);
  
  
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

   10/05/2018, TKr
*/
void con_awf1matc::matcap2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,ipp);
  
  
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

   01/05/2018, TKr
*/
void con_awf1matc::matcap3d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,c;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,ipp);
  
  
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
void con_awf1matc::rhs_u1 (matrix &d,long ri,long ci,long ipp)
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
void con_awf1matc::rhs1d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Cmu->ip[ipp].av[0];
  
  f = get_fu1(pw,ipp);
  
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
void con_awf1matc::rhs2d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
    
  pw = Cmu->ip[ipp].av[0];
  
  f = get_fu1(pw,ipp);
  
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
void con_awf1matc::rhs2d_ax_1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
    
  pw = Cmu->ip[ipp].av[0];
  
  f = get_fu1(pw,ipp);
  
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
void con_awf1matc::rhs3d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Cmu->ip[ipp].av[0];

  f = get_fu1(pw,ipp);
  
  fillm(0.0,d);
  d[0][0] = f*Cp->gr[0];
  d[1][0] = f*Cp->gr[1];
  d[2][0] = f*Cp->gr[2];
}





/**
   function computes degree of saturation(water retention curve)
   @param pw - water pressure
   @param ipp - number of integration point

   @retval sw - degree of saturation

   13/11/2021, TKr
*/
double con_awf1matc::get_sw(double pw, long ipp)
{
  double sw;
  sw = 0.0;
  
  switch (sr_type){
  case bazant_sr:{//Bazant
    sw = bazant_ret.sat(-pw,t0);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    sw = lewis_ret.sw(pw);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    sw = Tm->givenontransq(saturation_deg, ipp); //actual saturation degree
    break;
  }
  case table_sr:{//saturation degree and its derivative are obtained from table;
    //actual saturation degree    
    if (data.tfunc == tab)
      {
	if ((pw < data.tabf->x[0]) || (pw > data.tabf->x[data.tabf->asize-1]))
	  {
	    print_err("required value %le is out of table range <%le;%le> on ip=%ld\n", 
		      __FILE__, __LINE__, __func__, pw, data.tabf->x[0], data.tabf->x[data.tabf->asize-1], ipp);
	    abort();
	  }
      }
    sw = data.getval (pw);
    
    break;
  }
    
  /* case gardner_exponential_sr:{//partially saturated medium = Exponential model, Gardner 
     sw = gardner_ret.sw(pw);
     break;
     }
     case potts_log_linear_sr:{//partially saturated medium = Log linear model, Potts
     sw = potts_ret.sw(pw);
     break;
     }
     case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
     //sw = van_genuchten_ret.sw2(pw);
     break;
     }
  */
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    sw = van_genuchten_ret.sw(pw,t0);
    break;
  }
  case masin_sr:{//extended formulation from Brooks and Correy according to Masin
    double e=0.0,dpw=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }

    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    sw = masin_ret.sw(-pw,-dpw,e,t0);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    sw = febex_granit_ret.sw(-pw);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(sw);
}



/**
   function computes effective stress factor xi
   @param pw - water pressure
   @param ipp - number of integration point

   @retval xi - factor xi

   13/11/2023, TKr
*/
double con_awf1matc::get_xi(double pw, long ipp)
{
  double suc=0.0,sr=0.0,xi=0.0;
  
  switch (xi_type){
  case biot_xi:
    xi = get_sw(pw,ipp);
    break;
  case biot_reduced_xi:
    //xi = gamma*get_sw(pw,ipp);
    sr = get_sw(pw,ipp);
    //xi = pow(sr,(gamma/lambda0));
    xi = pow(sr,gamma);
    xi = (1-gamma)*xi;
    break;
  case biot_masin_xi:{//according to masin for testing
    sr = get_sw(pw,ipp);
    suc = -pw;
    if (suc>=s_entry)
      xi = pow((s_entry/suc),gamma);
    else
      xi = 1.0;
    if (suc>=s_entry)
      xi = (1-gamma)*xi;
    else
      xi = 1.0;
    break;
  }
  case masin_xi:{
    double e=0.0,dpw=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }
    
    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    xi = masin_ret.psi(-pw,-dpw,e,t0);//positive value of suction
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(xi);
}

/**
   function returns porosity

   @retval phi - porosity

   12/9/2008, TKr
*/
double con_awf1matc::get_phi(double /*pw*/,long ipp)
{
  double phi;

  switch (model_type){
  case lewis_and_schrefler_coup:{//Lewis and Schrefler's model
    phi = phi0;
    break;
  }
  case lewis_and_schrefler_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    phi = Tm->givenontransq(porosity, ipp); //actual porosity from models included in TRFEL
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
  return(phi);
}



/**
   function returns volume density of soil skeleton

   @retval rhos - volume density of concrete skeleton

   12/9/2008, TKr
*/
double con_awf1matc::get_rhos(double /*pw*/,long /*ipp*/)
{
  return(rhos0);
}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kuw - conductivity coefficient

   12/9/2008, TKr
*/
double con_awf1matc::get_kuw(double pw,long ipp)
{
  double chi,kuw;
  
  switch (model_type){
  case lewis_and_schrefler_coup:{//Lewis and Schrefler's model
    chi = get_xi(pw,ipp);
    kuw =  -chi;//effective stress parameter
    break;
  }
  case lewis_and_schrefler_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    kuw = 0.0;// pore pressure effect included in MEFEL
    if(pore_press_effectc == 1){
      // this below in not completed yet
      chi = 0.0;
      //chi = Tm->givenontransq(eff_stress_param, ipp);      //actual effective stress parameter from MEFEL
      kuw =  -chi;// pore pressure effect included in METR
    }

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
   function creates conductivity coefficient of the general material
   @param pw - water pressure

   @retval kwu - conductivity coefficient

   12/9/2008, TKr
*/
double con_awf1matc::get_kwu(double /*pw*/,long /*ipp*/)
{
  double kwu;

  kwu = 0.0;

  return(kwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capuw - capacity coefficient

   12/9/2008, TKr
*/
double con_awf1matc::get_capuw(double /*pw*/,long /*ipp*/)
{
  double capuw;
  
  capuw = 0.0;

  return(capuw);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure

   @retval capwu - capacity coefficient

   13/11/2023, TKr
*/
double con_awf1matc::get_capwu(double pw,long ipp)
{
  double capwu;
  double sw,sg,rhogw;
  double n,dsr_depsv;
  
  sw = get_sw(pw,ipp);
  rhogw = get_rhogw(pw);
  sg = 1.0 - sw;

  capwu = 0.0;  

  switch (model_type){
  case lewis_and_schrefler_coup:{//Lewis and Schrefler's model 
    capwu = (sg*rhogw + sw*rhow0)*alpha;//p. 394 or p. 95 multiplied by rhogw and plus vapour effect
    break;
  }
  case lewis_and_schrefler_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    capwu = 0.0;// volumetric strain effect included in TRFEL
    if(vol_strain_effectc == 1){
      sw = get_sw(pw,ipp);
      capwu = (sg*rhogw + sw*rhow0)*alpha;//volumetric strain effect included in METR
      if(wrc_vol_strain_effectc == 1){
	// this must be rewritten for METR elements
	n = get_phi(pw,ipp);
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	capwu = capwu + n*(rhow0 - rhogw)*dsr_depsv;//volumetric strain effect on retention curve
      }
    }
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
   function returns coefficient for righ-hand side of the general material 
   @param pw - water pressure
   @param pg - air pressure

   @retval fu1 - first part for right-hand side for balance equation

   13/11/2023, TKr
*/
double con_awf1matc::get_fu1(double pw,long ipp)
{
  double sw,n,fu1;  
  
  switch (model_type){
  case lewis_and_schrefler_coup:{//Lewis and Schrefler's model 
    sw = get_sw(pw,ipp);
    n = get_phi(pw,ipp);
    fu1 = rhos0*(1-n) + sw*n*rhow0;//p. 88
    break;
  }
  case lewis_and_schrefler_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    sw = get_sw(pw,ipp);
    n = get_phi(pw,ipp);
    fu1 = rhos0*(1-n) + sw*n*rhow0;//p. 88
    fu1 = 0.0;// temporarily no gravity acceleration effect is assumed
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
   function computes mass concentration of water vapour air in gas phase
   @param pw - pore water pressure

   @retval rhogw - mass concentration of water vapour air in gas phase

   13/11/2023, TKr
*/
double con_awf1matc::get_rhogw(double pw)
{
  double rhogw,pgw;

  pgw = get_pgw(pw);

  rhogw = pgw/t0/gasr*mw;

  return(rhogw);
}



/**
   function computes water vapour pressure = Kelvin equation
   @param pw - pore water pressure

   @retval pgw - water vapour pressure = Kelvin equation

   13/11/2023, TKr
*/
double con_awf1matc::get_pgw(double pw)
{
  double pgw,pgws,pc;

  pgws = get_pgws(t0);

  pc = -pw;
  
  pgw = pgws*exp(-1.0*pc*mw/rhow0/gasr/t0);

  check_math_err();

  return(pgw);
}



/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure

   13/11/2023, TKr
*/
double con_awf1matc::get_pgws(double t)
{
  double t1,t2,t3,pgws,psl;
  
  t1 = 1.0/t;
  t2 = t*t;
  t3 = t*t*t;

  psl = c8*t1 + c9 + c10*t + c11*t2 + c12*t3 + c13*log(t);
  pgws = exp(psl);

  return(pgws);
}

