/*
  File:             consol_hawf3.cpp
  Author:           Tomas Krejci, 24/11/2017
  Purpose:          computes conductivity and capacity matrices in a material point for consol_hawf3 porous media;
                    material model for saturated-nonsaturated air and water flow and heat transfer 
                    in a deforming porous medium (soils)
  unknowns:         number of unknowns=3, pw = pore water pressure, pg = pore gas(air) pressure, t = temperature
  sources:          
  
  FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
  --------------------------------------------------
  
  sources: 
  1. THE FINITE ELEMENT METHOD IN THE STATIC AND D?YNAMIC DEFORMATION AND CONSOLIDATION OF POROUS MEDIA
  R. W. Lewis, B.A. Schrefler, pp. 354-396
  
  2. NONLINEAR MODELLING OF CONCRETE AS CONSOL_HAWF3 POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
  Francesco Pesavento - doctoral thesis                      
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
#include "consol_hawf3c.h"
#include "intpoints.h"


con_hawf3matc::con_hawf3matc()
{
  compress = 0; //compressible grains: 0=no; 1=yes
  vol_strain_effectc = 0; //volumetric strain rate influence: 0=no; 1=yes
  wrc_vol_strain_effectc = 0; //volumetric strain rate influence on water retention curve: 0=no; 1=yes
  pore_press_effectc = 0; //pore pressure effect acording to effective stress concept:  0=no; 1=yes
  temper_effectc = 0;       //temperature effect on displacements
  sr_type = 1;           //retention curve calculation type
  xi_type = 1;           //effective stress parameter type
  betas_type = 0;        //thermal expansion calculation type

  //STATE VARIABLES
  
  mw = 18.01528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

  t0 = 273.15;
  p0 = 101325.0;

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  // from D.Gawin, F.Pesavento (PRVAP.f90)
  //Hyland and Wexler equation (1983) for the saturation water vapor pressure
  c8 = -5.8002206e+03;
  c9 = 1.3914993;
  c10 =-4.8640239e-02;
  c11 = 4.1764768e-05;
  c12 = -1.4452093e-08;
  c13 = 6.5459673;

  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // alphaw = 3.53e-8    (eq. 50)
  
  // PHYSICAL PROPERTIES OF WATER
  // from Dariusz Gawin (WATPROP.f90)
  //rhow0 = 999.84;//testing??!!
  rhow0 = 1000.0;//testing??!!
  tcr = 647.3;
  hvap0 =  2.7e+5;

  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // hvap0 = 2.672e+5    (eq. 49)

  a0 =  4.8863e-7; a1 = -1.6528e-9; a2 =  1.8621e-12;
  a3 =  2.4266e-13; a4 = -1.5996e-15; a5 =  3.3703e-18;
  b0 =  1.0213e3; b1 = -7.7377e-1; b2 =  8.7696e-3;

  // M. Starnoni 10-11-2010
  // Gaw-Pes-Sch "Modelling of hygro-thermal behaviour and damage of concrete at T above the cr point of water"
  // b0 = 1.02e-3    (eq. 35)
  
  b3 = -9.2118e-5; b4 =  3.3534e-7; b5 = -4.4034e-10;
  pr1 = 1.0e7; prif = 2.0e7;

  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha0 = 0.0;
  phi0 = 0.0;
  betas0 = 0.0;
  betas_dry = 0.0;
  betas_wet = 0.0;
  rhos0 = 0.0;
  emod = 1.3e6;
  nu = 0.25;

  gamma = 0.0;
  lambda0 = 0.0;
  s_entry = 0.0;
}

con_hawf3matc::~con_hawf3matc()
{}


/**
   function reads parameters
   
   @param in - input file

   29/10/2009, TKr
*/
void con_hawf3matc::read(XFILE *in)
{
  xfscanf (in,"%k%m","heatairwaterflowmechtype",&heatairwaterflowmechtype_kwdset, &model_type);
  xfscanf (in,"%lf %d %d %d  %d %lf %lf  %lf %lf %d %d %d", &alpha0, &vol_strain_effectc, &pore_press_effectc, &wrc_vol_strain_effectc, &temper_effectc, &phi0, &rhos0, &emod, &nu, &sr_type, &xi_type, &betas_type);
  
  switch (model_type){
  case lewis_and_schrefler3_coup:
  case lewis_and_schrefler3_mefel_coup: {//Lewis and Schrefler's book
    
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
    
    //thermal expansion calculation:
    switch (betas_type){
    case 0:{//constant
      xfscanf (in,"%le", &betas0);
      break;
    }
    case 1:{//dependent on moisture
      xfscanf (in,"%le %le", &betas_dry, &betas_wet);
      break;
    }
    default:{
      print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
    }
    } 
    
    break;
  }
    
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  } 
}


/**
   function reads parameters
   
   @param out - output file

   29/10/2009, TKr
*/
void con_hawf3matc::print(FILE *out)
{
  /*  fprintf (out,"\n %d ", int(model_type));
      switch (model_type){
      case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book
      fprintf (out,"\n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda0,emod,nu);
      break;
      }
      case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's book
      fprintf (out,"\n %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",mefel_units,alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda_dry,lambda_wet,sirr0,tau0);
      break;
      }
      case lewis_and_schrefler3_2_coup:{//Lewis and Schrefler's book p. 381
      fprintf (out,"\n %lf %lf %lf %lf %lf %lf %lf %lf \n",alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda0);
      break;
      }
      default:{
      fprintf (stderr,"\n unknown model type is required");
      fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
      }
      } 
  */
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
void con_hawf3matc::matcond_u (matrix &d,long ri,long ci,long ipp)
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
  
  if (ci == 2){
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
void con_hawf3matc::matcond1d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcond2d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcond2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
    
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcond3d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcap_u (matrix &d,long ri,long ci,long ipp)
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
void con_hawf3matc::matcap1d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcap2d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pw,pg,t,ipp);
  
  
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
void con_hawf3matc::matcap2d_ax_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pw,pg,t,ipp);
  
  
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
void con_hawf3matc::matcap3d_u (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  
  pw = Cmu->ip[ipp].av[0];
  pg = Cmu->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capuw(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pw,pg,t,ipp);
  
  
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
void con_hawf3matc::matcond_l (matrix &d,long ri,long ci,long ipp)
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
void con_hawf3matc::matcond1d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcond2d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcond2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
    
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcond3d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    k = get_kwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcap_l (matrix &d,long ri,long ci,long ipp)
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
void con_hawf3matc::matcap1d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pw,pg,t,ipp);
  
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
void con_hawf3matc::matcap2d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pw,pg,t,ipp);
  
  
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

   18/06/2018, TKr
*/
void con_hawf3matc::matcap2d_ax_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pw,pg,t,ipp);
  
  
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
void con_hawf3matc::matcap3d_l (matrix &d,long ri,long ci,long ipp)
{
  double pw,pg,c,t;
  c = 0.0;
  
  pw = Cml->ip[ipp].av[0];
  pg = Cml->ip[ipp].av[1];
  t = Cmu->ip[ipp].av[2];
  
  if((ri == 0) && (ci == 0))
    c = get_capwu(pw,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pw,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pw,pg,t,ipp);
  
  
  fillm(0.0,d);


  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = c;
  d[0][3] = 0.0;
  d[0][4] = 0.0;
  d[0][5] = 0.0;
}

//////////////////////////////////////////////////////////
////this below must be corrected:


/**
   function checks if gas pressure is greater than vapour pressure

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   TKr, 24/11/2017
*/
void con_hawf3matc::gaspress_check(double /*pw*/,double &/*pg*/,double /*t*/,long /*ipp*/)
{
  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book
    //this must be checked:
    //if (pg <= 0.0)
    //pg = 1.0;
    
    //storing into gauss point
    //Tm->ip[ipp].av[1] = pg;
    break;
  }
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    break;
  }
  case lewis_and_schrefler3_2_coup:{//Lewis and Schrefler's book p. 381
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function checks pore water pressure

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   TKr, 24/11/2017
*/
void con_hawf3matc::waterpress_check(double &/*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  switch (model_type){
  case lewis_and_schrefler3:{//Lewis and Schrefler's book
    //this must be checked:
    //if (pw >= 0.0) 
    //pw = -1.0;
    
    //storing into gauss point
    //Tm->ip[ipp].av[0] = pw;
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    break;
  }
  case lewis_and_schrefler3_2:{//Lewis and Schrefler's book p. 381
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function corrects values of variables
   
   @param nv - array with variables

   TKr, 24/11/2017
*/
void con_hawf3matc::values_correction (vector &nv)
{
  //  pore water pressure
  waterpress_check(nv[0],nv[1],nv[2],0);
  
  //  gas pressure
  gaspress_check(nv[0],nv[1],nv[2],0);
}



/**
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param pg - air pressure

   @retval sw - degree of saturation

   TKr 24/11/2017
*/
double con_hawf3matc::get_sw(double pw, double pg, double t, long ipp)
{
  double sw,pc;
  sw = 0.0;  
  pc = 0.0;
  
  switch (sr_type){
  case bazant_sr:{//Bazant
    pc = pg - pw;
    sw = bazant_ret.sat(pc,t);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;

    sw = lewis_ret.sw(-pc);
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
    
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    //pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    sw = van_genuchten_ret.sw(-pc,t);//suction = capillary pressure
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model type 2
    //sw = van_genuchten_ret.sw2(pw);
    break;
  }
  case masin_sr:{//extended formulation from Brooks and Correy according to Masin
    double e=0.0,dpw=0.0,dpg=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }

    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpg = Tm->ip[ipp].av[1]-Tm->ip[ipp].pv[1];
    dpc = dpg - dpw;
    
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;
    sw = masin_ret.sw(pc,dpc,e,t);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    //pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;

    sw = febex_granit_ret.sw(pc);
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
double con_hawf3matc::get_xi(double pw, double pg, double t, long ipp)
{
  double suc=0.0,sr=0.0,xi=0.0;
  
  switch (xi_type){
  case biot_xi:
    xi = get_sw(pw,pg,t,ipp);
    break;
  case biot_reduced_xi:
    //xi = gamma*get_sw(pw,ipp);
    sr = get_sw(pw,pg,t,ipp);
    //xi = pow(sr,(gamma/lambda0));
    xi = pow(sr,gamma);
    xi = (1-gamma)*xi;
    break;
  case biot_masin_xi:{//according to masin for testing
    sr = get_sw(pw,pg,t,ipp);
    suc = pg - p0 - pw;//pore gas pressure without atmospheric pressure
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
    double e=0.0,dpw=0.0,dpg=0.0,dpc=0.0,pc=0.0,por=0.0;
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
      e = por/(1-por);
    }
    else{
      e = phi0/(1-phi0);
    }
    
    dpw = Tm->ip[ipp].av[0]-Tm->ip[ipp].pv[0];
    dpg = Tm->ip[ipp].av[1]-Tm->ip[ipp].pv[1];
    dpc = dpg - dpw;

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;

    xi = masin_ret.psi(pc,dpc,e,t);//positive value of suction
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
   function computes Biot's constant
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   
   @retval alpha - Biot's constant

   TKr 24/11/2017
*/
double con_hawf3matc::get_alpha(double /*pw*/, double /*pg*/, double /*t*/,long /*ipp*/)
{
  double alpha;
  //double kt,ks;
  
  //kt = get_kt(pw,pg,t,ipp);
  //ks = get_ks(pw,pg,t,ipp);
  
  //alpha = 1.0 - kt/ks;  

  alpha = alpha0;
  
  return(alpha);
}

/**
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kuw - conductivity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_kuw(double pw, double pg, double t,long ipp)
{
  double chi,kuw;

  kuw = 0.0;

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book  
    chi = get_xi(pw,pg,t,ipp);
    kuw = -chi;//p. 362, but non-incremental eq. for displ.; effective stress param.
    break;
  }    
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
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
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kug - conductivity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_kug(double pw, double pg, double t,long ipp)
{
  double chi,kug;
  double sg;
  
  kug = 0.0;

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book; effective stress param.
    chi = get_xi(pw,pg,t,ipp);
    sg = 1.0 - chi;
    kug = -sg;//p. 362, but non-incremental eq. for displ.
    break;
  }    
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    kug = 0.0;// pore pressure effect included in MEFEL
    if(pore_press_effectc == 1){
      // this below in not completed yet
      chi = 0.0;
      //chi = Tm->givenontransq(eff_stress_param, ipp);      //actual effective stress parameter from MEFEL
      kug =  -chi;// pore pressure effect included in METR
    }
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   
  
  return(kug);
}


/**
   function returns conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kgu - conductivity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_kgu(double /*pw*/, double /*pg*/, double /*t*/,long /*ipp*/)
{
  double kgu;

  kgu = 0.0;//p. 362

  return(kgu);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kwu - conductivity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_kwu(double /*pw*/, double /*pg*/, double /*t*/,long /*ipp*/)
{
  double kwu;

  kwu = 0.0;//p. 362

  return(kwu);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capuw - capacity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_capuw(double /*pw*/, double /*pg*/, double /*t*/,long /*ipp*/)
{
  double capuw;
  
  capuw = 0.0;//p. 97 and p. 362, but non-incremental eq. for displ.
  
  return(capuw);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capwu - capacity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_capwu(double pw, double pg, double t,long ipp)
{
  double capwu;
  double sw,sg,alpha,rhow,rhogw;
  double n,dsr_depsv;
  
  capwu = 0.0;

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book  
    sw = get_sw(pw,pg,t,ipp);
    alpha = get_alpha(pw,pg,t,ipp);
    rhogw = get_rhogw(pw,pg,t);
    rhow = get_rhow(t);
    sg = 1.0 - sw;
    
    capwu = (sg*rhogw + sw*rhow)*alpha;//p. 394 or p. 95 multiplied by rhogw and plus vapour effect
    
    break;
  }    
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    capwu = 0.0;// volumetric strain effect included in TRFEL
    if(vol_strain_effectc == 1){
      sw = get_sw(pw,pg,t,ipp);
      alpha = get_alpha(pw,pg,t,ipp);
      rhogw = get_rhogw(pw,pg,t);
      rhow = get_rhow(t);
      sg = 1.0 - sw;
      
      capwu = (sg*rhogw + sw*rhow)*alpha;//p. 394 or p. 95 multiplied by rhogw and plus vapour effect
      
      if(wrc_vol_strain_effectc == 1){
	// this must be rewritten for METR elements
	n = get_porosity(pw,pg,t,ipp);
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	capwu = capwu + n*(rhow - rhogw)*dsr_depsv;//influence of volumetric strain on saturation degree
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
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capug - capacity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_capug(double /*pw*/, double /*pg*/, double /*t*/,long /*ipp*/)
{
  double capug;
  
  capug = 0.0;//p. 97 and p. 362, but non-incremental eq. for displ.
  
  return(capug);

}



/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capgu - capacity coefficient

   TKr 24/11/2017
*/
double con_hawf3matc::get_capgu(double pw, double pg, double t,long ipp)
{
  double capgu;
  double sg,sw,alpha,rhoga;
  double n,dsr_depsv;

  capgu = 0.0;

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book  
    
    sw = get_sw(pw,pg,t,ipp);
    alpha = get_alpha(pw,pg,t,ipp);
    rhoga = get_rhoga(pw,pg,t);
    
    sg = 1.0-sw;
    
    capgu = sg*rhoga*alpha;//p. 95 multiplied by rhoga or p. 395
    
    break;
  }    
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    capgu = 0.0;// this part is in TRFEL already included
    if(vol_strain_effectc == 1){
      
      sw = get_sw(pw,pg,t,ipp);
      alpha = get_alpha(pw,pg,t,ipp);
      rhoga = get_rhoga(pw,pg,t);
      
      sg = 1.0-sw;
      
      capgu = sg*rhoga*alpha;//p. 95 multiplied by rhoga or p. 395
      
      if(wrc_vol_strain_effectc == 1){
	// this must be rewritten for METR elements
	n = get_porosity(pw,pg,t,ipp);
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	capgu = capgu - n*rhoga*dsr_depsv;//influence of volumetric strain on saturation degree
      }
    }
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
   function creates conductivity coefficient of the general material
   @param pw - water capilary pressure
   @param pg - capilary gas pressure
   @param t - temperautre
   @param ipp - number of integration point

   TKr 24/11/2017
*/
double con_hawf3matc::get_kut(double /*pw*/,double /*pg*/,double /*t*/,long ipp)
{
  double kut;
  double betas;
  
  kut = 0.0;

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book  
    betas = get_betas(ipp);  
    kut = -1.0*betas/3.0;//p. 348
    break;
  }    
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    kut = 0.0;// this part is in MEFEL already included    
    if(temper_effectc == 1){
      
      betas = get_betas(ipp);  
      kut = -1.0*betas/3.0;//p. 348//temperature effect included in METR
    }
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
   function creates conductivity coefficient of the general material
   @param pw - water capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   TKr 24/11/2017
*/
double con_hawf3matc::get_ktu(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double ktu;

  ktu = 0.0;//p. 362
  
  return(ktu);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   TKr 24/11/2017
*/
double con_hawf3matc::get_caput(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double caput;
  
  caput = 0.0;//p. 362 but non-incremental eq. for displ.

  return(caput);
}

/**
   function creates capacity coefficient of the general material
   @param pw - water capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   TKr 24/11/2017
*/
double con_hawf3matc::get_captu(double pw,double pg,double t,long ipp)
{
  double captu;
  double dhvap,sw,rhow,alpha;
  double n,dsr_depsv;

  captu = 0.0;
  
  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book  
    dhvap = get_dhvap(t);
    sw = get_sw(pw,pg,t,ipp);
    rhow = get_rhow(t);
    alpha = get_alpha(pw,pg,t,ipp);
    
    captu = -1.0*dhvap*sw*rhow*alpha;//p. 395 but non-incremental eq. for displ.
    break;
  }    
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    captu = 0.0;// volumetric strain effect included in TRFEL
    if(vol_strain_effectc == 1){
      dhvap = get_dhvap(t);
      sw = get_sw(pw,pg,t,ipp);
      rhow = get_rhow(t);
      alpha = get_alpha(pw,pg,t,ipp);
      
      captu = -1.0*dhvap*sw*rhow*alpha;//p. 395 but non-incremental eq. for displ.
      
      if(wrc_vol_strain_effectc == 1){
	// this must be rewritten for METR elements
	n = get_porosity(pw,pg,t,ipp);
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	captu = captu - dhvap*rhow*n*dsr_depsv;//influence of volumetric strain on saturation degree
      }
    }
    break;
  }
  default:{
    fprintf (stderr,"\n unknown model type is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }   

  return(captu);
}


/**
   function computes right-hand-side volume matrix of the material
   in the required integration point for upper block
  
   @param d   - right-hand-side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   12/06/2018, TKr
*/
void con_hawf3matc::rhs_u1 (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.m;

  switch (m){
  case 1:{
    rhs1d_u1 (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    rhs2d_u1 (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    rhs2d_ax_u1 (d,ri,ci,ipp);//2D - axisymmetric
    break;
  }
  case 6:{
    rhs3d_u1 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   function creates right-hand-side matrix of the material for 1D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hawf3matc::rhs1d_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_fuw1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k*Cp->gr[0];
}

/**
   function creates right-hand-side matrix of the material for 2D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hawf3matc::rhs2d_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_fuw1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k*Cp->gr[0];
  d[1][0] = k*Cp->gr[1];
}

/**
   function creates right-hand-side matrix of the material for 2D axisymmetric problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hawf3matc::rhs2d_ax_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_fuw1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k*Cp->gr[0];
  d[1][0] = k*Cp->gr[1];
}

/**
   function creates right-hand-side matrix of the material for 3D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hawf3matc::rhs3d_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pw;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling  
  
  if((ri == 0) && (ci == 0))
    k = get_fuw1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pw,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k*Cp->gr[0];
  d[1][0] = k*Cp->gr[1];
  d[2][0] = k*Cp->gr[2];
}


/**
   function computes right-hand-side matrix of the material
   in the required integration point for upper block
  
   @param d   - right-hand-side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   12/06/2018, TKr
*/
void con_hawf3matc::rhs_u2 (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.m;
 
  switch (m){
  case 1:{
    rhs1d_u2 (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    rhs2d_u2 (d,ri,ci,ipp);//2D
    break;
  }
  case 4:{
    rhs2d_ax_u2 (d,ri,ci,ipp);//2D axisymmetric
    break;
  }
  case 6:{
    rhs3d_u2 (d,ri,ci,ipp);//3D      
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  if (ci == 2){
    matrix s(d.m,d.m);
    
    Cmu->matstiff (s,ipp);//only for temperature
    mxm(s,d,d);
    
    destrm (s);
  }
}


/**
   function creates right-hand-side volume matrix of the material for 1D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void con_hawf3matc::rhs1d_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pw;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scalin
			      
  k = get_fut2(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k;
}

/**
   function creates right-hand-side volume matrix of the material for 2D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void con_hawf3matc::rhs2d_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pw;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
			      
  k = get_fut2(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = 0.0;  
}


/**
   function creates right-hand-side volume matrix of the material for 2D axisymmetric problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   12/06/2018, TKr
*/
void con_hawf3matc::rhs2d_ax_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pw;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
			      
  k = get_fut2(pw,pg,t,ipp);
  
  fillm(0.0,d);
  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = k;
  d[3][0] = 0.0;
}

/**
   function creates right-hand-side volume matrix of the material for 3D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void con_hawf3matc::rhs3d_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pw,pg,t;
  k = 0.0;
  
  pw = Cmu->ip[ipp].av[0];//*scale_pw;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling  
  
  k = get_fut2(pw,pg,t,ipp);

  fillm(0.0,d);
  d[0][0]=k;  
  d[1][0]=k;   
  d[2][0]=k;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}


/**
   function creates first part right-hand side matrix of the general material
   @param pw - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval fuw1 - first part of right-hand side matrix fu
*/

double con_hawf3matc::get_fuw1(double pw,double pg,double t,long ipp)
{
  double fuw1;
  double phi,rhos,sw,rhow,rhog;

  fuw1 = 0.0;

  phi = get_porosity(pw,pg,t,ipp);
  rhos = get_rhos(t);
  sw = get_sw(pw,pg,t,ipp);
  rhow = get_rhow(t);
  rhog = get_rhog(pw,pg,t);

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book - Liakopoulos test
    fuw1 = -1.0*((1.0 - phi)*rhos + phi*sw*rhow + phi*(1.0 - sw)*rhog);
    break;
  }
  case lewis_and_schrefler3_2_coup:{//Lewis and Schrefler's book p. 381
    fuw1 = -1.0*((1.0 - phi)*rhos + phi*sw*rhow + phi*(1.0 - sw)*rhog);
    break;
  }
  case lewis_and_schrefler3_mefel_coup:{//saturation degree ind its derivative are obtained from mefel;
    fuw1 = -1.0*((1.0 - phi)*rhos + phi*sw*rhow + phi*(1.0 - sw)*rhog);
    fuw1 = 0.0;//temporarily
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
   function creates first part right-hand side matrix of the general material

   @param pw  - water capilary pressure
   @param pg  - capilary gas pressure
   @param t   - temperature
   @param ipp - number of integration point
   
   @retval fug1 - first part of right-hand side matrix fu
*/

double con_hawf3matc::get_fug1(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  return(0.0);
}



/**
   function creates first part right-hand side matrix of the general material

   @param pw  - water capilary pressure
   @param pg  - capilary gas pressure
   @param t   - temperature
   @param ipp - number of integration point
   
   @retval fut1 - first part of right-hand side matrix fu
*/

double con_hawf3matc::get_fut1(double /*pw*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  return(0.0);
}


/**
   function creates second part right-hand side matrix of the general material

   @param pw  - water capilary pressure
   @param pg  - capilary gas pressure
   @param t   - temperature
   @param ipp - number of integration point
   
   @retval fut1 - second part of right-hand side matrix fu
*/

double con_hawf3matc::get_fut2(double /*pw*/,double /*pg*/,double /*t*/,long ipp)
{
  double fut2;
  double betas;

  fut2 = 0.0;
  
  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book
    betas = get_betas(ipp);
    fut2 = -1.0*betas/3.0;
    break;
  }
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's model and moisture storage functions is included in MEFEL
    fut2 = 0.0;// this part is in MEFEL already included
    if(temper_effectc == 1){
      betas = get_betas(ipp);
      fut2 = -1.0*betas/3.0;//p. 348//temperature effect included in METR
    }
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
   function computes mass concentration of water vapour air in gas phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhogw - mass concentration of water vapour air in gas phase
*/
double con_hawf3matc::get_rhogw(double pw,double pg,double t)
{
  double rhogw,pgw;

  pgw = get_pgw(pw,pg,t);

  rhogw = pgw/t/gasr*mw;

  return(rhogw);
}

/**
   function computes water vapour pressure = Kelvin equation
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval pgw - water vapour pressure = Kelvin equation
*/
double con_hawf3matc::get_pgw(double pw,double pg,double t)
{
  double pgw,rhow,pgws,tt,pc;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  pc = pg - pw;
  
  
  pgw = pgws*exp(-1.0*pc*mw/rhow/gasr/tt);

  return(pgw);
}




/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure
*/
double con_hawf3matc::get_pgws(double t)
{
  double t1,t2,t3,pgws,psl;
  
  t1 = 1.0/t;
  t2 = t*t;
  t3 = t*t*t;

  //critical point of water check
  if (t < tcr){
    psl = c8*t1 + c9 + c10*t + c11*t2 + c12*t3 + c13*log(t);
    pgws = exp(psl);
  }
  else
    pgws = 21780137.37214;

  return(pgws);
}



/**
   function computes water density
   @param t - temperature

   @retval rhow - water density
*/
double con_hawf3matc::get_rhow(double t)
{
  double rhow,tem;
  if (t < tcr){
    tem = t - t0;
    rhow =  (b0+(b1+(b2+(b3+(b4+b5*tem)*tem)*tem)*tem)*tem) + (pr1-prif)*
      (a0+(a1+(a2+(a3+(a4+a5*tem)*tem)*tem)*tem)*tem);
  }
  else{
    tem = tcr - t0;
    rhow =  (b0+(b1+(b2+(b3+(b4+b5*tem)*tem)*tem)*tem)*tem) + (pr1-prif)*
      (a0+(a1+(a2+(a3+(a4+a5*tem)*tem)*tem)*tem)*tem);
  }
  
  rhow = rhow0;//testing??!!
  
  return(rhow);
}



/**
   function computes  volume density of concrete skeleton,
   changes of solid density, caused by dehydratation process

   @param t - temperature

   @retval rhos - volume density of soil skeleton
*/
double con_hawf3matc::get_rhos(double /*t*/)
{
  double rhos;
   
  rhos = rhos0;
  
  return(rhos);
}



/**
   function returns porosity

   @retval por - porosity

   12/06/2018, TKr
*/
double con_hawf3matc::get_porosity(double /*pw*/, double /*pg*/,double /*t*/,long ipp)
{
  double por;

  switch (model_type){
  case lewis_and_schrefler3_coup:{//Lewis and Schrefler's book
    por = phi0; //constant value
    break;
  }
  case lewis_and_schrefler3_mefel_coup:{//Lewis and Schrefler's book
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
    }
    else
      por = phi0; //constant value for testing
    break;
  }
  case lewis_and_schrefler3_2_coup:{//Lewis and Schrefler's book p. 381
    por = phi0;
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return por;
}



/**
   function computes mass concentration of dry air in gas phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhoga - mass concentration of dry air in gas phase
*/
double con_hawf3matc::get_rhoga(double pw,double pg,double t)
{
  double rhoga,pgw;

  //gas pressure check
  pgw = get_pgw(pw,pg,t);

  if (pgw <= pg)
    rhoga = (pg - pgw)*ma/gasr/t;
  else
    rhoga = 0.0;
    //rhoga = (pg - pgw)*ma/gasr/t;

  return(rhoga);
}



/**
   function computes gas phase density
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhog - gas phase density
*/
double con_hawf3matc::get_rhog(double pw,double pg,double t)
{
  double rhog,pgw;

  pgw = get_pgw(pw,pg,t);

  rhog = (pg*ma + (mw - ma)*pgw)/gasr/t;
  
  return(rhog);
}


/**
   function computes cubic thermal expansion coefficient of solid (K-1)

   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
*/
double con_hawf3matc::get_betas(long /*ipp*/)
{
  double betas;

  betas = betas0;

  switch (betas_type){
  case 0:{ //constant
    betas = betas0;
    break;
  }
  case 1:{ //moisture dependent - not finished
    break;
  }
 default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(betas);
}


/**
   function returns Young's modulus of solid [Pa]

   @retval e -  Young's modulus of solid [Pa]
*/
double con_hawf3matc::get_e(long /*ipp*/)
{
  return(emod); //provisionally
}



/**
   function returns Poisson's ratio [-]

   @retval nu -  Poisson's ratio [-]
*/
double con_hawf3matc::get_nu(long /*ipp*/)
{
  return(nu); //provisionally
}



/**
   function computes enthalpy of evaporation (latent heat of vaporization)
   @param t - temperature

   @retval - enthalpy of evaporation (latent heat of vaporization)
*/
double con_hawf3matc::get_dhvap(double t)
{
  double dhvap,tem;

  tem = tcr - t;
  if (t < tcr)
    dhvap = hvap0*pow(tem,0.38);
  else
    dhvap = 0.0;

  // M. Starnoni 10-11-2010
  // Gaw-Maj-Sch "Numerical analysis of hygro thermal behaviour and damage of concrete at high T"
  // dhvap = 2.672e+5 * pow((t-tcr),0.38)    (eq. 49)

  //dhvap = dhvap*0.0;//debug??!!

  return(dhvap);
}
