#include "cebfip78.h"
#include "mechmat.h"
#include "tablefunct.h"
#include "iotools.h"
#include "matrix.h"
#include "global.h"
#include <stdlib.h>
#include <math.h>


/** 
  The constructor defines the variables and fills them nought.

  created 1.12.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/ 
cebfip78::cebfip78()
{
  t_end_curing=0.0;
  t_loading=0.0;
  humidity=0.0;
  cs_thickness=0.0;
  fcyl28=0.0;
  p6=1.0;
}


cebfip78::~cebfip78()
{
}

/** 
  The function reads input characteristics (time end of curing, time at loading, humidity, cross-sectional thickness, compressive strength of concrete at 28 days) from the file.
  Parameters:
  @param in[in] is the name of input file  

  created 1.12.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/ 
void cebfip78::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf",&t_end_curing,&t_loading,&humidity,&cs_thickness,&fcyl28);
}


/** 
  The function computes compliance function. 
  Parameters:
  @param t_current[in]   - the concrete age [days]
  @param fi_t_t_dash[in] - the concrete compliance at time t since loading to t_loading [1/kPa]
  @param fcyl_t_dash[in] - the concrete strength in compression at time of loading [kPa]
  @param eps_shr_t[in]   - the shrinkage at time %t (from time_end_curing)

  created 1.12.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/ 
void cebfip78::compliance (double t_current,double &fi_t_t_dash,double &fcyl_t_dash,double &eps_shr_t)
  //  input & output units:
  //  kN, m, day, kPa, C
  //
  //  internal units:
  //  kN, cm, day, kPa, C
  //
{
  double thick,d, eps_sh1, eps_sh2, eps_sh0, k1, k2, k3, k4, beta_sh_t, beta_sh_t0, e_28;
  double beta_i, beta_a, beta_d, phi_d, phi_f1, phi_f2, phi_f, beta_f_t, beta_f_t_dash;
  double beta_c_28, beta_c_t_dash, E_t_dash_SI, phi_t_t_dash;
  double iv=0;
  tablefunct tb;
  tb.itype=piecewiselin;
  tb.asize=4;
  tb.x=new double[4];
  tb.y=new double[4];
  tb.x[0]=0.4;
  tb.x[1]=0.7;
  tb.x[2]=0.9;
  tb.x[3]=1.0;
  tb.y[0]=1.0;
  tb.y[1]=1.5;
  tb.y[2]=5.0;
  tb.y[3]=30.0;
  
  //Units conversion from input to internal units
  //*********************************************
  thick=cs_thickness*100.0;                  //conversion [m] -> [cm]
  
  if (humidity>=0.4 && humidity<=1.0) iv=tb.getval(humidity);
  else{
    perror("\ncebfip78::compliance\nRelative humidity of environment is out of range (cebfip78.cpp)");
    exit(2);
  }
  d=2.0*thick*iv;
  
  //Shrinkage
  //*********
  if (t_current>t_end_curing){
    eps_sh1   =(775.0*humidity*humidity*humidity-1565.0*humidity*humidity+1103.25*humidity-303.25)*0.00001;//(139)
    eps_sh2   =exp(0.00174*d-0.32/d-log(pow(d,0.251)/1.9));//(139)
    eps_sh0   =eps_sh1*eps_sh2;//(138)
    k3        =11.8*d+16.0;
    k4        =exp(-0.00257*d+0.32/d+log(0.22*pow(d,0.4)));
    beta_sh_t =pow(t_current/(t_current+k3),k4);//(140)
    beta_sh_t0=pow(t_end_curing/(t_end_curing+k3),k4);//(140)
    eps_shr_t =eps_sh0*(beta_sh_t-beta_sh_t0)*p6;//(137)
  }
  else{
    perror("\ncebfip78::compliance\nCurrent time is less than time and of curing (cebfip78.cpp)");
    exit(2);
  }
  
  //Strength
  //********
  beta_c_28    =pow(28.0/(28.0+47.0),1./2.45);
  beta_c_t_dash=pow(t_loading/(t_loading+47.0),1./2.45);
  e_28         =1.25*9500*pow(fcyl28/1000.0,1./3.);
  fcyl_t_dash  =fcyl28*beta_c_t_dash/beta_c_28;
  
  //Creep
   //*****
  if (t_current>t_loading){
    beta_i       =0.875*pow((t_loading+47.0)/t_loading,1.0/7.35);//(130)
    beta_a       =0.8*pow(1.0-t_loading/(t_loading+47.0),1.0/2.45);//(131)
    beta_d       =pow((t_current-t_loading)/(t_current-t_loading+328.0),1.0/4.2);//(132)
    phi_d        =0.4;//(132)
    phi_f1       =4.45-3.5*humidity;
    phi_f2       =exp(4.4e-5*d-0.357/d-log(pow(d,0.1667)/2.6));//(134)
    phi_f        =phi_f1*phi_f2;//(133)
    k1           =exp(5.02/d+log(6.95*pow(d,1.25)));//(136)
    k2           =exp(0.00144*d-1.1/d-log(1.005*pow(d,0.2954)));//(136)
    beta_f_t     =pow(t_current/(t_current+k1),k2);//(135)
    beta_f_t_dash=pow(t_loading/(t_loading+k1),k2);//(135)
    fi_t_t_dash  =(beta_i+beta_a+phi_d*beta_d+phi_f*(beta_f_t-beta_f_t_dash))/e_28/1000.0;
    E_t_dash_SI  =e_28/(beta_i+beta_a)*1000.0;  //Modulus of elasticity [kPa]
  }
  else{
    beta_i     =0.875*pow((t_loading+47.0)/t_loading,1.0/7.35);//(130)
    beta_a     =0.8*pow(1.0-t_loading/(t_loading+47.0),1.0/2.45);//(131)
    E_t_dash_SI=e_28/(beta_i+beta_a)*1000.0;  //Modulus of elasticity [kPa]
    fi_t_t_dash=1.0/e_28/1000.0;
  }
  
  //Final results
  //*************
  fcyl_t_dash =fcyl_t_dash;
  fi_t_t_dash =fi_t_t_dash;
  phi_t_t_dash=E_t_dash_SI*fi_t_t_dash-1.0;
  eps_shr_t   =eps_shr_t;
  
  return;
}


/** 
The function returns tangent stiffness matrix of material 
Parameters:
  @param d[out]  - elastic stiffness matrix
  @param ipp[in] - integration point pointer
  @param ido[in] - index of internal variables for given material in the ipp other array

created 1.12.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/
void cebfip78::matstiff (matrix &d, long ipp, long ido)
{
  Mm->elmatstiff (d, ipp, ido);
}
