#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aci78.h"
#include "tablefunct.h"
#include "matrix.h"
#include "global.h"
#include "mechmat.h"

/** 
  The constructor defines the variables and fills them nought except that parameter p6=1 (real value without the correction).

  created 11.9.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/ 
aci78::aci78()
{
  t_end_curing=0.0;
  slump=0.0;
  density=0.0;
  ratio_ac=0.0;
  ratio_wc=0.0;
  ratio_as=0.0;
  humidity=0.0;
  cs_thickness=0.0;
  air_content=0.0;
  p6=1.0;
  fcyl28=0.0;
  curing=0;
 }



aci78::~aci78()
{}



/** 
  The function reads input characteristics (time end of curing, time at loading, slump of fresh concrete, concrete density, 
  aggregate-cement ratio, water-cement ratio, aggregate-send ratio, humidity, cross-sectional thickness, air content, 
  compressive strength of concrete at 28 days, type of curing, concrete type) from the file.
  Parameters:
  @param in[in] - the name of input file  

  created 9.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/ 
void aci78::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %ld %ld",
	  &t_end_curing,&t_loading,&slump,&density,&ratio_ac,&ratio_wc,
	  &ratio_as,&humidity,&cs_thickness,&air_content,&fcyl28,&curing,&concrete_type);
}



/** 
  The function computes compliance function. 
  Parameters:
  @param t_current[in] - the concrete age [days]
  @param fi_t_t_dash[in] - the concrete compliance at time t since loading to t_loading [1/kPa]
  @param fcyl_t_dash[in] - the concrete strength in compression at time of loading [kPa]
  @param eps_shr_t[in] - the shrinkage at time %t (from time_end_curing)

  created 9.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/ 
void aci78::compliance (double t_current,double &fi_t_t_dash,double &fcyl_t_dash,double &eps_shr_t)
  //  input & output units:
  //  kN, m, day, kPa, C
  //
  //  internal units:
  //  kN, mm, day, kPa, C
{
  //declaration
  //***********
  //declaration of internal parameters
  double c_SI, acon=0, bcon=0;
  double eps_sh_inf, k_1, k_2, k_3, k_4=0, k_6, k_7;
  double q_1=0, q_3, q_4, q_5, q_6=0, q_7, q_8, phi_inf;
  double phi_t_t_dash,E_t_dash_SI;
  tablefunct tb;
  
  //Units conversion from input to internal
  //************************************
  double thick_mm   = cs_thickness*1000.0;     //conversion [m] -> [mm]
  double d_av_mm    = 4.0*thick_mm;
  double slump_mm   = slump*1000.;
  
  //Calculated parameters
  //*********************
  c_SI=density/(1.0+ratio_ac+ratio_wc);  //mass of cement [kg/m3]
  
  //Shrinkage
  //*********
  if (t_current>t_end_curing){
    //coefficient q_1 (119)
    if (humidity>=0.4 && humidity<=0.8) q_1=1.4-humidity;
    else if (humidity>0.8 && humidity<=1.0)  q_1=3.0-3.0*humidity;
    else{
      perror("\naci78::compliance\nRelative humidity of environment is out of range");
      exit(2);
    }
    
    //coefficient q_3 (120)
    q_3=0.89+0.00264*slump_mm;
    
    //coefficient q_4 (124-126)
    if (d_av_mm <=50.0) q_4=1.35;
    else if (d_av_mm<=150.0){ 
      tb.itype=lagrange;
      tb.asize=5;
      tb.x=new double[5];
      tb.y=new double[5];
      tb.x[0]=50.0;
      tb.x[1]=75.0;
      tb.x[2]=100.0;
      tb.x[3]=125.0;
      tb.x[4]=150.0;
      tb.y[0]=1.35;
      tb.y[1]=1.25;
      tb.y[2]=1.17;
      tb.y[3]=1.08;
      tb.y[4]=1.00;
      q_4=tb.getval(d_av_mm);
    }
    else if (d_av_mm<=380.0){
      if (t_current-t_end_curing<=365.0) q_4=1.23-0.0015*d_av_mm;
      else q_4=1.17-0.0015*d_av_mm;
    }
    else q_4=1.2*exp(-0.00473*thick_mm);
    
    //coefficient q_5 (127)
    q_5=1.0;
    if(curing!=STEAM_CURING){
      //if (t_end_curing<=1.0)  q_5=1.0;
      if (t_end_curing<=90.0){ 
	tb.itype=lagrange;
	tb.asize=6;
	tb.x=new double[6];
	tb.y=new double[6];
	tb.x[0]=1.0;
	tb.x[1]=3.0;
	tb.x[2]=7.0;
	tb.x[3]=14.0;
	tb.x[4]=28.0;
	tb.x[5]=90.0;
	tb.y[0]=1.2;
	tb.y[1]=1.1;
	tb.y[2]=1.0;
	tb.y[3]=0.93;
	tb.y[4]=0.86;
	tb.y[5]=0.75;
	q_5=tb.getval(t_end_curing);
      }
      else q_5=0.75;
    }
    
    //coefficient q_6 (121)
    if(1.0/ratio_as<=0.5) q_6=0.3+1.4/ratio_as;
    else q_6=0.9+0.2/ratio_as;
    
    //coefficient q_7 (122)
    q_7=0.95+0.008*air_content;
    
    //coefficient q_8 (123)
    q_8=0.75+0.00061*c_SI;
    
    eps_sh_inf=780.0*1.0e-6*q_1*q_3*q_4*q_5*q_6*q_7*q_8;//(118)
    
    //eps_shr_t (117)
    if (curing!=STEAM_CURING) eps_shr_t=(t_current-t_end_curing)/(35.0+t_current-t_end_curing)*eps_sh_inf*p6;//(117)
    else eps_shr_t=(t_current-t_end_curing)/(55.0+t_current-t_end_curing)*eps_sh_inf*p6;//(117)
  }
  else{
    perror("\naci78::compliance\nCurrent time is less than time and of curing");
    exit(2); 
  }
  
  
  //Creep
  
  if (t_current>t_loading){
    k_1=1.27-0.67*humidity;//(106)
    if(curing==STEAM_CURING) k_2=1.13*pow(t_loading,-0.095);//(107)
    else k_2=1.25*pow(t_loading,-0.118);
    
    k_3=0.82+0.00264*slump_mm;//(108)
    
    //coefficient k_4
    if (d_av_mm <=50.0) k_4=1.3;
    else if (d_av_mm<=150.0){ 
      tb.itype=lagrange;
      tb.asize=5;
      tb.x=new double[5];
      tb.y=new double[5];
      tb.x[0]=50.0;
      tb.x[1]=75.0;
      tb.x[2]=100.0;
      tb.x[3]=125.0;
      tb.x[4]=150.0;
      tb.y[0]=1.30;
      tb.y[1]=1.17;
      tb.y[2]=1.11;
      tb.y[3]=1.04;
      tb.y[4]=1.10;
      k_4=tb.getval(d_av_mm);
    }
    else if (d_av_mm<=380.0){
      if (t_current-t_loading<=365.0) k_4=1.14-0.91e-3*d_av_mm;//(112)
      else k_4=1.14-0.67e-3*d_av_mm;//(112)
    }
    else k_4=2.0/3.0*(1.0+1.13*exp(-0.0212*thick_mm));//(113)
    
    
    k_6=0.88+0.24/ratio_as;//(109)
    k_7=0.46+0.09*air_content;  //(110)
    if (k_7<1.0) k_7=1.0;
    
    phi_inf =2.35*k_1*k_2*k_3*k_4*k_6*k_7;                            //(105)
    phi_t_t_dash=pow(t_current-t_loading,0.6)/(10.0+pow(t_current-t_loading,0.6))*phi_inf;      //(104)
  }
  else{
    fi_t_t_dash = 0.0;
    phi_t_t_dash  = 0.0;
  }
  
  //acon, bcon - constants for calculation E, fcyl, according to ACI78, PhD thesis of Libor, Table 8.3 pg. 286
  if (concrete_type==1){
    if (curing != STEAM_CURING){
      acon=4.00;
      bcon=0.85;
    }
    else{
      acon=1.00;
      bcon=0.95;
    }
  }
  else if (concrete_type==3){
    if (curing != STEAM_CURING){
      acon=2.30;
      bcon=0.92;
    }
    else{
      acon=0.70;
      bcon=0.98;
    }
  }
  else{
    perror("\naci78::compliance\nUnsupported concrete type");
    exit(2);      
  }
  
  //Final results
  //*************
  fcyl_t_dash=fcyl28*t_loading/(acon+bcon*t_loading);  //concrete compression strength [kPa], (116)
  E_t_dash_SI   =42.8*pow(density*density*density*fcyl_t_dash*0.001,0.5); //(115) [kPa] //Modulus of elasticity [GPa] (1.16)
  fi_t_t_dash =(1.0+phi_t_t_dash)/E_t_dash_SI;       //(114) [1/kPa]
  eps_shr_t=-eps_shr_t;
}


/** 
  The function returns tangent stiffness matrix of material 
  Parameters:

  @param d[out]  - elastic stiffness matrix
  @param ipp[in] - integration point pointer
  @param ido[in] - index of internal variables for given material in the ipp other array

  created 9.11.2001 by Martin Kulhavy, mkulhavy@cml.fsv.cvut.cz
*/
void aci78::matstiff (matrix &d, long ipp, long ido)
{
  Mm->elmatstiff (d, ipp, ido);
}

