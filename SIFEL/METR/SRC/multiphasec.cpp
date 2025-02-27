/*
    File:             multiphasec.cpp
    Author:           Tomas Krejci
    Purpose:          computes coupled conductivity and capacity matrices in a material point for multiphase porous medium
                      (solid-water-vapour) for rigid prous materials

    FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
    --------------------------------------------------

    sources: 
    1. FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS BY AN ALGEBRAIC MULTIGRID METHOD
    Wang Xicheng, B.A. Schrefler (there are few mistakes)
    
    2. NUMERICAL ANALYSIS OF HYGRO-THERMAL BEHAVIOUR AND DAMAGE OF CONCRETE AT HIGH TEMPERATURE
    D. Gawin, C.E. Majorana, B.A. Schrefler

    3. NONLINEAR MODELLING OF CONCRETE AS MULTIPHASE POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
    Francesco Pesavento - doctoral thesis    
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "multiphasec.h"
#include "constrel.h"
#include "constrelcl.h"
#include "constrelcu.h"
#include "globalt.h"
#include "globalc.h"

multiphc::multiphc()
{
  scale_pc = Tp->scale[0];
  scale_pg = Tp->scale[1];
  scale_t = Tp->scale[2];
}
multiphc::~multiphc()
{}


/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function computes conductivity matrix of the material
   in the required integration point for upper block
  
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void multiphc::matcond_u (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.m;

  switch (m){
  case 1:{
    if (ci == 2){
      matrix s(d.m,d.m),sd(d.m,d.n);
      
      matcond1d_u (sd,ri,ci,ipp);//1D
      
      Cmu->matstiff (s,ipp);//only for temperature
      mxm(s,sd,d);
      
      destrm (s); destrm (sd);
    }
    else
      matcond1d_u (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    if (ci == 2){
      matrix s(d.m,d.m),sd(d.m,d.n);
      
      matcond2d_u (sd,ri,ci,ipp);//2D
      
      Cmu->matstiff (s,ipp);//only for temperature
      mxm(s,sd,d);
      
      destrm (s); destrm (sd);
    }
    else
      matcond2d_u (d,ri,ci,ipp);//2D
    break;
  }
  case 6:{
    if (ci == 2){
      matrix s(d.m,d.m),sd(d.m,d.n);
      
      matcond3d_u (sd,ri,ci,ipp);//3D
      
      Cmu->matstiff (s,ipp);//only for temperature
      mxm(s,sd,d);
      
      destrm (s); destrm (sd);
    }
    else
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
   function computes capacity matrix of the material
   in the required integration point for upper block
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void multiphc::matcap_u (matrix &d,long ri,long ci,long ipp)
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
   function creates conductivity matrix of the material for 1D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcond1d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kuc(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pc,pg,t,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcond2d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kuc(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = 0.0;  
}

/**
   function creates conductivity matrix of the material for 3D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcond3d_u (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kuc(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kug(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_kut(pc,pg,t,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[1][0]=k;   
  d[2][0]=k;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}


/**
   function creates capacity matrix of the material for 1D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcap1d_u (matrix &d,long ri,long ci,long ipp)
{
  double c;
  double pc,pg,t;
  c = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capuc(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = c;
}


/**
   function creates capacity matrix of the material for 2D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcap2d_u (matrix &d,long ri,long ci,long ipp)
{
  double c;
  double pc,pg,t;
  c = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capuc(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = 0.0;
}


/**
   function creates capacity matrix of the material for 3D problems
   for upper block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcap3d_u (matrix &d,long ri,long ci,long ipp)
{
  double c;
  double pc,pg,t;
  c = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capuc(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capug(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    c = get_caput(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = c;
  d[1][0] = c;
  d[2][0] = c;
  d[3][0] = 0.0;
  d[4][0] = 0.0;
  d[5][0] = 0.0;
}



/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function computes conductivity matrix of the material
   in the required integration point for lower block
  
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void multiphc::matcond_l (matrix &d,long ri,long ci,long ipp)
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
   function computes capacity matrix of the material
   in the required integration point for lower block
   
   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void multiphc::matcap_l (matrix &d,long ri,long ci,long ipp)
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
   function creates conductivity matrix of the material for 1D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcond1d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cml->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cml->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cml->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kcu(pc,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pc,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pc,pg,t,ipp);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcond2d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cml->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cml->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cml->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kcu(pc,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pc,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pc,pg,t,ipp);
  
  fillm(0.0,d);


  d[0][0] = k;
  d[0][1] = k;
  d[0][2] = 0.0;  
}

/**
   function creates conductivity matrix of the material for 3D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcond3d_l (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cml->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cml->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cml->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_kcu(pc,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgu(pc,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    k = get_ktu(pc,pg,t,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[0][1]=k;   
  d[0][2]=k;
  d[0][3]=0.0;   
  d[0][4]=0.0;  
  d[0][5]=0.0;
}


/**
   function creates capacity matrix of the material for 1D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcap1d_l (matrix &d,long ri,long ci,long ipp)
{
  double c;
  double pc,pg,t;
  c = 0.0;
  
  pc = Cml->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cml->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cml->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capcu(pc,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pc,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = c;
}


/**
   function creates capacity matrix of the material for 2D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcap2d_l (matrix &d,long ri,long ci,long ipp)
{
  double c;
  double pc,pg,t;
  c = 0.0;
  
  pc = Cml->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cml->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cml->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capcu(pc,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pc,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = 0.0;
}


/**
   function creates capacity matrix of the material for 3D problems
   for lower block

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::matcap3d_l (matrix &d,long ri,long ci,long ipp)
{
  double c;
  double pc,pg,t;
  c = 0.0;
  
  pc = Cml->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cml->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cml->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    c = get_capcu(pc,pg,t,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgu(pc,pg,t,ipp);
  if((ri == 2) && (ci == 0))
    c = get_captu(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = c;
  d[0][1] = c;
  d[0][2] = c;
  d[0][3] = 0.0;
  d[0][4] = 0.0;
  d[0][5] = 0.0;
}


/*---------------------------------------------------------------------------------------------------------------------*/

/**
   function checks if gass pressure is greater than vapour pressure
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point


   TKr, 18.5.2005
*/
void multiphc::gaspress_check(double pc,double &pg,double t,long ipp)
{
  double pgw;
  state_eq tt;
  
  //general
  //gas pressure check
  pgw = tt.get_pgw(pc,t);
  if (pgw > pg) 
    pg = pgw;
  
  //general
  //fully saturated state
  if(pc < 100.0)
    pg = 101325.0;
  
  //storing into gauss point
  Tm->ip[ipp].av[1] = pg;
}


/**
   function checks if capillary pressure is higher than maximum
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void multiphc::cappress_check(double &/*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  /*  double pgw,pcmax,rhow;
      state_eq tt;
      
      //only boundary nodes
      //maximum capillary pressure check
      pgw = tt.get_pgw(pc,t);
      rhow = tt.get_rhow(t);
      pcmax = 0.0;
      if (pg < pgw && pg > 0.0)
      pcmax = -1.0*rhow*gasr*t/mw*log(pg/pgw);
      else
      pcmax = 0.0;
      
      if (pg < pgw){
      if (pcmax > 0.0){
      pc = pcmax;
      }
      }
      
      
      if (pc > pcmax)
      pc = pcmax;
      
      //storing into gauss point
      Tm->ip[ipp].av[0] = pc;
  */
}


/**
   function fixed capillary pressure if on one element is temperature higher than critical temperature for water
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   
   TKr, 18.5.2005
*/
void multiphc::cappress_stop(double &/*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  //storing into gauss point
  //Tm->ip[ipp].av[0] = pc;
}


void multiphc::values_correction (vector &nv)
{
  //  capillary pressure
  cappress_check(nv[0],nv[1],nv[2],0);

  //  gas pressure
  gaspress_check(nv[0],nv[1],nv[2],0);
}

/**
   function creates conductivity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval kcu - conductivity coefficient of the general material
*/
double multiphc::get_kcu(double pc,double pg,double t,long ipp)
{
  double rhow,kintr,krw,muw,kcu;
  state_eq tt;
  state_eqcl cl;

  //gas pressure check
  //gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  //cappress_check(pc,pg,t,ipp);
  
  rhow = tt.get_rhow(t);
  kintr = cl.get_kintr(pc,pg,t,ipp);
  krw = cl.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  
  kcu = -1.0*(rhow*kintr*krw/muw*rhow);
  
  return(kcu);
}


/**
   function creates capacity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval capcu - capacity coefficient of the general material
*/
double multiphc::get_capcu(double pc,double pg,double t,long ipp)
{
  double capcu,s,rhogw,rhow,alpha;
  state_eq tt;
  state_eqcl cl;
  
  s = cl.get_s(pc,pg,t,ipp);
  rhogw = tt.get_rhogw(pc,t);
  rhow = tt.get_rhow(t);
  alpha = cl.get_alpha(pc,pg,t,ipp);
    
  capcu = ((1.0-s)*rhogw + s*rhow)*alpha;

  return(capcu);  
}

/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function creates conductivity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_kgu(double /*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double kgu;

  kgu = 0.0;

  return(kgu);

}

/**
   function creates capacity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_capgu(double pc,double pg,double t,long ipp)
{
  double capgu,s,rhoga,alpha;
  state_eq tt;
  state_eqcl cl;
  
  s = cl.get_s(pc,pg,t,ipp);
  rhoga = tt.get_rhoga(pc,pg,t);
  alpha = cl.get_alpha(pc,pg,t,ipp);
  
  capgu = (1.0-s)*rhoga*alpha;

  return(capgu); 
}

/*---------------------------------------------------------------------------------------------------------------------*/

/**
   function creates conductivity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_ktu(double pc,double pg,double t,long ipp)
{
  double dhvap,rhow,kintr,krw,muw,ktu;
  state_eq tt;
  state_eqcl cl;
  
  //gas pressure check
  //gaspress_check(pc,pg,t,ipp);
  //capillary pressure check
  //cappress_check(pc,pg,t,ipp);
  
  dhvap = tt.get_dhvap(t);
  rhow = tt.get_rhow(t);
  kintr = cl.get_kintr(pc,pg,t,ipp);
  krw = cl.get_krw(pc,t,ipp);
  muw = tt.get_muw(t);
  
  ktu = dhvap*(rhow*kintr*krw/muw*rhow);
  
  return(ktu);
}

/**
   function creates capacity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_captu(double pc,double pg,double t,long ipp)
{
  double captu,dhvap,s,rhow,alpha;
  state_eq tt;
  state_eqcl cl;
  
  dhvap = tt.get_dhvap(t);
  s = cl.get_s(pc,pg,t,ipp);
  rhow = tt.get_rhow(t);
  alpha = cl.get_alpha(pc,pg,t,ipp);
  
  captu = -1.0*dhvap*s*rhow*alpha;

  return(captu);
}


/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function creates conductivity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_kuc(double pc,double pg,double t,long ipp)
{
  double alpha,s,kuc;
  state_eqcu cu;
  
  alpha = cu.get_alpha(pc,pg,t,ipp);
  s = cu.get_s(pc,pg,t,ipp);

  kuc = alpha*s;

  return(kuc);

}

/**
   function creates capacity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_capuc(double /*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double capuc;
  
  capuc = 0.0;

  return(capuc);
}


/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function creates conductivity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_kug(double pc,double pg,double t,long ipp)
{
  double alpha,kug;
  state_eqcu cu;

  alpha = cu.get_alpha(pc,pg,t,ipp);

  kug = -1.0*alpha;

  return(kug);

}


/**
   function creates capacity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_capug(double /*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double capug;
  
  capug = 0.0;

  return(capug);
}


/*---------------------------------------------------------------------------------------------------------------------*/



/**
   function creates conductivity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperautre
   @param ipp - number of integration point
*/
double multiphc::get_kut(double pc,double pg,double t,long ipp)
{
  double betas,kut;
  state_eqcu cu;
  
  betas = cu.get_betas(pc,pg,t,ipp);
  
  kut = -1.0*betas/3.0;

  return(kut);
}


/**
   function creates capacity coefficient of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
*/
double multiphc::get_caput(double /*pc*/,double /*pg*/,double /*t*/,long /*ipp*/)
{
  double caput;
  
  caput = 0.0;

  return(caput);
}

/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function computes right-hand-side matrix of the material
   in the required integration point for upper block
  
   @param d   - right-hand-side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void multiphc::rhs_u1 (matrix &d,long ri,long ci,long ipp)
{
  long m;

  m = d.m;

  switch (m){
  case 1:{
    if (ci == 2){
      matrix s(d.m,d.m),sd(d.m,d.n);
      
      rhs1d_u1 (sd,ri,ci,ipp);//1D
      
      Cmu->matstiff (s,ipp);//only for temperature
      mxm(s,sd,d);
      
      destrm (s); destrm (sd);
    }
    else
      rhs1d_u1 (d,ri,ci,ipp);//1D
    break;
  }
  case 3:{
    if (ci == 2){
      matrix s(d.m,d.m),sd(d.m,d.n);
      
      rhs2d_u1 (sd,ri,ci,ipp);//2D
      
      Cmu->matstiff (s,ipp);//only for temperature
      mxm(s,sd,d);
      
      destrm (s); destrm (sd);
    }
    else
      rhs2d_u1 (d,ri,ci,ipp);//2D
    break;
  }
  case 6:{
    if (ci == 2){
      matrix s(d.m,d.m),sd(d.m,d.n);
      
      rhs3d_u1 (sd,ri,ci,ipp);//3D
      
      Cmu->matstiff (s,ipp);//only for temperature
      mxm(s,sd,d);
      
      destrm (s); destrm (sd);
    }
    else
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

   17.07.2005
*/
void multiphc::rhs1d_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_fuc1(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pc,pg,t,ipp);
  
  d[0][0] = k;
}

/**
   function creates right-hand-side matrix of the material for 2D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::rhs2d_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
  
  if((ri == 0) && (ci == 0))
    k = get_fuc1(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pc,pg,t,ipp);
  
  fillm(0.0,d);

  d[0][0] = k;
  d[1][0] = k;
  d[2][0] = 0.0;  
}

/**
   function creates right-hand-side matrix of the material for 3D problems
   for upper block

   @param d   - right-hand-side %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

   17.07.2005
*/
void multiphc::rhs3d_u1 (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pc,pg,t;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling  
  
  if((ri == 0) && (ci == 0))
    k = get_fuc1(pc,pg,t,ipp);
  if((ri == 0) && (ci == 1))
    k = get_fug1(pc,pg,t,ipp);
  if((ri == 0) && (ci == 2))
    k = get_fut1(pc,pg,t,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;  
  d[1][0]=k;   
  d[2][0]=k;
  d[3][0]=0.0;   
  d[4][0]=0.0;  
  d[5][0]=0.0;
}

/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function creates first part right-hand side matrix of the general material

   @param pc  - capilary pressure
   @param pg  - capilary gas pressure
   @param t   - temperature
   @param ipp - number of integration point
   
   @retval fug1 - first part of right-hand side matrix fu
*/

double multiphc::get_fuc1(double pc,double pg,double t,long ipp)
{
  double alpha,s,fuc1;
  state_eqcu cu;
  
  alpha = cu.get_alpha(pc,pg,t,ipp);
  s = cu.get_s(pc,pg,t,ipp);
  
  fuc1 = alpha*s;
  
  return(fuc1);
}

/**
   function creates first part right-hand side matrix of the general material

   @param pc  - capilary pressure
   @param pg  - capilary gas pressure
   @param t   - temperature
   @param ipp - number of integration point
   
   @retval fug1 - first part of right-hand side matrix fu
*/

double multiphc::get_fug1(double pc,double pg,double t,long ipp)
{
  double alpha,fug1;
  state_eqcu cu;
  
  alpha = cu.get_alpha(pc,pg,t,ipp);
  
  fug1 = -1.0*alpha;

  return(fug1);
}


/**
   function creates first part right-hand side matrix of the general material

   @param pc  - capilary pressure
   @param pg  - capilary gas pressure
   @param t   - temperature
   @param ipp - number of integration point
   
   @retval fut1 - first part of right-hand side matrix fu
*/

double multiphc::get_fut1(double pc,double pg,double t,long ipp)
{
  double fut1,betas;
  state_eqcu cu;
  
  betas = cu.get_betas(pc,pg,t,ipp);

  fut1 = -1.0*betas/3.0;
  
  return(fut1);
}


/*---------------------------------------------------------------------------------------------------------------------*/


/**
   function computes right-hand-side volume matrix of the material
   in the required integration point for upper block
  
   @param d   - right-hand-side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
   
   17.07.2005
*/
void multiphc::rhs_u2 (matrix &d,long ri,long ci,long ipp)
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
  case 6:{
    rhs3d_u2 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
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
void multiphc::rhs1d_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pc,pg,t;
  double g;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scalin
			      
  k = get_fu2(pc,pg,t,ipp);
  
  //g = Gt->g;//complete vector of gravity acceleration
  g = 0.0;//temp
  
  fillm(0.0,d);
  d[0][0] = k*g;
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
void multiphc::rhs2d_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pc,pg,t;
  double *g;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling
			      
  k = get_fu2(pc,pg,t,ipp);
  
  g = new double [2];
  //g = Gt->g;//complete vector of gravity acceleration
  g[0] = 0.0;//temp
  g[1] = 0.0;//temp
  
  fillm(0.0,d);
  d[0][0] = k*g[0];
  d[1][0] = k*g[1];
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
void multiphc::rhs3d_u2 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double k;
  double pc,pg,t;
  double *g;
  k = 0.0;
  
  pc = Cmu->ip[ipp].av[0];//*scale_pc;//scaling
  pg = Cmu->ip[ipp].av[1];//*scale_pg;//scaling
  t  = Cmu->ip[ipp].av[2];//*scale_t;//scaling  
  
  k = get_fu2(pc,pg,t,ipp);

  g = new double [3];
  //g = Gt->g;//complete vector of gravity acceleration
  g[0] = 0.0;//temp
  g[1] = 0.0;//temp
  g[2] = 0.0;//temp
  
  fillm(0.0,d);
  d[0][0] = k*g[0];
  d[1][0] = k*g[1];
  d[2][0] = k*g[2];
}


/**
   function creates first part right-hand side matrix of the general material
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature
   @param ipp - number of integration point
   
   @retval fu2 - first part of right-hand side matrix fu
*/

double multiphc::get_fu2(double pc,double pg,double t,long ipp)
{
  double fu2,phi,rhos,s,rhow,rhog;
  state_eq tt;
  state_eqcu cu;

  phi = tt.get_phi(t,ipp);
  rhos = cu.get_rhos(pc,pg,t,ipp);
  s = cu.get_s(pc,pg,t,ipp);
  rhow = tt.get_rhow(t);
  rhog = tt.get_rhog(pc,pg,t);

  fu2 = -1.0*((1.0 - phi)*rhos + phi*s*rhow + phi*(1.0 - s)*rhog);

  fu2 = -1.0*fu2;

  return(fu2);
}

/*---------------------------------------------------------------------------------------------------------------------*/
