/*
    File:             ped.cpp
    Author:           Tomas Krejci, 1/12/2002, corrected 29/10/2007
    Purpose:          computes conductivity and capacity matrices 
                      for coupled heat and moisture transfer (two media - w ... water content [kg/kg], t ... temperature [K])
    Source:           my Doctoral Thesis; Bazant and Najjar, 1972; Pedersen, 1990;
    Assumptions:      water vapor is the only driving mechanism; relative humidity is up to capilary saturation
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pedersen.h"
#include "globalt.h"

pedmat::pedmat()
{
  por = 0.0;
  rho = 0.0; ceff = 0.0; chieff = 0.0;        //thermal coefficients
  a_0 = 0.0; nn = 0.0; phi_c = 0.0;           //water vapour permeability 
  delta_dry = 0.0; delta_wet = 0.0;           //water vapour permeability 
  w_h_sorp = 0.0; n_sorp = 0.0; a_sorp = 0.0; //sorption isotherm

  //new added 20.3.2004
  mw = 18.01528e-3; //molar mass of water kg.mol-1
  gasr = 8.31441;   //universal gas constant J.mol-1.K-1
  rhow = 998.0;     //kg/m^3 = water density
  awet = 0.0; bwet = 0.0;
  k_wg = 0.0; ak = 0.0; bk = 0.0; nk = 0.0; //parameters for hydralic conductivity (S-shape function)
  w_cr = 0.0; w_cap = 0.0; w_vac = 0.0; w_98 = 0.0;
  dhvap = 2.7e+5;   //enthalpy of evaporation (latent heat of vaporization)


  //sorption isothem data
  //  general function
  gf=NULL;
  s_type = 0;
}

pedmat::~pedmat()
{}


/**
   JK, 11.4.2019
*/
double pedmat::give_hum (long nn)
{
  long k;
  double h;
  
  k=Gtt->give_dof(nn,0);
  if (k>0)   {h = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {h = 0.0;}
  if (k<0)   {h = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return h;
}

/**
   JK, 11.4.2019
*/
double pedmat::give_temp (long nn)
{
  long k;
  double t;
  
  k=Gtt->give_dof(nn,1);
  if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
  if (k==0)  {t = 0.0;}
  if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval(Tp->time);}
  
  return t;
}


/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void pedmat::matcond (matrix &d,long ri,long ci,long ipp)
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
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates conductivity matrix of the material for 1D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void pedmat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double w,t;
  k = 0.0;
  
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = perm_ww(w,t);
  if((ri == 0) && (ci == 1))
    k = perm_wt(w,t);
  if((ri == 1) && (ci == 0))
    k = perm_tw(w,t);
  if((ri == 1) && (ci == 1))
    k = perm_tt(w,t);
  
  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void pedmat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double w,t;
  k = 0.0;
  
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = perm_ww(w,t);
  if((ri == 0) && (ci == 1))
    k = perm_wt(w,t);
  if((ri == 1) && (ci == 0))
    k = perm_tw(w,t);
  if((ri == 1) && (ci == 1))
    k = perm_tt(w,t);
  
  fillm(0.0,d);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
}

/**
   function creates conductivity matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/

void pedmat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double w,t;
  k = 0.0;
  
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = perm_ww(w,t);
  if((ri == 0) && (ci == 1))
    k = perm_wt(w,t);
  if((ri == 1) && (ci == 0))
    k = perm_tw(w,t);
  if((ri == 1) && (ci == 1))
    k = perm_tt(w,t);
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}




/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void pedmat::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;

  switch (n){
  case 1:{
    matcond1d_2 (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d_2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    matcond3d_2 (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}




/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void pedmat::matcond1d_2 (matrix &d,long ri,long ci,long ipp)
{
  double w,t,phi,delta_gw,kwg,pgw,rhogw;
    
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];

  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  kwg = get_kw_g(w);
  pgw = phi*exp(23.5771 - 4042.9/(t-37.58));
  rhogw = pgw/t/gasr*mw;

   
  fillm(0.0,d);
  
  if((ri == 0) && (ci == 0)){
    //water + vapour
    d[0][0] = -1.0*kwg*rhow*Tp->gr[0] - delta_gw*rhogw*Tp->gr[0];
  }
}


/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void pedmat::matcond2d_2 (matrix &d,long ri,long ci,long ipp)
{
  double w,t,phi,delta_gw,kwg,pgw,rhogw;
    
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];

  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  kwg = get_kw_g(w);
  pgw = phi*exp(23.5771 - 4042.9/(t-37.58));
  rhogw = pgw/t/gasr*mw;
  
  fillm(0.0,d);
  
  if((ri == 0) && (ci == 0)){    
    d[0][0] = -1.0*kwg*rhow*Tp->gr[0] - delta_gw*rhogw*Tp->gr[0];
    d[0][1] = -1.0*kwg*rhow*Tp->gr[1] - delta_gw*rhogw*Tp->gr[1];
  }
}



/**
   function creates conductivity matrix of the material for 1D problems (convective term)

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void pedmat::matcond3d_2 (matrix &d,long ri,long ci,long ipp)
{
  double w,t,phi,delta_gw,kwg,pgw,rhogw;
    
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];

  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  kwg = get_kw_g(w);
  pgw = phi*exp(23.5771 - 4042.9/(t-37.58));
   rhogw = pgw/t/gasr*mw;
  
  fillm(0.0,d);
  
  if((ri == 0) && (ci == 0)){
    d[0][0] = -1.0*kwg*rhow*Tp->gr[0] - delta_gw*rhogw*Tp->gr[0];
    d[0][1] = -1.0*kwg*rhow*Tp->gr[1] - delta_gw*rhogw*Tp->gr[1];
    d[0][2] = -1.0*kwg*rhow*Tp->gr[2] - delta_gw*rhogw*Tp->gr[2];
  }
}





/**
   function creates capacity matrix of the material

   @param c   - capacity matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void pedmat::matcap (double &c,long ri,long ci,long ipp)
{
  double w,t;
  c = 0.0;
  
  w = Tm->ip[ipp].av[0];
  t = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = c_ww(w,t);
  if((ri == 0) && (ci == 1))
    c = c_wt(w,t);
  if((ri == 1) && (ci == 0))
    c = c_tw(w,t);
  if((ri == 1) && (ci == 1))
    c = c_tt(w,t);
}



/**
   Function calculates permability water content-water content (k_ww)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval k_ww - conductivity
*/
double pedmat::perm_ww(double w,double t)
{
  double k_ww,phi,delta_gw,dphi_dw,p_gws;
  double kwg;
  double help1,help2,help3;

  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  dphi_dw = get_dphi_dw(w);
  p_gws = get_p_gws(t);

  k_ww = delta_gw*p_gws*dphi_dw;

  //added term for water flux 20.3.2004
  //-------------------------------------------------------------
  //old:
  //kwg = get_kw_g(w);
  //dpc_dw = get_dpc_dw(w,t);
  //k_ww = k_ww + kwg*(p_gws*dphi_dw - dpc_dw);
  //-------------------------------------------------------------
  //new:
  //suction curve function 7.8.2007 TKr:
  if (w < w_cap){//moisture content control (transient region II)
    
    kwg = get_kw_g(w);
    help1 = pow((w_cap-w)/awet,(1.0/bwet));
    help2 = exp(help1);
    help3 = bwet/(w-w_cap);

    //expression for suction pressure
    //k_ww = k_ww - kwg*pow((w_cap-w)/awet,(1.0/bwet))*exp(pow((w_cap-w)/awet,(1.0/bwet)))/bwet/(w-w_cap);
    //-------------------------------------------------------------
    k_ww = k_ww - kwg*help1*help2/help3;
  }
  else{
    w = w_cap;
  }
  
  return(k_ww);
}

/**
   Function calculates permability water content-temperature (k_wt)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval k_wt - conductivity
*/
double pedmat::perm_wt(double w,double t)
{
  double k_wt,phi,delta_gw,dpgw_dt;

  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  dpgw_dt = get_dpgw_dt(t,phi);

  k_wt = delta_gw*dpgw_dt;

  //added term for water flux 20.3.2004
  //-------------------------------------------------------------
  //old:
  //kwg = get_kw_g(w);
  //dpc_dt = get_dpc_dt(w,t);
  //k_wt = k_wt + kwg*(dpgw_dt - dpc_dt);
  //-------------------------------------------------------------
  //new:
  //suction curve function 7.8.2007 TKr:
  if (w < w_cap){//moisture content control (transient region II)
    //suction stress is set to be function of water content only
  }
  else{
    w = w_cap;
  }
  
  return(k_wt);
}

/**
   Function calculates permability water temperature-temperature (k_wt)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval k_tw - conductivity
*/
double pedmat::perm_tw(double w,double t)
{
  double k_tw,phi,delta_gw,dphi_dw,p_gws;

  //added term for vaporization 20.3.2004
  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  dphi_dw = get_dphi_dw(w);
  p_gws = get_p_gws(t);
  
  k_tw = dhvap*delta_gw*p_gws*dphi_dw;
  
  return(k_tw);
}

/**
   Function calculates permability temperature-temperature (k_tt)
   Function calculates effective thermal conductivity

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval k_tt - conductivity
*/
double pedmat::perm_tt(double w,double t)
{
  double k_tt,phi,delta_gw,dpgw_dt;

  //added term for vaporization 20.3.2004
  phi = inverse_sorption_isotherm(w);
  delta_gw = get_delta_gw(phi,w);
  dpgw_dt = get_dpgw_dt(t,phi);
    
  k_tt = chieff + dhvap*delta_gw*dpgw_dt;
  
  return(k_tt);
}

/**
   Function calculates capacity water content-water content (c_ww)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval c_ww - capacity coefficient 
*/
double pedmat::c_ww (double /*w*/,double /*t*/)
{
  double c;
  
  c=1.0*rho;

  return(c);
}

/**
   Function calculates capacity water content-temperature (c_wt)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval c_wt - capacity coefficient 
*/
double pedmat::c_wt (double /*w*/,double /*t*/)
{
  double c;
  
  c=0.0;//simplification

  return(c);
}

/**
   Function calculates capacity temperature-water content (c_tw)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval c_tw - capacity coefficient 
*/
double pedmat::c_tw (double w,double t)
{
  double phi,dphi_dw,pgw,s,rhogw,ds_dw,c;
  
  //new:
  //added term for vaporization 20.3.2004
  phi = inverse_sorption_isotherm(w);
  dphi_dw = get_dphi_dw(w);
  pgw = phi*exp(23.5771 - 4042.9/(t-37.58));
  rhogw = pgw/t/gasr*mw;
  
  s = ((1.0 - por)*rho*w - por*rhogw)/(por*rhow - por*rhogw);
  ds_dw = ((1.0 - por)*rho)/(por*rhow - por*rhogw);
  c = dhvap*por*rhogw*(phi*ds_dw - (1.0-s)*dphi_dw);
  //simplification
  //c = 0.0;

  return(c);
}

/**
   Function calculates effective thermal capacity (c_tt)

   @param w - water content [kg/kg]
   @param t - temperature [K]

   @retval c_tt - thermal capacity
*/
double pedmat::c_tt (double /*w*/,double /*t*/)
{
  double c;
  
  c = ceff * rho;

  return(c);
}


void pedmat::values_correction (vector &nv)
{
  //  moisture content
  moisture_check(nv[0],nv[1],0);
}


/**
   function checks if relative humidity is grater than 1.0
   @param w  - moisture content kg/kg
   @param t  - temperature
   @param ipp - number of integration point


   TKr, 8.8.2007
*/
void pedmat::moisture_check(double &w,double /*t*/,long ipp)
{
  //suction curve function 7.8.2007 TKr:
  if (w > w_cap){//moisture content control (transient region II)
    w = w_cap;
    //fprintf (stderr,"\n\n Moisture content w = %lf is out of region II (over capillary saturation)!!! (%s, line %d).\n",w,__FILE__,__LINE__);
  }

  //storing into gauss point
  Tm->ip[ipp].av[0] = w;
}



/**
   Function calculates the water vapor permeability delta_gw (S-shape function)
   delta_gw ... by Z. P. Bazant and L. J. Najjar (1972), Nonlinear water diffusion in nonsaturated concrete,
   MATERIAUX ET CONSTRUCTIONS, Vol. 5, No. 25, pp. 3 -20.
   phi ... relative humidity
   a_0, nn, phi_c, delta_wet ... constants obtained from experiments
   
   @param phi - relative humidity

   @retval delta_gw - water vapor permeability delta_gw
*/
double pedmat::get_delta_gw(double phi,double w)
{
  double delta_gw;

  // water vapor permeability
  //delta_gw = delta_wet * (a_0 + (1.0 - a_0)/(1.0 + pow((1.0-phi)/(1.0-phi_c),nn)));

  //according to Pedersen:
  if(phi <= 0.60){
    delta_gw = delta_dry;
  }
  if(phi > 0.60 && phi < 0.98){    		
    delta_gw = delta_dry + ((phi-0.60)/(0.98-phi))*(delta_wet-delta_dry);
  }
  if(phi >= 0.98 && w <= w_vac){
    delta_gw = delta_wet * ((w_cap - w)/(w_vac - w_98));
  }
  if(w > w_vac){
    delta_gw = 0.0;
  }
  
  delta_gw = fabs(delta_gw);

  return(delta_gw);
}

/**
   Function computes hydraulic conductivity

   @param w - water content

   @retval kw_g - hydraulic conductivity
*/
double pedmat::get_kw_g(double w)
{
  double kw_g,a0_k,wc_k;
  
  //according to Pedersen:
  //if(w <= w_cr)
  //kw_g = 0.0;
  //if((w_cr < w) && (w < w_cap))
  //kw_g = ak*exp(bk*w);
  //if(w_cap <= w)
  //kw_g = ak*exp(bk*w_cap);
  
  //S-shape function:
  a0_k = ak;// a0_k = K_min/Kmax; k_wg = Kmax
  wc_k = bk;// moisture content for average K = (Kmin + Kmax)/2
  kw_g = k_wg*(a0_k + (1.0 - a0_k)/(1 + pow((wc_k/w),nk)));

  return(kw_g);
}


/**
   Function calculates water content from relative humidity.
*/

double pedmat::sorption_isotherm(double phi)
{
  double w;

  if (s_type == 0)
    w = sorption_isotherm_data_get_w(phi);
  if (s_type == 1)
    w = sorption_isotherm_hansen(w_h_sorp,a_sorp,n_sorp,phi);

  return(w);
}


/**
   Function calculates relative humidity from water content.
*/

double pedmat::inverse_sorption_isotherm(double w)
{
  double phi;

  if (s_type == 0)
    phi = inverse_sorption_isotherm_data_get_phi(w);
  if (s_type == 1)
    phi = inverse_sorption_isotherm_hansen(w_h_sorp,a_sorp,n_sorp,w);

  return(phi);
}


/**
   Function calculates derivative of relative humidity with respect to water content.
*/

double pedmat::get_dphi_dw(double w)
{
  double dphi_dw;

  if (s_type == 0)
    dphi_dw = inverse_sorption_isotherm_data_get_derivative(w);
  if (s_type == 1)
    dphi_dw = get_dphi_dw_hansen(w);

  return(dphi_dw);
}


/**
   Function calculates water content from relative humidity.
   Sorption isotherm by C. R. Pedersen (1990) according to Hansen, Combined heat and moisture transfer in building constructions, 
   PhD-thesis, Technical University of Denmark, Lingby.
   w (kg/kg) ... water content
   phi ... relative humidity
   w_h, n, a ... constants obtained from experiments
   
   @param phi - relative humidity

   @retval w - water content [kg/kg]
*/
double pedmat::sorption_isotherm_hansen(double w_h_s, double a_s, double n_s, double phi)
{
  double w;
  
  // water content
  w = w_h_s*pow((1.0-log(phi)/a_s),(-1.0/n_s));

  return(w);
}

double pedmat::inverse_sorption_isotherm_hansen(double w)
{
  double phi;
  
  if (w > w_cap)//moisture content control (transient region II)
    w = w_h_sorp;

  // relative humidity
  phi = inverse_sorption_isotherm_hansen(w_h_sorp,a_sorp,n_sorp,w);

  return(phi);
}


double pedmat::suction_curve(double s)
{
  double w;

  if (s <= 0)//moisture content control (transient region II)
    w = w_cap;
  else
    w = w_cap - awet*pow(log(s),bwet);

  return(w);
}


double pedmat::inverse_suction_curve(double w,double /*t*/)
{
  double s,help;

  if (w > w_cap)//moisture content control (transient region II)
    w = w_cap;

  help = pow(((w_cap-w)/awet),(1.0/bwet));
  s = exp(help);

  return(s);
}

/**
   Function calculates relative humidity from water content (inverse relation form sorption isotherm).
   Sorption isotherm by C. R. Pedersen (1990) according to Hansen, Combined heat and moisture transfer in building constructions, 
   PhD-thesis, Technical University of Denmark, Lingby.
   w (kg/kg) ... water content
   phi ... relative humidity
   w_h, n, a ... constants obtained from experiments
   
   @param w - water content

   @retval phi - relative humidity
*/
double pedmat::inverse_sorption_isotherm_hansen(double w_h_s, double a_s, double n_s, double w)
{
  double phi;  
  
  if (w > w_cap)//moisture content control (transient region II)
    w = w_h_sorp;

  // relative humidity
  phi = exp(a_s*(1.0-pow((w_h_s/w),(n_s))));

  return(phi);
}

/**
   Function calculates derivative of relative humidity with respect to water content (inverse relation form sorption isotherm).
   Sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions, 
   PhD-thesis, Technical University of Denmark, Lingby.
   w (kg/kg) ... water content
   phi ... relative humidity
   w_h, n, a ... constants obtained from experiments
   
   @param w - water content

   @retval dphi_dw ...
*/
double pedmat::get_dphi_dw_hansen(double w)
{

  double  dphi_dw;
  
  if (w > w_cap)//moisture content control (transient region II)
    w = w_h_sorp;
  
  dphi_dw = exp(a_sorp*(1.0-pow((w_h_sorp/w),n_sorp)))*a_sorp*n_sorp*pow(w_h_sorp,n_sorp)*pow(w,(-1.0-n_sorp)); 

  return(dphi_dw);
}  

/**
   Function calculates saturation water vapor pressure
   saturation water vapor pressure by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions, 
   PhD-thesis, Technical University of Denmark, Lingby.
   t ... temperature [K]
   phi ... relative humidity
   
   @param t - temperature [K]

   @retval p_gws - saturation water vapor pressure
*/
double pedmat::get_p_gws(double t)
{
  double p_gws;
  
  //saturation water vapor pressure ... analytical expression
  p_gws = exp(23.5771 - 4042.9/(t-37.58));
  
 return(p_gws);
}

/**
   Function calculates derivative of water vapor pressure with respect to temperature
   saturation water vapor pressure by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions, 
   PhD-thesis, Technical University of Denmark, Lingby.
   t ... temperature [K]
   phi ... relative humidity
   
   @param t - temperature [K]
   @param phi - relative humidity

   @retval dpgw_dt
*/
double pedmat::get_dpgw_dt(double t, double phi)
{
  double dp_gw_dt,dp_gws_dt;
  
  dp_gws_dt = exp(23.5771 - 4042.9/(t-37.58))*4042.9/(t-37.58)/(t-37.58);
  
  dp_gw_dt = phi*dp_gws_dt;
  
 return(dp_gw_dt);
}

/**
   Function computes derivative of capillary presure with respect to water content
   
   @param t - temperature [K]
   @param w - relative humidity

   @retval dpc_dw
*/
double pedmat::get_dpc_dw(double w, double t)
{
  double dphi_dw,rh,dpc_drh,dpc_dw;
  
  dphi_dw = get_dphi_dw(w);
  rh = inverse_sorption_isotherm(w);

  dpc_drh = -1.0*rhow/gasr/t*mw/rh;

  dpc_dw = dpc_drh*dphi_dw;
  
 return(dpc_dw);
}

/**
   Function computes derivative of capillary presure with respect to temperature
   
   @param t - temperature [K]
   @param w - relative humidity

   @retval dpc_dt
*/
double pedmat::get_dpc_dt(double w, double /*t*/)
{
  double rh,dpc_dt;
  
  rh = inverse_sorption_isotherm(w);

  dpc_dt = -1.0*log(rh)/mw*rhow*gasr;
  
 return(dpc_dt);
}




/**
   function computes new transmission coefficient (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double pedmat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc,w,t;
  new_trc = 0.0;
  
  w = give_hum (nn);
  t = give_temp (nn);

  if((ri == 0) && (ci == 0))
    new_trc = get_transmission_transcoeff_ww(w,t,bc,ipp);
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;
  
  if((ri == 1) && (ci == 0))
    new_trc = 0.0;
  if((ri == 1) && (ci == 1))
    new_trc = get_transmission_transcoeff_tt(w,t,bc,ipp);

  new_trc = new_trc*trc;

  return (new_trc);
}

/**
   function computes new nodal value (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double pedmat::transmission_nodval(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  double new_nodval,w,t;
  new_nodval = 0.0;
  
  w = give_hum (nn);
  t = give_temp (nn);
  
  if((ri == 0) && (ci == 0))
    new_nodval = get_transmission_nodval_ww(nodval,w,t,bc,ipp);
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;
  
  if((ri == 1) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 1))
    new_nodval = get_transmission_nodval_tt(nodval,w,t,bc,ipp);

  return (new_nodval);
}


/**
   function computes flux (for transmission_vector) 
   for boundary condition (third kind of boundary condition)

   @param nodval     - prescribed nodal value
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double pedmat::transmission_flux(double nodval,double /*trc2*/,long ri,long ci,long nn,long bc,long ipp)
{
  double flux,w,t;
  flux = 0.0;

  w = give_hum (nn);
  t = give_temp (nn);
  
  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_ww(nodval,w,t,bc,ipp);
  if((ri == 0) && (ci == 1))
    flux = 0.0;
  
  if((ri == 1) && (ci == 0))
    flux = 0.0;
  if((ri == 1) && (ci == 1))
    flux = get_transmission_flux_tt(nodval,w,t,bc,ipp);

  return (flux);
}


/**
   function creates transfer coefficient on the boundary for prescribed condition (water vapour pressure)

   @param w   - water content
   @param t   - temperature
   @param bc  - type of boundary condition
*/
double pedmat::get_transmission_transcoeff_ww(double w,double t,long bc,long /*ipp*/)
{
  double trc;
  double pgws;

  switch (bc){//type of prescribed variable
  case 30:{// relative humidity
    //pgw -> w approximated by Taylor's serie in point w
    pgws = get_p_gws(t);
    
    trc = pgws*get_dphi_dw(w);
    break;
  }
  case 31:{//water vapour pressure pgw
    //pgw -> w approximated by Taylor's serie in point w
    pgws = get_p_gws(t);
    
    trc = pgws*get_dphi_dw(w);
    break;
  }
  case 32:{//water vapour pressure pgw
    trc = 0.0;
    break;
  }
  case 33:{//water vapour pressure pgw
    trc = 0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  
  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (water vapour pressure)

   @param bv - prescribed value near the boundary
   @param w - actual water conctent on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double pedmat::get_transmission_nodval_ww(double bv,double w,double t,long bc,long /*ipp*/)
{  
  double nodval;
  double pgws,rh;

  switch (bc){//type of prescribed variable
  case 30:{// relative humidity -> pgw (approximation by Taylor's serie in point w)
    
    pgws = get_p_gws(t);
    nodval = bv*pgws;
    nodval = nodval - pgws*inverse_sorption_isotherm(w) + pgws*w*get_dphi_dw(w);//adding of 0th member and part of 1st member of serie
    break;
  }
  case 31:{//water vapour pressure pgw
    pgws = get_p_gws(t);
    nodval = bv;
    nodval = nodval - pgws*inverse_sorption_isotherm(w) + pgws*w*get_dphi_dw(w);//adding of 0th member and part of 1st member of serie
    break;
  }
  case 32:{//relative humidity
    //w -> pgw
    pgws = get_p_gws(t);
    rh = inverse_sorption_isotherm(w);

    nodval = pgws*rh;
    bv = pgws*bv;
    nodval = bv - nodval;
    break;
  }
  case 33:{//water vapour pressure pgw
    //w -> pgw
    pgws = get_p_gws(t);
    rh = inverse_sorption_isotherm(w);
    
    nodval = pgws*rh;
    nodval = bv - nodval;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates flux on the boundary (convective mass transfer)

   @param bv - prescribed value near the boundary
   @param w - water content
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double pedmat::get_transmission_flux_ww(double bv,double w,double t,long bc,long /*ipp*/)
{
  double flux;
  double pgws,rh;

  switch (bc){//type of prescribed variable
  case 30:{//relative humidity
    //w -> pgw
    pgws = get_p_gws(t);
    rh = inverse_sorption_isotherm(w);

    flux = pgws*rh;
    bv = pgws*bv;
    flux = bv - flux;
    break;
  }
  case 31:{//water vapour pressure pgw
    //w -> pgw
    pgws = get_p_gws(t);
    rh = inverse_sorption_isotherm(w);
    
    flux = pgws*rh;
    flux = bv - flux;
    break;
  }
  case 32:{//relative humidity
    //w -> pgw
    pgws = get_p_gws(t);
    rh = inverse_sorption_isotherm(w);

    flux = pgws*rh;
    bv = pgws*bv;
    flux = bv - flux;
    break;
  }
  case 33:{//water vapour pressure pgw
    //w -> pgw
    pgws = get_p_gws(t);
    rh = inverse_sorption_isotherm(w);
    
    flux = pgws*rh;
    flux = bv - flux;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function creates transfer coefficient on the boundary for prescribed condition (temperature)

   @param w   - water content
   @param t   - temperature
   @param bc  - type of boundary condition
*/
double pedmat::get_transmission_transcoeff_tt(double /*w*/,double /*t*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission 
    trc = 1.0;
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)    
    trc = 0.0;
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  
  return(trc);
}

/**
   function computes the right prescribed value on the boundary for prescribed condition (temperature)

   @param bv - prescribed value near the boundary
   @param w - actual water conctent on the boundary
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double pedmat::get_transmission_nodval_tt(double bv,double /*w*/,double t,long bc,long /*ipp*/)
{  
  double nodval;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    nodval = bv;
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    nodval = (bv - t);
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(nodval);
}


/**
   function creates heat flux on the boundary

   @param bv - prescribed value near the boundary
   @param w - water content
   @param t - actual temperature on the boundary
   @param bc - type of boundary condition
*/
double pedmat::get_transmission_flux_tt(double bv,double /*w*/,double t,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//heat transmission
    flux = (bv - t);
    break;
  }
  case 31:{//heat transmission for testing (and for boundary flux)
    flux = (bv - t);
    break;
  }
  default:{
    fprintf (stderr,"\n\n No real boundary condition is prescribed (%s, line %d).\n",__FILE__,__LINE__);
    exit(0);
  }
  }
  return(flux);
}


/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param pc - capillary pressure on actual node
   @param pg - gas pressure on actual node
   @param t - temperature on actual node
*/

double pedmat::get_othervalue(long compother,double w,double t)
{
  double other;

  switch (compother){
  case 0:{//moisture content
    other = w;
      break;
  }
  case 1:{//relative humidity
    other = inverse_sorption_isotherm(w);
    break;
  }
  case 2:{//temperature
    other = t;
      break;
  }
  case 3:{//water vapour pressure pgw
    other = inverse_sorption_isotherm(w) * get_p_gws(t);
    break;
  }
  case 4:{//water pressure pw
    double pc,pgw,rh;

    pgw = inverse_sorption_isotherm(w) * get_p_gws(t);
    rh = inverse_sorption_isotherm(w);
    pc = -rhow*gasr*t/mw*log(rh);
    other = pgw - pc;
    break;
  }
  case 5:{//capillary pressure pc (from Kelvin-Laplace law)
    double rh;
    
    rh = inverse_sorption_isotherm(w);
    other = -rhow*gasr*t/mw*log(rh);
    break;
  }
  case 6:{//suction pressure s = pc (from suction curve)
    other = inverse_suction_curve(w,t);
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return (other);

}

/**
     function prints names of all variables in nodes
     @param out - output file
     @param compother - number of other components
*/
void pedmat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//moisture content
    fprintf (out,"Moisture content u (kg/kg)");
    break;
  }
  case 1:{//relative humidity
    fprintf (out,"Relative humidity ()      ");
    break;
  }
  case 2:{//temperature
    fprintf (out,"Temperature (K)           ");
    break;
  }
  case 3:{//water vapour pressure pgw
    fprintf (out,"Water vapour pressure (Pa)");
    break;
  }
  case 4:{//water pressure pw
    fprintf (out,"Water pressure (Pa)       ");
    break;
  }
  case 5:{//capillary pressure pc (from Kelvin-Laplace law)
    fprintf (out,"Capillary pressure (Pa)   ");
    break;
  }
  case 6:{//suction pressure s = pc (from suction curve)
    fprintf (out,"Suction pressure (Pa)     ");
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function reads parameters
   
   @param in - input file
*/
void pedmat::read(XFILE *in)
{
  //old
  //xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &rho, &ceff, &chieff, &a_0, &nn, &phi_c, &delta_wet, &w_h, &n, &a, &awet, &bwet, &w_cr, &w_cap, &w_vac, &w_98);
  
  //new 7.8.2007
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %d", &por, &rho, &ceff, &chieff, &delta_dry, &delta_wet, &s_type);

  if(s_type == 0)
    sorption_isotherm_data_read(in);
  else
    xfscanf (in,"%lf %lf %lf", &w_h_sorp, &a_sorp, &n_sorp);

  xfscanf (in," %lf %lf %lf %lf %lf %lf %lf %lf %lf", &awet, &bwet, &w_cr, &w_cap, &w_vac, &k_wg, &ak, &bk, &nk);
  
  w_98 = sorption_isotherm(0.98);
}


/**
   function prints parameters
   
   @param out - output file
*/
void pedmat::print(FILE *out)
{
  //old
  //fprintf (out,"  %e %e %e %e %e %e %e %e %e %e %e %e %e %e", rho, ceff, chieff, a_0, nn, phi_c, delta_wet, w_h, n, a, awet, bwet, w_cr, w_cap, w_vac, w_98);

  //new 7.8.2007
  fprintf (out," %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e   %e ", por,rho, ceff, chieff, delta_dry, delta_wet, w_h_sorp, a_sorp, n_sorp, awet, bwet, w_cr, w_cap, w_vac, k_wg, ak, bk, nk);

}




/* section for sorption isotherm from measured data */

/**
   function deals with sorption isother from experiments


*/
double pedmat::sorption_isotherm_data_get_w(double phi)
{
  double help,w;

  w = gf->getval3(phi,help); //isotherma je zadana obracene!!

  return w;
}


/**
   function deals with sorption isother from experiments


*/
double pedmat::inverse_sorption_isotherm_data_get_phi(double w)
{
  double phi;

  phi = gf->getval (w); //isotherma je zadana obracene!!

  return phi;
}


/**
   function deals with sorption isother from experiments

   returns dphi/dw

*/
double pedmat::inverse_sorption_isotherm_data_get_derivative(double w)
{
  double dphi_dw,phi;

  phi = gf->getval2(w, dphi_dw); //isotherma je zadana obracene!!

  return dphi_dw;
}


/**
   function deals with sorption isother from experiments


*/
void pedmat::sorption_isotherm_data_read(XFILE *in)
{
  gf = new tablefunct;
  gf->read(in);
}



/**
   function deals with sorption isother from experiments


*/
void pedmat::sorption_isotherm_data_print(FILE *out)
{
  if (gf != NULL)
    gf->print (out);
}

