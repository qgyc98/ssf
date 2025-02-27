/*
    File:            consol_awf1.cpp; This model is for testing.
    Author:          Tomas Krejci, 12/09/2008, revised 03/10/2023
    Purpose:         material model for saturated-nonsaturated one-phase flow (water) in deforming medium plus vapour diffusion
    sources:         Lewis and Schrefler's book
    unknowns:        number of unknowns=1, pw = liqud water pressure
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "constrel.h"
#include "consol_awf1.h"
#include "globalt.h"
#include "globmatt.h"

con_awf1mat::con_awf1mat()
{
  por_type = 0;          //porosity calculation type
  kintr_type = 0;        //intrinsic permability calculation type
  krw_type = 0;          //relative permeability calculation type
  deff_type = 0;         // diffusion calculation type
  sr_type = 1;           //retention curve calculation type

  vol_strain_effect = 0; //volumetric strain rate influence: 0=no; 1=yes 
  mefel_units = 1.0;//basic units for pressures = Pa (Pascals)

  //STATE VARIABLES
  mw = 18.01528; //molar mass of water kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

  t0 = 273.15;

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  dv0 = 2.58e-5; //effective diffusion coefficient of vapour m^2/s at reference temperature 273.15 
  c8 = -5.8002206e+03;
  c9 = 1.3914993;
  c10 =-4.8640239e-02;
  c11 = 4.1764768e-05;
  c12 = -1.4452093e-08;
  c13 = 6.5459673;

  // PHYSICAL PROPERTIES OF WATER
  rhow = 1000.0;//water density [kg/m^3]
  muw0 = 1.0e-3;//water viscosity [Pa.s]
  kw = 2.0e9;//bulk modulus of water [Pa]

  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha0 = 1.0;//Biot's constant [-] alpha = 1 - kt/ks
  ks = .167e9;//bulk modulus of solid phase (grains) [Pa]
  phi0 = 0.297;//initial porosity [-]
  kintr0 = 4.5e-13;//intrinsic permeability [m^2]
  rhos0 = 2500.0; //solid density
  tau0 = 0.0;
  sirr = 0.0;
  ssat = 1.0;

  bb1 = 0.0;
  phi01 = 0.0;

  pw_bc = 0.0; //free boundary pressure 
}

con_awf1mat::~con_awf1mat()
{}


/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void con_awf1mat::matcond (matrix &d,long ri,long ci,long ipp)
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
void con_awf1mat::matcond1d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  kk = get_kww(pw,ipp);
  
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
void con_awf1mat::matcond2d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  kk = get_kww(pw,ipp);
  
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

void con_awf1mat::matcond3d (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double kk;
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  kk = get_kww(pw,ipp);
  
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
void con_awf1mat::matcap (double &cc,long /*ri*/,long /*ci*/,long ipp)
{
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  cc = get_capww(pw,ipp);
}

/**
   function reads parameters
   
   @param in - input file

   03/10/2023, TKr
*/
void con_awf1mat::read(XFILE *in)
{
  xfscanf (in,"%k%m","waterflowtype",&waterflowtype_kwdset, &model_type); //water flow model type
  xfscanf (in,"%d", &compress);                                           //compressibility of grains
  // common material parameters
  xfscanf (in,"%le %le %le %le %le %le %le %d %d %d %d %d %d", &alpha0, &ks, &rhos0, &pw_bc, &t0, &tau0, &kintr0, &por_type, &kintr_type, &krw_type, &deff_type, &sr_type, &xi_type);

  switch (model_type){
  case lewis_and_schrefler:{//Lewis and Schrefler's model approach
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's model approach coupled with MEFEL
    xfscanf (in,"%le %d", &mefel_units, &vol_strain_effect);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  //porosity calculation:
  switch (por_type){
  case 0:{//constant
    xfscanf (in,"%le",&phi0);
    break;
  }
  case 1:{//dependent on mefel calculation
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 


  //intrinsic permebility calculation:
  switch (kintr_type){
  case 0:{//constant
    break;
  }
  case 1:{//dependent on porosity
    xfscanf (in,"%le %le", &b1, &phi01);
    break;
  }
  case 2:{//dependent on porosity - cubic and quadratic
    xfscanf (in,"%le", &phi01);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  //relative permebility calculation:
  switch (krw_type){
  case 0:{//constant
    xfscanf (in,"%le", &krw0);
    break;
  }
  case 1:
  case 2:{//dependent on saturation degree:
    xfscanf (in,"%le %le", &sirr, &ssat);
    break;
  }
  case 3:{//exponential
    xfscanf (in,"%le", &lambda_krw);
    break;
  }
  case 4:{//liakopoulos
    break;
  }
  case 5:{//double exponential
    xfscanf (in,"%le %le", &sirr, &ssat);
    xfscanf (in,"%le", &beta_krw);
    break;
  }
  case 6:{//FEBEX granit
    //xfscanf (in,"%le %le %le %le", &m1, &m2, &m3, &m4);
    break;
  }
  case 7:{//van_genuchten
    break;
  }
  case 9:{//bazant
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

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
    /* case gardner_exponential_sr:{//dependent on saturation degree:
       gardner_ret.read(in);
       break;
       }
       case potts_log_linear_sr:{//exponential
       potts_ret.read(in);
       break;
       }
       
       case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
       //van_genuchten_ret.read2(in);
       break;
       }
    */
  case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
    van_genuchten_ret.read(in);
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
}


/**
   function reads parameters
   
   @param in - input file

   13/11/2023, TKr
*/
void con_awf1mat::print(FILE *out)
{
}





/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf1mat::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
    print_err ("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}


/**
   function computes volume part 2 of right-hand side matrix
   in the required integration point
   
   @param cc - right-hand side coefficient of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf1mat::rhs_volume2 (double &cc,long /*ri*/,long /*ci*/,long ipp)
{
  double pw;
  
  pw = Tm->ip[ipp].av[0];
  cc = get_fwu(pw,ipp);
}



/**
   function creates volume right-hand side matrix of the material for 1D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf1mat::rhs1d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Tm->ip[ipp].av[0];
  f = get_fw1(pw,ipp);
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
void con_awf1mat::rhs2d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Tm->ip[ipp].av[0];
  
  f = get_fw1(pw,ipp);
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
void con_awf1mat::rhs3d1 (matrix &d,long /*ri*/,long /*ci*/,long ipp)
{
  double f,pw;
  
  pw = Tm->ip[ipp].av[0];
  
  f = get_fw1(pw,ipp);
  fillm(0.0,d);

  d[0][0] = f*Tp->gr[0];
  d[1][0] = f*Tp->gr[1];
  d[2][0] = f*Tp->gr[2];
}



/**
   function computes effective stress factor xi
   @param pw - water pressure
   @param ipp - number of integration point

   @retval xi - factor xi

   03/10/2023, TKr
*/
double con_awf1mat::get_xi(double pw, long ipp)
{
  double suc=0.0,sr=0.0,xi=0.0,temp=294.0;
  
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
    xi = masin_ret.psi(-pw,-dpw,e,temp);//positive value of suction
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
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param ipp - number of integration point

   @retval sw - degree of saturation

   03/10/2023, TKr
*/
double con_awf1mat::get_sw(double pw, long ipp)
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
    
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    sw = van_genuchten_ret.sw(pw,t0);
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model type 2
    //sw = van_genuchten_ret.sw2(pw);
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
   function computes partial derivative of degree of saturation with respect to pw, specific wter content
   @param pw - water pressure
   @param ipp - number of integration point

   @retval dsw_dpw - partial derivative of degree of saturation with respect to pw

   03/10/2023, TKr
*/
double con_awf1mat::get_dsw_dpw(double pw, long ipp)
{
  double dsw_dpw;
  dsw_dpw = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    dsw_dpw = -bazant_ret.dsat_dpc(-pw,t0);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    dsw_dpw = lewis_ret.dsw_dpw(pw);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    dsw_dpw = Tm->givenontransq(der_saturation_deg, ipp); //actual derivative of saturation degree
    dsw_dpw = dsw_dpw/mefel_units; //basic units = Pa
    break;
  }

  case table_sr:{//saturation degree and its derivative are obtained from table;
    dsw_dpw = data.getderiv (pw); //actual derivative of saturation degree
    break;
  }
    /*
      case gardner_exponential_sr:{//partially saturated medium = Exponential model, Gardner 
      dsw_dpw = gardner_ret.dsw_dpw(pw);
      break;
      }
      case potts_log_linear_sr:{//partially saturated medium = Log linear model, Potts
      dsw_dpw = potts_ret.dsw_dpw(pw);
      break;
      }
      case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
      //dsw_dpw = van_genuchten_ret.dsw_dpw2(pw);
      break;
      }
    */
  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    dsw_dpw = van_genuchten_ret.dsw_dpw(pw,t0);
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
    dsw_dpw = -masin_ret.dsw_dpw(-pw,-dpw,e,t0);//positive value of suction
    //dsw_dpw = dsw_dpw/mefel_units; //basic units = Pa
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    dsw_dpw = -febex_granit_ret.dsw_ds(-pw);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  return(dsw_dpw);
}



/**
   function computes intrinsic permeability
   @param pw - water pressure
   @param ipp - number of integration point
   
   @retval kintr - intrinsic permeability

   03/10/2023, TKr
*/
double con_awf1mat::get_kintr(double pw, long ipp)
{
  double kintr,phi;

  switch (kintr_type){
  case 0:{//constant
    kintr = kintr0;
    break;
  }
  case 1:{//dependent on porosity
    phi = get_porosity(ipp);
    kintr = kintr0*exp(bb1*(phi - phi01));
    break;
  }
  case 2:{//dependent on porosity - cubic and quadratic
    phi = get_porosity(ipp);
    kintr = kintr0*phi*phi*phi/(1 - phi01)*(1 - phi01)*(1 - phi)*(1 - phi)/(phi01*phi01*phi01);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(kintr);
}

/**
   function computes water relative permeability
   @param pw - water pressure
   @param ipp - number of integration point
   
   @retval krw - water relative permeability

   03/10/2023, TKr
*/
double con_awf1mat::get_krw(double pw, long ipp)
{
  double sw,krw,sef=0.0;
  
  switch (krw_type){
  case 0:{//constant
    krw = 1.0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,ipp);
    krw = (sw-sirr)/(ssat-sirr);
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,ipp);
    krw = (sw-sirr)*(sw-sirr)*(sw-sirr)/((ssat-sirr)*(ssat-sirr)*(ssat-sirr));
    break;
  }
  case 3:{//exponential
    sw = get_sw(pw,ipp);
    krw = pow(sw,lambda_krw);
    break;
  }
  case 4:{//liakopoulos
    sw = get_sw(pw,ipp);
    krw = 1.0-2.207*pow((1.0-sw),1.0121);
    break;
  }
  case 5:{//double exponential
    sw = get_sw(pw,ipp);
    sef = (sw-sirr)/(ssat-sirr); //effective saturation degree
    krw = pow(sef,(1.0/beta_krw));
    krw = pow((1.0-krw),beta_krw);
    krw = pow((1.0-krw),2.0);
    krw = pow(sef,0.5)*krw;
    break;
  }
  case 6:{//FEBEX granit
    sw = get_sw(pw,ipp);
    krw = febex_granit_ret.get_krw(sw);
    break;
  }
  case 7:{
    krw = van_genuchten_ret.get_krw(pw,t0);
  }
  case 9:{//bazant
    sw = get_sw(pw,ipp);
    krw = exp(10.0*log(sw));
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(krw);
}


/**
   function returns porosity
   @param ipp - number of integration point

   @retval phi - porosity

   03/04/2023, TKr
*/
double con_awf1mat::get_porosity(long ipp)
{
  double por;
  
  //porosity calculation:
  switch (por_type){
  case 0:{//constant
    por = phi0;
    break;
  }
  case 1:{//dependent on mefel calculation
    if(Tm->nontransq != NULL){
      por = Tm->givenontransq(porosity, ipp);// from mefel
    }
    else
      por = phi0;//for testing
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  
  return(por);
}


/**
   function returns Biot's constant

   @retval alpha - Biot's constant

   12/9/2008, TKr
*/
double con_awf1mat::get_alpha()
{
  return(alpha0);
}


/**
   function returns kw bulk modulus of water [Pa]

   @retval kw - bulk modulus of water [Pa]

   12/9/2008, TKr
*/
double con_awf1mat::get_kw()
{
  return(kw);
}


/**
   function returns muw water viscosity [Pa.s]

   @retval muw - water viscosity [Pa.s]

   12/9/2008, TKr
*/
double con_awf1mat::get_muw()
{
  return(muw0);
}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param ipp - number of integration point

   @retval kww - conductivity coefficient

   03/10/2023, TKr
*/
double con_awf1mat::get_kww(double pw, long ipp)
{
  double krw,kww,muw,kintr,dv,drhogw_dpw;
  
  //water pressure check
  waterpress_check(pw,ipp);
  
  krw = get_krw(pw,ipp);
  muw = get_muw();
  kintr = get_kintr(pw,ipp);
  drhogw_dpw = get_drhogw_dpw(pw);
  dv = get_dv(pw,ipp);
  
  kww = rhow*krw*kintr/muw0;//liquid water
  kww = kww + dv*drhogw_dpw;//water vapour

  return(kww);
}



/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param ipp - number of integration point

   @retval capww - capacity coefficient

   03/10/2023, TKr
*/
double con_awf1mat::get_capww(double pw, long ipp)
{
  double sw,ds_dpc,n,capww,rhogw,sg,dpgw_dpc,alpha;  

  //water pressure check
  waterpress_check(pw,ipp);

  
  n = get_porosity(ipp);
  sw = get_sw(pw,ipp);
  rhogw = get_rhogw(pw);
  sg = 1.0 - sw;
  dpgw_dpc = get_dpgw_dpc(pw);
  ds_dpc = -get_dsw_dpw(pw,ipp);//this is equal to minus derivative with respect to suction
  alpha = get_alpha();


  if(compress == 1){
    //compressible grains:
    capww = (alpha - n)/ks*sw*(rhogw*sg + rhow*sw) + rhow*sw*n/kw;
    capww = capww - sg*n*mw/t0/gasr*dpgw_dpc;
    capww = capww - ((alpha - n)/ks*(rhogw*sg*(pw) + rhow*sw*(pw)) + n*(rhow - rhogw))*ds_dpc;
  }
  else{
    //incompressible grains:
    capww = -sg*n*mw/t0/gasr*dpgw_dpc;
    capww = capww - n*(rhow - rhogw)*ds_dpc;    
  }

  return(capww);
}




/**
   function computes mass concentration of water vapour air in gas phase
   @param pw - pore water pressure

   @retval rhogw - mass concentration of water vapour air in gas phase

   03/10/2023, TKr
*/
double con_awf1mat::get_rhogw(double pw)
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

   03/10/2023, TKr
*/
double con_awf1mat::get_pgw(double pw)
{
  double pgw,pgws,pc;

  pgws = get_pgws(t0);

  pc = -pw;
  
  pgw = pgws*exp(-1.0*pc*mw/rhow/gasr/t0);

  check_math_err();

  return(pgw);
}




/**
   function computes partial derivative of rhogw with respect to pw (Kelvin equation)
   @param pw - pore water pressure

   @retval drhogw_dpw - partial derivative of pgw with respect to pw (Kelvin equation)

   03/10/2023, TKr
*/
double con_awf1mat::get_drhogw_dpw(double pw)
{
  double drhogw_dpw,rhogw;

  rhogw = get_rhogw(pw);
  
  drhogw_dpw = rhogw*mw/rhow/gasr/t0;

  return(drhogw_dpw);
}



/**
   function computes partial derivative of pgw with respect to pc (Kelvin equation)
   @param pw - pore water pressure

   @retval dpgw_dpc - partial derivative of pgw with respect to pc (Kelvin equation)

   03/10/2023, TKr
*/
double con_awf1mat::get_dpgw_dpc(double pw)
{
  double dpgw_dpc,pgws,pc;

  pgws = get_pgws(t0);
  
  pc = -pw;
  
  dpgw_dpc = -pgws*exp(-pc*mw/rhow/gasr/t0)*mw/rhow/gasr/t0;

  check_math_err();
  
  return(dpgw_dpc);
}



/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure

   03/10/2023, TKr
*/
double con_awf1mat::get_pgws(double t)
{
  double t1,t2,t3,pgws,psl;
  
  t1 = 1.0/t;
  t2 = t*t;
  t3 = t*t*t;

  psl = c8*t1 + c9 + c10*t + c11*t2 + c12*t3 + c13*log(t);
  pgws = exp(psl);

  return(pgws);
}




/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pw - pore water pressure
   @param ipp - number of integration point

   @retval deff - effective diffusion coefficient of vapour inside pores
*/
double con_awf1mat::get_dv(double pw,long ipp)
{
  double dv;
  double phi,sw,tau;

  dv = dv0;

  switch (deff_type){
  case 0:{//constant
    dv = dv0;
    break;
  }
  case 1:{//moisture dependent
    phi = get_porosity(ipp);
    sw = get_sw(pw,ipp);
    tau = tau0;
    dv = phi*(1.0 - sw)*tau*dv0*pow((t0/273.15),1.8);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(dv);
}




/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param ipp - number of integration point

   @retval fw - first part for right-hand side for continutiy equation


   03/10/2023, TKr
*/
double con_awf1mat::get_fw1(double pw, long ipp)
{
  double krw,muw,kintr,fw1;
  
  fw1 = 0.0;
  
  krw = get_krw(pw,ipp);
  muw = get_muw();
  kintr = get_kintr(pw,ipp);

  krw = get_krw(pw,ipp);
  fw1 = rhow*krw*kintr/muw0*rhow;//p. 89

  return(fw1);
}




/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param ipp - number of integration point

   @retval fwu - volumetric strain rate effect on pore water pressure; only for partially coupled METR version

   03/10/2023, TKr
*/
double con_awf1mat::get_fwu(double pw, long ipp)
{
  double depsv_r,sw,fwu,dsr_depsv,n,sg,rhogw,alpha;
  
  fwu = 0.0;
  
  switch (model_type){
  case lewis_and_schrefler2:{
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book
    if(vol_strain_effect == 1){
      n = get_porosity(ipp);
      depsv_r = Tm->givenontransq(strain_vol_rate, ipp);      //actual rate of volumetric strain from MEFEL
      Tm->ip[ipp].eqother[0] = depsv_r; //this is not necessary
      
      
      /*
	if (sr_type == mefel_sr){
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	Tm->ip[ipp].eqother[1] = dsr_depsv; //this is not necessary
	}
      */
      
      switch (sr_type){
      case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
	dsr_depsv = 0.0;
	break;
      }
      case gardner_exponential_sr:{//dependent on saturation degree:
	dsr_depsv = 0.0;
	break;
      }
      case potts_log_linear_sr:{//exponential
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case van_genuchten2_sr:{//partially saturated medium =Van Genuchten model
	dsr_depsv = 0.0;
	break;
      }
      case mefel_sr:{//from MEFEL
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	break;
      }
      case table_sr:{//from table
	dsr_depsv = 0.0;
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
	dsr_depsv = masin_ret.dsw_depsv(-pw,-dpw,e,t0);//positive value of suction
	break;
      }
      default:{
	print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
	abort();
      }
      } 

      
      sw = get_sw(pw,ipp);
      sg = 1.0 - sw;
      rhogw = get_rhogw(pw);
      alpha = get_alpha();
      
      fwu = -1.0*(sg*rhogw + sw*rhow)*alpha*depsv_r;//volumetric strain effect        
      
      fwu = fwu - n*(rhow - rhogw)*dsr_depsv*depsv_r;//volumetric strain effect
    }
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }    
  return(fwu);
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
double con_awf1mat::transmission_transcoeff(double trc,long /*ri*/,long /*ci*/,long nn,long bc,long ipp)
{
  double c,pw;
  c = 0.0;
  
  pw = nodalval (nn,0);
  
  //  if((ri == 0) && (ci == 0))
  c = get_transmission_transcoeff_ww(pw,bc,ipp);// *scale_pw;//scaling

  c = c*trc;

  return (c);
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
double con_awf1mat::transmission_transcoeff(double trc,long /*ri*/,long /*ci*/,long nn,long bc,long ipp, int flag)
{
  double c,pw;
  c = 0.0;
  
  pw = nodalval (nn,0);

  //  if((ri == 0) && (ci == 0))
  c = get_transmission_transcoeff_ww(pw,bc,ipp,flag);// *scale_pw;//scaling
  
  c = c*trc;
  
  return (c);
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
double con_awf1mat::transmission_nodval(double nodval,double /*trc2*/,long /*ri*/,long /*ci*/,long nn,long bc,long ipp)
{
  double c,pw;
  c = 0.0;
  
  pw = nodalval (nn,0);

  //if((ri == 0) && (ci == 0))
  c = get_transmission_nodval_ww(nodval,pw,bc,ipp);// *scale_pw;//scaling

  return (c);
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
double con_awf1mat::transmission_flux(double nodval,double /*trc2*/,long /*ri*/,long /*ci*/,long nn,long bc,long ipp)
{
  double c,pw;
  c = 0.0;
  
  pw = nodalval (nn,0);

  //if((ri == 0) && (ci == 0))
  c = get_transmission_flux_ww(nodval,pw,bc,ipp);// *scale_pw;//scaling

  return (c);
}



/**
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf1mat::get_transmission_transcoeff_ww(double pw,long bc,long ipp)
{
  double trc;
  
  //water pressure check
  waterpress_check(pw,ipp);

  //other conditions types will be added according to the neeed
  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure

    trc = 1.0;
    
    break;
  }
  case 40:{//simulation of free soil surface
    if(pw < pw_bc)
      trc = 0.0;
    else
      trc = 1.0e+03;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(trc);
}



/**
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf1mat::get_transmission_transcoeff_ww(double pw,long bc,long ipp,int flag)
{
  double trc;
  
  //water pressure check
  waterpress_check(pw,ipp);

  //other conditions types will be added according to the neeed
  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure
    
    trc = 1.0;
    
    break;
  }
  case 33:{//pore water pressure - flux
    if(flag == 1)
      trc=0.0;//into right hand side (matrix)
    else{
      trc=1.0;//into left hand side (flux)
    }
    break;
  }
  case 40:{//simulation of free soil surface
    if(pw < pw_bc)
      trc = 0.0;
    else
      trc = 1.0e+03;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(trc);
}

/**
   function creates correct new nodal value on the boundary (transmission) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf1mat::get_transmission_nodval_ww(double bv,double pw,long bc,long ipp)
{
  double new_nodval;
  
  //water pressure check
  waterpress_check(pw,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure
    new_nodval = bv;
    break;
  }
  case 33:{//pore water pressure - flux
    new_nodval = bv - pw;//minus sign
    break;
  }
  case 40:{//simulation of free soil surface
    new_nodval = bv;
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(new_nodval);
}



/**
   function creates flux on the boundary (transmission - convective mass transfer) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf1mat::get_transmission_flux_ww(double bv,double pw,long bc,long ipp)
{
  double flux,trc;
  
  //water pressure check
  waterpress_check(pw,ipp);

  switch (bc){//type of prescribed variable
  case 30:{//pore water pressure
    flux = bv - pw;//minus sign
    break;
  }
  case 33:{//pore water pressure - flux
    flux = bv - pw;//minus sign
    break;
  }
  case 40:{//simulation of free soil surface
    if(pw < pw_bc)
      trc = 0.0;
    else
      trc = 1.0e+03;
    flux = trc*(bv - pw);//minus sign
    break;
  }
  default:{
    print_err ("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }

  return(flux);
}




/**
   function computes all other variables at nodes
   @param compother - number of other components
   @param pw - water capillary pressure on actual node
   @param ipp - first integration point on element

   @retval other - other variable

   03/03/2011, TKr
*/

double con_awf1mat::get_othervalue(long compother,double pw, long ipp)
{
  double other;
  state_eq tt;

  switch (compother){
  case 0:{//capillary pressure
    other = -pw;
    break;
  }
  case 1:{//saturation
    other = get_sw(pw,ipp);
    break;
  }
  case 2:{//specific moisture content
    other = get_w(pw,ipp);
    break;
  }    
  default:{
    print_err("unknown type of component is required in function", __FILE__, __LINE__, __func__);
    abort();
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
void con_awf1mat::print_othervalue_name(FILE *out,long compother)
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
  case 2:{//moisture content
    fprintf (out,"Moisture (water) content (kg/kg)      ");
    //fprintf (out,"Moisture content (kg/kg)      ");
    break;
  }    
  default:{
    print_err ("unknown type of component is required in function", __FILE__, __LINE__, __func__);
    abort();
  }
  }
}



/**
   function checks if computed unknowns are physically reasonable
   @param nv  - vector of unknowns
   @param ipp - number of integration point


   28/02/2011, TKr
*/
void con_awf1mat::values_correction (vector &nv, long ipp)
{
  //  pore water pressure control
  waterpress_check(nv[0],ipp);
}


/**
   function checks if water pressure is non-positive
   @param pw - pore water pressure
   @param ipp - number of integration point


   28/02/2011, TKr
*/
void con_awf1mat::waterpress_check(double &/*pw*/,long /*ipp*/)
{
}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.

   24/08/2012, TKr
*/
void con_awf1mat::updateval (long /*ipp*/)
{
}


/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number

   24/08/2012, TKr
*/
void con_awf1mat::initval(long /*ipp*/)
{
}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Krejci according to Tomas Koudelka, 28/04/2014
*/
void con_awf1mat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_water_press;
}




/**
   function returns effective pore pressure (pressure has negative sign - mechanical convention)

   @param ipp - integration point number

   @retval pw - pore water pressure

   03/10/2023, TKr
*/
double con_awf1mat::give_effective_pore_pressure(long ipp)
{
  double pw,xi=1.0;

  pw = Tm->ip[ipp].av[0];
  
  switch (model_type){
  case lewis_and_schrefler:{//Lewis and Schrefler's book
    xi =  get_xi(pw,ipp);
    pw = -xi*pw;//returns effective pore pressure
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    xi = get_xi(pw,ipp); //effective stress factor
    pw = -xi*pw/mefel_units;//returns effective pore pressure //corrected units for mefel //basic units = Pa  
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }   
  return(pw);
}


/**
   function returns water pressure

   @param ipp - integration point number

   @retval pw - pore water pressure
*/
double con_awf1mat::give_water_pressure(long ipp)
{
  double pw;

  pw = Tm->ip[ipp].av[0];

  return(pw);
}


/**
   function returns pore pressure (pressure has negative sign - mechanical convention)

   @param ipp - integration point number

   @retval pw - pore water pressure
*/
double con_awf1mat::give_pore_pressure(long ipp)
{
  double pw;

  pw = -Tm->ip[ipp].av[0];

  return(pw);
}

/**
   function computes suction stress s = -pw = pc; for MEFEL stress notation s = pw = -pc

   @param ipp - integration point number

   @retval suction - suction stress [Pa]

   03/10/2023, TKr
*/
double con_awf1mat::give_suction(long ipp)
{
  double pw,suction;
  
  pw = Tm->ip[ipp].av[0];

  switch (model_type){
  case lewis_and_schrefler:{//Lewis and Schrefler's book
    suction = pw;//this is correct for Lewis and Schrefler's notation
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    suction = pw;//this is correct, because capillary water pressure is negative
    
    suction = suction/mefel_units;//corrected units for mefel //basic units = Pa    
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(suction);
}


/**
   function returns the degree of saturation
   @param ipp - integration point number

   @retval saturation degree [-]

   25/9/2012, TKr
*/
double con_awf1mat::give_saturation_degree(long ipp)
{
  double pw,s;
  
  pw = Tm->ip[ipp].av[0];
  s = get_sw(pw,ipp);
  
  return(s);
}

/**
   function computes water content w [kg/kg]
   @param pw - water capillary pressure

   @retval w - computes water content w [kg/kg]
*/
double con_awf1mat::get_w(double pw,long ipp)
{
  double w,n,s;

  n = get_porosity(ipp);
  s = get_sw(pw,ipp);

  w = (n*s*rhow)/(1.0 - n)/rhos0;

  return(w);
}


/**
  The funtion marks required non-transport quantities in the array antq.

  @param antq - array with flags for used material types
                antq[i] = 1 => quantity type nontransquant(i+1) is required
                antq[i] = 0 => quantity type nontransquant(i+1) is not required

  @return The function does not return anything, but it may change content of antq array.
*/
void con_awf1mat::give_reqntq(long *antq)
{
  switch (model_type){
  case lewis_and_schrefler:{//Lewis and Schrefler's book
    break;
  }
  case lewis_and_schrefler_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    antq[porosity-1] = 1;
    if (sr_type == mefel_sr){
      antq[saturation_deg-1] = 1;
      antq[der_saturation_deg-1] = 1;
      antq[der_saturation_deg_depsv-1] = 1;
      antq[der_saturation_deg_dtemp-1] = 1;
    }
    if(vol_strain_effect == 1)
      antq[strain_vol_rate-1] = 1;
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
}
