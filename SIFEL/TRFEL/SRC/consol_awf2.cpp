/*
    File:             consol_awf2.cpp
    Author:           Tomas Krejci, 29/10/2009, revised 03/10/2023
    Purpose:          material model for saturated-nonsaturated air and water flow in a deforming porous medium - advection of liquid water, vapour diffusion, and advection of gas (moist air)
    sources:          Lewis and Schrefler pp. 93-97, material parameters are set for benchmark on page n. 168
    unknowns:         number of unknowns=2, pw = liquid water pressure, pg = gas(air) pressure
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "constrel.h"
#include "consol_awf2.h"
#include "globalt.h"
#include "globmatt.h"

con_awf2mat::con_awf2mat()
{
  compress = 0;          //compressible grains: 0=no; 1=yes
  por_type = 0;          //porosity calculation type
  kintr_type = 0;        //intrinsic permability calculation type
  krw_type = 0;          //relative permeability calculation type
  krg_type = 0;          //relative permability calculation type
  deff_type = 0;         // diffusion calculation type
  sr_type = 1;           //retention curve calculation type

  vol_strain_effect = 0; //volumetric strain rate influence: 0=no; 1=yes 
  mefel_units = 1.0;//basic units for pressures = Pa (Pascals)

  //STATE VARIABLES
  mw = 18.01528; //molar mass of water kg.mol-1
  ma = 28.9645;   //molar mass of dry air kg.mol-1
  gasr = 8314.41; //universal gas constant J.mol-1.K-1

  t0 = 273.15;
  p0 = 101325.0;
  p_atm = 101325.0; //atmospheric pressure [Pa]
  p_atm_kpa = 101.325; //athmospheric pressure in kPa

  // PHYSICAL PROPERTIES OF DRY AIR
  muga0 = 17.17e-6;//[Pa.s]
  
  // PHYSICAL PROPERTIES OF WATER VAPOUR
  dv0 = 2.58e-5; //effective diffusion coefficient of vapour m^2/s at reference temperature 273.15 
  bv0 = 1.667;
  mug0 = 1.8e-5;//air viscosity [Pa.s]
  mugw0 = 8.85e-6;//[Pa.s]
  c8 = -5.8002206e+03;
  c9 = 1.3914993;
  c10 =-4.8640239e-02;
  c11 = 4.1764768e-05;
  c12 = -1.4452093e-08;
  c13 = 6.5459673;

  // PHYSICAL PROPERTIES OF WATER
  tcr = 647.3; //critical point of water [K]
  rhow0 = 1000.0;//water density [kg/m^3]
  muw0 = 1.0e-3;//water viscosity [Pa.s]
  kw0 = 2.0e9;//bulk modulus of water [Pa]

  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha0 = 1.0;//Biot's constant [-] alpha = 1 - kt/ks
  ks0 = .167e9;//bulk modulus of solid phase (grains) [Pa]
  kt0 = 0.0;
  phi0 = 0.297;//initial porosity [-]
  kintr0 = 4.5e-13;//intrinsic permeability [m^2]
  deff0 = 0.0;
  rhos0 = 2500.0; //solid density
  cdiff0 = 25.4e-6;//initial water vapour diffusivity in air
  tau0 = 0.0;
  sirr = 0.0;
  ssat = 1.0;

  bb1 = 0.0;
  phi01 = 0.0;

  pw_bc = 0.0; //free boundary pressure 
  rel_gas_press = 0;
}

con_awf2mat::~con_awf2mat()
{}


/**
   function reads parameters
   
   @param in - input file

   29/10/2009, TKr
*/
void con_awf2mat::read(XFILE *in)
{
  xfscanf (in,"%k%m","airwaterflowtype",&airwaterflowtype_kwdset, &model_type);
  xfscanf (in,"%d", &compress);

  // common material parameters
  xfscanf (in,"%le %le %le %le %le %le %le %d %d %d %d %d %d %d", &alpha0, &ks0, &rhos0, &pw_bc, &t0, &tau0, &kintr0, &por_type, &kintr_type, &krw_type, &krg_type, &deff_type, &sr_type, &xi_type);
  
  switch (model_type){
  case lewis_and_schrefler2:{//Lewis and Schrefler's model approach
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's approach coupled with mechanics
    xfscanf (in,"%le %d", &mefel_units, &vol_strain_effect);
    break;
  }
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
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
    xfscanf (in,"%le %le", &bb1, &phi01);
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
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  //gas relative permebility calculation:
  switch (krg_type){
  case 0:{//constant
    xfscanf (in,"%le", &krg0);
    break;
  }
  case 1:{//dependent on saturation degree:
    xfscanf (in,"%le %le", &s_crit, &ag);
    break;
  }
  case 2:{//dependent on saturation degree:
    xfscanf (in,"%le %le", &s_crit, &ag);
    break;
  }
  case 3:{//exponential dependence on void ratio and saturation for Kg
    xfscanf (in,"%le %le", &kg0, &kgn); //not finished
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
}


/**
   function prints parameters
   
   @param out - output file

   29/10/2009, TKr
*/
void con_awf2mat::print(FILE */*out*/)
{
  /*  fprintf (out,"\n %d ", int(model_type));
      fprintf (out,"\n %d ", compress);
      
      //opravit??!!!
      
      switch (model_type){
      case lewis_and_schrefler2:{//Lewis and Schrefler's book
      fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le\n",alpha0,ks0,phi0,kw0,rhow0,muw0,mug0,kintr0,rhos0,t0);
      lewis_ret.print(out);
      break;
      }
      case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book
      fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le %le %le\n",mefel_units,alpha0,ks0,phi0,kw0,rhow0,muw0,mug0,kintr0,rhos0,t0,pw_bc);
      break;
      }
      case van_genuchten2:{//partially saturated medium =Van Genuchten model
      fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le %le %le\n", mefel_units,alpha0,ks0,phi0,kw0,rhow0,muw0,mug0,kintr0,rhos0,t0,pw_bc);
      van_genuchten_ret.print(out);
      break;
      }
      default:{
      print_err("unknown model type is required", __FILE__, __LINE__, __func__);
      abort();
      }
      }
  */
}




/**
   function computes conductivity matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point

*/
void con_awf2mat::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
    abort();
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
void con_awf2mat::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,ipp);
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,ipp);

  d[0][0] = k;
}

/**
   function creates conductivity matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2mat::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,ipp);
  check_math_errel(0);      
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,ipp);
  check_math_errel(0);      
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,ipp);
  check_math_errel(0);      
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,ipp);
  check_math_errel(0);      
  
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
void con_awf2mat::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  double pw,pg;
  k = 0.0;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    k = get_kww(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    k = get_kwg(pw,pg,ipp);
  if((ri == 1) && (ci == 0))
    k = get_kgw(pw,pg,ipp);
  if((ri == 1) && (ci == 1))
    k = get_kgg(pw,pg,ipp);
  
  fillm(0.0,d);
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}


/**
   function creates capacity %matrix of the material

   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2mat::matcap (double &c,long ri,long ci,long ipp)
{
  double pw,pg;
  c = 0.0;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if((ri == 0) && (ci == 0))
    c = get_capww(pw,pg,ipp);
  if((ri == 0) && (ci == 1))
    c = get_capwg(pw,pg,ipp);
  if((ri == 1) && (ci == 0))
    c = get_capgw(pw,pg,ipp);
  if((ri == 1) && (ci == 1))
    c = get_capgg(pw,pg,ipp);

}


/**
   function computes volume part of right-hand side matrix
   in the required integration point
   
   @param d - right-hand side %matrix of material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2mat::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
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
void con_awf2mat::rhs_volume2 (double &cc,long ri,long /*ci*/,long ipp)
{
  double pw,pg;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];

  if(ri == 0){
    cc = get_fwu(pw,pg,ipp);
  }
  if(ri == 1){
    cc = get_fgu(pw,pg,ipp);
  }
}




/**
   function creates volume right-hand side matrix of the material for 1D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2mat::rhs1d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if(ri == 0){
    f = get_fw1(pw,pg,ipp);
    
    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
  if(ri == 1){
    f = get_fg1(pw,pg,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
  }
}

/**
   function creates volume right-hand side matrix of the material for 2D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2mat::rhs2d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if(ri == 0){
    f = get_fw1(pw,pg,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
  if(ri == 1){
    f = get_fg1(pw,pg,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
  }
}

/**
   function creates volume right-hand side matrix of the material for 3D problems
   @param d - right-hand %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - number of integration point
*/
void con_awf2mat::rhs3d1 (matrix &d,long ri,long /*ci*/,long ipp)
{
  double f,pw,pg;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  
  if(ri == 0){
    f = get_fw1(pw,pg,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
  if(ri == 1){
    f = get_fg1(pw,pg,ipp);

    fillm(0.0,d);
    d[0][0] = f*Tp->gr[0];
    d[1][0] = f*Tp->gr[1];
    d[2][0] = f*Tp->gr[2];
  }
}


/**
   function computes degree of saturation(sorption curve)
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval sw - degree of saturation

   03/10/2023, TKr
*/
double con_awf2mat::get_sw(double pw, double /*pg*/, long ipp)
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
   function computes specific moisture content (partial derivative of degree of saturation with respect to pc)
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval dsw_dpc - specific moisture content = partial derivative of degree of saturation with respect to pc


   03/10/2023, TKr
*/
double con_awf2mat::get_ds_dpc(double pw, double pg, long ipp)
{
  double dsw_dpc,pc;
  dsw_dpc = 0.0;

  switch (sr_type){
  case bazant_sr:{//Bazant
    pc = pg - pw;
    dsw_dpc = bazant_ret.dsat_dpc(pc,t0);
    break;
  }
  case lewis_and_schrefler_sr:{//Lewis and Schrefler's book
    dsw_dpc = -lewis_ret.dsw_dpw(pw);
    break;
  }
  case mefel_sr:{//saturation degree ind its derivative are obtained from mefel;
    dsw_dpc = -Tm->givenontransq(der_saturation_deg, ipp); //actual derivative of saturation degree
    dsw_dpc = dsw_dpc/mefel_units; //basic units = Pa
    break;
  }

  case table_sr:{//saturation degree and its derivative are obtained from table;
    dsw_dpc = -data.getderiv (pw); //actual derivative of saturation degree
    break;
  }

  case van_genuchten_sr:{//partially saturated medium = Van Genuchten model
    dsw_dpc = -van_genuchten_ret.dsw_dpw(pw,t0);
    break;
  }
  case van_genuchten2_sr:{//partially saturated medium = Van Genuchten model
    //dsw_dpc = -van_genuchten_ret.dsw_dpw2(pw);
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
    
    dsw_dpc = masin_ret.dsw_dpw(pc,dpc,e,t0);//positive value of suction
    break;
  }
  case febex_granit_sr:{//FEBEX granit
    //pg = pg - p0;//pore gas pressure without atmospheric pressure//debug??
    pc = pg - pw;
    
    dsw_dpc = febex_granit_ret.dsw_ds(pc);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
  return(dsw_dpc);
}




/**
   function computes effective stress factor xi
   @param pw - water pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval xi - factor xi

   03/10/2023, TKr
*/
double con_awf2mat::get_xi(double pw, double pg, long ipp)
{
  double suc=0.0,sr=0.0,xi=0.0;
  
  switch (xi_type){
  case biot_xi:
    xi = get_sw(pw,pg,ipp);
    break;
  case biot_reduced_xi:
    //xi = gamma*get_sw(pw,ipp);
    sr = get_sw(pw,pg,ipp);
    //xi = pow(sr,(gamma/lambda0));
    xi = pow(sr,gamma);
    xi = (1-gamma)*xi;
    break;
  case biot_masin_xi:{//according to masin for testing
    sr = get_sw(pw,pg,ipp);
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

    xi = masin_ret.psi(pc,dpc,e,t0);//positive value of suction
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
   function computes water relative permeability
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval krw - water relative permeability

   03/10/2023, TKr
*/
double con_awf2mat::get_krw(double pw, double pg, long ipp)
{
  double pc,sw,sef,krw;
  krw = 1.0;

  switch (krw_type){
  case 0:{//constant
    krw = krw0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,pg,ipp);
    krw = (sw-sirr)/(ssat-sirr);
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,pg,ipp);
    krw = (sw-sirr)*(sw-sirr)*(sw-sirr)/((ssat-sirr)*(ssat-sirr)*(ssat-sirr));
    break;
  }
  case 3:{//exponential
    sw = get_sw(pw,pg,ipp);
    krw = pow(sw,lambda_krw);
    break;
  }
  case 4:{//liakopoulos
    sw = get_sw(pw,pg,ipp);
    krw = 1.0-2.207*pow((1.0-sw),1.0121);
    break;
  }
  case 5:{//double exponential
    sw = get_sw(pw,pg,ipp);
    sef = (sw-sirr)/(ssat-sirr); //effective saturation degree
    krw = pow(sef,(1.0/beta_krw));
    krw = pow((1.0-krw),beta_krw);
    krw = pow((1.0-krw),2.0);
    krw = pow(sef,0.5)*krw;
    break;
  }
  case 6:{//FEBEX granit
    sw = get_sw(pw,pg,ipp);
    krw = febex_granit_ret.get_krw(sw);
    break;
  }
  case 7:{
    pg = pg - p0;//pore gas pressure without atmospheric pressure
    pc = pg - pw;
    
    krw = van_genuchten_ret.get_krw(-pc,t0);
  }
  case 9:{//bazant
    sw = get_sw(pw,pg,ipp);
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
double con_awf2mat::get_porosity(long ipp)
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
   function computes air relative permeability
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval krg - gas relative permeability

   03/10/2023, TKr
*/
double con_awf2mat::get_krg(double pw, double pg, long ipp)
{
  double krg,sw,n,kg,mug,g,rhog,kintr,e;
  krg=1.0;

  switch (krg_type){
  case 0:{//constant
    krg = krg0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,pg,ipp);
    krg = 1.0 - pow((sw/s_crit),ag);
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,pg,ipp);
    krg = pow((1.0-sw),ag);
    break;
  }
  case 3:{//exponential dependence of Kg on void ration and saturation degree - FEBEX bentonite
    sw = get_sw(pw,pg,ipp);
    n = get_porosity(ipp);
    e = n/(1.0-n);
    mug = get_mug(pw,pg,t0);
    g = 9.81;
    rhog = get_rhog(pw,pg,t0);
    kintr = get_kintr(pw,pg,ipp);
    kg = kg0*pow((e*(1.0-sw)),kgn);
    krg = kg*mug/(g*rhog*kintr);
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(krg);
}


/**
   function computes intrinsic permeability
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval kintr - intrinsic permeability

   03/10/2023, TKr
*/
double con_awf2mat::get_kintr(double pw, double pg, long ipp)
{
  double kintr,phi;
  kintr = kintr0;

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
    //Kozeny's approach for bentonite:
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
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval kww - conductivity coefficient

   03/10/2023, TKr
*/
double con_awf2mat::get_kww(double pw, double pg, long ipp)
{
  double kww;
  double rhow,krw,muw,rhog,mg,dg,dpgw_dpc,kintr;

  rhow = get_rhow(t0);
  krw = get_krw(pw,pg,ipp);
  muw = get_muw(t0);
  rhog =get_rhog(pw,pg,t0);
  mg = get_mg(pw,pg,t0);
  dg = get_dg(pw,pg,t0,ipp);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  kintr = get_kintr(pw,pg,ipp);

  kww = rhow*kintr*krw/muw;//water

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  kww = kww - rhog*ma*mw/mg/mg*dg/pg*dpgw_dpc;//water vapour

  return(kww);
}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval kwg - conductivity coefficient

   03/10/2023, TKr
*/
double con_awf2mat::get_kwg(double pw, double pg, long ipp)
{
  double kwg;
  double rhogw,krg,mug,dg,pgw,rhow,krw,muw,rhog,mg,dpgw_dpc,kintr;
  
  rhogw = get_rhogw(pw,pg,t0);
  krg = get_krg(pw,pg,ipp);
  mug = get_mug(pw,pg,t0);
  dg = get_dg(pw,pg,t0,ipp);
  pgw = get_pgw(pw,pg,t0);
  rhow = get_rhow(t0);
  krw = get_krw(pw,pg,ipp);
  muw = get_muw(t0);
  rhog = get_rhog(pw,pg,t0);
  mg = get_mg(pw,pg,t0);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  kintr = get_kintr(pw,pg,ipp);

  kwg = rhogw*kintr*krg/mug;//air

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  kwg = kwg + rhog*ma*mw/mg/mg*dg*(dpgw_dpc/pg - pgw/pg/pg);//water vapour

  return(kwg);
}

/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval kgw - conductivity coefficient

   29/10/2009, TKr
*/
double con_awf2mat::get_kgw(double pw, double pg, long ipp)
{
  double kgw,rhog,dg,mg,dpgw_dpc;
  
  rhog = get_rhog(pw,pg,t0);
  dg = get_dg(pw,pg,t0,ipp);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  mg = get_mg(pw,pg,t0);
  
  if(rel_gas_press == 1)
    pg = pg + p_atm;
  kgw = rhog*ma*mw/mg/mg*dg/pg*dpgw_dpc;//water vapour

  return(kgw);
}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval kgg - conductivity coefficient

   29/10/2009, TKr
*/
double con_awf2mat::get_kgg(double pw, double pg, long ipp)
{
  double kgg,rhoga,krg,mug,rhog,dg,mg,dpgw_dpc,pgw,kintr;
  
  rhoga = get_rhoga(pw,pg,t0);
  krg = get_krg(pw,pg,ipp);
  mug = get_mug(pw,pg,t0);
  rhog = get_rhog(pw,pg,t0);
  dg = get_dg(pw,pg,t0,ipp);
  mg = get_mg(pw,pg,t0);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  pgw = get_pgw(pw,pg,t0);
  kintr = get_kintr(pw,pg,ipp);

  kgg = rhoga*kintr*krg/mug;//air

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  kgg = kgg - rhog*ma*mw/mg/mg*dg*(dpgw_dpc/pg-pgw/pg/pg);//water vapour

  return(kgg);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval capww - capacity coefficient

   03/10/2023, TKr
*/
double con_awf2mat::get_capww(double pw, double pg, long ipp)
{
  double capww;
  double alpha,n,ks,sw,rhogw,sg,rhow,kw,dpgw_dpc,ds_dpc;

  alpha = get_alpha(pw,pg,ipp);
  n = get_porosity(ipp);
  ks = get_ks(pw,pg,ipp);
  sw = get_sw(pw,pg,ipp);
  rhogw = get_rhogw(pw,pg,t0);
  sg = 1.0 - sw;
  rhow = get_rhow(t0);
  kw = get_kw(pw,pg,ipp);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  ds_dpc = get_ds_dpc(pw,pg,ipp);//this is equal to minus derivative with respect to suction
  
  if(rel_gas_press == 1)
    pg = pg + p_atm;

  if(compress == 1){
    //compressible grains:
    capww = (alpha - n)/ks*sw*(rhogw*sg + rhow*sw) + rhow*sw*n/kw;
    capww = capww - sg*n*mw/t0/gasr*dpgw_dpc;
    capww = capww - ((alpha - n)/ks*(rhogw*sg*(pw-pg) + rhow*sw*(pw-pg)) + n*(rhow - rhogw))*ds_dpc;
  }
  else{
    //incompressible grains:
    capww = -sg*n*mw/t0/gasr*dpgw_dpc;
    capww = capww - n*(rhow - rhogw)*ds_dpc;    
  }
  return(capww);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval capwg - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2mat::get_capwg(double pw, double pg, long ipp)
{
  double capwg;
  double alpha,n,ks,sw,rhogw,sg,rhow,kw,dpgw_dpc,ds_dpc;

  alpha = get_alpha(pw,pg,ipp);
  n = get_porosity(ipp);
  ks = get_ks(pw,pg,ipp);
  sw = get_sw(pw,pg,ipp);
  rhogw = get_rhogw(pw,pg,t0);
  sg = 1.0 - sw;
  rhow = get_rhow(t0);
  kw = get_kw(pw,pg,ipp);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  ds_dpc = get_ds_dpc(pw,pg,ipp);//this is equal to minus derivative with respect to suction
  
  if(rel_gas_press == 1)
    pg = pg + p_atm;

  if(compress == 1){
    //compressible grains:
    capwg = (alpha - n)/ks*sg*(rhogw*sg + rhow*sw);
    capwg = capwg + sg*n*mw/t0/gasr*dpgw_dpc;
    capwg = capwg + ((alpha - n)/ks*(rhogw*sg*(pw-pg) + rhow*sw*(pw-pg)) + n*(rhow - rhogw))*ds_dpc;
  }    
  else{
    //incompressible grains:
    capwg = sg*n*mw/t0/gasr*dpgw_dpc;
    capwg = capwg + n*(rhow - rhogw)*ds_dpc;    
  }

  return(capwg);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval capgw - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2mat::get_capgw(double pw, double pg, long ipp)
{
  double capgw,alpha,n,ks,sw,sg,pc,rhoga,ds_dpc,dpgw_dpc;
  
  alpha = get_alpha(pw,pg,ipp);
  n = get_porosity(ipp);
  ks = get_ks(pw,pg,ipp);
  sw = get_sw(pw,pg,ipp);
  sg = 1.0 - sw;
  pc = pg - pw;
  rhoga = get_rhoga(pw,pg,t0);
  ds_dpc = get_ds_dpc(pw,pg,ipp);
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);//this is equal to minus derivative with respect to suction
  
  if(rel_gas_press == 1)
    pg = pg + p_atm;

  if(compress == 1){
    //compressible grains:
    capgw = (alpha - n)/ks*sw*sg*rhoga;
    capgw = capgw + ((alpha - n)/ks*sg*pc + n)*rhoga*ds_dpc;
    capgw = capgw + sg*n*mw/t0/gasr*dpgw_dpc;
  }
  else{
    //incompressible grains:
    capgw =  n*rhoga*ds_dpc;
    capgw = capgw + sg*n*mw/t0/gasr*dpgw_dpc;
  }

  return(capgw);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval capgg - capacity coefficient

   29/10/2009, TKr
*/
double con_awf2mat::get_capgg(double pw, double pg, long ipp)
{
  double capgg,alpha,n,ks,sw,sg,pc,rhoga,ds_dpc,dpgw_dpc;
  
  alpha = get_alpha(pw,pg,ipp);
  n = get_porosity(ipp);
  ks = get_ks(pw,pg,ipp);
  sw = get_sw(pw,pg,ipp);
  sg = 1.0 - sw;
  rhoga = get_rhoga(pw,pg,t0);
  pc = pg - pw;
  ds_dpc = get_ds_dpc(pw,pg,ipp);//this is equal to minus derivative with respect to suction
  dpgw_dpc = get_dpgw_dpc(pw,pg,t0);
  
  if(rel_gas_press == 1)
    pg = pg + p_atm;

  if(compress == 1){
    //compressible grains:
    capgg = (alpha - n)/ks*sg*sg*rhoga;
    capgg = capgg - ((alpha - n)/ks*sg*pc + n)*rhoga*ds_dpc;
    capgg = capgg + sg*n*ma/t0/gasr;
    capgg = capgg - sg*n*mw/t0/gasr*dpgw_dpc;
  }
  else{
    //incompressible grains:
    capgg = -n*rhoga*ds_dpc;
    capgg = capgg + sg*n*ma/t0/gasr;
    capgg = capgg - sg*n*mw/t0/gasr*dpgw_dpc;    
  }
  return(capgg);
}

/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval fw - first part for right-hand side coeficient for continutiy equation

   29/10/2009, TKr
*/
double con_awf2mat::get_fw1(double pw, double pg, long ipp)
{
  double fc1;
  double rhow,rhogw,rhog,krw,muw,krg,mug,kintr;
  
  rhow = get_rhow(t0);
  rhogw = get_rhogw(pw,pg,t0);
  rhog = get_rhog(pw,pg,t0);  
  krw = get_krw(pw,pg,ipp);
  muw = get_muw(t0);
  krg = get_krg(pw,pg,ipp);
  mug = get_mug(pw,pg,t0);
  kintr = get_kintr(pw,pg,ipp);

  fc1 = rhogw*kintr*krg*rhog/mug + rhow*kintr*krw*rhow/muw;
  
  return(fc1);
}


/**
   function returns coefficient for righ-hand side of the general material 
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval fg - first part for right-hand side for continutiy equation

   29/10/2009, TKr
*/
double con_awf2mat::get_fg1(double pw, double pg, long ipp)
{
  double fg1;
  double rhoga,rhog,krg,mug,kintr;
  
  rhoga = get_rhoga(pw,pg,t0);
  rhog = get_rhog(pw,pg,t0);
  krg = get_krg(pw,pg,ipp);
  mug = get_mug(pw,pg,t0);
  kintr = get_kintr(pw,pg,ipp);
  
  fg1 = rhoga*kintr*krg*rhog/mug;
  fg1 = 0.0;//no gravity force is included for the moist air

  return(fg1);
}


/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval fwu - volumetric strain rate effect on pore water pressure; only for partially coupled METR version

   03/10/2023, TKr
*/
double con_awf2mat::get_fwu(double pw, double pg, long ipp)
{
  double depsv_r,sw,fwu,dsr_depsv,n,sg,rhow,rhogw,alpha;
  
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

      
      sw = get_sw(pw,pg,ipp);
      sg = 1.0 - sw;
      rhow = get_rhow(t0);
      rhogw = get_rhogw(pw,pg,t0);
      alpha = get_alpha(pw,pg,ipp);
      
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
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure
   @param ipp - number of integration point

   @retval fgu - volumetric strain rate effect on pore gas pressure; only for partially coupled METR version

   16/05/2018, TKr
*/
double con_awf2mat::get_fgu(double pw, double pg, long ipp)
{
  double depsv_r,sw,fgu,dsr_depsv,n,sg,rhoga,alpha;
  
  fgu = 0.0;
  
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

      sw = get_sw(pw,pg,ipp);
      sg = 1.0 - sw;
      rhoga = get_rhoga(pw,pg,t0);
      alpha = get_alpha(pw,pg,ipp);
      fgu = -alpha*sg*rhoga*depsv_r;
      
      fgu = fgu + n*rhoga*dsr_depsv*depsv_r;//volumetric strain rate effect
    }
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }    
  return(fgu);
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
double con_awf2mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,pg,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_transcoeff_gg(pw,pg,bc,ipp);// *scale_pg;//scaling

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
double con_awf2mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp, int flag)
{
  double c,pw,pg;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_transcoeff_ww(pw,pg,bc,ipp,flag);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_transcoeff_gg(pw,pg,bc,ipp);// *scale_pg;//scaling
  
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
double con_awf2mat::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_nodval_ww(nodval,pw,pg,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_nodval_gg(nodval,trc2,pw,pg,bc,ipp);// *scale_pg;//scaling

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
double con_awf2mat::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double c,pw,pg;
  c = 0.0;
  
  pw = nodalval (nn,0);
  pg = nodalval (nn,1);

  if((ri == 0) && (ci == 0))
    c = get_transmission_flux_ww(nodval,pw,pg,bc,ipp);// *scale_pw;//scaling
  if((ri == 0) && (ci == 1))
    c = 0.0;
  
  if((ri == 1) && (ci == 0))
    c = 0.0;
  if((ri == 1) && (ci == 1))
    c = get_transmission_flux_gg(nodval,trc2,pw,pg,bc,ipp);// *scale_pg;//scaling

  return (c);
}



/**
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure on the boundary
   @param pg - pore gas pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_transcoeff_ww(double pw,double /*pg*/,long bc,long /*ipp*/)
{
  double trc;
  
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
    print_err("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(trc);
}



/**
   function creates correct transfer coefficient on the boundary (transmission) for pw

   @param pw - pore water pressure on the boundary
   @param pg - pore gas pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_transcoeff_ww(double pw,double /*pg*/,long bc,long /*ipp*/,int flag)
{
  double trc;
  
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
   @param pg - actual pore gas pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_nodval_ww(double bv,double pw,double /*pg*/,long bc,long /*ipp*/)
{
  double new_nodval;
  
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
    print_err("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }
  return(new_nodval);
}



/**
   function creates flux on the boundary (transmission - convective mass transfer) for pw

   @param bv - prescribed value near the boundary
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_flux_ww(double bv,double pw,double /*pg*/,long bc,long /*ipp*/)
{
  double flux,trc;
  
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
    print_err("no real boundary condition is prescribed", __FILE__, __LINE__, __func__);
    exit(0);
  }
  }

  return(flux);
}


/**
   function creates correct transfer coefficient on the boundary (transmission) for t medium

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_transcoeff_gg(double /*pw*/,double /*pg*/,long bc,long /*ipp*/)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//gas transmission
    trc = 1.0;
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
   function creates correct new nodal value on the boundary (transmission) for t medium

   @param bv - value of prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_nodval_gg(double bv,double /*trr*/,double /*pw*/,double /*pg*/,long bc,long /*ipp*/)
{
  double new_nodval;

  switch (bc){//type of prescribed variable
  case 30:{//gas transmission
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
   function creates flux on the boundary (transmission - convective mass transfer) for c medium

   @param bv - prescribed value near the boundary
   @param trr - trr coefficient
   @param pw - actual pore water pressure on the boundary
   @param pg - actual pore gas pressure on the boundary
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
*/
double con_awf2mat::get_transmission_flux_gg(double bv,double /*trr*/,double /*pw*/,double pg,long bc,long /*ipp*/)
{
  double flux;

  switch (bc){//type of prescribed variable
  case 30:{//gas transmission - boundary flux
    flux = (bv - pg);//minus sign
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
   function computes all other variables at nodes
   @param compother - number of other components
   @param pw - water capillary pressure on actual node
   @param pg - gas(air) pressure on actual node
   @param ipp - first integration point on element

   @retval other - other variable

   03/03/2011, TKr
*/

double con_awf2mat::get_othervalue(long compother,double pw, double pg, long ipp)
{
  double other;
  state_eq tt;

  switch (compother){
  case 0:{//capillary pressure
    other = pg - pw;
      break;
  }
  case 1:{//gas pressure
    other = pg;
    break;
  }
  case 2:{//saturation
    other = get_sw(pw,pg,ipp);
    break;
  }
  case 3:{//liquid water pressure
    other = pw;
    break;
  }
  case 4:{//moisture content
    other = get_w(pw,pg,ipp);
    break;
  }    
  case 5:{//suction
    other = give_suction(ipp);
    break;
  }    
  case 6:{//rel.hum
    other = get_pgw(pw,pg,t0)/get_pgws(t0);
    break;
  }    
  case 7:{//dsw_dpw
    other = -get_ds_dpc(pw,pg,ipp);
    break;
  }    
  case 8:{//depsv
    other = Tm->ip[ipp].eqother[1];
    break;
  }   

  default:{
    print_err (" unknown type of component is required in function (%s, line %d).\n",__FILE__,__LINE__,__func__);
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
void con_awf2mat::print_othervalue_name(FILE *out,long compother)
{
  switch (compother){
  case 0:{//capillary pressure
    fprintf (out,"Capillary pressure (Pa)");
    break;
  }
  case 1:{//gas pressure
    fprintf (out,"Gas pressure (Pa)             ");
    break;
  }
  case 2:{//saturation
    fprintf (out,"Degree of saturation ()       ");
    break;
  }
  case 3:{//liquid water pressure
    fprintf (out,"Pore water pressure (Pa)      ");
    break;
  }
  case 4:{//moisture content
    fprintf (out,"Moisture (water) content (kg/kg)      ");
    //fprintf (out,"Moisture content (kg/kg)      ");
    break;
  }    
  case 5:{//capillary pressure
    fprintf (out,"Suction (Pa)");
    break;
  }
  case 6:{//rel. hum
    fprintf (out,"Relative humidity (-) ");
    break;
  }    
  case 7:{//dS_dpw
    fprintf (out,"dS_dpw Derivative of Saturation degree with respect to pore water pressure      ");
    break;
  }    
  case 8:{//depsv
    fprintf (out,"Volumetric strain rate (-)       ");
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


   03/03/2011, TKr
*/
void con_awf2mat::values_correction (vector &nv, long ipp)
{
  //  pore water pressure control
  waterpress_check(nv[0],nv[1],ipp);
  
  //  gas pressure
  gaspress_check(nv[0],nv[1],ipp);
}


/**
   function checks if gas pressure is greater than vapour pressure

   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param ipp - number of integration point

   26/07/2017 TKr
*/
void con_awf2mat::gaspress_check(double pw,double &pg,long ipp)
{
  double pgw;
  
  //general
  //gas pressure check
  pgw = get_pgw(pw,pg,t0);
  if (pgw > pg){
    pg = pgw;
    //fprintf(Outt,"gas pressure was modified in integration point No. %ld",ipp);
  }

  //general
  //fully saturated state
  if(pgw < 100.0){
    pg = 101325.0;
  }

  // M. Starnoni 24-11-2010
  if (pg < 0.0){
    print_err("gas pressure was modified in integration point No. %ld",__FILE__,__LINE__,__func__,ipp);
    pg = 0.0;
  }
}


/**
   function checks if water pressure is non-positive
   @param pw - pore water pressure
   @param ipp - number of integration point


   03/03/2011, TKr
*/
void con_awf2mat::waterpress_check(double &/*pw*/,double /*pg*/,long /*ipp*/)
{
}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.

   23/05/2016, TKr
*/
void con_awf2mat::updateval (long /*ipp*/)
{
}


/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number

   23/05/2016, TKr
*/
void con_awf2mat::initval(long /*ipp*/)
{
}



/**
  The function returns ordered dof names of primary unknowns 
  required by the model.

  @param dofname   - array of uknown name for particular nodal dofs (output)
                     dofname[i] = name of i-th nodal unknown (for names see aliast.h - enum namevart)
  @param ntm       - number of transported media = number of nodal dof = length of array dofname

  Created by Tomas Krejci according to Tomas Koudelka, 23/05/2016
*/
void con_awf2mat::give_dof_names(namevart *dofname, long ntm)
{
  if (ntm < 1)
  {
    print_err("the model defines %ld unknowns but number of transported media is %ld", 
              __FILE__, __LINE__, __func__, 1, ntm);
    abort();
  }
  dofname[0] = trf_water_press;
  dofname[1] = trf_gas_press;
}



/**
   function returns effective pore pressure (pressure has negative sign - mechanical convention)

   @param ipp - integration point number

   @retval pw - pore water pressure

   03/10/2023, TKr
*/
double con_awf2mat::give_effective_pore_pressure(long ipp)
{
  double ps,pw,pg,t,xi=1.0;

  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];
  t = Tm->ip[ipp].av[2];


  //needs to be corrected??!!
  switch (model_type){
  case lewis_and_schrefler3:{//Lewis and Schrefler's book
    xi =  get_xi(pw,pg,ipp);//effective stress factor

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    ps = xi*pw + (1.0 - xi)*pg;
    ps = -ps;//returns effective pore pressure
    break;
  }
  case lewis_and_schrefler3_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;

    xi = get_xi(pw,pg,ipp); //effective stress factor

    pg = pg - p0;//pore gas pressure without atmospheric pressure
    ps = xi*pw + (1.0 - xi)*pg;
    ps = -ps/mefel_units;//returns effective pore pressure //corrected units for mefel //basic units = Pa  
    
    //debug???!!!
    //pw = -pw/mefel_units;//returns effective pore pressure //corrected units for mefel //basic units = Pa  
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }   
  return(ps);
}



/**
   function returns water pressure

   @param ipp - integration point number

   @retval pw - pore water pressure

   27/05/2016, TKr
*/
double con_awf2mat::give_water_pressure(long ipp)
{
  double pw;

  pw = Tm->ip[ipp].av[0];

  return(pw);
}

/**
   function returns pore pressure

   @param ipp - integration point number

   @retval pp - pore pressure

   27/05/2016, TKr
*/
double con_awf2mat::give_pore_pressure(long ipp)
{
  double pp,pg,pw;

  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];

  pp = pg - pw;

  return(pp);
}

/**
   function returns gas pressure

   @param ipp - integration point number

   27/05/2016, TKr
*/
double con_awf2mat::give_gas_pressure(long ipp)
{
  double pg;

  pg = Tm->ip[ipp].av[1];

  switch (model_type){
  case lewis_and_schrefler2:{//Lewis and Schrefler's book
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel
    break;
  }
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(pg);
}


/**
   function computes suction stress s = -(pg - pw) = -pc;
   @param ipp - integration point number

   @retval suction - suction stress [Pa]

   03/10/2023, TKr
*/
double con_awf2mat::give_suction(long ipp)
{
  double pw,pg,suction;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];

  if(rel_gas_press == 0)
    pg = pg - p_atm;
  
  switch (model_type){
  case lewis_and_schrefler2:{//Lewis and Schrefler's book
    suction = -1.0*(pg - pw);//this is correct, because capillary water pressure is negative
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel
    
    suction = -1.0*(pg - pw);//this is correct, because capillary water pressure is negative
    suction = suction/mefel_units;//corrected units for mefel //basic units = Pa    
    break;
  }
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 
  return(suction);
}


/**
   function returns the degree of saturation
   @param ipp - integration point number

   @retval saturation degree [-]

   27/05/2016, TKr
*/
double con_awf2mat::give_saturation_degree(long ipp)
{
  double pw,pg,s;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];

  s = get_sw(pw,pg,ipp);
  
  return(s);
}
/**
   function computes water content w [kg/kg]
   @param pw - water capillary pressure
   @param pg - gas pore pressure
   

   @retval w - computes water content w [kg/kg]
*/
double con_awf2mat::get_w(double pw,double pg,long ipp)
{
  double w,n,s;

  n = get_porosity(ipp);
  s = get_sw(pw,pg,ipp);

  w = (n*s*rhow0)/(1.0 - n)/rhos0;

  return(w);
}




/**
  The funtion marks required non-transport quantities in the array antq.

  @param antq - array with flags for used material types
                antq[i] = 1 => quantity type nontransquant(i+1) is required
                antq[i] = 0 => quantity type nontransquant(i+1) is not required

  @return The function does not return anything, but it may change content of antq array.
  
  03/10/2023, TKr
*/
void con_awf2mat::give_reqntq(long *antq)
{
  switch (model_type){
  case lewis_and_schrefler2:{//Lewis and Schrefler's book
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
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
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
  }
  } 
}



/****************************************************/
/** The basic functions and constitutive relations **/
/****************************************************/


/**
   function computes  volume density of concrete skeleton,
   changes of solid density, caused by dehydratation process

   @param t - temperature

   @retval rhos - volume density of soil skeleton
*/
double con_awf2mat::get_rhos(double /*t*/)
{
  double rhos;
   
  rhos = rhos0;
  
  return(rhos);
}

/**
   function computes water vapour diffusivity in air
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval cdiff - water vapour diffusivity in air

   03/10/2023, TKr
*/
double con_awf2mat::get_cdiff(double /*pw*/,double pg,double t)
{
  double cdiff;
  
  if(rel_gas_press == 1)
    pg = pg + p_atm;

  //conversion pow(a,b) -> exp(b*log(a))
  //cdiff = dv0*p0/pg*pow((t/t0),bv0);
  if (t < tcr)
    cdiff = dv0*p0/pg*exp(bv0*log(t/t0));
  else
    cdiff = dv0*p0/pg*exp(bv0*log(tcr/t0));

  //this equation below is from Patek's thesis:
  //dv0 = 5.9e-6;
  //cdiff = dv0*pow(t,2.3)/pg;

  return(cdiff);
}

/**
   function computes effective diffusion coefficient of vapour inside pores
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature
   @param ipp - number of integration point

   @retval deff - effective diffusion coefficient of vapour inside pores
   
   03/10/2023, TKr
*/
double con_awf2mat::get_dg(double pw,double pg,double t,long ipp)
{
  double deff;
  double phi,cdiff,sw,tau;

  deff = deff0;

  switch (deff_type){
  case 0:{//constant
    deff = cdiff0;
    break;
  }
  case 1:{//moisture dependent
    phi = get_porosity(ipp);
    cdiff = get_cdiff(pw,pg,t);
    sw = get_sw(pw,pg,ipp);
    tau = tau0;
    deff = phi*(1.0 - sw)*tau*cdiff;
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(deff);
}


/**
   function computes Biot's constant
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param ipp - number of integration point

   @retval alpha - Biot's constant
*/
double con_awf2mat::get_alpha(double pw,double pg,long ipp)
{
  double alpha;
  //double kt,ks;
  
  //kt = get_kt(pw,pg,ipp);
  //ks = get_ks(pw,pg,ipp);
  
  //alpha = 1.0 - kt/ks;  

  alpha = alpha0;
  
  return(alpha);
}

/**
   function computes bulk modulus of solid phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param ipp - number of integration point

   @retval ks - bulk modulus of solid phase
*/
double con_awf2mat::get_ks(double /*pw*/,double /*pg*/,long /*ipp*/)
{
  double ks;

  ks = ks0;

  return(ks);
}


/**
   function computes bulk modulus of porous medium
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param ipp - number of integration point

   @retval kt - bulk modulus of porous medium
*/
double con_awf2mat::get_kt(double /*pw*/,double /*pg*/,long /*ipp*/)
{
  double kt;
  
  kt = kt0;
  
  return(kt);
}


/**
   function computes mass concentration of water vapour air in gas phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhogw - mass concentration of water vapour air in gas phase
*/
double con_awf2mat::get_rhogw(double pw,double pg,double t)
{
  double rhogw,pgw;

  pgw = get_pgw(pw,pg,t0);

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
double con_awf2mat::get_pgw(double pw,double pg,double t)
{
  double pgw,rhow,pgws,tt,pc;

  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);

  if(rel_gas_press == 1)
    pg = pg + p0;
  
  pc = pg - pw;
  
  pgw = pgws*exp(-1.0*pc*mw/rhow/gasr/tt);

  return(pgw);
}



/**
   function computes partial derivative of pgw with respect to pc (Kelvin equation)
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval dpgw_dpc - partial derivative of pgw with respect to pc (Kelvin equation)
*/
double con_awf2mat::get_dpgw_dpc(double pw,double pg,double t)
{
  double dpgw_dpc,pgws,rhow,tt,pc;
  
  //critical point of water check
  if(t >= tcr)
    tt = tcr;
  else
    tt = t;

  pgws = get_pgws(tt);
  rhow = get_rhow(tt);
  
  if(rel_gas_press == 1)
    pg = pg + p0;

  pc = pg - pw;
    
  dpgw_dpc = -pgws*exp(-pc*mw/rhow/gasr/tt)*mw/rhow/gasr/tt;
  
  return(dpgw_dpc);
}

/**
   function computes water vapour saturation pressure
   @param t - temperature

   @retval pgws - water vapour saturation pressure
*/
double con_awf2mat::get_pgws(double t)
{
  /* //Clausius-Clapeyron equation
     double pgws,pgws0,t0,dhvap;
     
     dhvap = get_dhvap(pc,pg,t0);
     
     pgws = pgws0*exp(-mw*dhvap/gasr*(1.0/t - 1.0/t0));
  */

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
double con_awf2mat::get_rhow(double /*t*/)
{
  double rhow;

  rhow = rhow0;
  
  return(rhow);
}


/**
   function computes compresibility coefficient of water
   
   @retval kw - compresibility coefficient of water

*/
double con_awf2mat::get_kw(double /*pw*/,double /*pg*/,long /*ipp*/)
{
  double kw;
  
  kw = kw0;

  return(kw);
}

/**
   function computes dynamic viscosity of water = 1000e-6 Pa*s at 20 C
   @param t - temperature

   @retval muw dynamic viscosity of water = 1000e-6 Pa*s at 20 deg. C
*/
double con_awf2mat::get_muw(double /*t*/)
{
  double muw;
  
  muw = muw0;

  return(muw);
}

/**
   function computes molar mass of moist air
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval mg - molar mass of moist air
*/
double con_awf2mat::get_mg(double pw,double pg,double /*t*/)
{
  double mg,pgw;

  pgw = get_pgw(pw,pg,t0);

  if(rel_gas_press == 1)
    pg = pg + p_atm;

  mg = ma + (mw-ma)*pgw/pg;
  
  return(mg);
}


/**
   function computes gas phase density
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhog - gas phase density
*/
double con_awf2mat::get_rhog(double pw,double pg,double t)
{
  double rhog,pgw;

  pgw = get_pgw(pw,pg,t0);

  if(rel_gas_press == 1)
    pg = pg + p0;

  rhog = (pg*ma + (mw - ma)*pgw)/gasr/t;
  
  return(rhog);
}



/**
   function computes dynamic viscosity of moist air
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval mug - dynamic viscosity of moist air
*/
double con_awf2mat::get_mug(double pw,double pg,double t)
{
  double mug,mugw,muga,pga,pgw;

  pgw = get_pgw(pw,pg,t0);

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  
  pga = pg - pgw;
  muga = get_muga(t);
  mugw = get_mugw(t);

  //gas pressure check
  if(pgw <= pg)
    mug = mugw + (muga - mugw)*pow((1.0 - pgw/pg),0.6083);
  else
    mug = mugw;

  return(mug);
}


/**
   function computes dynamic viscosity of dry air
   @param t - temperature

   @retval muga - dynamic viscosity of dry air
*/
double con_awf2mat::get_muga(double /*t*/)
{
  double muga;

  muga = muga0;

  return(muga);
}

/**
   function computes dynamic viscosity of water vapour
   @param t - temperature

   @retval mugw - dynamic viscosity of water vapour
*/
double con_awf2mat::get_mugw(double /*t*/)
{
  double mugw;
 
  mugw = mugw0;

  return(mugw);
}



/**
   function computes mass concentration of dry air in gas phase
   @param pw - pore water pressure
   @param pg - pore gas pressure
   @param t - temperature

   @retval rhoga - mass concentration of dry air in gas phase
*/
double con_awf2mat::get_rhoga(double pw,double pg,double t)
{
  double rhoga,pgw;

  //gas pressure check
  pgw = get_pgw(pw,pg,t0);

  if(rel_gas_press == 1)
    pg = pg + p0;

  if (pgw <= pg)
    rhoga = (pg - pgw)*ma/gasr/t;
  else
    rhoga = 0.0;
    //rhoga = (pg - pgw)*ma/gasr/t;

  return(rhoga);
}


/***********************************************************/
/** The end of basic functions and constitutive relations **/
/***********************************************************/


