/*
    File:             consol_wf2.cpp
    Author:           Tomas Krejci, 22/05/2019
    Purpose:          material model for saturated-nonsaturated water flow in a deforming porous medium
    sources:          Lewis and Schrefler pp. 93-97, material parameters are set for benchmark on page n. 168
    unknowns:         number of unknowns=2, pw = liqud water pressure, pg = gas(air) pressure
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "constrel.h"
#include "consol_wf2.h"
#include "globalt.h"
#include "globmatt.h"

con_wf2mat::con_wf2mat()
{
  compress = 0;          //compressible grains: 0=no; 1=yes
  por_type = 0;          //porosity calculation type
  kintr_type = 0;        //intrinsic permability calculation type
  krw_type = 0;          //relative permeability calculation type
  krg_type = 0;          //relative permability calculation type
  sr_type = 1;           //retention curve calculation type

  vol_strain_effect = 0; //volumetric strain rate influence: 0=no; 1=yes 
  mefel_units = 1.0;//basic units for pressures = Pa (Pascals)

  p_atm = 101325.0; //atmospheric pressure [Pa]
  p_atm_kpa = 101.325; //athmospheric pressure in kPa

  // PHYSICAL PROPERTIES OF WATER
  rhow0 = 1000.0;//water density [kg/m^3]
  muw0 = 1.0e-3;//water viscosity [Pa.s]
  kw = 2.0e9;//bulk modulus of water [Pa]
  
  mug0 = 1.8e-5;//air viscosity [Pa.s]

  //PHYSICAL PROPERTIES OF SOIL set to zero
  alpha = 1.0;//Biot's constant [-] alpha = 1 - kt/ks
  ks = 2.167e9;//bulk modulus of solid phase (grains) [Pa]
  kt0 = 0.0;
  phi0 = 0.297;//initial porosity [-]
  kintr0 = 4.5e-13;//intrinsic permeability [m^2]
  rhos0 = 2500.0; //solid density
  sirr = 0.0;
  ssat = 1.0;

  bb1 = 0.0;
  phi01 = 0.0;

  pw_bc = 0.0; //free boundary pressure 
  rel_gas_press = 0;

  t0 = 293.15; // reference temperature
}

con_wf2mat::~con_wf2mat()
{}


/**
   function reads parameters
   
   @param in - input file

   29/10/2009, TKr
*/
void con_wf2mat::read(XFILE *in)
{
  xfscanf (in,"%k%m","airwaterflowtype",&airwaterflowtype_kwdset, &model_type);
  xfscanf (in,"%d", &compress);

  // common material parameters
  xfscanf (in,"%le %le %le %le %le %d %d %d %d %d", &alpha, &ks, &rhos0, &pw_bc, &kintr0, &por_type, &kintr_type, &krw_type, &krg_type, &sr_type);
  
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
    break;
  }
  case 1:
  case 2:{//dependent on saturation degree:
    //xfscanf (in,"%le %le", &sirr, &ssat); //not finished
    break;
  }
  case 3:{//exponential
    //xfscanf (in,"%le", &lambda_krw); //not finished
    break;
  }
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  //retention curve type:
  switch (sr_type){
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
void con_wf2mat::print(FILE */*out*/)
{
  /*   fprintf (out,"\n %d ", int(model_type));
       fprintf (out,"\n %d ", compress);
       
       switch (model_type){
       case lewis_and_schrefler2:{//Lewis and Schrefler's book
       fprintf (out,"\n %le %le %le %le %le %le %le %le %le \n",alpha,ks,phi0,kw,rhow,muw0,mug0,kintr0,rhos0);
       lewis_ret.print(out);
       break;
       }
       case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book
       fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le %le %d \n",mefel_units,alpha,ks,phi0,kw,rhow,muw0,mug0,kintr0,rhos0, pw_bc, vol_strain_effect);
       break;
       }
       case van_genuchten2:{//partially saturated medium =Van Genuchten model
       fprintf (out,"\n %le %le %le %le %le %le %le %le %le %le %le %d \n", mefel_units,alpha,ks,phi0,kw,rhow,muw0,mug0,kintr0,rhos0, pw_bc, vol_strain_effect);
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
void con_wf2mat::matcond (matrix &d,long ri,long ci,long ipp)
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
void con_wf2mat::matcond1d (matrix &d,long ri,long ci,long ipp)
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
void con_wf2mat::matcond2d (matrix &d,long ri,long ci,long ipp)
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
void con_wf2mat::matcond3d (matrix &d,long ri,long ci,long ipp)
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
void con_wf2mat::matcap (double &c,long ri,long ci,long ipp)
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
void con_wf2mat::rhs_volume (matrix &d,long ri,long ci,long ipp)
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
void con_wf2mat::rhs_volume2 (double &cc,long ri,long /*ci*/,long ipp)
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
void con_wf2mat::rhs1d1 (matrix &d,long ri,long /*ci*/,long ipp)
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
void con_wf2mat::rhs2d1 (matrix &d,long ri,long /*ci*/,long ipp)
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
void con_wf2mat::rhs3d1 (matrix &d,long ri,long /*ci*/,long ipp)
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

   @retval sw - degree of saturation

   29/10/2009, TKr
*/
double con_wf2mat::get_sw(double pw, double /*pg*/, long ipp)
{
  double sw;
  sw = 0.0;

  switch (sr_type){
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

   @retval cs - specific moisture content = partial derivative of degree of saturation with respect to pc


   29/10/2009, TKr
*/
double con_wf2mat::get_cs(double pw, double /*pg*/, long ipp)
{
  double dsw_dpc;
  dsw_dpc = 0.0;

  switch (sr_type){
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
  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

  return(dsw_dpc);
}


/**
   function computes water relative permeability
   @param pw - water pressure
   @param pg - air pressure
   
   @retval krw - water relative permeability

   12/9/2008, TKr
*/
double con_wf2mat::get_krw(double pw, double pg, long ipp)
{
  double krw,sw;
  krw = 1.0;

  switch (krw_type){
  case 0:{//constant
    krw = 1.0;
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

  default:{
    print_err ("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  } 

  return(krw);
}



/**
   function computes air relative permeability
   @param pw - water pressure
   @param pg - air pressure
   
   @retval krg - gas relative permeability

   12/9/2008, TKr
*/
double con_wf2mat::get_krg(double pw, double pg, long ipp)
{
  double krg,sw;
  krg=1.0;

  
  switch (krg_type){
  case 0:{//constant
    krg = 1.0;
    break;
  }
  case 1:{//dependent on saturation degree
    sw = get_sw(pw,pg,ipp);
    krg = 1.0; //not finished
    break;
  }
  case 2:{//dependent on saturation degree
    sw = get_sw(pw,pg,ipp);
    krg = 1.0; //not finished
    break;
  }
  case 3:{//exponential
    sw = get_sw(pw,pg,ipp);
    krg = 1.0; //not finished
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
   @param t - temperature
   
   @retval kintr - intrinsic permeability

   13/11/2018, TKr
*/
double con_wf2mat::get_kintr(double /*pw*/, double /*pg*/, long ipp)
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
  case 2:{//not finished
    //Kozeny's approach for bentonite:
    //phi = get_porosity(ipp);
    //kintr = kintr0*(phi*phi*phi)/(1 - phi)/(1 - phi)*(1 - phi0)*(1 - phi0)/phi0/phi0/phi0;

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

   @retval kww - conductivity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_kww(double pw, double pg, long ipp)
{
  double krw,kww,kintr;

  kintr = get_kintr(pw,pg,ipp);
  krw = get_krw(pw,pg,ipp);
  kww = krw*kintr/muw0;//p. 97

  return(kww);
}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kwg - conductivity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_kwg(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  double kwg;
  
  kwg = 0.0;

  return(kwg);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kgw - conductivity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_kgw(double /*pw*/, double /*pg*/, long /*ipp*/)
{
  double kgw;
  
  kgw = 0.0;

  return(kgw);

}


/**
   function creates conductivity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval kgg - conductivity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_kgg(double pw, double pg, long ipp)
{
  double krg,kgg,kintr;

  kintr = get_kintr(pw,pg,ipp);
  krg = get_krg(pw,pg,ipp);
  kgg = krg*kintr/mug0;//p. 97

  return(kgg);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capww - capacity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_capww(double pw, double pg, long ipp)
{
  double sw,cs,n,capww;  
  
  sw = get_sw(pw,pg,ipp);
  cs = get_cs(pw,pg,ipp);
  n = get_porosity(ipp);

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  
  if(compress == 1){
    //compressible grains:
    // note: cs = ds_dpc;
    capww = (alpha-n)/ks*sw*(sw - pw*cs + pg*cs) + n*sw/kw - n*cs;//this equation is correct
  }
  else{
    //incompressible grains:
    capww = -n*cs;
  }  
  check_math_errel(0);

  return(capww);

}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capwg - capacity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_capwg(double pw, double pg, long ipp)
{
  double sw,sg,cs,n,capwg;  
  
  sw = get_sw(pw,pg,ipp);
  sg = 1.0 - sw;
  cs = get_cs(pw,pg,ipp);
  n = get_porosity(ipp);

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  
  if(compress == 1){
    //compressible grains:
    // note: cs = ds_dpc;
    capwg = (alpha-n)/ks*sw*(sg + pw*cs - pg*cs) + n*cs;//this equation is correct
  }
  else{
    //incompressible grains:
    capwg = n*cs;
  }

  return(capwg);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capgw - capacity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_capgw(double pw, double pg, long ipp)
{
  double sw,sg,cs,n,capgw;  
  
  sw = get_sw(pw,pg,ipp);
  sg = 1.0 - sw;
  cs = get_cs(pw,pg,ipp);
  n = get_porosity(ipp);

  if(rel_gas_press == 1)
    pg = pg + p_atm;

  if(compress == 1){
    //compressible grains:
    // note: cs = ds_dpc;
    capgw = (alpha-n)/ks*sg*(sw + pg*cs - pw*cs) + n*cs;//this equation is correct
  }
  else{
    //incompressible grains:
    capgw = n*cs;
  }

  return(capgw);
}


/**
   function creates capacity coefficient of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval capgg - capacity coefficient

   29/10/2009, TKr
*/
double con_wf2mat::get_capgg(double pw, double pg, long ipp)
{
  double sw,sg,cs,n,capgg;  
  
  sw = get_sw(pw,pg,ipp);
  sg = 1.0 - sw;
  cs = get_cs(pw,pg,ipp);
  n = get_porosity(ipp);

  if(rel_gas_press == 1)
    pg = pg + p_atm;
  
  if(compress == 1){
    //compressible grains:
    // note: cs = ds_dpc;
    capgg = (alpha-n)/ks*sg*(sg - pg*cs + pw*cs) + n*sg/pg - n*cs;//this equation is correct
  }
  else{
    //incompressible grains:
    capgg = n*sg/pg - n*cs;//this equation is correct
  }  
  return(capgg);

}


/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval fw - first part for right-hand side coeficient for continutiy equation

   29/10/2009, TKr
*/
double con_wf2mat::get_fw1(double pw, double pg, long ipp)
{
  double fw1,krw,kintr;

  krw = get_krw(pw,pg,ipp);
  kintr = get_kintr(pw,pg,ipp);
  fw1 = krw*kintr/muw0*rhow0;//p. 97

  return(fw1);

}


/**
   function returns coefficient for righ-hand side of the general material 
   @param pw - water pressure
   @param pg - air pressure

   @retval fg - first part for right-hand side for continutiy equation

   29/10/2009, TKr
*/
double con_wf2mat::get_fg1(double pw, double pg, long ipp)
{
  double fg1,krg,rhog,kintr;

  krg = get_krg(pw,pg,ipp);
  rhog = 1.25;//temporarilly
  kintr = get_kintr(pw,pg,ipp);
  fg1 = krg*kintr/mug0*rhog;//p. 97
  fg1 = 0.0;//no gravity force is included for the air
  return(fg1);

}



/**
   function returns coefficient for righ-hand side of the general material
   @param pw - water pressure
   @param pg - air pressure

   @retval fwu - volumetric strain rate effect on pore water pressure; only for partially coupled METR version

   16/05/2018, TKr
*/
double con_wf2mat::get_fwu(double pw, double pg, long ipp)
{
  double depsv_r,sw,fwu,dsr_depsv,n;
  
  fwu = 0.0;
  
  switch (model_type){
  case lewis_and_schrefler2:{
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book
    if(vol_strain_effect == 1){
      n = get_porosity(ipp);
      depsv_r = Tm->givenontransq(strain_vol_rate, ipp);      //actual rate of volumetric strain from MEFEL
      Tm->ip[ipp].eqother[0] = depsv_r;
      if (sr_type == mefel_sr){
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
	Tm->ip[ipp].eqother[1] = dsr_depsv; //this is not necessary
      }
      sw = get_sw(pw,pg,ipp);
      fwu = -alpha*sw*depsv_r;//volumetric strain effect
      
      if (sr_type == mefel_sr){
	fwu = fwu - n*dsr_depsv*depsv_r;//volumetric strain effect
      }
    }
    break;
  }
  case van_genuchten2:{
    fwu = 0.0;
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

   @retval fgu - volumetric strain rate effect on pore gas pressure; only for partially coupled METR version

   16/05/2018, TKr
*/
double con_wf2mat::get_fgu(double pw, double pg, long ipp)
{
  double depsv_r,sw,fgu,dsr_depsv,n;
  
  fgu = 0.0;
  
  switch (model_type){
  case lewis_and_schrefler2:{
    break;
  }
  case lewis_and_schrefler2_mefel:{
    if(vol_strain_effect == 1){
      n = get_porosity(ipp);
      depsv_r = Tm->givenontransq(strain_vol_rate, ipp);      //actual rate of volumetric strain from MEFEL
      
      if (sr_type == mefel_sr){
	dsr_depsv = Tm->givenontransq(der_saturation_deg_depsv, ipp); //actual derivative of saturation degree with respect to volumetric strain
      }
      
      sw = get_sw(pw,pg,ipp);
      fgu = -alpha*(1.0 - sw)*depsv_r;
      
      if (sr_type == mefel_sr){
	fgu = fgu + n*dsr_depsv*depsv_r;//volumetric strain rate effect
      }
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
double con_wf2mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
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
double con_wf2mat::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp, int flag)
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
double con_wf2mat::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
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
double con_wf2mat::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
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
double con_wf2mat::get_transmission_transcoeff_ww(double pw,double /*pg*/,long bc,long /*ipp*/)
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
double con_wf2mat::get_transmission_transcoeff_ww(double pw,double /*pg*/,long bc,long /*ipp*/,int flag)
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
double con_wf2mat::get_transmission_nodval_ww(double bv,double pw,double /*pg*/,long bc,long /*ipp*/)
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
double con_wf2mat::get_transmission_flux_ww(double bv,double pw,double /*pg*/,long bc,long /*ipp*/)
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
double con_wf2mat::get_transmission_transcoeff_gg(double /*pw*/,double /*pg*/,long bc,long /*ipp*/)
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
double con_wf2mat::get_transmission_nodval_gg(double bv,double /*trr*/,double /*pw*/,double /*pg*/,long bc,long /*ipp*/)
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
double con_wf2mat::get_transmission_flux_gg(double bv,double /*trr*/,double /*pw*/,double pg,long bc,long /*ipp*/)
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

double con_wf2mat::get_othervalue(long compother,double pw, double pg, long ipp)
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
    other = get_suction(pw,pg,ipp);
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
void con_wf2mat::print_othervalue_name(FILE *out,long compother)
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
void con_wf2mat::values_correction (vector &nv, long ipp)
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
void con_wf2mat::gaspress_check(double /*pw*/,double &/*pg*/,long /*ipp*/)
{
  switch (model_type){
  case lewis_and_schrefler2://Lewis and Schrefler's book
  case van_genuchten2:
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    break;
  }
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }
  
}


/**
   function checks if water pressure is non-positive
   @param pw - pore water pressure
   @param ipp - number of integration point


   03/03/2011, TKr
*/
void con_wf2mat::waterpress_check(double &/*pw*/,double /*pg*/,long /*ipp*/)
{
  switch (model_type){
  case lewis_and_schrefler2:{//Lewis and Schrefler's book
    break;
  }
  case van_genuchten2:
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    break;
  }
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
    abort();
  }
  }

}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.

   23/05/2016, TKr
*/
void con_wf2mat::updateval (long /*ipp*/)
{
}


/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number

   23/05/2016, TKr
*/
void con_wf2mat::initval(long /*ipp*/)
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
void con_wf2mat::give_dof_names(namevart *dofname, long ntm)
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
   function returns water pressure

   @param ipp - integration point number

   @retval pw - pore water pressure

   27/05/2016, TKr
*/
double con_wf2mat::give_water_pressure(long ipp)
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
double con_wf2mat::give_pore_pressure(long ipp)
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
double con_wf2mat::give_gas_pressure(long ipp)
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

   27/05/2016, TKr
*/
double con_wf2mat::give_suction(long ipp)
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
   function computes suction stress s = -(pg - pw) = -pc;
   @param ipp - integration point number

   @retval suction - suction stress [Pa]

   @param pw - pore water pressure
   @param pw - pore gas pressure

   01/05/2018, TKr
*/
double con_wf2mat::get_suction(double pw, double pg, long /*ipp*/)
{
  double suction;

  if(rel_gas_press == 0)
    pg = pg - p_atm;
  
  switch (model_type){
  case van_genuchten2:
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
double con_wf2mat::give_saturation_degree(long ipp)
{
  double pw,pg,s;
  
  pw = Tm->ip[ipp].av[0];
  pg = Tm->ip[ipp].av[1];

  s = get_sw(pw,pg,ipp);
  
  return(s);
}


/**
   function returns porosity

   @retval phi - porosity

   15/04/2021, TKr
*/
double con_wf2mat::get_porosity(long ipp)
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
   function computes water content w [kg/kg]
   @param pw - water capillary pressure
   @param pg - gas pore pressure
   

   @retval w - computes water content w [kg/kg]
*/
double con_wf2mat::get_w(double pw,double pg,long ipp)
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
  
  27/05/2016, TKr
*/
void con_wf2mat::give_reqntq(long *antq)
{
  switch (model_type){
  case lewis_and_schrefler2:{//Lewis and Schrefler's book
    break;
  }
  case lewis_and_schrefler2_mefel:{//Lewis and Schrefler's book, saturation degree ind its derivative are obtained from mefel;
    antq[saturation_deg-1] = 1;
    antq[der_saturation_deg-1] = 1;
    antq[der_saturation_deg_depsv-1] = 1;
    antq[porosity-1] = 1;
    antq[strain_vol_rate-1] = 1;
    break;
  }  
  default:{
    print_err("unknown model type is required", __FILE__, __LINE__, __func__);
  }
  } 
}
