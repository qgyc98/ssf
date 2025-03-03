#include "globalt.h"
#include "saltmat4.h"
#include "globmatt.h"

/**********************************************
*                                             *
* General material model for 3 media transfer *
*      9 independent parameters               *
********************************************* */

/*
   class describes material model which deals with simultaneous
   transport of heat, moisture, salt and salt crystals
   
   
   ordering of unknowns:
   w - the volumetric moisture content (m^3/m^3)
   C_f - the concentration of free salts in water (kg/m^3 of solution)
   C_c - the amount of crystallized salt (kg/m^3 of sample)
   T - temperature (K)
   
   components in CORD:
   2 - density
   3 - porosity
   4 - faktor difusniho odporu
   5 - kappa
   6 - sorption izoterm
   7 - saturated moisture
   8 - none
   9 - Cecko
   10 - Lambda - the thermal conductivity (W/m/K)
   11 - not used
   12 - not used
   13 - not used
   14 - Dcoef - the salt diffusion coefficient (m$^2$/s)
   15 - binding isotherm
   16 - cfmax - the saturated free salt concentration (kg/m$^3$ of solution)
   17 - ws
   18 - not used
   19 - not used
   
   ordering of values in eqother array in integration points
   eqother[0] - the water vapour diffusion permeability (s)
   eqother[1] - relative humidity
   eqother[2] - derivative of the phi with respect to w
   eqother[3] - saturated volumetric moisture content
   eqother[4] - maximum
   
*/
FILE *in1a,*in2a,*in3a,*in4a;

saltmat4::saltmat4 ()
{
  //  molar mass of water kg.mol-1
  mw = 0.01801528;
  //  molar mass of dry air kg.mol-1
  ma = 28.9645;
  //  universal gas constant J.mol-1.K-1
  gasr = 8.31441;
  
  //  influence of damage on permeability
  daminfl=off;
  //  parameter for the generalized Heaviside function
  eps=1.0;
}

saltmat4::~saltmat4 ()
{
}


/**
   function reads material characteristics
   
   @param in - input file

   JM, 12.10.2011
*/
void saltmat4::read (XFILE *in)
{
  //  density
  rho.read (in);
  //  porosity
  por.read (in);
  //  water vapour diffusion resistance factor
  mu.read (in);
  //  moisture diffusivity
  kappa.read (in);
  //  sorption isotherm
  sorpiso.read (in);
  //  saturated moisture
  sm.read (in);
  //  specific heat capacity
  c.read (in);
  //  thermal conductivity
  lambda.read (in);
  //  Dcoef
  dcoef.read (in);
  //  binding isotherm
  bindiso.read (in);
  //  cfmax
  cfmax.read (in);
  //  ws
  ws.read (in);
  
  //  damage influence
  xfscanf (in,"%m",&flagsw_kwdset,&daminfl);
 
}


/**
   function prints material characteristics
   
   @param out - output file

   JM, 12.10.2011
*/
void saltmat4::print (FILE *out)
{
  //  density
  rho.print (out);
  //  porosity
  por.print (out);
  //  water vapour diffusion resistance factor
  mu.print (out);
  //  moisture diffusivity
  kappa.print (out);
  //  sorption isotherm
  sorpiso.print (out);
  //  saturated moisture
  sm.print (out);
  //  specific heat capacity
  c.print (out);
  //  thermal conductivity
  lambda.print (out);
  //  Dcoef
  dcoef.print (out);
  //  binding isotherm
  bindiso.print (out);
  //  cfmax
  cfmax.print (out);
  //  ws
  ws.print (out);

  //  damage influence
  fprintf (out," %d\n",daminfl);
}





/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index
   @param ipp - integration point id
   
*/
void saltmat4::matcond (matrix &d,long ri,long ci,long ipp)
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
    print_err("unknown number of space dimension is required",__FILE__,__LINE__,__func__);
  }
  }
  
  if (daminfl == on){
    //  conductivity matrix is modified due to damage
    damper.matcond (d,ipp);
  }
  
}

/**
   function computes conductivity %matrix of the material
   in the required integration point
   
   @param d - conductivity %matrix of material
   @param ri  - row index
   @param ci  - column index
   @param ipp - integration point id
   
*/
void saltmat4::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long n;
  n = d.n;
  
  switch (n){
  case 1:{
    //matcond1d (d,ri,ci,ipp);//1D
    break;
  }
  case 2:{
    matcond2d2 (d,ri,ci,ipp);//2D
    break;
  }
  case 3:{
    //matcond3d (d,ri,ci,ipp);//3D
    break;
  }
  default:{
    print_err("unknown number of space dimension is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function assembles conductivity %matrix of the material for 1D problems
   
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index
   @param ipp - integration point id
*/
void saltmat4::matcond1d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  
  if((ri == 0) && (ci == 0))
    k = k11 (ipp);
  if((ri == 0) && (ci == 1))
    k = k12 ();
  if((ri == 0) && (ci == 2))
    k = k13 ();
  if((ri == 0) && (ci == 3))
    k = k14 (ipp);
  
  if((ri == 1) && (ci == 0))
    k = k21 (ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (ipp);
  if((ri == 1) && (ci == 2))
    k = k23 ();
  if((ri == 1) && (ci == 3))
    k = k24 ();
  
  if((ri == 2) && (ci == 0))
    k = k31 ();
  if((ri == 2) && (ci == 1))
    k = k32 ();
  if((ri == 2) && (ci == 2))
    k = k33 ();
  if((ri == 2) && (ci == 3))
    k = k34 ();
  
  if((ri == 3) && (ci == 0))
    k = k41 (ipp);
  if((ri == 3) && (ci == 1))
    k = k42 ();
  if((ri == 3) && (ci == 2))
    k = k43 ();
  if((ri == 3) && (ci == 3))
    k = k44 (ipp);
  
  d[0][0] = k;
}

/**
   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index
   @param ipp - integration point id
*/
void saltmat4::matcond2d2 (matrix &d,long /*ri*/,long /*ci*/,long /*ipp*/)
{
  d[0][0] = 0.0;  d[0][1] = 0.0;
}

/**
   function creates conductivity %matrix of the material for 2D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index
   @param ipp - integration point id
*/
void saltmat4::matcond2d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  
  if((ri == 0) && (ci == 0))
    k = k11 (ipp);
  if((ri == 0) && (ci == 1))
    k = k12 ();
  if((ri == 0) && (ci == 2))
    k = k13 ();
  if((ri == 0) && (ci == 3))
    k = k14 (ipp);

  if((ri == 1) && (ci == 0))
    k = k21 (ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (ipp);
  if((ri == 1) && (ci == 2))
    k = k23 ();
  if((ri == 1) && (ci == 3))
    k = k24 ();

  if((ri == 2) && (ci == 0))
    k = k31 ();
  if((ri == 2) && (ci == 1))
    k = k32 ();
  if((ri == 2) && (ci == 2))
    k = k33 ();
  if((ri == 2) && (ci == 3))
    k = k34 ();

  if((ri == 3) && (ci == 0))
    k = k41 (ipp);
  if((ri == 3) && (ci == 1))
    k = k42 ();
  if((ri == 3) && (ci == 2))
    k = k43 ();
  if((ri == 3) && (ci == 3))
    k = k44 (ipp);
  
  d[0][0] = k;   d[0][1] = 0.0;
  d[1][0] = 0.0; d[1][1] = k;
  
}

/**
   function assembles conductivity %matrix of the material for 3D problems

   @param d   - conductivity %matrix of the material
   @param ri  - row index
   @param ci  - column index
   @param ipp - integration point id
*/
void saltmat4::matcond3d (matrix &d,long ri,long ci,long ipp)
{
  double k;
  
  if((ri == 0) && (ci == 0))
    k = k11 (ipp);
  if((ri == 0) && (ci == 1))
    k = k12 ();
  if((ri == 0) && (ci == 2))
    k = k13 ();
  if((ri == 0) && (ci == 3))
    k = k14 (ipp);
  
  if((ri == 1) && (ci == 0))
    k = k21 (ipp);
  if((ri == 1) && (ci == 1))
    k = k22 (ipp);
  if((ri == 1) && (ci == 2))
    k = k23 ();
  if((ri == 1) && (ci == 3))
    k = k24 ();

  if((ri == 2) && (ci == 0))
    k = k31 ();
  if((ri == 2) && (ci == 1))
    k = k32 ();
  if((ri == 2) && (ci == 2))
    k = k33 ();
  if((ri == 2) && (ci == 3))
    k = k34 ();
  
  if((ri == 3) && (ci == 0))
    k = k41 (ipp);
  if((ri == 3) && (ci == 1))
    k = k42 ();
  if((ri == 3) && (ci == 2))
    k = k43 ();
  if((ri == 3) && (ci == 3))
    k = k44 (ipp);
  
  
  d[0][0]=k;   d[0][1]=0.0; d[0][2]=0.0;
  d[1][0]=0.0; d[1][1]=k;   d[1][2]=0.0;
  d[2][0]=0.0; d[2][1]=0.0; d[2][2]=k;
}

/**
   function assembles capacity %matrix of the material
   
   @param c   - capacity %matrix of the material
   @param ri  - row index
   @param ci  - column index   
   @param ipp - integration point id
*/
void saltmat4::matcap (double &c,long ri,long ci,long ipp)
{
  if((ri == 0) && (ci == 0))
    c = c11 (ipp);
  if((ri == 0) && (ci == 1))
    c = c12 ();
  if((ri == 0) && (ci == 2))
    c = c13 ();
  if((ri == 0) && (ci == 3))
    c = c14 ();
  
  if((ri == 1) && (ci == 0))
    c = c21 (ipp);
  if((ri == 1) && (ci == 1))
    c = c22 (ipp);
  if((ri == 1) && (ci == 2))
    c = c23 ();
  if((ri == 1) && (ci == 3))
    c = c24 ();
  
  if((ri == 2) && (ci == 0))
    c = c31 (ipp);
  if((ri == 2) && (ci == 1))
    c = c32 (ipp);
  if((ri == 2) && (ci == 2))
    c = c33 ();
  if((ri == 2) && (ci == 3))
    c = c34 ();
  
  if((ri == 3) && (ci == 0))
    c = c41 ();
  if((ri == 3) && (ci == 1))
    c = c42 ();
  if((ri == 3) && (ci == 2))
    c = c43 ();
  if((ri == 3) && (ci == 3))
    c = c44 (ipp);
}



/**
   @param w - volumetric moisture content
   @param cf - concentration of free salts in water
   @param cc - amount of crystallized salt in sample
   @param t - temperature
   
   JM, JK, 3. 10. 2013
*/



/******************************************************
Conductivity and capacity terms for C and K matrices 
*******************************************************/

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k11 (long ipp)
{
  double k11,w,cf,t;
  double kapak, ps, delta, dfdw, rhow, wsat, cfmaxi,hf;
  
  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  t  = Tm->ip[ipp].av[3];

  //  saturated volumetric water content
  //wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  //  the moisture diffusivity (m$^2$/s)
  kapak = kappa.getval (w);

  //  partial pressure of saturated water vapour in the air
  ps = pgws(t);
  
  //  the water vapour diffusion permeability (s)
  delta = Tm->ip[ipp].eqother[0];
  //delta = permeabilitavodnipary(w,t);
  
  //  derivative of relative humidity with respect to moisture content
  //dfdw = sorpiso.derivative_inverse_isotherm_value (w);
  dfdw = Tm->ip[ipp].eqother[2];
  
  //  water density
  rhow = water_density ();
  
  k11 = rhow*kapak + ps*delta*dfdw;
  
  //  the saturated free salt concentration (kg/m$^3$ of solution)
  if (cf<0.0)
    cf=0.0;

  //cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];
  hf = genheaviside (cfmaxi-cf,eps);
  
  return k11*hf;
}


/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k12 ()
{
  double k12;
  
  k12 = 0.0;
  
  return k12;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k13 ()
{
  double k13;

  k13 = 0.0;

  return k13;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k14 (long ipp)
{
  double k14,delta,dpvsdt,phi,wsat,cfmaxi,hf,w,cf,t;
  
  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  t  = Tm->ip[ipp].av[3];
  
  //  saturated volumetric water content
  //wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  
  //   the water vapour diffusion permeability
  //delta = permeabilitavodnipary (w,t);
  delta = Tm->ip[ipp].eqother[0];
  
  //  relative_humidity
  //phi = sorpiso.inverse_isotherm_value (w);
  phi = Tm->ip[ipp].eqother[1];
  
  dpvsdt =  derivative_saturation_water_vapour_pressure_temperature (t);
  
  k14 = delta*phi*dpvsdt;
  
  //  the saturated free salt concentration (kg/m$^3$ of solution)
  if (cf<0.0)
    cf=0.0;
  //  cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];
  hf = genheaviside (cfmaxi-cf,eps);

  return k14*hf;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k21 (long ipp)
{
  double k21, kapak,wsat,w,cf;
  
  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  
  //  saturated volumetric water content
  //  wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;


  //  the moisture diffusivity (m$^2$/s)
  kapak = kappa.getval (w);

  k21 = kapak*cf;
  
  return k21;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k22 (long ipp)
{
  double k22;
  double dcoeff,wsat,w,cf;
  
  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];

  //  saturated volumetric water content
  //wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];
  
  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  dcoeff = dcoef.getval (cf);
  
  k22 = dcoeff*w;
  
  return k22;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k23 ()
{
  double k23;

  k23 = 0.0;
  
  return k23;
}

/**

   JM, JK, 3. 10. 2013
*/
double saltmat4::k24 ()
{
  double k24;

  k24 = 0.0;

  return k24;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k31 ()
{
  double k31;

  k31 = 0.0;

  return k31;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k32 ()
{
  double k32;

  k32 = 0.0;

  return k32;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k33 ()
{
  double k33;

  k33 = 0.0;

  return k33;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k34 ()
{
  double k34;

  k34 = 0.0;

  return k34;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k41 (long ipp)
{
  double k41;
  double wsat,lv1,ps,delta,dfdw,w,t;
  
  w  = Tm->ip[ipp].av[0];
  t  = Tm->ip[ipp].av[3];

  //  saturated volumetric water content
  //  wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  //  latent heat of evaporation of water (J/kg)
  lv1 = latent_heat_of_evaporation_of_water (t);

  //  partial pressure of saturated water vapour in the air
  ps = pgws (t);

  //  the water vapour diffusion permeability (s)
  //delta = permeabilitavodnipary(w,t);
  delta = Tm->ip[ipp].eqother[0];
  
  //  derivative of relative humidity with respect to the volumetric moisture content
  //dfdw = sorpiso.derivative_inverse_isotherm_value (w);
  dfdw =  Tm->ip[ipp].eqother[2];
  
  k41= lv1*ps*delta*dfdw;
  
  return k41;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k42 ()
{
  double k42;

  k42 = 0.0;

  return k42;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k43 ()
{
  double k43;

  k43 = 0.0;

  return k43;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::k44 (long ipp)
{
  double k44;
  double lv1,lamb,dpvsdt,delta,wsat,w,t;
  
  w  = Tm->ip[ipp].av[0];
  t  = Tm->ip[ipp].av[3];

  //  saturated volumetric water content
  // wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if(w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;


  //  latent heat of evaporation of water (J/kg)
  lv1 = latent_heat_of_evaporation_of_water (t);

  //  the water vapour diffusion permeability (s)
  //delta = permeabilitavodnipary(w,t);
  delta = Tm->ip[ipp].eqother[0];

  //  the thermal conductivity (W/m/K)
  lamb = lambda.getval (w);

  //  derivative of saturation water vapour pressure with respect to temperature
  dpvsdt =  derivative_saturation_water_vapour_pressure_temperature (t);
  
  k44= lamb + lv1 * delta * dpvsdt;
  
  return k44;
}


/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c11 (long ipp)
{
  double c11;
  double ps,phi,dfdw,wsoli,poro,rhow,wsat,hf,cfmaxi,w,cf,t;

  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];
  t  = Tm->ip[ipp].av[3];

  //  saturated volumetric water content
  //  wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  
  //  partial pressure of saturated water vapour in the air
  ps = pgws (t);

  //  porosity
  poro = por.getval (w);

  //  relative_humidity
  //phi = sorpiso.inverse_isotherm_value (w);
  phi = Tm->ip[ipp].eqother[1];

  //  derivative of relative humidity with respect to moisture content
  //dfdw = sorpiso.derivative_inverse_isotherm_value (w);
  dfdw = Tm->ip[ipp].eqother[2];
  
  //  water density
  rhow = water_density ();
  
  wsoli = ws.getval (w);
  
  c11 = rhow + mw/(gasr*294.15)*ps*(poro-w-wsoli)*dfdw-mw*ps*phi/(gasr*294.15);

  //  the saturated free salt concentration (kg/m$^3$ of solution)
  if (cf<0.0)
    cf=0.0;
  //  cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];
  hf = genheaviside (cfmaxi-cf,eps);
  
  return c11*hf;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c12 ()
{
  double c12;

  c12 = 0.0;

  return c12;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c13 ()
{
  double c13;

  c13 = 0.0;

  return c13;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c14 ()
{
  double c14;

  c14 = 0.0;

  return c14;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c21 (long ipp)
{
  double c21,cfmaxi,hf,cf;
  
  cf = Tm->ip[ipp].av[1];

  //  the saturated free salt concentration (kg/m$^3$ of solution)
  if (cf<0.0)
    cf=0.0;
  //  cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];

  hf = genheaviside (cfmaxi-cf,eps);
  
  c21 = cf*hf;
  
  return c21;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c22 (long ipp)
{
  double c22;
  double dcbdcf,cfmaxi,wsat,hf,w,cf;

  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];

  //  the saturated free salt concentration (kg/m$^3$ of solution)
  if (cf<0.0)
    cf=0.0;
  //  cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];
  
  //  saturated volumetric water content
  //wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  
  //  derivative of the binding isotherm with respect to the concentration of free salts in water
  dcbdcf=bindiso.derivative_isotherm_value (cf);
  hf = genheaviside (cfmaxi-cf,eps);
  
  c22 = w*hf + dcbdcf;
  
  return c22;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c23 ()
{
  double c23;

  c23 = 1.0;

  return c23;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c24 ()
{
  double c24;

  c24 = 0.0;

  return c24;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c31 (long ipp)
{
  double c31;
  double hf,cfmaxi,cf;

  cf = Tm->ip[ipp].av[1];

  //  the saturated free salt concentration (kg/m$^3$ of solution)
  //  cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];
  if (cf<0.0)
    cf=0.0;
  
  
  hf = genheaviside (cf-cfmaxi,eps);
  
  c31 = (cfmaxi-cf)*hf;
  
  return c31;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c32 (long ipp)
{
  double c32;
  double hf,cfmaxi,wsat,w,cf;
  
  w  = Tm->ip[ipp].av[0];
  cf = Tm->ip[ipp].av[1];

  //  the saturated free salt concentration (kg/m$^3$ of solution)
  //  cfmaxi = cfmax.getval (cf);
  cfmaxi = Tm->ip[ipp].eqother[4];
  if (cf<0.0)
    cf=0.0;

  //  saturated volumetric water content
  //  wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];

  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  
  hf = genheaviside (cf-cfmaxi,eps);
  
  c32 = -1.0*w*hf;
  
  return c32;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c33 ()
{
  double c33;

  c33 = 1.0;

  return c33;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c34 ()
{
  double c34;

  c34 = 0.0;

  return c34;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c41 ()
{
  double c41;

  c41 = 0.0;

  return c41;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c42 ()
{
  double c42;

  c42 = 0.0;

  return c42;
}

/**
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c43 ()
{
  double c43;

  c43 = 0.0;

  return c43;
}

/**
   @param ipp - integration point id
   
   JM, JK, 3. 10. 2013
*/
double saltmat4::c44 (long ipp)
{
  double c44;
  double rom,cap,wsat,w;
  
  w  = Tm->ip[ipp].av[0];

  //  saturated volumetric water content
  //  wsat = sm.getval (w);
  wsat = Tm->ip[ipp].eqother[3];
  
  if (w > wsat){
    //  the volumetric moisture content is greater than the saturated volumetric moisture content
    w = wsat;
  }
  if (w<0.0)
    w=0.0;
  
  //  density of the material
  rom = rho.getval (0.0);
  //  cpecific heat capacity
  cap = c.getval (w);
  c44 = rom * cap;
  
  return c44;
}









/*********************
  Boundary conditions
*********************/
	
/**
   function determines transmission coefficient
   in some cases, the boundary conditions are prescribed in different
   variables than variables used in the problem
   for example, moisture content is used in the problem but boundary
   condition is prescribed with the help of pressures

   @param trc - prescribed transmission coefficient on the boundary
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double saltmat4::transmission_transcoeff (double trc,long ri,long ci,long nid,long /*bc*/)
{
  double new_trc,w,cf,cc,t;
  new_trc = 0.0;
  
  //  moisture content
  w  = nodalval(nid, 0);
  //  concentration of free salts in water
  cf = nodalval(nid, 1);
  //  amount of crystallized salt
  cc = nodalval(nid, 2);
  //  temperature
  t  = nodalval(nid, 3);
  
  
  if((ri == 0) && (ci == 0))
    //new_trc = get_transmission_transcoeff_11 (w,cf,cc,t,bc,ipp);
    new_trc = trc;
  if((ri == 0) && (ci == 1))
    new_trc = 0.0;
  if((ri == 0) && (ci == 2))
    new_trc = 0.0;
  if((ri == 0) && (ci == 3))
    new_trc = 0.0;
  

  if((ri == 1) && (ci == 0))
    new_trc = 0.0;
  if((ri == 1) && (ci == 1))
    new_trc = trc;
  if((ri == 1) && (ci == 2))
    new_trc = 0.0;
  if((ri == 1) && (ci == 2))
    new_trc = 0.0;
  if((ri == 1) && (ci == 3))
    new_trc = 0.0;
  

  if((ri == 2) && (ci == 0))
    new_trc = 0.0;
  if((ri == 2) && (ci == 1))
    new_trc = 0.0;
  if((ri == 2) && (ci == 2))
    new_trc = 0.0;
  if((ri == 2) && (ci == 3))
    new_trc = 0.0;
  
  if((ri == 3) && (ci == 0))
    new_trc = 0.0;
  if((ri == 3) && (ci == 1))
    new_trc = 0.0;
  if((ri == 3) && (ci == 2))
    new_trc = 0.0;
  if((ri == 3) && (ci == 3))
    new_trc = trc;

  return new_trc;
}

/**
   function creates correct transfer coefficient on the boundary (transmission) for 1st medium
   @param f11        - correct transfer coefficient
   @param w ... cc  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
/*
double saltmat4::get_transmission_transcoeff_11 (double w,double cf,double cc,double t,long bc,long ipp)
{
  double trc;

  switch (bc){//type of prescribed variable
  case 30:{//transmission
    trc=1.0;//should be changed
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }

  return(trc);
}
*/


/**
   function determines nodal value
   in some cases, the boundary conditions are prescribed in different
   variables than variables used in the problem
   for example, moisture content is used in the problem but boundary
   condition is prescribed with the help of pressures

   @param nodval - prescribed nodal value
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double saltmat4::transmission_nodval (double nodval,long ri,long ci,long nid,long /*bc*/)
{
  double new_nodval,w,cf,cc,t;
  new_nodval = 0.0;
  
  //  moisture content
  w  = nodalval(nid, 0);
  //  concentration of free salts in water
  cf = nodalval(nid, 1);
  //  amount of crystallized salt
  cc = nodalval(nid, 2);
  //  temperature
  t  = nodalval(nid, 3);
 
  if((ri == 0) && (ci == 0))
    //new_nodval = get_transmission_nodval_11(nodval,w,cf,cc,t,bc,ipp);
    new_nodval = nodval;
  if((ri == 0) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 0) && (ci == 2))
    new_nodval = 0.0;
if((ri == 0) && (ci == 3))
    new_nodval = 0.0;
  
  if((ri == 1) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 1))
    new_nodval = nodval;
  if((ri == 1) && (ci == 2))
    new_nodval = 0.0;
  if((ri == 1) && (ci == 3))
    new_nodval = 0.0;
 
  if((ri == 2) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 2) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 2) && (ci == 2))
    new_nodval = 0.0;
  if((ri == 2) && (ci == 3))
    new_nodval = 0.0;

  if((ri == 3) && (ci == 0))
    new_nodval = 0.0;
  if((ri == 3) && (ci == 1))
    new_nodval = 0.0;
  if((ri == 3) && (ci == 2))
    new_nodval = 0.0;
  if((ri == 3) && (ci == 3))
    new_nodval = nodval;


  return new_nodval;
}


/**
   function creates correct new nodal value on the boundary (transmission) for 1st medium
   @param new_nodval - new prescribed value near the boundary
   @param bv         - value of prescribed value near the boundary
   @param w ... cc  - actual unknowns on the boundary
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
/*
double saltmat4::get_transmission_nodval_11 (double bv,double w,double cf,double cc,double t,long bc,long ipp)
{
  double nodval;

  switch (bc){//type of prescribed variable
  case 30:{//transmission
    nodval = bv;//should be changed
    break;
  } 
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return(nodval);
}
*/

/**
   function computes flux through boundary
   
   @param nodval - prescribed nodal value
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double saltmat4::transmission_flux (double nodval,long ri,long ci,long nid,long bc)
{
  double flux,w,cf,cc,t;
  flux = 0.0;
  
  //  moisture content
  w  = nodalval(nid, 0);
  //  concentration of free salts in water
  cf = nodalval(nid, 1);
  //  amount of crystallized salt
  cc = nodalval(nid, 2);
  //  temperature
  t  = nodalval(nid, 3);
  
  
  if((ri == 0) && (ci == 0))
    flux = get_transmission_flux_ww (nodval,w,bc);
  if((ri == 0) && (ci == 1))
    flux = 0.0;
  if((ri == 0) && (ci == 2))
    flux = 0.0;
  if((ri == 0) && (ci == 3))
    flux = 0.0;
  
  if((ri == 1) && (ci == 0))
    flux = 0.0;
  if((ri == 1) && (ci == 1))
    flux = 0.0;
  if((ri == 1) && (ci == 2))
    flux = 0.0;
  if((ri == 1) && (ci == 3))
    flux = 0.0;
  
  
  if((ri == 2) && (ci == 0))
    flux = 0.0;
  if((ri == 2) && (ci == 1))
    flux = 0.0;
  if((ri == 2) && (ci == 2))
    flux = 0.0;
  if((ri == 2) && (ci == 3))
    flux = 0.0;
  
  
  if((ri == 3) && (ci == 0))
    flux = 0.0;
  if((ri == 3) && (ci == 1))
    flux = 0.0;
  if((ri == 3) && (ci == 2))
    flux = 0.0;
  if((ri == 3) && (ci == 3))
    flux = 0.0;
  
  return (flux);
}


/**
   function computes flux through the boundary (transmission - convective mass transfer) for the first medium
   
   @param bv - prescribed value on the boundary
   @param w - actual moisture content on the boundary
   @param bc - type of boundary condition
*/
double saltmat4::get_transmission_flux_ww (double bv,double w,long bc)
{
  double flux;
  
  switch (bc){//type of prescribed variable
  case 30:{//transmission - boundary flux
    flux = (bv - w);//should be changed
    break;
  }
  default:{
    print_err("no acceptable boundary condition is prescribed",__FILE__,__LINE__,__func__);
    exit(0);
  }
  }
  
  return(flux);
}




/**
   function computes all variables in nodes
   @param compother  - number of other components
   @param ipp        - first integration point on element
   @param w ... cc  - actual unknowns on the boundary
   
*/
double saltmat4::get_othervalue(long compother,long /*ipp*/, double w,double cf,double cc,double /*t*/)
{
  double other,cb=0.0;
  
  switch (compother){
  case 1:{
    //CorD(15,1000,0,cf,cb,b,a3);
    other= w*cf+cb+fabs(cc);
  }
    break;
  case 2://other  = Tm->ip[ipp].eqother[9];
    other = 0.0;
    break;
  case 3://other  = Tm->ip[ipp].eqother[10];
    other = 0.0;
    break;
  case 4://other  = Tm->ip[ipp].eqother[11];
    other = 0.0;
    break;
  case 5:
    //get_rel_hum2 (w, other, ppp); // Tm->ip[ipp].eqother[5];
    break;
  case 6:other  = 0.0; //Tm->ip[ipp].eqother[5];
    break;
  case 7: //CorD(15,kod,1000,0,cf,other,b,a3);
    other = 0.0;
    break;
  case 0:{//first unknown
    
	  // nefunguje.. ja to mam v uzlech vypocitane, ale tady to jde pres IPP.....
    /*	  if(ipp >11) ipp = 10;
      if( Tt->nodes[ipp].eqother[9] == w)
      {
      other = Tt->nodes[ipp].eqother[8];
      }
      else
      {
      other = -10000.0;
      }*/
    //CorD(15,1000,0,cf,cb,b,a3);
    other= w*cf+cb+cc;///pozor..nevim, zda to tam ma byt...
      //other = 0.0;
      break; 
  }
  default:{
    print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
  return (other);
}



/**
     function prints names of all variables in nodes
     @param out       - output file
     @param compother - number of other components
*/
void saltmat4::print_othervalue_name(FILE *out,long compother)
{
  
  switch (compother){
  case 0:{//Relative humidity
    fprintf (out,"Moisture content w (m3/m3)       ");
    break;
  }
  case 1:{//Temperature
    fprintf (out,"Temperature (K)             ");
    break;
  }
  case 2:{//Moisture content w
    fprintf (out,"Moisture content w (m3/m3)   ");
    break;
  }
  case 3:{//water vapour pressure pgw
    fprintf (out,"Water vapour pressure (Pa)  ");
    break;
  }
  case 4:{//Moisture content u 
    fprintf (out,"Moisture content u (kg/kg)   ");
    break;
  }
  case 5:{//Moisture content u 
    fprintf (out,"nevim (kg/kg)   ");
    break;
  }
  case 6:{//Moisture content u 
    fprintf (out,"nevim (kg/kg)   ");
    break;
  }
  case 7:{//Moisture content u 
    fprintf (out,"nevim (kg/kg)   ");
    break;
  }
  case 8:{//Moisture content u 
    fprintf (out,"nevim u (kg/kg)   ");
    break;
  }
  default:{
    print_err("unknown type of component is required",__FILE__,__LINE__,__func__);
  }
  }
}


/****************************************************
*****************************************************
    auxiliary functions
*****************************************************
****************************************************/

/**
   water density
   
   3. 10. 2013
*/
double saltmat4::water_density ()
{
  double rhow;
  
  rhow = 1000.0;
  
  return rhow;
}

/**
   pressure of saturated water vapour
   
   @param t - temperature (K)
*/
double saltmat4::pgws (double t)
{
  double pgw;
  
  pgw = exp(23.5771 - 4042.9/(t - 37.58));
  
  return pgw;
}

/**
   water vapour diffusion permeability (s)
   
   @param w - the volumetric moisture content (m^3/m^3)
   @param t - temperature (K)
   
*/
double saltmat4::permeabilitavodnipary (double w, double t)
{
  double Da, Dp,rv,pa, mi, p;
  
  mi = mu.getval (w);

  rv = 461.5;
  pa = 101325.0;
  p = 101325.0;
  
  Da = (2.306e-5 * pa)/(rv * t * p)*pow((t/294.15),1.81);
  Dp = Da/mi;
  
  return Dp;
}


/**
   derivative of saturation water vapour pressure with respect to temperature
   
   @param t - temperature (K)
   
   2. 10. 2013
*/
double saltmat4::derivative_saturation_water_vapour_pressure_temperature (double t)
{
  return ((4042.9 * exp(23.5771-4042.9/(t - 37.58)))/((t-37.58)*(t-37.58)));
}

/**
   latent heat of evaporation of water (J/kg)
   
   param t - temperature (K)
   
   2. 10. 2013
*/
double saltmat4::latent_heat_of_evaporation_of_water (double t)
{
  return (2.5008e6)*pow((273.15/t),(0.167+t*3.67e-4));
}

















void saltmat4::der_value_hyst (int /*matchar*/,int kod, double /*pv*/, double & /*outvalue*/,double & /*outvalue2*/, long /*ipp*/)
{
  
  switch (kod){
	case 0:
		break;
	case 1:{
		break;
		   }
	case 2:{
		break;
		   }
	case 30:{
		break;
			}
	case 31:{
		break;
			}
	case 40:{

	  //outvalue = derivation_dy_dx(matchar,pv,0,1);
	  //outvalue2 = derivation_dy_dx(matchar,pv,2,3);

		}
		break;
	}

}






long saltmat4::cycle_detection (double *r,double *pr, double *ppr)
{
	double t,pvt, ppvt;

	t = r[3];
	pvt = pr[3];
	ppvt = ppr[3];

	if (ppvt > 273.15 && pvt > t && t < 273.15 && ppvt > pvt && pvt < 273.15)
	{
		//fprintf(Outt, "\nTime je %lf  t je %lf a pvt je %lf a ppvt je %lf",Tp->time,t, pvt, ppvt);
		//printf ("\nTime je %lf  t je %lf a pvt je %lf a ppvt je %lf",Tp->time,t, pvt, ppvt);
		return(1);
	}
	else {
		return(0);
	}

}





void saltmat4::read_Sourcet(XFILE *in)
{
double x, y;

	
	xfscanf(in, "%ld",&pocet_radku);
    source = new double *[pocet_radku+1];
	  for( long j = 0; j< pocet_radku+1; j++)
	  {
		  source [j] = new double [2];
	  }

//source[0][1] = pocet_radku;
	 
 for(int i=1; i< pocet_radku+1; i++)
	{
		xfscanf(in, "%lf %lf",&x, &y);
		source[i][0] = x;
		source[i][1] = y;
	}
source[0][0] = x;
source[0][1] = y;

}

double saltmat4::getval_source (double t)
{
double y;

while(t>=31536000)
	{
	t = t-31536000;
	}	

long i=1;

if (source[i][0] > t)
	{
		return(0.0);
	}


if(source[0][0] > t)
	{
	while (source[i][0] <= t){
		i= i+1;
		}
	   
		if (source[i-1][0]+3600 > t)
			{
		   //	y = source[i-1][1] +((0.0 - source[i-1][1])/((source[i-1][0]+3600) - source[i-1][0])*(t - source[i-1][0]));
			y = source[i-1][1];
			}
		else{
		  //	y = source[i-1][1] +((source[i][1] - source[i-1][1])/(source[i][0] - source[i-1][0])*(t - source[i-1][0]));
			return(0.0);
			}
}
else{
		i = long(source[0][1]);
		y = source[i][1] +((0.0 - source[i][1])/((source[i][0]+3600) - source[i][0])*(t - source[i][0]));
}
		


if (y<0.0) return(0.0);

double fdt;
fdt = Tp->timecont.forwarddt;
if(fdt < 60)
{
 fprintf(Outt, "%le  %le \n",  Tp->time, fdt);// kapa a decko

}

return(fdt*y/3600);

}

/**
   function compare obtained values with limits
   
   @param nv - array of values
   
*/
void saltmat4::values_correction (vector &nv,long ipp)
{
  double cfmax,cf,wsat,w;
  
  //  volumetric moisture content
  w = nv[0];
  //  saturated volumetric moisture content
  wsat = Tm->ip[ipp].eqother[3];
  
  if (w<0.0)
    w=0.0;
  if (w>wsat)
    w=wsat;
  
  
  //  concentration of free salts in water
  cf = nv[1];
  //  maximum concentration
  cfmax=Tm->ip[ipp].eqother[4];
  
  if (cf<0.0)
    cf=0.0;
  if (cf>cfmax)
    cf=cfmax;
  
  
  //fprintf (Outt,"\n ipp %3ld   w  %le  wsat  %le   cf %le   cfmax %le",ipp,w,wsat,cf,cfmax);

  nv[0]=w;
  nv[1]=cf;

  // je to jen pro urcity material... nutno opravit... 25.ledna 2008
  //if (nv[1] > cfmax){
  //	fprintf(Outt, "%le  %le  %ld\n",  Tp->time, nv[1],ipp);// kapa a decko
  //nv[1] = 190;
  //}
  //
  //if (nv[1] < 0.0){
  //	fprintf(Outt, "%le  %le  %ld\n",  Tp->time, nv[1],ipp);// kapa a decko
  //nv[1] = 0.0;
  //}
  //if (nv[1] < 1e-8){
  //	fprintf(Outt, "%le  %le  %ld\n",  Tp->time, nv[1],ipp);// kapa a decko
  //nv[1] = 0.0;
  //}
  //if(nv[0] < 0){
  //  nv[0] = 1e-6;
  //}
  //if(nv[2] < 0){
  //nv[2] = 1.0e-12;
  //}
  //
  //if(nv[0] >= wsat){
  //nv[0] = wsat;
  //}
  //
  //if(nv[3] <= 240){
  //nv[3] = 240;
  //}
  //if(nv[3] >= 350){
  //nv[3] = 350;
  // }
  //if(nv[2] < 0){
  //fprintf(Outt, " Cfkryst %le  %le  %ld\n",  Tp->time, nv[2],ipp);// kapa a decko
  //nv[2] = 0.0;
  //}

}





/**
   function selects auxiliary values
   
   @param ipp - integration point id
   @param av - array of actual values
   @param pv - array of previous values
   @param eq - array of values stored in eqother array
   
   JK, 7.1.2008
*/
void saltmat4::give_values (long ipp,double *av, double *pv,double *eq)
{
  //  actual values
  av[0] = Tm->ip[ipp].av[0];
  av[1] = Tm->ip[ipp].av[1];
  av[2] = Tm->ip[ipp].av[2];
  av[3] = Tm->ip[ipp].av[3];
  
  //  values from the previous time step
  pv[0] = Tm->ip[ipp].pv[0];
  pv[1] = Tm->ip[ipp].pv[1];
  pv[2] = Tm->ip[ipp].pv[2];
  pv[3] = Tm->ip[ipp].pv[3];
  
  // ********************************
  //  values stored in eqother array
  // ********************************
  //  water vapour diffusion permeability
  eq[0] = Tm->ip[ipp].eqother[0];
  //  relative humidity
  eq[1] = Tm->ip[ipp].eqother[1];
  //  derivative of the relative humidity with repsect to the moisture content
  eq[2] = Tm->ip[ipp].eqother[2];
  //  saturated volumetric moisture content
  eq[3] = Tm->ip[ipp].eqother[3];
  //  maximum concentration
  eq[4] = Tm->ip[ipp].eqother[4];
  //  total salt concentration
  eq[5] = Tm->ip[ipp].eqother[5];
}

/**
   function computes auxiliary values which are necessary for future computation
   
   input and output values have to be sent via function argument because this
   function is called in integration points as well as in nodes, therefore,
   the input and output values cannot be obtained directly in the function
   
   @param ipp - integration point id
   @param in - array with actual values (read the text above)
   @param inp - array with values from the previous time step
   @param ine - array with components of eqother array
   @param out - array with computed values (read the text above)
   
   JM, 29.5.2007, revision 16. 10. 2013
*/
void saltmat4::aux_values (long /*ipp*/,double *in,double */*inp*/,double */*ine*/,double *out)
{
  double time,w,cf,cb,cc,t,ctot,delta,phi,dphidw,smc,cfmaxi;
  
  //  actual time
  time = Tp->time;
  
  //  the volumetric moisture content
  w  = in[0];
  //  the concentration of free salts in water
  cf = in[1];
  //  the amount of crystallized salt
  cc = in[2];
  //  the temperature
  t  = in[3];
  
  //  the water vapour diffusion permeability (s)
  delta = permeabilitavodnipary(w,t);
  
  //  relative humidity
  phi = sorpiso.inverse_isotherm_value (w);
  
  //  derivative of the phi with respect to w
  dphidw = sorpiso.derivative_inverse_isotherm_value (w);
  
  //  saturated volumetric moisture content
  smc = sm.getval (w);
  
  //  maximum 
  cfmaxi = cfmax.getval (cf);
  
  //  concentration of bonded salt in the whole porous body
  cb = bindiso.isotherm_value (cf);
  
  //  total salt concentration
  ctot = cf*w+cc+cb;

  out[0]=delta;
  out[1]=phi;
  out[2]=dphidw;
  out[3]=smc;
  out[4]=cfmaxi;
  out[5]=ctot;
  out[6]=cf;
  out[7]=cc;
  out[8]=cb;
}



/**
   function saves auxiliary values
   
   @param ipp - integration point id
   @param out - array with auxiliary values
   
   JK, 7.1.2008
*/
void saltmat4::save_values (long ipp,double *out)
{
  //  water vapour diffusion permeability
  Tm->ip[ipp].eqother[0]=out[0];
  //  relative humidity
  Tm->ip[ipp].eqother[1]=out[1];
  //  derivative of the relative humidity with repsect to the moisture content
  Tm->ip[ipp].eqother[2]=out[2];
  //  saturated volumetric moisture content
  Tm->ip[ipp].eqother[3]=out[3];
  //  maximum concentration
  Tm->ip[ipp].eqother[4]=out[4];
  //  total salt concentration
  Tm->ip[ipp].eqother[5]=out[5];
  //  salt concentration
  Tm->ip[ipp].eqother[6]=out[6];
  //  crystallized salt amount
  Tm->ip[ipp].eqother[7]=out[7];
  //  salt bonded in the body
  Tm->ip[ipp].eqother[8]=out[8];
}




/**
   The funtion marks required non-transport quantities in the array antq.
   
   @param antq - array with flags for used material types
                 antq[i] = 1 => quantity type nontransquant(i+1) is required
                 antq[i] = 0 => quantity type nontransquant(i+1) is not required

   @return The function does not return anything, but it may change content of antq array.
   
   24. 10. 2013
*/
void saltmat4::give_reqntq(long *antq)
{
  if (daminfl == on){
    //  damage parameter
    antq[scal_iso_damage-1] = 1;
    //  process zone length
    antq[proc_zone_length-1] = 1;
    //  crack width
    antq[crack_width-1] = 1;
  }
}



/**
   Function returns temperature in the given integration point.
   
   @param ipp - integration point id

   @return Funtion returns value of temperature stored in the integartion point.

   24. 10. 2013, JK
*/
double saltmat4::give_temperature (long ipp)
{
  return Tm->ip[ipp].av[3];
}



/**  
  Function returns initial temperature in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of initial temperature stored in the integartion point.

  Created by TKo, 21.11.2013
*/
double saltmat4::give_inittemperature (long /*ipp*/)
{
  print_err("Not yet implemented", __FILE__, __LINE__, __func__);
  abort();
}



/**  
  Function returns relative humidity in the given integration point.
 
  @param ipp - integration point id

  @return Funtion returns value of relative humidity stored in the integartion point.

  Created by TKo, 21.11.2013
*/
double saltmat4::give_rel_hum (long ipp)
{
  return Tm->ip[ipp].eqother[1];
}



/**
  Function returns pore pressure in the given integration point.
   
  @param ipp - integration point id
   
  @return Funtion returns value of pore pressure stored in the integartion point.

  Created by JK+TKo,  1. 11. 2013
*/
double saltmat4::give_pore_pressure (long /*ipp*/)
{
  return -1.0e5;
}
