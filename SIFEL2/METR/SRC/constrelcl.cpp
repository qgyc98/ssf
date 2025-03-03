/*
    File:             constrelcl.cpp
    Author:           Tomas Krejci
    Purpose:          constitutive relations for thermo-hydro-mechanical problem


    FEM FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS
    --------------------------------------------------

    sources: 
    1. FULLY COUPLED THERMO-HYDRO-MECHANICAL ANALYSIS BY AN ALGEBRAIC MULTIGRID METHOD
    Wang Xicheng, B.A. Schrefler
    
    2. NUMERICAL ANALYSIS OF HYGRO-THERMAL BEHAVIOUR AND DAMAGE OF CONCRETE AT HIGH TEMPERATURE
    D. Gawin, C.E. Majorana, B.A. Schrefler

    3. NONLINEAR MODELLING OF CONCRETE AS MULTIPHASE POROUS MATERIAL IN HIGH TEMPERATURE CONDITIONS
    Francesco Pesavento - doctoral thesis
    
    LIST OF USED VARIABLES:
    -----------------------
    
    cp    ... effective specific heat of porous medium (J.kg-1.K-1)    
    cpg   ... effective specific heat of gaz mixture (J.kg-1.K-1)    
    cpw   ... effective specific heat of liquid (water) (J.kg-1.K-1)    
    cps   ... effective specific heat of solid matrix (J.kg-1.K-1)    
    ct    ... tangent matrix (Pa)
    deff  ... effective diffusivity of gas mixture (m2.s-1)
    g     ... acceleration of gravity (m.s-2)
    kh    ... hydraulic conductivity
    kap   ... absolute permeability (m2)
    krg   ... relative permeability of gas phase
    krw   ... relative permeability of liquid phase
    ks    ... bulk modulus of solid phase
    kt    ... bulk modulus of porous medium
    m     ... molar mass of gas mixture (moist air) (kg.kmol-1)
    ma    ... molar mass of dry air (kg.kmol-1)
    mw    ... molar mass of water vapour (kg.kmol-1)
    gasr  ... universal gas constant 8.3144e3 J.mol-1.K-1
    pav   ... average pressure of mixture (Pa)
    patm  ... atmospheric pressure (Pa)
    pc    ... capillary pressure (Pa)
    pg    ... presure of gas phase (Pa)
    pgw   ... water vapour partial presure (Pa)
    pgws  ... water vapour saturation presure (Pa)
    pw    ... liquid water presure (Pa)
    s     ... liquid phase volumic saturation (=Sw = liquid volume/pore volume)
    t     ... temperature (K)

    alpha ... Biot's constant
    alphac... convective heat transfer coefficient on the boundary
    betac ... mass transfer coefficient on the boundary
    betas ... cubic thermal expansion coefficient of solid (K-1)
    betaswg ... volume thermal dilatation of solid - water vapor for incompressible solid grains
    betasw ... volume thermal dilatation of solid - water for incompressible solid grains
    betaw ... volume thermal dilatation of water
    dhvap ... enthalpy of vaporization per unit mass (J.kg-1)
    eps   ... strain tensor
    epsv  ... volumetric strain
    epst  ... thermoelastic strain
    lambdaeff ... effective thermal conductivity (W.m-1.K-1)
    rho   ... effective density of porous medium (kg.m-3)
    rhog  ... gas phase density (kg.m-3)
    rhoga ... mass concentration of dry air in gas phase (kg.m-3)
    rhogw ... mass concentration of water vapour in gas phase (kg.m-3)
    rhow  ... liquid phase density (kg.m-3)
    phi   ... porosity (pore volume/total volume)
    mug   ... dynamic viscosity of gas phase (kg.m-1.s-1)
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliasc.h"
#include "globalc.h"
#include "constrelcl.h"
#include "constrel.h"

state_eqcl::state_eqcl()
{}
state_eqcl::~state_eqcl()
{}

//DEGREE OF SATURATION
/**
   function computes degree of saturation
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature

   @retval s - degree of saturation
*/
double state_eqcl::get_s(double pc,double /*pg*/,double t,long ipp)
{
  long i;
  double s;
  
  i = Cml->ip[ipp].idm;

  switch (Cml->ip[ipp].tm){//type of material
  case concreteBc:{
    s = Cml->concretec[i].concreteB_sw(pc,t);//function from concreteht.cpp
    break;
  }
  case baroghelBc:{
    s = Cml->baroghelc[i].baroghel_sw(pc,t);//function from baroghel.cpp
    break;
  }
  case C60baroghelBc:{
    s = Cml->C60baroghelc[i].sat(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    s = Cml->C30baroghelc[i].sat(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    s = Cml->o30bazantc[i].sat(pc,t);//function from o30bazant.cpp
    break;
  }
  case C60bazantBc:{
    s = Cml->C60bazantc[i].sat(pc,t);//function from o30bazant.cpp
    break;
  }
    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(s);
}

// BIOT'S CONSTANT
/**
   function computes Biot's constant
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature

   @retval alpha - Biot's constant
*/
double state_eqcl::get_alpha(double /*pc*/, double /*pg*/, double /*t*/,long ipp)
{
  long i;
  double alpha;
  
  i = Cml->ip[ipp].idm;
  
  switch (Cml->ip[ipp].tm){//type of material
    /* case concreteBc:{
       alpha = Cml->concretec[i].concreteB_alpha();//function from concreteB.cpp
       break;
       }
    */
  case baroghelBc:{
    alpha = Cml->baroghelc[i].baroghel_alpha();//function from concreteB.cpp
    break;
  }    
  case C60baroghelBc:{
    alpha = Cml->C60baroghelc[i].C60bar_alpha();//value from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    alpha = Cml->C30baroghelc[i].C30bar_alpha();//value from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    alpha = Cml->o30bazantc[i].o30baz_alpha();//value from o30bazant.cpp
    break;
  }
  case C60bazantBc:{
    alpha = Cml->C60bazantc[i].C60baz_alpha();//value from o30bazant.cpp
    break;
  }
  default:{
    alpha = 1.0;//provisionally
    //fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  
  return(alpha);
}


//AVERAGED DENSITY

/**
   function computes density of solid phase
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature

   @retval rho - averaged density of solid phase
*/
double state_eqcl::get_rhos(double /*pc*/,double /*pg*/,double t,long ipp)
{
  long i;
  double rhos;

  i = Cml->ip[ipp].idm;

  switch (Cml->ip[ipp].tm){//type of material
  case concreteBc:{
    rhos = Cml->concretec[i].concreteB_rhos(t);//function from concreteB.cpp
    break;
  }
  case baroghelBc:{
    rhos = Cml->baroghelc[i].baroghel_rhos();//function from baroghel.cpp
    break;
  }
  case C60baroghelBc:{
    rhos = Cml->C60baroghelc[i].C60bar_rhos();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    rhos = Cml->C30baroghelc[i].C30bar_rhos();//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    rhos = Cml->o30bazantc[i].o30baz_rhos(t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantBc:{
    rhos = Cml->C60bazantc[i].C60baz_rhos();//function from o30bazant.cpp
    break;
  }    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
  return(rhos);
}


/**
   function computes cubic thermal expansion coefficient of solid (K-1)
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature

   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
*/
double state_eqcl::get_betas(double /*pc*/,double /*pg*/,double /*t*/,long ipp)
{
  long i;
  double betas;

  i = Cml->ip[ipp].idm;

  switch (Cml->ip[ipp].tm){//type of material
  case concreteBc:{
    betas = Cml->concretec[i].concreteB_betas();//function from concreteB.cpp
    break;
  }
  case baroghelBc:{
    betas = Cml->baroghelc[i].baroghel_betas();//function from baroghel.cpp
    break;
  }
  case C60baroghelBc:{
    betas = Cml->C60baroghelc[i].C60bar_betas();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    betas = Cml->C30baroghelc[i].C30bar_betas();//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    betas = Cml->o30bazantc[i].o30baz_betas();//function from o30bazant.cpp
    break;
  }                
  case C60bazantBc:{
    betas = Cml->C60bazantc[i].C60baz_betas();//function from o30bazant.cpp
    break;
  }                
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(betas);
}


/**
   function computes water relative permeability
   @param pg - capillary gas pressure
   @param t - temperature

   @retval krw - water relative permeability
*/
double state_eqcl::get_krw(double pc,double t,long ipp)
{
  long i;
  double krw,s,rh;
  state_eq tt;

  i = Cml->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//type of material
  case concreteBc:{
    s = Cml->concretec[i].concreteB_sw(pc,t);//function from concreteB.cpp
    rh = tt.get_rh(pc,t);
    krw = Cml->concretec[i].concreteB_krw(s,rh);//function from concreteB.cpp
    break;
  }
  case baroghelBc:{
    krw = Cml->baroghelc[i].baroghel_krw(pc,t);//function from baroghelB.cpp
    break;
  }
  case C60baroghelBc:{
    rh = tt.get_rh(pc,t);
    krw = Cml->C60baroghelc[i].C60bar_krw(pc,t,rh);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    rh = tt.get_rh(pc,t);
    krw = Cml->C30baroghelc[i].C30bar_krw(pc,t,rh);//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    rh = tt.get_rh(pc,t);
    krw = Cml->o30bazantc[i].o30baz_krw(pc,t,rh);//function from o30bazant.cpp
    break;
  }        
  case C60bazantBc:{
    krw = Cml->C60bazantc[i].C60baz_krw(pc,t);//function from C60bazant.cpp
    break;
  }        
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(krw);
}

/**
   function computes intrinsic permeability
   @param pc - capillary pressure
   @param pg - capillary gas pressure
   @param t - temperature

   @retval kintr - intrinsic permeability
*/
double state_eqcl::get_kintr(double pc,double pg,double t,long ipp)
{
  long i;
  double kintr,dam=0.0;

  i = Cml->ip[ipp].idm;
  //dam = Tm->damage[ipp];//??!!

  switch (Tm->ip[ipp].tm){//type of material
  case concreteBc:{
    //kintr = Cml->concrete[i].concreteB_kintr(pg,t,dam);//function from concreteB.cpp
    break;
  }
  case baroghelBc:{
    kintr = Cml->baroghelc[i].baroghel_kintr();//function from baroghelB.cpp
    break;
  }
  case C60baroghelBc:{
	//  dam = Mm->givemechq(scal_iso_damage, ipp);
    kintr = Cml->C60baroghelc[i].C60bar_kintr(pc,pg,t,dam);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
	  //  dam = Mm->givemechq(scal_iso_damage, ipp);
	  kintr = Cml->C30baroghelc[i].C30bar_kintr(pc,pg,t,dam);//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
	  //  dam = Mm->givemechq(scal_iso_damage, ipp);
	  kintr = Cml->o30bazantc[i].o30baz_kintr(pc,pg,t,dam);//function from o30bazant.cpp
    break;
  }            
  case C60bazantBc:{
    kintr = Cml->C60bazantc[i].C60baz_kintr();//function from C60bazant.cpp
    break;
  }            
    
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(kintr);
}
