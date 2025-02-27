/*
    File:             constrelcu.cpp
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
#include "constrelcu.h"
#include "glasgowmatc.h"
#include "consol_hwf2c.h"

state_eqcu::state_eqcu()
{
  scale_pc = 1.0;//Tp->scale[0];
  scale_pg = 1.0;//Tp->scale[1];
  scale_t = 1.0;//Tp->scale[2];
  scale_u = 1.0;//Mp->scale;
}
state_eqcu::~state_eqcu()
{}

//DEGREE OF SATURATION
/**
   function computes degree of saturation
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature

   @retval s - degree of saturation
*/
double state_eqcu::get_s(double pc,double /*pg*/,double t,long ipp)
{
  long i;
  double s;
  
  i = Cmu->ip[ipp].idm;

  switch (Cmu->ip[ipp].tm){//type of material
  case concreteBc:{
    s = Cmu->concretec[i].concreteB_sw(pc,t);//function from concreteht.cpp
    break;
  }
  case baroghelBc:{
    s = Cmu->baroghelc[i].baroghel_sw(pc,t);//function from baroghel.cpp
    break;
  }
  case C60baroghelBc:{
    s = Cmu->C60baroghelc[i].sat(pc,t);//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    s = Cmu->C30baroghelc[i].sat(pc,t);//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    s = Cmu->o30bazantc[i].sat(pc,t);//function from o30bazant.cpp
    break;
  }
  case C60bazantBc:{
    s = Cmu->C60bazantc[i].sat(pc,t);//function from o30bazant.cpp
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
double state_eqcu::get_alpha(double /*pc*/, double /*pg*/, double /*t*/,long ipp)
{
  long i;
  double alpha;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//type of material
    /* case concreteBc:{
       alpha = Cmu->concretec[i].concreteB_alpha();//function from concreteB.cpp
       break;
       }
    */
  case baroghelBc:{
    alpha = Cmu->baroghelc[i].baroghel_alpha();//function from concreteB.cpp
    break;
  }    
  case C60baroghelBc:{
    alpha = Cmu->C60baroghelc[i].C60bar_alpha();//value from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    alpha = Cmu->C30baroghelc[i].C30bar_alpha();//value from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    alpha = Cmu->o30bazantc[i].o30baz_alpha();//value from o30bazant.cpp
    break;
  }
  case C60bazantBc:{
    alpha = Cmu->C60bazantc[i].C60baz_alpha();//value from o30bazant.cpp
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
double state_eqcu::get_rhos(double /*pc*/,double /*pg*/,double t,long ipp)
{
  long i;
  double rhos;

  i = Cmu->ip[ipp].idm;

  switch (Cmu->ip[ipp].tm){//type of material
  case concreteBc:{
    rhos = Cmu->concretec[i].concreteB_rhos(t);//function from concreteB.cpp
    break;
  }
  case baroghelBc:{
    rhos = Cmu->baroghelc[i].baroghel_rhos();//function from baroghel.cpp
    break;
  }
  case C60baroghelBc:{
    rhos = Cmu->C60baroghelc[i].C60bar_rhos();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    rhos = Cmu->C30baroghelc[i].C30bar_rhos();//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    rhos = Cmu->o30bazantc[i].o30baz_rhos(t);//function from o30bazant.cpp
    break;
  }    
  case C60bazantBc:{
    rhos = Cmu->C60bazantc[i].C60baz_rhos();//function from o30bazant.cpp
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
double state_eqcu::get_betas(double /*pc*/,double /*pg*/,double /*t*/,long ipp)
{
  long i;
  double betas;

  i = Cmu->ip[ipp].idm;

  switch (Cmu->ip[ipp].tm){//type of material
  case concreteBc:{
    betas = Cmu->concretec[i].concreteB_betas();//function from concreteB.cpp
    break;
  }
  case baroghelBc:{
    betas = Cmu->baroghelc[i].baroghel_betas();//function from baroghel.cpp
    break;
  }
  case C60baroghelBc:{
    betas = Cmu->C60baroghelc[i].C60bar_betas();//function from C60baroghel.cpp
    break;
  }
  case C30baroghelBc:{
    betas = Cmu->C30baroghelc[i].C30bar_betas();//function from C30baroghel.cpp
    break;
  }
  case o30bazantBc:{
    betas = Cmu->o30bazantc[i].o30baz_betas();//function from o30bazant.cpp
    break;
  }                
  case C60bazantBc:{
    betas = Cmu->C60bazantc[i].C60baz_betas();//function from o30bazant.cpp
    break;
  }                
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(betas);
}


/**
   function creates elastic stiffness matrix for coupled problems
   
   @param d     - stiffness matrix of the material
   @param mssst - stress-strain state
   @param ipp   - index of integration point
   
   @return The function returns material elastic stiffness %matrix in the parameter d.
   
   11/09/2001
*/
void state_eqcu::matstiff (matrix &d,strastrestate mssst, long ipp)
{
  switch (mssst){
  case bar:{
    matstiff_bar (d,ipp);
    break;
  }
  case planestrain:{
    matstiff_plstress (d,ipp);
    break;
  }
  case planestress:{
    matstiff_plstress (d,ipp);
    break;
  }
  case spacestress:{
    matstiff_plstress (d,ipp);
    break;
  }
  case axisymm:{
    matstiff_axi(d,ipp);
    break;
  }
  default:{
    fprintf (stderr,"\n unknown number of components of stress tensor is required");
    fprintf (stderr,"\n in function elastisomat::matstiff (%s, line %d).\n",__FILE__,__LINE__);
  }
  }
}


/**
   function creates elastic stiffness matrix for 1D coupled problems
   
   @param d   - stiffness matrix of the material

   @param ipp - index of integration point
   
   @return The function returns elastic material stiffness %matrix in the parameter d.

   11/09/2001
*/
void state_eqcu::matstiff_bar (matrix &d,long ipp)
{
  d[0][0] = give_e(ipp);
}



/**
   function creates elastic stiffness matrix for 2D coupled problems (plane stress)
   
   @param d   - stiffness matrix of the material

   @param ipp - index of integration point
   
   @return The function returns elastic material stiffness %matrix in the parameter d.

   05/05/2010, TKr
*/
void state_eqcu::matstiff_plstress (matrix &d,long ipp)
{
  double c,e,nu;
  
  fillm(0.0,d);
  e = give_e(ipp);
  nu = give_nu(ipp);
  c = e/(1.0-nu*nu);
  
  d[0][0] = c;     d[0][1] = c*nu;  d[0][2] = 0.0;
  d[1][0] = c*nu;  d[1][1] = c;     d[1][2] = 0.0;
  d[2][0] = 0.0;   d[2][1] = 0.0;   d[2][2] = e/2.0/(1.0+nu);
}


/**
   Function creates elastic stiffness %matrix for 2D coupled problems (plane strain)
   
   @param d   - stiffness %matrix of the material

   @param ipp - index of integration point
   
   @return The function returns elastic material stiffness %matrix in the parameter d.
   
   05/05/2010, TKr
*/
void state_eqcu::matstiff_plstrain (matrix &d,long ipp)
{
  double c,e,nu;
  
  fillm(0.0,d);
  e = give_e(ipp);
  nu = give_nu(ipp);

  c = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0] = c*(1.0-nu);   d[0][1] = c*nu;         d[0][2] = 0.0;
  d[1][0] = c*nu;         d[1][1] = c*(1.0-nu);   d[1][2] = 0.0;
  d[2][0] = 0.0;          d[2][1] = 0.0;          d[2][2] = e/2.0/(1.0+nu);

  if (d.m > 3)
    {
      d[0][3] = d[0][1]; d[1][3] = d[1][0];
      d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
    }
}


/**
   Function creates elastic stiffness %matrix for 3D coupled problems
   
   @param d   - stiffness %matrix of the material

   @param ipp - index of integration point
   
   @return The function returns elastic material stiffness %matrix in the parameter d.
   
   05/05/2010, TKr
*/
void state_eqcu::matstiff_spacestress (matrix &d,long ipp)
{
  double e,nu;
  double g,s;
  
  fillm(0.0,d);
  e = give_e(ipp);
  nu = give_nu(ipp);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=s*nu;
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;         d[4][4]=g;        d[5][5]=g;
}




/**
  Function creates elastic stiffness %matrix of the elastic
  isotropic material for 2D axisymmetric problems.
   
  @param d   - stiffness %matrix of the material

  @param ipp - index of integration point
  
  @return The function returns elastic material stiffness %matrix in the parameter d.

  17/07/2018, TKr
*/
void state_eqcu::matstiff_axi (matrix &d,long ipp)
{
  double e,nu;
  double g,s;
  
  fillm(0.0,d);
  e = give_e(ipp);
  nu = give_nu(ipp);
  
  g = e/2.0/(1.0+nu);
  s = e/(1.0+nu)/(1.0-2.0*nu);
  
  d[0][0]=s*(1-nu);  d[0][1]=s*nu;     d[0][2]=d[0][1];
  d[1][0]=d[0][1];   d[1][1]=d[0][0];  d[1][2]=d[0][1];
  d[2][0]=d[0][1];   d[2][1]=d[0][1];  d[2][2]=d[0][0];

  d[3][3]=g;

}




/**
   function computes Young's modulus (Pa)

   @param ipp - number of integration point

   @retval emod -  Young's modulus (Pa)
*/
double state_eqcu::give_e(long ipp)
{
  long i;
  double emod;

  i = Cmu->ip[ipp].idm;

  switch (Cmu->ip[ipp].tm){//type of material
   case isotransmatc:{
    emod = Cmu->itrmc[i].get_e();//from isotrmat.cpp
    break;
  }
  case consolhwf2c:{
    emod = Cmu->consol_hwf2c[i].get_e(ipp);//provisionally
    break; 
  }
  case consolhawf3c:{
    emod = Cmu->consol_hawf3c[i].get_e(ipp);//provisionally
    break; 
  }
  case concreteBc:{
    emod = Cmu->concretec[i].concreteB_emod();//from concreteBc.cpp
    break;
  }
  case baroghelBc:{
    emod = Cmu->baroghelc[i].baroghel_emod();//function from baroghelBc.cpp
    break;
  }
  case C60baroghelBc:{
    //emod = Cmu->C60baroghelc[i].C60bar_emod(pc,pg,t);//function from C60baroghelc.cpp
    break;
  }
  case C30baroghelBc:{
    //emod = Cmu->C30baroghelc[i].C30bar_emod(pc,pg,t);//function from C30baroghelc.cpp
    break;
  }
  case o30bazantBc:{
    emod = Cmu->o30bazantc[i].o30baz_emod();//function from o30bazantc.cpp
    break;
  }                
  case C60bazantBc:{
    //emod = Cmu->C60bazantc[i].C60baz_emod(pc,pg,t);//function from o30bazantc.cpp
    break;
  }                
  case glasgowc:{
    emod = Cmu->tenchc[i].emod();//function from glasgowc.cpp//temporarily??!!
    break;
  }                
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(emod);
}

/**
   function computes cubic thermal expansion coefficient of solid (K-1)
   @param pc - capilary pressure
   @param pg - capilary gas pressure
   @param t - temperature

   @retval betas -  cubic thermal expansion coefficient of solid (K-1)
*/
double state_eqcu::give_nu(long ipp)
{
  long i;
  double nu;

  i = Cmu->ip[ipp].idm;

  switch (Cmu->ip[ipp].tm){//type of material
   case isotransmatc:{
    nu = Cmu->itrmc[i].get_nu();//from isotrmat.cpp
    break;
  }
  case consolhwf2c:{
    nu = Cmu->consol_hwf2c[i].get_nu(ipp);//provisionally
    break; 
  }
  case consolhawf3c:{
    nu = Cmu->consol_hawf3c[i].get_nu(ipp);//provisionally
    break; 
  }
  case concreteBc:{
    nu = Cmu->concretec[i].concreteB_nu();//from concreteBc.cpp
    break;
  }
  case baroghelBc:{
    nu = Cmu->baroghelc[i].baroghel_nu();//function from baroghelBc.cpp
    break;
  }
  case C60baroghelBc:{
    nu = Cmu->C60baroghelc[i].C60bar_nu();//from C60baroghelc.cpp
    break;
  }
  case C30baroghelBc:{
    nu = Cmu->C30baroghelc[i].C30bar_nu();//from C30baroghelc.cpp
    break;
  }
  case o30bazantBc:{
    nu = Cmu->o30bazantc[i].o30baz_nu();//from o30bazantc.cpp
    break;
  }                
  case C60bazantBc:{
    nu = Cmu->C60bazantc[i].C60baz_nu();//from o30bazantc.cpp
    break;
  }                
    //case glasgowc:{
    //nu = Cmu->tenchc[i].nu();//function from glasgowc.cpp//temporarily??!!
    //break;
    //}                
  default:{
    fprintf (stderr,"\n\n unknown material type is required in function (%s, line %d).\n",__FILE__,__LINE__);
  }
  }

  return(nu);
}
