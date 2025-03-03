/*
    File:             threemedia.cpp
    Author:           Tomas Krejci, 20.12.2002
    Purpose:          computes conductivity and capacity matrices in a material point for coupled three media transport
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "threemedia.h"
#include "multiphase.h"
#include "gmultiphase.h"
#include "constrel.h"
#include "globalt.h"
#include "glasgowmat.h"
//#include "meshtransfert.h"

med3::med3()
{}
med3::~med3()
{}

/**
   function computes conductivity %matrix D in a material point for three media transfer
   
   @param d - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   revision, JK, 24.11.2008
*/
void med3::matcond (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    mtph.matcond(d,ri,ci,ipp);
    break;
  }    
  case glasgow:{
    Tm->tench[i].matcond(d,ri,ci,ipp);
    break;
  }
  case soilmat1:{
    gmultiph gmtph;
    gmtph.matcond(d,ri,ci,ipp);
    break;
  }    
  case lincoupledmat:{
    Tm->lcmat[i].matcond (d,ri,ci,ipp);
    break;
  }
    
  case salt2mat:{
    Tm->salt2[i].matcond(d,ri,ci,ipp);
    break;
  }
  case salt3mat:{
    Tm->salt3[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolhawf3:{
    Tm->consol_hawf3[i].matcond(d,ri,ci,ipp);
    break; 
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }  

}

/**
   function computes second type of conductivity %matrix D in a material point for three media transfer
   
   @param d - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   revision, JK, 24.11.2008
*/
void med3::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    mtph.matcond2(d,ri,ci,ipp);
    break;
  }    
  case glasgow:{
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    gmtph.matcond2(d,ri,ci,ipp);
    break;
  }    
  case salt2mat:{
    Tm->salt2[i].matcond2 (d,ri,ci,ipp);
    break;
  }
  case salt3mat:{
    Tm->salt3[i].matcond2 (d,ri,ci,ipp);
    break;
  }
  case consolhawf3:{
    Tm->consol_hawf3[i].matcond2(d,ri,ci,ipp);
    break; 
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}
 
/**
   function computes capacity %matrix C in a material point for three media transfer
   
   @param c - capacity coefficient
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   revision, JK, 24.11.2008
*/
void med3::matcap (double &c,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    mtph.matcap(c,ri,ci,ipp);
    break;
  }    
  case glasgow:{
    Tm->tench[i].matcap(c,ri,ci,ipp);
    break;
  }
  case soilmat1:{
    gmultiph gmtph;
    gmtph.matcap(c,ri,ci,ipp);
    break;
  }    
  case lincoupledmat:{
    Tm->lcmat[i].matcap (c,ri,ci,ipp);
    break;
  }
  case salt2mat:{
    Tm->salt2[i].matcap(c,ri,ci,ipp);
    break;
  }
  case salt3mat:{
    Tm->salt3[i].matcap(c,ri,ci,ipp);
    break;
  }
  case consolhawf3:{
    Tm->consol_hawf3[i].matcap(c,ri,ci,ipp);
    break; 
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function computes volume part of right-hand
   in the required integration point
   
   @param d - right-hand side %matrix of a material
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void med3::rhs_volume (matrix &d,long ri, long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    
    mtph.rhs_volume(d,ri,ci,ipp);
    break;
  }    
  case glasgow:{
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    gmtph.rhs_volume(d,ri,ci,ipp);
    break;
  }    
    
  case salt2mat:{
    //  Tm->salt2[i].matcond(d,ri,ci,ipp);
    break;
  }
  case salt3mat:{
    //  Tm->salt3[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolhawf3:{
    Tm->consol_hawf3[i].rhs_volume(d,ri,ci,ipp);
    break; 
  }
    
  default:{
    print_err ("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }

}


/**
   function computes volume part 2 of right-hand
   in the required integration point
   
   @param d - right-hand side %matrix of a material
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/

void med3::rhs_volume2(double &c,long ri,long ci,long ipp)
{
  long i;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    //multiph mtph;

    //mtph.rhs_volume2(d,ri,ci,ipp);//must be completed
    break;
  }
  case consolhawf3:{
    Tm->consol_hawf3[i].rhs_volume2(c,ri,ci,ipp);
    break; 
  }
  default:{
    print_err("\n unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}




/**
   function computes new transmission coefficient for transmission on the boundary
   (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double med3::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp)
{
  double new_trc;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    
    new_trc = mtph.transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }    
  case glasgow:{
    new_trc = Tm->tench[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    new_trc = gmtph.transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }    
  case salt2mat:{
    new_trc = Tm->salt2[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }
  case salt3mat:{
    new_trc = Tm->salt3[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }
  case consolhawf3:{
    new_trc = Tm->consol_hawf3[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break; 
  }

  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(new_trc);
}



/**
   function computes new transmission coefficient for transmission on the boundary
   (third kind of boundary condition)

   @param trc     - prescribed transmission coefficient on the boundary
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
   @param flag - coefficient is computing for what 0=matrix,1=loading vector
*/
double med3::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp, int flag)
{
  double new_trc;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    
    new_trc = mtph.transmission_transcoeff(trc,ri,ci,nn,bc,ipp,flag);
    break;
  }    
  case glasgow:{
    new_trc = Tm->tench[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    new_trc = gmtph.transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }    
  case salt2mat:{
    new_trc = Tm->salt2[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }
  case salt3mat:{
    new_trc = Tm->salt3[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break;
  }
  case consolhawf3:{
    new_trc = Tm->consol_hawf3[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break; 
  }    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(new_trc);
}


/**
   function computes new nodal value for transmission on the boundary
   (third kind of boundary condition)

   @param nodval     - prescribed value on the boundary
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - number of node
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double med3::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double new_nodval;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    
    new_nodval = mtph.transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }    
  case glasgow:{
    new_nodval = Tm->tench[i].transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    new_nodval = gmtph.transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }    
  case salt2mat:{
    new_nodval = Tm->salt2[i].transmission_nodval(nodval,ri,ci,nn,bc,ipp);
    break;
  }
  case salt3mat:{
    new_nodval = Tm->salt3[i].transmission_nodval(nodval,ri,ci,nn,bc,ipp);
    break;
  }
  case consolhawf3:{
    new_nodval = Tm->consol_hawf3[i].transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break; 
  }
    
  default:  {
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(new_nodval);
}


/**
   function computes flux on the boundary for transmission on the boundary
   (third kind of boundary condition)

   @param nodval  - prescribed value on the boundary
   @param trc2    - second prescribed transmission coefficient on the boundary, 
                    if is needed (for example heat radiation coef.)
   @param ri      - row index
   @param ci      - column index
   @param nn      - number of node
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double med3::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  double flux;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    
    flux = mtph.transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }    
  case glasgow:{
    flux = Tm->tench[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    flux = gmtph.transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }    
  case salt2mat:{
    flux = Tm->salt2[i].transmission_flux(nodval,ri,ci,nn,bc,ipp);
    break;
  }
  case salt3mat:{
    flux = Tm->salt3[i].transmission_flux(nodval,ri,ci,nn,bc,ipp);
    break;
  }
  case consolhawf3:{
    flux = Tm->consol_hawf3[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break; 
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(flux);
}



//this part below will be corrected and changed!!??

/**
   function computes all variables in nodes

   @param compother - number of other components
   @param ipp       - first integration point on element
   @param r         - %vector of unknowns on actual node
*/

double med3::compute_othervalues (long compother,long ipp,double *r)
{
  double other;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){
  case glasgow:{
    other = Tm->tench[Tm->ip[ipp].idm].get_othervalue(compother,ipp,r);
    break;
  }
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:
    {
      multiph mtph;
      
      other = mtph.get_othervalue(compother,ipp,r);
      break;
    }
    
  case soilmat1:{
    gmultiph gmtph;
    other = gmtph.get_othervalue(compother,ipp,r);
    break;
  }    
    
  case salt2mat:{
    // new_trc = Tm->salt2[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    other=0;
    break;
  }
  case salt3mat:{
    // new_trc = Tm->salt3[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    other=0;
    break;
  }
  case consolhawf3:{
    other=Tm->consol_hawf3[i].get_othervalue(compother,ipp,r);
    break; 
  }

  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return(other);
}


/**
   function prints names of all variables in nodes

   @param out       - output file
   @param compother - number of other components
   @param ipp       - first integration point on element
*/

void med3::print_othervaluesnames (FILE *out,long ipp,long compother)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){
  case glasgow:{
    Tm->tench[Tm->ip[ipp].idm].print_othervalue_name(out,compother);
    break;
  }
  case concreteB:
  case baroghelB:
  case C60baroghelB:
  case C30baroghelB:
  case o30bazantB:
  case C30bazantB:
  case C60bazantB:{
    multiph mtph;
    
    mtph.print_othervalue_name(out,compother);
    break;
  }
    
  case soilmat1:{
    gmultiph gmtph;
    gmtph.print_othervalue_name(out,compother);
    break;
  }    
  case consolhawf3:{
    Tm->consol_hawf3[i].print_othervalue_name(out,compother);
    break; 
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}


