/*
    File:             twomedia.cpp
    Author:           Tomas Krejci,  20.12.2002
    Purpose:          computes conductivity and capacity matrices in a material point for coupled two media trasport
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "twomedia.h"
#include "globalt.h"

med2::med2()
{  
  //not completed
  scale_w = 1.0;
  scale_t = 1.0;
}
med2::~med2()
{}

/**
   function computes conductivity %matrix D in a material point for two media transfer
   
   @param d - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   revision, JK, 24.11.2008
*/
void med2::matcond (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case isotransmat:{
    Tm->itrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case nlisotransmat:{
    Tm->nlitrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case homomat:{
    Tm->hommat[i].matcond(d,ri,ci,ipp);
    break;
  }
  case tdisotransmat:{
    Tm->tditrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case damisotransmat:{
    Tm->damitrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case bazantpedersen:{
    Tm->bazped[i].matcond(d,ri,ci,ipp);
    break;
  }
  case pedersen:{
    Tm->ped[i].matcond(d,ri,ci,ipp);
    break;
  }
  case kunzel:{
    Tm->kun[i].matcond(d,ri,ci,ipp);
    break;
  }
  case kunzel2:{
    Tm->kun2[i].matcond(d,ri,ci,ipp);
    break;
  }
  case moistheat:{
    Tm->moisth[i].matcond(d,ri,ci,ipp);
    break;
  }
  case grunewald:{
    Tm->grunw[i].matcond(d,ri,ci,ipp);
    break;
  }
  case simplediscmat:{
    Tm->sdmat[i].matcond(d,ri,ci,ipp);
    break;
  }
  case devries:{
    Tm->dvries[i].matcond(d,ri,ci,ipp);
    break;
  }
 case milly:{
    Tm->mill[i].matcond(d,ri,ci,ipp);
    break;
  }
  case lincoupledmat:{
    Tm->lcmat[i].matcond (d,ri,ci,ipp);
    break;
  }
  case salt1mat:{
    Tm->salt1[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolwf2:{
    Tm->consol_wf2[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolawf2:{
    Tm->consol_awf2[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolhwf2:{
    Tm->consol_hwf2[i].matcond(d,ri,ci,ipp);
    break;
  }
  case interfacem:{
    Tm->ifacemat[i].matcond(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes conductivity %matrix D in a material point for two media transfer
   
   @param d  - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   revision, JK, 24.11.2008
*/
void med2::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case isotransmat:{
    break;
  }
  case tdisotransmat:{
    break;
  }
  case nlisotransmat:{
    break;
  }
  case homomat:{
    break;
  }
  case damisotransmat:{
    break;
  }
  case bazantpedersen:{
    Tm->bazped[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case pedersen:{
    //Tm->ped[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case kunzel:{
    //Tm->kun[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case kunzel2:{
    //Tm->kun2[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case moistheat:{
    //Tm->moisth[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case grunewald:{
    //Tm->grunw[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case simplediscmat:{
    //Tm->sdmat[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case devries:{
    //Tm->dvries[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case milly:{
    //Tm->dvries[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case lincoupledmat:{
    break;
  }
  case salt1mat:{
    //Tm->salt1[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case consolwf2:{
    //Tm->consol_wf2[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case consolawf2:{
    //Tm->consol_awf2[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case consolhwf2:{
    Tm->consol_hwf2[i].matcond2(d,ri,ci,ipp);
    break;
  }
  case richardsmat:{
    Tm->richar[i].matcond2 (d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes capacity %matrix C in a material point for two media transfer
   
   @param c - capacity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   revision, JK, 24.11.2008
*/
void med2::matcap (double &c,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case isotransmat:{
    Tm->itrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case nlisotransmat:{
    Tm->nlitrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case homomat:{
    Tm->hommat[i].matcap(c,ri,ci,ipp);
    break;
  }
  case tdisotransmat:{
    Tm->tditrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case damisotransmat:{
    Tm->damitrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case bazantpedersen:{
    Tm->bazped[i].matcap(c,ri,ci,ipp);
    break;
  }
  case pedersen:{
    Tm->ped[i].matcap(c,ri,ci,ipp);
    break;
  }
  case kunzel:{
    Tm->kun[i].matcap(c,ri,ci,ipp);
    break;
  }
  case kunzel2:{
    Tm->kun2[i].matcap(c,ri,ci,ipp);
    break;
  }
  case moistheat:{
    Tm->moisth[i].matcap(c,ri,ci,ipp);
    break;
  }
  case grunewald:{
    Tm->grunw[i].matcap(c,ri,ci,ipp);
    break;
  }
  case simplediscmat:{
    Tm->sdmat[i].matcap(c,ri,ci,ipp);
    break;
  }
  case devries:{
    Tm->dvries[i].matcap(c,ri,ci,ipp);
    break;
  }
  case milly:{
    Tm->mill[i].matcap(c,ri,ci,ipp);
    break;
  }
  case lincoupledmat:{
    Tm->lcmat[i].matcap (c,ri,ci,ipp);
    break;
  }
  case salt1mat:{
    Tm->salt1[i].matcap(c,ri,ci,ipp);
    break;
  }
  case consolwf2:{
    Tm->consol_wf2[i].matcap(c,ri,ci,ipp);
    break;
  }
  case consolawf2:{
    Tm->consol_awf2[i].matcap(c,ri,ci,ipp);
    break;
  }
  case consolhwf2:{
    Tm->consol_hwf2[i].matcap(c,ri,ci,ipp);
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
void med2::rhs_volume (matrix &d,long ri, long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  
  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:
  case nlisotransmat:
  case homomat:
  case tdisotransmat:
  case damisotransmat:
  case bazantpedersen:
  case pedersen:
  case kunzel:
  case kunzel2:
  case moistheat:
  case grunewald:
  case simplediscmat:
  case devries:
  case milly:
  case lincoupledmat:
  case salt1mat:
  case interfacem:{
    break;
  }
  case consolwf2:{
    Tm->consol_wf2[i].rhs_volume(d,ri,ci,ipp);
    break;
  }
  case consolawf2:{
    Tm->consol_awf2[i].rhs_volume(d,ri,ci,ipp);
    break;
  }
  case consolhwf2:{
    Tm->consol_hwf2[i].rhs_volume(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n unknown material type is required in ", __FILE__, __LINE__, __func__);
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

void med2::rhs_volume2(double &c,long ri,long ci,long ipp)
{
  long i;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:
  case nlisotransmat:
  case homomat:
  case tdisotransmat:
  case damisotransmat:
  case bazantpedersen:
  case pedersen:
  case kunzel:
  case kunzel2:
  case moistheat:
  case grunewald:
  case simplediscmat:
  case devries:
  case milly:
  case lincoupledmat:
  case salt1mat:
  case interfacem:{
    break;
  }
  case consolwf2:{
    Tm->consol_wf2[i].rhs_volume2(c,ri,ci,ipp);
    break;
  }
  case consolawf2:{
    Tm->consol_awf2[i].rhs_volume2(c,ri,ci,ipp);
    break;
  }
  case consolhwf2:{
    Tm->consol_hwf2[i].rhs_volume2(c,ri,ci,ipp);
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

   @param trc - prescribed transmission coefficient on the boundary
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
   
   revision, JK, 24.11.2008
*/
double med2::transmission_transcoeff(double trc,long ri,long ci,long nid,long bc,long ipp)
{
  double new_trc;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case bazantpedersen:{
    new_trc = Tm->bazped[i].transmission_transcoeff (trc,ri,ci,nid,bc);
    break;
  }    
  case pedersen:{
    new_trc = Tm->ped[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  }
  case homomat:{
    new_trc = Tm->hommat[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  }
  case kunzel:{
    new_trc = Tm->kun[i].transmission_transcoeff (trc,ri,ci,nid,bc);
    break;
  }
  case kunzel2:{
    new_trc = Tm->kun2[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  } 
  case moistheat:{
    new_trc = Tm->moisth[i].transmission_transcoeff(trc,ri,ci,nid,bc);
    break;
  } 
  case grunewald:{
    new_trc = Tm->grunw[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  }      
  case simplediscmat:{
    new_trc = Tm->sdmat[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  }      
  case devries:{
    new_trc = Tm->dvries[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  }      
  case milly:{
    new_trc = Tm->mill[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break;
  }
  case lincoupledmat:{
    new_trc = trc;
    break;
  }
  case consolwf2:{
    new_trc = Tm->consol_wf2[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break; 
  }
  case consolawf2:{
    new_trc = Tm->consol_awf2[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
    break; 
  }
  case consolhwf2:{
    new_trc = Tm->consol_hwf2[i].transmission_transcoeff(trc,ri,ci,nid,bc,ipp);
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

   @param nodval     - prescribed transmission nodal value on the boundary
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if is needed (for example heat radiation coef.)
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element

   revision, JK, 24.11.2008
*/
double med2::transmission_nodval(double nodval,double trc2,long ri,long ci,long nid,long bc,long ipp)
{
  double new_nodval;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case bazantpedersen:{
    new_nodval = Tm->bazped[i].transmission_nodval (nodval,trc2,ri,ci,nid,bc);
    break;
  }
  case pedersen:{
    new_nodval = Tm->ped[i].transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;
  }
  case homomat:{
    new_nodval = Tm->hommat[i].transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;
  }
  case kunzel:{
    new_nodval = Tm->kun[i].transmission_nodval (nodval,ri,ci,nid,bc);
    break;
  }
  case kunzel2:{
    new_nodval = Tm->kun2[i].transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break;
  }
  case moistheat:{
    new_nodval = Tm->moisth[i].transmission_nodval(nodval,ri,ci,nid,bc);
    break;
  }
  case grunewald:{
    new_nodval = Tm->grunw[i].transmission_nodval(nodval,ri,ci,nid,bc,ipp);
    break;
  }
  case simplediscmat:{
    new_nodval = Tm->sdmat[i].transmission_nodval(nodval,ri,ci,nid,bc,ipp);
    break;
  }
  case devries:{
    new_nodval = Tm->dvries[i].transmission_nodval(nodval,ri,ci,nid,bc,ipp);
    break;
  }
  case milly:{
    new_nodval = Tm->mill[i].transmission_nodval(nodval,ri,ci,nid,bc,ipp);
    break;
  }
  case lincoupledmat:{
    new_nodval = nodval;
    break;
  }
  case consolwf2:{
    new_nodval = Tm->consol_wf2[i].transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break; 
  }
  case consolawf2:{
    new_nodval = Tm->consol_awf2[i].transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break; 
  }
  case consolhwf2:{
    new_nodval = Tm->consol_hwf2[i].transmission_nodval(nodval,trc2,ri,ci,nid,bc,ipp);
    break; 
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }  

  return new_nodval;
}


/**
   function computes flux on the boundary for transmission on the boundary
   (third kind of boundary condition)

   @param nodval  - prescribed transmission nodal value on the boundary
   @param trc2    - second prescribed transmission coefficient on the boundary, 
                    if is needed (for example heat radiation coef.)
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
   @param ipp - number of first integration point on element
   
   revision, JK, 24.11.2008
*/
double med2::transmission_flux (double nodval,double trc2,long ri,long ci,long nid,long bc,long ipp)
{
  double flux;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case bazantpedersen:{      
    flux = Tm->bazped[i].transmission_flux(nodval,trc2,ri,ci,nid,bc);
    break;
  }    
  case pedersen:{
    flux = Tm->bazped[i].transmission_flux(nodval,trc2,ri,ci,nid,bc);
    break;
  }
  case homomat:{
    flux = Tm->hommat[i].transmission_flux(nodval,trc2,ri,ci,nid,bc,ipp);
    break;
  }  
  case kunzel:{
    flux = Tm->kun[i].transmission_flux (nodval,ri,ci,nid,bc);
    break;
  }  
  case kunzel2:{
    flux = Tm->kun2[i].transmission_flux(nodval,trc2,ri,ci,nid,bc,ipp);
    break;
  }
  case grunewald:{
    flux = Tm->grunw[i].transmission_flux(nodval,ri,ci,nid,bc,ipp);
    break;
  }       
  case simplediscmat:{
    flux = Tm->sdmat[i].transmission_flux(nodval,ri,ci,nid,bc,ipp);
    break;
  }       
  case devries:{
    flux = Tm->dvries[i].transmission_flux(nodval,ri,ci,nid,bc,ipp);
    break;
  }       
  case milly:{
    flux = Tm->mill[i].transmission_flux(nodval,ri,ci,nid,bc,ipp);
    break;
  }
  case moistheat:{
    flux = Tm->moisth[i].transmission_flux (nodval,ri,ci,nid,bc);
    break;
  }
  case consolwf2:{
    flux = Tm->consol_wf2[i].transmission_flux(nodval,trc2,ri,ci,nid,bc,ipp);
    break; 
  }
  case consolawf2:{
    flux = Tm->consol_awf2[i].transmission_flux(nodval,trc2,ri,ci,nid,bc,ipp);
    break; 
  }
  case consolhwf2:{
    flux = Tm->consol_hwf2[i].transmission_flux(nodval,trc2,ri,ci,nid,bc,ipp);
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
   @param ipp - first integration point on element
   @param r - %vector of unknowns on actual node
*/

double med2::compute_othervalues (long compother,long ipp,double *r)
{
  long i;
  double other;

  //  material id
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:
    break;
  case consolwf2:{
    double pw,pg;
    
    pw = r[0];
    pg = r[1];
    
    other = Tm->consol_wf2[i].get_othervalue(compother,pw,pg,ipp);
    break;
  }
  case consolawf2:{
    double pw,pg;
    
    pw = r[0];
    pg = r[1];
    
    other = Tm->consol_awf2[i].get_othervalue(compother,pw,pg,ipp);
    break;
  }
  case consolhwf2:{
    double pw,t;
    
    pw = r[0];
    t = r[1];
    
    other = Tm->consol_hwf2[i].get_othervalue(compother,ipp,r);
    break;
  }
  case nlisotransmat:{
    break;
  }
  case bazantpedersen:{
    double w,t;
    
    w = r[0];
    t = r[1];
    
    other = Tm->bazped[i].get_othervalue(compother,w,t);
    break;
  }
  case pedersen:{
    double w,t;
    
    w = r[0];
    t = r[1];
    
    other = Tm->ped[i].get_othervalue(compother,w,t);
    break;
  }
  case homomat:{
    double rh,t;
    
    rh = r[0];
    t = r[1];
    
    other = Tm->hommat[i].get_othervalue(compother,rh,t,ipp);
    break;
  }
  case kunzel:{
    double rh,t;
    
    rh = r[0];
    t = r[1];
    
    other = Tm->kun[i].get_othervalue(compother,rh,t,ipp);
    break;
  }
  case kunzel2:{
    double rh,t;
    
    rh = r[0];
    t = r[1];
    
    other = Tm->kun2[i].get_othervalue(compother,rh,t,ipp);
    break;
  }
  case moistheat:{
    double rh,t;
    
    rh = r[0];
    t = r[1];
    
    other = Tm->moisth[i].get_othervalue(compother,rh,t,ipp);
    break;
  }
  case grunewald:{
    double w,t;
    
    w = r[0];
    t = r[1];
    
    other = Tm->grunw[i].get_othervalue(compother,ipp,w,t,0);
    break;
  }
  case simplediscmat:{
    double w,t;
    
    w = r[0];
    t = r[1];
    
    other = Tm->sdmat[i].get_othervalue(compother,ipp,w,t);
    break;
  }
  case devries:{
    double w,t;
    
    w = r[0];
    t = r[1];
    
    other = Tm->dvries[i].get_othervalue(compother,ipp,w,t,0);
    break;
  }
  case milly:{
    double psi,t;
    
    psi = r[0];
    t = r[1];
    
    other = Tm->mill[i].get_othervalue(compother,ipp,psi,t);
    break;
  }
  case salt1mat:{
    double w,cf;
    
    w = r[0];
    cf = r[1];
    
    //other = Tm->bazped[i].get_othervalue(compother,w,t);
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
   @param out - output file
   @param compother - number of other components
   @param ipp - first integration point on element
*/

void med2::print_othervaluesnames (FILE *out,long ipp,long compother)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case nlisotransmat:{
    break;
  }
  case bazantpedersen:{
    Tm->bazped[i].print_othervalue_name(out,compother);
    break;
  }
  case consolwf2:{
    Tm->consol_wf2[i].print_othervalue_name(out,compother);
    break;
  }
  case consolawf2:{
    Tm->consol_awf2[i].print_othervalue_name(out,compother);
    break;
  }
  case consolhwf2:{
    Tm->consol_hwf2[i].print_othervalue_name(out,compother);
    break;
  }
  case pedersen:{
    Tm->ped[i].print_othervalue_name(out,compother);
    break;
  }
  case homomat:{
    Tm->hommat[i].print_othervalue_name(out,compother);
    break;
  }
  case kunzel:{
    Tm->kun[i].print_othervalue_name(out,compother);
    break;
  }
  case moistheat:{
    Tm->moisth[i].print_othervalue_name(out,compother);
    break;
  }
  case kunzel2:{
    Tm->kun2[i].print_othervalue_name(out,compother);
    break;
  }
  case grunewald:{
    Tm->grunw[i].print_othervalue_name(out,compother);
    break;
  }
  case simplediscmat:{
    Tm->sdmat[i].print_othervalue_name(out,compother);
    break;
  }
  case devries:{
    Tm->dvries[i].print_othervalue_name(out,compother);
    break;
  }
  case milly:{
    Tm->mill[i].print_othervalue_name(out,compother);
    break;
  }
  case salt1mat:{
    Tm->salt1[i].print_othervalue_name(out,compother);
    break;
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }

}

