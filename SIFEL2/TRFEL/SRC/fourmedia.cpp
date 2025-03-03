/*
    File:             fourmedia.cpp
    Author:           Tomas Krejci, 20.12.2002, Jaroslav Kruis, 23.7.2007
    Purpose:          computes conductivity and capacity matrices in a material point for coupled four media transport
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "fourmedia.h"
#include "multiphase.h"
#include "gmultiphase.h"
#include "constrel.h"
#include "globalt.h"
#include "glasgowmat.h"
//#include "meshtransfert.h"

med4::med4()
{
  //not completed
  scale_w = 1.0;
  scale_t = 1.0;
}

med4::~med4()
{}

/**
   function computes conductivity %matrix D in a material point for four media transfer

   @param d - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void med4::matcond (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    Tm->salt4[i].matcond(d,ri,ci,ipp);
    break;
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }  

}

/**
   function computes second type of conductivity %matrix D in a material point for four media transfer
   
   @param d - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void med4::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    Tm->salt4[i].matcond2 (d,ri,ci,ipp);
    break;
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }
}
 
/**
   function computes capacity %matrix C in a material point for four media transfer

   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void med4::matcap (double &c,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;

  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    Tm->salt4[i].matcap(c,ri,ci,ipp);
    break;
  }
    
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
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
*/
double med4::transmission_transcoeff(double trc,long ri,long ci,long nid,long bc,long ipp)
{
  double new_trc;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    new_trc = Tm->salt4[i].transmission_transcoeff (trc,ri,ci,nid,bc);
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
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double med4::transmission_nodval(double nodval,double /*trc2*/,long ri,long ci,long nid,long bc,long ipp)
{
  double new_nodval;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    new_nodval = Tm->salt4[i].transmission_nodval (nodval,ri,ci,nid,bc);
    break;
  }
  default:{
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
   @param ri - row index
   @param ci - column index
   @param nid - node id
   @param bc - type of boundary condition
*/
double med4::transmission_flux(double nodval,double /*trc2*/,long ri,long ci,long nid,long bc,long ipp)
{
  double flux;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    flux = Tm->salt4[i].transmission_flux (nodval,ri,ci,nid,bc);
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

double med4::compute_othervalues (long compother,long ipp,double *r)
{
  double other;
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case salt4mat:{
    
    double w,cf,cfmax,t;
    w = r[0];
    cf = r[1];
    cfmax = r[2];
    t = r[3];
    
    w = w*scale_w;//scaling
    t = t*scale_w;//scaling
    cf = cf*scale_w;//scaling
    cfmax = cfmax*scale_w;//scaling
    
    
    other = Tm->salt4[i].get_othervalue(compother,ipp,w,cf,cfmax,t);
    
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

void med4::print_othervaluesnames (FILE */*out*/,long /*ipp*/,long /*compother*/)
{
/*  switch (Tp->mednam){//names of transported media
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }*/
}

