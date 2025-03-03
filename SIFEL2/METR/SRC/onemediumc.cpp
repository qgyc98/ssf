/*
    File:             onemediumc.cpp
    Author:           Tomas Krejci
    Purpose:          computes conductivity and capacity matrices in a material point for coupled three media trasport with mechanics
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "onemediumc.h"
#include "coupmatu.h"
#include "coupmatl.h"
#include "globalc.h"
#include "global.h"
#include "globalt.h"
#include "isotrmatc.h"
#include "intpoints.h"
#include "consol_awf1.h"
#include "consol_wf1.h"

medc1::medc1()
{  
  /* scale_t = Tp->scale[0];
     scale_u = Mp->scale;
  */
  //not completed
  scale_t = 1.0;
  scale_u = 1.0;
}
medc1::~medc1()
{}


/**
   function computes conductivity matrix D in a material point for one medium transfer
   for upper block
   
   @param d   - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

*/
void medc1::matcond_u (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//material type
  case isotransmatc:{
    Cmu->itrmc[i].matcond_u(d,ri,ci,ipp);
    break;
  }
  case consolawf1c:{
    Cmu->consol_awf1c[i].matcond_u(d,ri,ci,ipp);
    break;
  }
  case consolwf1c:{
    Cmu->consol_wf1c[i].matcond_u(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}




/**
   function computes capacity matrix D in a material point for one medium transfer
   for upper block
   
   @param d   - cpapcity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

*/
void medc1::matcap_u (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//material type
  case isotransmatc:{
    Cmu->itrmc[i].matcap_u(d,ri,ci,ipp);
    break;
  }
  case consolawf1c:{
    Cmu->consol_awf1c[i].matcap_u(d,ri,ci,ipp);
    break;
  }
   case consolwf1c:{
    Cmu->consol_wf1c[i].matcap_u(d,ri,ci,ipp);
    break;
  }
 default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function computes conductivity matrix D in a material point for one medium transfer
   for lower block
   
   @param d   - conductivity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

*/
void medc1::matcond_l (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//material type
  case isotransmatc:{
    Cml->itrmc[i].matcond_l(d,ri,ci,ipp);
    break;
  }
  case consolawf1c:{
    Cml->consol_awf1c[i].matcond_l(d,ri,ci,ipp);
    break;
  }
  case consolwf1c:{
    Cml->consol_wf1c[i].matcond_l(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}




/**
   function computes capacity matrix D in a material point for one medium transfer
   for lower block
   
   @param d   - cpapcity %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

*/
void medc1::matcap_l (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//material type
  case isotransmatc:{
    Cml->itrmc[i].matcap_l(d,ri,ci,ipp);
    break;
  }
  case consolawf1c:{
    Cml->consol_awf1c[i].matcap_l(d,ri,ci,ipp);
    break;
  }
  case consolwf1c:{
    Cml->consol_wf1c[i].matcap_l(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   function computes right-hand-side matrix D in a material point for one medium transfer - influence of gravity accel.
   for upper block
   
   @param d   - right-hand-side %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

*/
void medc1::rhs_u1 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//material type
  case isotransmatc:{
    break;
  }
  case consolawf1c:{
    Cmu->consol_awf1c[i].rhs_u1(d,ri,ci,ipp);
    break;
  }
  case consolwf1c:{
    Cmu->consol_wf1c[i].rhs_u1(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   function computes right-hand-side matrix D in a material point for one medium transfer - influence of initial values (temperature)
   for upper block
   
   @param d   - right-hand-side %matrix
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point

*/
void medc1::rhs_u2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  i = Cmu->ip[ipp].idm;
  
  switch (Cmu->ip[ipp].tm){//material type
  case isotransmatc:{
    Cmu->itrmc[i].rhs_u2(d,ri,ci,ipp);
    break;
  }
  case consolawf1c:{
    break;
  }
  case consolwf1c:{
    break;
  }
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}
