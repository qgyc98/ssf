/*
    File:             onemedium.cpp
    Author:           Tomas Krejci,  20.12.2002
    Purpose:          computes conductivity and capacity matrices in a material point for one medium trasport
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aliast.h"
#include "onemedium.h"
#include "globalt.h"

med1::med1()
{  
  //not completed
  scale = 1.0;
}
med1::~med1()
{}


/**
   function computes conductivity %matrix D in a material point for one medium transfer
   
   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
*/
void med1::matcond (matrix &d,long ri,long ci,long ipp)
{
  long i;

  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case bazantpedersen:{//material type
    Tm->bazped[i].matcond(d,ri,ci,ipp);
    break;
  }
  case pedersen:{
    Tm->ped[i].matcond(d,ri,ci,ipp);
    break;
  }
    
  case carb1mat:{
    Tm->carb1[i].matcond(d,ri,ci,ipp);
    break;
  }  
  case sejtkr:{
    Tm->sejtkrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolawf1:{
    Tm->consol_awf1[i].matcond(d,ri,ci,ipp);
    break;
  }
  case consolwf1:{
    Tm->consol_wf1[i].matcond(d,ri,ci,ipp);
    break;
  }
  case richardsmat:{
    Tm->richar[i].matcond (d,0,0,ipp);
    break;
  } 
  case isotransmat:{
    Tm->itrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case nlisotransmat:{
    Tm->nlitrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case discontisotrmat:{
    Tm->ditrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case tdisotransmat:{
    Tm->tditrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case interfacem:{
    Tm->ifacemat[i].matcond(d,ri,ci,ipp);
    break;
  }
  case cernyconcrete:{
    Tm->cernym[i].matcond(d,ri,ci,ipp);
    break;
  }
  case cementhydrmat:{
    Tm->cemhydr[i].matcond(d,ri,ci,ipp);
    break;
  }
  case lincoupledmat:{
    Tm->lcmat[i].matcond(d,ri,ci,ipp);
    break;
  }
  case damisotransmat:{
    Tm->damitrm[i].matcond(d,ri,ci,ipp);
    break;
  }
  case radiationmater:{
    Tm->radmat[i].matcond(d,ri,ci,ipp);
    break;
  }
  case homomat:{
    Tm->hommat[i].matcond(d,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n unknown material type is required in ", __FILE__, __LINE__, __func__);
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
void med1::matcond2 (matrix &d,long ri,long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case isotransmat:{
    //Tm->itrm[i].matcond(d,ri,ci,ipp);//zde pro cylindr??!!
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
  case tdisotransmat:{
    break;
  }
  case sejtkr:{
    break;
  }
  case consolawf1:{
    break;
  }
  case consolwf1:{
    break;
  }
  case bazantpedersen:{
    //Tm->bazped[i].matcond2(d,ri,ci,ipp);
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
  case consolawf2:{
    //Tm->consol_awf2[i].matcond2(d,ri,ci,ipp);
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
   function computes capacity %matrix C in a material point for one medium transfer

   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
*/
void med1::matcap (double &c,long ri,long ci,long ipp)
{
  long i;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case bazantpedersen:{
    Tm->bazped[i].matcap(c,ri,ci,ipp);    
    break;
  }
  case pedersen:{
    Tm->ped[i].matcap(c,ri,ci,ipp);
    break;
  }
  case carb1mat:{
    Tm->carb1[i].matcap(c,ri,ci,ipp);
    break;
  }
  case richardsmat:{
    Tm->richar[i].matcap (c,ri,ci,ipp);
    break;
  }
  case consolawf1:{
    Tm->consol_awf1[i].matcap(c,ri,ci,ipp);
    break;
  }
  case consolwf1:{
    Tm->consol_wf1[i].matcap(c,ri,ci,ipp);
    break;
  }
  case sejtkr:{
    Tm->sejtkrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case isotransmat:{
    Tm->itrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case nlisotransmat:{
    Tm->nlitrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case damisotransmat:{
    Tm->damitrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case tdisotransmat:{
    Tm->tditrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case discontisotrmat:{
    Tm->ditrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case cernyconcrete:{
    Tm->cernym[i].matcap(c,ri,ci,ipp);
    break;
  }
  case cementhydrmat:{
    Tm->cemhydr[i].matcap(c,ri,ci,ipp);
    break;
  }
  case lincoupledmat:{
    Tm->lcmat[i].matcap(c,ri,ci,ipp);
    break;
  }
  case radiationmater:{
    Tm->radmat[i].matcap(c,ri,ci,ipp);
    break;
  }
  case homomat:{
    Tm->hommat[i].matcap(c,ri,ci,ipp);
    break;
  }
  default:{
    print_err("\n unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}

/**
   function computes reaction coefficient R in a material point for one medium transfer

   @param ri  - row index
   @param ci  - column index
   @param ipp - number of integration point
   
   JK, 12. 6. 2019
*/
void med1::matreact (double &r,long ri,long ci,long ipp)
{
  long i;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:{
    Tm->itrm[i].matreact (r,ri,ci,ipp);
    break;
  }
    /*
  case bazantpedersen:{
    Tm->bazped[i].matcap(c,ri,ci,ipp);    
    break;
  }
  case pedersen:{
    Tm->ped[i].matcap(c,ri,ci,ipp);
    break;
  }
  case carb1mat:{
    Tm->carb1[i].matcap(c,ri,ci,ipp);
    break;
  }
  case richardsmat:{
    Tm->richar[i].matcap (c,ri,ci,ipp);
    break;
  }
  case consolawf1:{
    Tm->consol_awf1[i].matcap(c,ri,ci,ipp);
    break;
  }
  case sejtkr:{
    Tm->sejtkrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case nlisotransmat:{
    Tm->nlitrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case damisotransmat:{
    Tm->damitrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case tdisotransmat:{
    Tm->tditrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case discontisotrmat:{
    Tm->ditrm[i].matcap(c,ri,ci,ipp);
    break;
  }
  case cernyconcrete:{
    Tm->cernym[i].matcap(c,ri,ci,ipp);
    break;
  }
  case cementhydrmat:{
    Tm->cemhydr[i].matcap(c,ri,ci,ipp);
    break;
  }
  case lincoupledmat:{
    Tm->lcmat[i].matcap(c,ri,ci,ipp);
    break;
  }
  case radiationmater:{
    Tm->radmat[i].matcap(c,ri,ci,ipp);
    break;
  }
  case homomat:{
    Tm->hommat[i].matcap(c,ri,ci,ipp);
    break;
  }
    */
  default:{
    print_err("\n unknown material type is required in ", __FILE__, __LINE__, __func__);
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
void med1::rhs_volume (matrix &d,long ri, long ci,long ipp)
{
  long i;
  
  //  material id
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){
  case sejtkr:
    Tm->sejtkrm[i].rhs_volume(d,ri,ci,ipp);
    break;
  case consolawf1:
    Tm->consol_awf1[i].rhs_volume(d,ri,ci,ipp);
    break;
  case consolwf1:
    Tm->consol_wf1[i].rhs_volume(d,ri,ci,ipp);
    break;
  case isotransmat:
  case tdisotransmat:
  case nlisotransmat:
  case discontisotrmat:
  case cernyconcrete:
  case damisotransmat:
  case cementhydrmat:
  case lincoupledmat:
  case bazantpedersen:
  case pedersen:
  case carb1mat:
  case richardsmat:
  case homomat:
    break;
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

void med1::rhs_volume2(double &c,long ri,long ci,long ipp)
{
  long i;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:
  case tdisotransmat:
  case nlisotransmat:
  case discontisotrmat:
  case cernyconcrete:
  case damisotransmat:
  case cementhydrmat:
  case lincoupledmat:
  case bazantpedersen:
  case pedersen:
  case carb1mat:
  case richardsmat:
  case homomat:
    break;
  case consolawf1:{
    Tm->consol_awf1[i].rhs_volume2(c,ri,ci,ipp);
    break;
  }
  case consolwf1:{
    Tm->consol_wf1[i].rhs_volume2(c,ri,ci,ipp);
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
   @param nn      - node id
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
   @param flag       - coefficient is computing for what 0=matrix,1=loading vector
*/
double med1::transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag)
{
  long i;
  double new_trc;
  
  //new:
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case isotransmat:{      
    new_trc = Tm->itrm[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp,flag);
    break;
  }
  case consolawf1:{
    new_trc = Tm->consol_awf1[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break; 
  }
  case consolwf1:{
    new_trc = Tm->consol_wf1[i].transmission_transcoeff(trc,ri,ci,nn,bc,ipp);
    break; 
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }  

  //old  
  /* 
     switch (Tp->mednam){//names of transported media
     case moisture:{
     
     switch (bc){//type of prescribed variable
     case 30:{//mass transmission
     new_trc=1.0;
     break;
     }
     default:{
     print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
     exit(0);
     }
     }
     
     break;
     }
     case heat:{
     
     switch (bc){//type of prescribed variable
     case 11:{
     new_trc=trc;
     }
     case 30:{//heat transmission
     new_trc=1.0;
     break;
     }
     case 31:{//heat transmission for testing
     new_trc=0.0;
     break;
     }
     case 90:{//radiation
     if(flag == 1)//for rhs vector
     new_trc=1.0;
     else
     new_trc=0.0;//for matrix
     break;
     }
     default:{
     print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
     exit(0);
     }
     }
     
     break;
     }
     default:{
     print_err("\n Unknown media name is required in ", __FILE__, __LINE__, __func__);
     }
     } 
     new_trc = new_trc*trc;
  
  */
 
  return(new_trc);
}



/**
   function computes new nodal value for transmission on the boundary
   (third kind of boundary condition)

   @param nodval     - nodal value defined on boundary
   @param trc2       - second prescribed transmission coefficient on the boundary, 
                       if it is needed (for example heat radiation coef.)
   @param ri         - row index
   @param ci         - column index
   @param nn         - node id
   @param bc         - type of boundary condition
   @param ipp        - number of first integration point on element
*/
double med1::transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  long i;
  double new_nodval;
  //long k;
  //double t,w;
  
  //new:
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case isotransmat:{      
    new_nodval = Tm->itrm[i].transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }
  case consolawf1:{
    new_nodval = Tm->consol_awf1[i].transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break; 
  }
  case consolwf1:{
    new_nodval = Tm->consol_wf1[i].transmission_nodval(nodval,trc2,ri,ci,nn,bc,ipp);
    break; 
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }  

  //old  
  /* switch (Tp->mednam){//names of transported media
     case moisture:{
     
     k=Gtt->give_dof(nn,0);
     if (k>0)   {w = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
     if (k==0)  {w = 0.0;}
     if (k<0)   {w = Tb->lc[0].pv[0-k-1].getval();}
     
     switch (bc){//type of prescribed variable
     case 30:{//mass transmission
     new_nodval = nodval;
     break;
     }
     default:{
     print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
     exit(0);
     }
     }
     
     break;
     }
     case heat:{
     
     k=Gtt->give_dof(nn,0);
     if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
     if (k==0)  {t = 0.0;}
     if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval();}
     
     switch (bc){//type of prescribed variable
     case 11:{
     new_nodval = nodval;
     break;
     }
     case 30:{//heat transmission
     new_nodval = nodval;
     break;
     }
     case 31:{//heat transmission for testing (and for boundary flux)
     new_nodval = (nodval - t);
     break;
     }
     case 90:{//radiation
     //new_nodval = (nodval - t) + trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
     new_nodval = trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
     break;
     }
     default:{
     print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
     exit(0);
     }
     }
     
     break;
     }
     default:{
     print_err("\n Unknown media name is required in ", __FILE__, __LINE__, __func__);
     }
     }   
  */
 
  return(new_nodval);
}



/**
   function computes flux on the boundary for transmission on the boundary
   (third kind of boundary condition)

   @param nodval  - prescribed nodal value on boundary
   @param trc2    - second prescribed transmission coefficient on the boundary, 
                    if it is needed (for example heat radiation coef.)
   @param ri      - row index
   @param ci      - column index
   @param nn      - node id
   @param bc      - type of boundary condition
   @param ipp     - number of first integration point on element
*/
double med1::transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp)
{
  long i;
  double flux;
  //long k;
  //double t,w;
  
  //new:
  //  material id
  i = Tm->ip[ipp].idm;
        
  switch (Tm->ip[ipp].tm){
  case isotransmat:{      
    flux = Tm->itrm[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break;
  }
  case consolawf1:{
    flux = Tm->consol_awf1[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break; 
  }
  case consolwf1:{
    flux = Tm->consol_wf1[i].transmission_flux(nodval,trc2,ri,ci,nn,bc,ipp);
    break; 
  }
  default:{
    print_err("unknown material type is required",__FILE__,__LINE__,__func__);
  }
  }  
 

  //old:
  /* 
     switch (Tp->mednam){//names of transported media
     case moisture:{
     
     k=Gtt->give_dof(nn,0);
     if (k>0)   {w = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
     if (k==0)  {w = 0.0;}
     if (k<0)   {w = Tb->lc[0].pv[0-k-1].getval();}
     
     switch (bc){//type of prescribed variable
     case 30:{//mass transmission
     flux = (nodval - w);
     break;
     }
     default:{
     print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
     exit(0);
     }
     }
     
     break;
     }
     case heat:{
     
     k=Gtt->give_dof(nn,0);
     if (k>0)   {t = Lsrst->lhs[k-1]+Lsrst->lhsi[k-1];}
     if (k==0)  {t = 0.0;}
     if (k<0)   {t = Tb->lc[0].pv[0-k-1].getval();}
     
     switch (bc){//type of prescribed variable
     case 30:{//heat transmission
     flux = (nodval - t);
     break;
     }
     case 31:{//heat transmission for testing (and for boundary flux)
     flux = (nodval - t);
     break;
     }
     case 90:{//radiation
     //flux = (nodval - t) + trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
     flux = trc2*(nodval*nodval*nodval*nodval - t*t*t*t);//inverse eq. (trc2 = e_sigma0/trc)
     break;
     }
     default:{
     print_err("\n No real boundary condition is prescribed in ", __FILE__, __LINE__, __func__);
     exit(0);
     }
     }
     
     break;
     }
     default:{
     print_err("\n Unknown media name is required in ", __FILE__, __LINE__, __func__);
     }
     }    
  */

  return(flux);

}



/**
   function computes all variables in nodes
   @param compother - number of other components
   @param ipp - first integration point on element
   @param r - %vector of unknowns on actual node
*/

double med1::compute_othervalues (long compother,long ipp,double *r)
{
  long i;
  double other;

  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:{
    double t;

    t = r[0];

    other = Tm->itrm[i].get_othervalue(compother,t,ipp);
    break;

  }
  case sejtkr:{
    double pw;
    
    pw = r[0];
    
    other = Tm->sejtkrm[i].get_othervalue(compother,pw,ipp);
    break;
  }
  case consolawf1:{
    double pw;
    
    pw = r[0];
    
    other = Tm->consol_awf1[i].get_othervalue(compother,pw,ipp);
    break;
  }
  case consolwf1:{
    double pw;
    
    pw = r[0];
    
    other = Tm->consol_wf1[i].get_othervalue(compother,pw,ipp);
    break;
  }
  case bazantpedersen:{
    double w,t;
    
    w = r[0];
    t = 0.0;
    
    other = Tm->bazped[i].get_othervalue(compother,w,t);
    break;
  }
  case pedersen:{
    double w,t;
    
    w = r[0];
    t = 0.0;
    
    other = Tm->ped[i].get_othervalue(compother,w,t);
    break;
  }
    
  case richardsmat:{
    break;
  }
    
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
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

void med1::print_othervaluesnames (FILE *out,long ipp,long compother)
{
  long i;
  
  i = Tm->ip[ipp].idm;
  
  switch (Tm->ip[ipp].tm){//material type
  case isotransmat:{
    Tm->itrm[i].print_othervalue_name(out,compother);
    break;
  }
  case bazantpedersen:{
    Tm->bazped[i].print_othervalue_name(out,compother);
    break;
  }
  case pedersen:{
    Tm->ped[i].print_othervalue_name(out,compother);
    break;
  }
  case sejtkr:{
    Tm->sejtkrm[i].print_othervalue_name(out,compother);
    break;
  }
  case consolawf1:{
    Tm->consol_awf1[i].print_othervalue_name(out,compother);
    break;
  }
  case consolwf1:{
    Tm->consol_wf1[i].print_othervalue_name(out,compother);
    break;
  }
  case damisotransmat:
    Tm->damitrm[i].print_othervalue_name(out,compother);
    break;
  default:{
    print_err("\n Unknown material type is required in ", __FILE__, __LINE__, __func__);
  }
  }
}
