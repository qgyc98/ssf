#include "elemswitch.h"
#include "elemswitcht.h"
#include "elemswitchc.h"
#include "globalc.h"
#include "mechtop.h"


/**
   function computes coupling matrix of particular elements
   
   @param eid - element id
   @param lcid - load case id
   @param vm - coupling matrix of one element
   
   JK, 28.7.2001
*/
void upper_cond_coupl_mat (long eid,long /*lcid*/,matrix &vm)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case coupbar:{
    Cbar->res_upper_cond_coup_matrix (eid,vm);
    break;
  }
  case coupquad:{
    Cquad->res_upper_cond_coup_matrix (eid,vm);
    break;
  }
  case coupaxiquad:{
    Caxiq->res_upper_cond_coup_matrix (eid,vm);
    break;
  }
  case couphex:{
    Chex->res_upper_cond_coup_matrix (eid,vm);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function computes coupling matrix of particular elements
   
   @param eid - element id
   @param lcid - load case id
   @param vm - coupling matrix of one element

   JK, 28.7.2001
*/
void lower_cond_coupl_mat (long eid,long /*lcid*/,matrix &vm)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case coupbar:{
    Cbar->res_lower_cond_coup_matrix (eid,vm);
    break;
  }
  case coupquad:{
    Cquad->res_lower_cond_coup_matrix (eid,vm);
    break;
  }
  case coupaxiquad:{
    Caxiq->res_lower_cond_coup_matrix (eid,vm);
    break;
  }
  case couphex:{
    Chex->res_lower_cond_coup_matrix (eid,vm);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function computes coupling matrix of particular elements
   
   @param eid - element id
   @param lcid - load case id
   @param vm - coupling matrix of one element
   
   JK, 28.7.2001
*/
void upper_cap_coupl_mat (long eid,long /*lcid*/,matrix &vm)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case coupbar:{
    Cbar->res_upper_cap_coup_matrix (eid,vm);
    break;
  }
  case coupquad:{
    Cquad->res_upper_cap_coup_matrix (eid,vm);
    break;
  }
  case coupaxiquad:{
    Caxiq->res_upper_cap_coup_matrix (eid,vm);
    break;
  }
  case couphex:{
    Chex->res_upper_cap_coup_matrix (eid,vm);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function computes coupling matrix of particular elements
   
   @param eid - element id
   @param lcid - load case id
   @param vm - coupling matrix of one element
   
   JK, 28.7.2001
*/
void lower_cap_coupl_mat (long eid,long /*lcid*/,matrix &vm)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case coupbar:{
    Cbar->res_lower_cap_coup_matrix (eid,vm);
    break;
  }
  case coupquad:{
    Cquad->res_lower_cap_coup_matrix (eid,vm);
    break;
  }
  case coupaxiquad:{
    Caxiq->res_lower_cap_coup_matrix (eid,vm);
    break;
  }
  case couphex:{
    Chex->res_lower_cap_coup_matrix (eid,vm);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function assembles zero order %matrix of a single element
   
   @param eid - element id
   @param km - zero order %matrix of an element
   
   JK, 11.4.2019
*/
void zero_order_matrix (long eid,matrix &km)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case axisymfc:{
    Caxifc->zero_order_matrix (eid,km);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function assembles first order %matrix of a single element
   
   @param eid - element id
   @param cm - first order %matrix of an element
   
   JK, 11.4.2019
*/
void first_order_matrix (long eid,matrix &cm)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case axisymfc:{
    Caxifc->first_order_matrix (eid,cm);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function computes values on element int. points on coupling mesh
   
   JK, 28.7.2001
*/
void intpointvaluesc (long eid)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);

  switch (te){
  case coupbar:{
    Cbar->intpointval (eid);
    break;
  }
  case coupquad:{
    Cquad->intpointval (eid);
    break;
  }
  case coupaxiquad:{
    Caxiq->intpointval (eid);
    break;
  }
  case couphex:{
    Chex->intpointval (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function computes strains on element int. points on coupling mesh
   
   JK, 28.7.2001
*/
void intpointstrainsc (long eid)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);

  switch (te){
  case coupbar:{
    Cbar->res_mainip_strains (eid);
    break;
  }
  case coupquad:{
    Cquad->res_mainip_strains (eid);
    break;
  }
  case coupaxiquad:{
    Caxiq->res_mainip_strains (eid);
    break;
  }
  case couphex:{
    Chex->res_mainip_strains (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function computes gradients on element int. points on coupling mesh
   
   JK, 28.7.2001
*/
void intpointgradientsc (long eid)
{
  elemtypec te;
  
  te=Ct->give_elem_type (eid);

  switch (te){
  case coupbar:{
    Cbar->intpointgrad (eid);
    break;
  }
  case coupquad:{
    Cquad->intpointgrad (eid);
    break;
  }
  case coupaxiquad:{
    Caxiq->intpointgrad (eid);
    break;
  }
  case couphex:{
    Chex->intpointgrad (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function assembles vector of coupling right hand side
   
   @param lcid - load case id
   @param rhs - right hand side
   
   21.3.2004, JK+TKr
*/
void assemble_coup (long /*lcid*/,double *rhs,long n)
{
  long i,mndofe,*cn,ne;
  vector lv;
  elemtypec te;
  
  ne = Ct->ne;
  
  nullv (rhs,n);
  
  //  contributions from conduction
  for (i=0;i<ne;i++){
    te=Ct->give_elem_type (i);
    mndofe=Gtm->give_ndofe (i);
    allocv (mndofe,lv);
    fillv (0.0,lv);

    switch (te){
    case coupbar:{
      Cbar->res_upper_cond_coup_vector (lv,i);
      break;
    }
    case coupquad:{
      Cquad->res_upper_cond_coup_vector (lv,i);
      break;
    }
    case coupaxiquad:{
      Caxiq->res_upper_cond_coup_vector (lv,i);
      break;
    }
   case couphex:{
      Chex->res_upper_cond_coup_vector (lv,i);
      break;
    }
    default:{
      print_err("unknown element type is required",__FILE__,__LINE__,__func__);
    }
    }
    
    cn = new long [mndofe];
    Mt->give_code_numbers (i,cn);
    locglob (rhs,lv.a,cn,mndofe);
    destrv (lv);
    delete [] cn;
  }  
}



/**
   function computes volume integral
   
   @param lv - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   
   09/03/2011, TKr
*/
void volume_rhs_vectorc (vector &/*lv*/,long /*lcid*/,long eid)
{
  elemtypec te;
  //  element type
  te=Ct->give_elem_type (eid);
  
  switch (te){
  case coupbar:{
    //Cbar->res_volume_rhs_vector (lv,eid,lcid); //has to be corrected in barelc.cpp
    break;
  }
  case coupquad:{
    //Cquad->res_volume_rhs_vector (lv,eid,lcid); //not finished yet in quadrilatc.cpp
    break;
  }
  case coupaxiquad:{
    //Caxiq->res_volume_rhs_vector (lv,eid,lcid); //not finished yet in axiquadc.cpp
    break;
  }
  case couphex:{
    break;
  }  
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  } 
}



/**
  Function computes global coordinates of integration point (ipp) on
  element (eid).
   
  Function is used in adjacip and in mesh transer values procedures.
   
  @param eid - element id
  @param ipp - integration point pointer
  @param ri,ci - row and column indices
  @param ipcoord - %vector containing coordinates of integration point (output)
   
  @return The function returns coordinates in the parameter ipcoord.

  Created by JK, 11.4.2019
*/
/*
void ipcoord (long eid,long ipp,long ri,long ci,vector &ipcoord)
{
  switch (Ct->give_elem_type(eid)){
  case axisymmlq:
    Asymlq->ipcoord(eid, ipp, ri, ci, ipcoord);
    break;
  case axisymmqq: 
    Asymqq->ipcoord(eid, ipp, ri, ci, ipcoord);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}
*/


/**
   function assembles %vector of coupling nodal internal fluxes and forces

   @param lcid - load case id
   @param intflux - array containing internal nodal fluxes and forces
   
   22/09/2021, TKr
*/
void internal_coup_fluxes (long /*lcid*/,double *intflux)
{
  long i,ne,ndofe,*cn;
  elemtypec te;
  vector lif;
  
  ne = Ct->ne;

  for (i=0;i<ne;i++){
    te = Ct->give_elem_type (i);
    ndofe = Gtt->give_ndofe (i);
    reallocv (ndofe,lif);
    fillv (0.0,lif);

    cn = new long [ndofe];
    Gtt->give_code_numbers (i,cn);
    
    switch (te){
    case coupbar:{
      //Cbar->res_lower_internal_fluxes (i,lif); //not finished yet
      break;
    }
    case coupquad:{
      //Cquad->res_lower_internal_fluxes (i,lif); //not finished yet
      break;
    }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      abort ();
    }
    }
    
    locglob (intflux,lif.a,cn,ndofe);
    delete [] cn;
    
  }
  //fluxes on boundary
  //bound_fluxes (0,intflux,n);//not completed yet
}



/**
   function assembles %vector of coupling nodal internal forces

   @param lcid - load case id
   @param intflux - array containing internal nodal fluxes
   
   22/09/2021, TKr
*/
void internal_coup_forces (long /*lcid*/,double *intflux)
{
  long i,ne,ndofe,*cn;
  elemtypec te;
  vector lif;
  
  ne = Ct->ne;

  for (i=0;i<ne;i++){
    te = Ct->give_elem_type (i);
    ndofe = Gtm->give_ndofe (i);//mefel ndofe
    reallocv (ndofe,lif);
    fillv (0.0,lif);

    cn = new long [ndofe];
    Gtm->give_code_numbers (i,cn);//mefel cn
    
    switch (te){
    case coupbar:{
      Cbar->res_upper_internal_forces (i,lif);
      break;
    }
    case coupquad:{
      Cquad->res_upper_internal_forces (i,lif);
      break;
    }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
      abort ();
    }
    }
    
    locglob (intflux,lif.a,cn,ndofe);
    delete [] cn;
    
  }

  //forces on boundary
  //bound_forces (0,intflux,n);//not completed yet
}

