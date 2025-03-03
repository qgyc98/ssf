#include "elemswitcht.h"
#include "globalt.h"
#include "globmatt.h"
#include "ipmap.h"

/**
   function computes conductivity %matrix of one element
   
   @param eid - element id
   @param lcid - load case id
   @param km - conductivity %matrix of one element

*/
void conductmat (long eid,long lcid,matrix &km)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case barlint:{
    Lbt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barlint3d:{
    Lbt3d->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barlintax:{
    Lbat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barquadt:{
    Qbt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case barquadtax:{
    Qbat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case trlint:{
    Ltt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case trlaxisym:{
    Ltat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadlint:{
    Lqt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadquadt:{
    Qqt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadquadtax:{
    Qqat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadlaxisym:{
    Lqat->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case ifacequadel:{
    Ifcquadt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case lineartett:{
    Ltett->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case linearhext:{
    Lht->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case quadratichext:{
    Qht->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  case linearwedget:{
    Lwt->res_conductivity_matrix (eid,lcid,km);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
   function computes capacity %matrix of one element

   @param eid - element id
   @param lcid - load case id
   @param cm - capacity %matrix of one element

*/
void capacmat (long eid,long /*lcid*/,matrix &cm)
{
  elemtypet te;

  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->res_capacity_matrix (eid,cm);
    break;
  }
  case barlint3d:{
    Lbt3d->res_capacity_matrix (eid,cm);
    break;
  }
  case barlintax:{
    Lbat->res_capacity_matrix (eid,cm);
    break;
  }
  case barquadt:{
    Qbt->res_capacity_matrix (eid,cm);
    break;
  }
  case barquadtax:{
    Qbat->res_capacity_matrix (eid,cm);
    break;
  }
  case trlint:{
    Ltt->res_capacity_matrix (eid,cm);
    break;
  }
  case trlaxisym:{
    Ltat->res_capacity_matrix (eid,cm);
    break;
  }
  case quadlint:{
    Lqt->res_capacity_matrix (eid,cm);
    break;
  }
  case quadquadt:{
    Qqt->res_capacity_matrix (eid,cm);
    break;
  }
  case quadquadtax:{
    Qqat->res_capacity_matrix (eid,cm);
    break;
  }
  case quadlaxisym:{
    Lqat->res_capacity_matrix (eid,cm);
    break;
  }
  case ifacequadel:{
    Ifcquadt->res_capacity_matrix (eid,cm);
    break;
  }
  case lineartett:{
    Ltett->res_capacity_matrix (eid,cm);
    break;
  }
  case linearhext:{
    Lht->res_capacity_matrix (eid,cm);
    break;
  }
  case quadratichext:{
    Qht->res_capacity_matrix (eid,cm);
    break;
  }
  case linearwedget:{
    Lwt->res_capacity_matrix (eid,cm);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}

/**
   function computes reaction %matrix of one element

   @param eid - element id
   @param lcid - load case id
   @param rm - reaction %matrix of one element

*/
void reactmat (long eid,long /*lcid*/,matrix &rm)
{
  elemtypet te;

  te=Tt->give_elem_type (eid);

  switch (te){
    /*
  case barlint:{
    Lbt->res_capacity_matrix (eid,cm);
    break;
  }
  case barlint3d:{
    Lbt3d->res_capacity_matrix (eid,cm);
    break;
  }
  case barlintax:{
    Lbat->res_capacity_matrix (eid,cm);
    break;
  }
  case barquadt:{
    Qbt->res_capacity_matrix (eid,cm);
    break;
  }
  case barquadtax:{
    Qbat->res_capacity_matrix (eid,cm);
    break;
  }
  case trlint:{
    Ltt->res_capacity_matrix (eid,cm);
    break;
  }
  case trlaxisym:{
    Ltat->res_capacity_matrix (eid,cm);
    break;
  }
  */
  case quadlint:{
    Lqt->res_reaction_matrix (eid,rm);
    break;
  }
    /*
  case quadquadt:{
    Qqt->res_capacity_matrix (eid,cm);
    break;
  }
  case quadquadtax:{
    Qqat->res_capacity_matrix (eid,cm);
    break;
  }
  case quadlaxisym:{
    Lqat->res_capacity_matrix (eid,cm);
    break;
  }
  case ifacequadel:{
    Ifcquadt->res_capacity_matrix (eid,cm);
    break;
  }
  case lineartett:{
    Ltett->res_capacity_matrix (eid,cm);
    break;
  }
  case linearhext:{
    Lht->res_capacity_matrix (eid,cm);
    break;
  }
  case quadratichext:{
    Qht->res_capacity_matrix (eid,cm);
    break;
  }
  */
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   The function computes gradient %matrix of the given element element at point with
   natural coordinates xi, eta, zeta.

   @param eid - element id [in]
   @param x, y ,z - vectors of nodal coordinates of the given element [in]
   @param xi, eta, zeta  - natural coordinates of the point for which the %matrix will be calculated
   @param gm - resulting gradient %matrix [out]

  @return The function returns resulting stiffness matrix in the argument gm.
 
  Created by Tomas Koudelka, 19.12.2017
*/
void grad_matrix (long eid, vector &x, vector &y, vector &z, double xi, double eta, double zeta, matrix &gm)
{
  elemtypet te;
  double jac, det;
  vector b, c, d;

  te=Tt->give_elem_type (eid);

  switch (te){
    case barlint:{
      Lbt->grad_matrix(gm, x, xi, jac);
      break;
    }
    case barlint3d:{
      Lbt3d->grad_matrix(gm, x, y, z, xi, jac);
      break;
    }
    case barlintax:{
      Lbat->grad_matrix(gm, x, xi, jac);
      break;
     }
    case barquadt:{
      Qbt->grad_matrix(gm, x, xi, jac);
      break;
    }
    case barquadtax:{
      Qbat->grad_matrix(gm, x, xi, jac);
      break;
    }
    case trlint:{
      det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
      reallocv(RSTCKVEC(3, b));
      reallocv(RSTCKVEC(3, c));
      plsb(b.a, y.a, det);
      plsc(c.a, x.a, det);
      Ltt->grad_matrix(gm, b, c);
      break;
    }
    case trlaxisym:{
      det = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
      reallocv(RSTCKVEC(3, b));
      reallocv(RSTCKVEC(3, c));
      plsb(b.a, y.a, det);
      plsc(c.a, x.a, det);
      Ltat->grad_matrix(gm, b, c);
      break;
    }
    case quadlint:{
      Lqt->grad_matrix(gm, x, y, xi, eta, jac);
      break;
    }
    case quadquadt:{
      Qqt->grad_matrix(gm, x, y, xi, eta, jac);
      break;
    }
    case quadquadtax:{
      Qqat->grad_matrix(gm, x, y, xi, eta, jac);
     break;
    }
    case quadlaxisym:{
      Lqat->grad_matrix(gm, x, y, xi, eta, jac);
      break;
    }
    /*
    case ifacequadel:{
      //  virtual width of the element
      h = Tm->ifacemat[Tm->ip[app].idm].h;
      Ifcquadt->grad_matrix(gm, x, y, l, h);
      break;
    }
    */
    case lineartett:{
      reallocv(RSTCKVEC(3, b));
      reallocv(RSTCKVEC(3, c));
      reallocv(RSTCKVEC(3, d));
      det = det3d (x.a, y.a, z.a);
      volb_3d (b.a, y.a, z.a, det);
      volc_3d (c.a, x.a, z.a, det);
      vold_3d (d.a, x.a, y.a, det);
      Ltett->grad_matrix(gm, b, c, d);
      break;
    }
    case linearhext:{
      Lht->grad_matrix(gm, x, y, z, xi, eta, zeta, jac);
      break;
    }
    case quadratichext:{
      Qht->grad_matrix(gm, x, y, z, xi, eta, zeta, jac);
      break;
    }
    case linearwedget:{
      Lwt->grad_matrix(gm, x, y, z, xi, eta, zeta, jac);
      break;
    }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
}



/**
   function computes advection %matrix of one element
   
   @param eid - element id
   @param lcid - load case id
   @param hm - advection %matrix of one element

*/
void advectmat (long eid,long lcid,matrix &hm)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case barlint:{
    Lbt->res_advection_matrix (eid,lcid,hm);
    break;
  }
  case barlint3d:{
    Lbt3d->res_advection_matrix (eid,lcid,hm);
    break;
  }
  case trlint:{
    Ltt->res_advection_matrix (eid,lcid,hm);
    break;
  }
  case quadlint:{
    Lqt->res_advection_matrix (eid,lcid,hm);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function approximates nodal values to integration points
   
   @param eid - element id
   
   JK
*/
void intpointvalues (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointval (eid);
    break;
  }
  case barlint3d:{
    Lbt3d->intpointval (eid);
    break;
  }
  case barlintax:{
    Lbat->intpointval (eid);
    break;
  }
  case barquadt:{
    Qbt->intpointval (eid);
    break;
  }
  case barquadtax:{
    Qbat->intpointval (eid);
    break;
  }
  case trlint:{
    Ltt->intpointval (eid);
    break;
  }
  case trlaxisym:{
    Ltat->intpointval (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointval (eid);
    break;
  }
  case quadquadt:{
    Qqt->intpointval (eid);
    break;
  }
  case quadquadtax:{
    Qqat->intpointval (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointval (eid);
    break;
  }
  case ifacequadel:{
    //Ifcquadt->intpointval (eid);
    break;
  }
  case lineartett:{
    Ltett->intpointval (eid);
    break;
  }
  case linearhext:{
    Lht->intpointval (eid);
    break;
  }
  case quadratichext:{
    Qht->intpointval (eid);
    break;
  }
  case linearwedget:{
    Lwt->intpointval (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function approximates nodal values of actual unknown to integration points of PUC
   
   @param eid - element id
   
   TKr, 08/09/2010 
*/
void intpointvalues_puc (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case trlint:{
    Ltt->intpointval_puc (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointval_puc (eid);
    break;
  }
  case lineartett:{
    Ltett->intpointval_puc (eid);
    break;
  }
  case linearhext:{
    Lht->intpointval_puc (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   function approximates initial nodal values to integration points
   
   @param eid - element id
   
   TKo, 4.7.2018
*/
void initintpointvalues (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->initintpointval (eid);
    break;
  }
  case barlint3d:{
    Lbt3d->initintpointval (eid);
    break;
  }
    //case trlint:{
    //Ltt->initintpointval (eid);
    //break;
    //}
  case trlaxisym:{
    Ltat->initintpointval (eid);
    break;
  }
  case quadlint:{
    Lqt->initintpointval (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->initintpointval (eid);
    break;
  }
  case quadquadtax:{
    Qqat->initintpointval (eid);
    break;
  }
  case lineartett:{
    Ltett->initintpointval (eid);
    break;
  }
  case linearhext:{
    Lht->initintpointval (eid);
    break;
  }
  case linearwedget:{
    Lwt->initintpointval (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function computes gradients on integration points
   
   @param eid - element id
   
   JK, 10.7.2008
*/
void intpointgradients (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointgrad (eid);
    break;
  }
  case barlint3d:{
    Lbt3d->intpointgrad (eid);
    break;
  }
  case barlintax:{
    Lbat->intpointgrad (eid);
    break;
  }
  case barquadt:{
    Qbt->intpointgrad (eid);
    break;
  }
  case barquadtax:{
    Qbat->intpointgrad (eid);
    break;
  }
  case trlint:{
    Ltt->intpointgrad (eid);
    break;
  }
  case trlaxisym:{
    Ltat->intpointgrad (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointgrad (eid);
    break;
  }
  case quadquadt:{
    Qqt->intpointgrad (eid);
    break;
  }
  case quadquadtax:{
    Qqat->intpointgrad (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointgrad (eid);
    break;
  }
    //case ifacequadel:{
    //Ifcquadt->intpointgrad (eid);
    //break;
    //}
  case lineartett:{
    Ltett->intpointgrad (eid);
    break;
  }
  case linearhext:{
    Lht->intpointgrad (eid);
    break;
  }
  case quadratichext:{
    Qht->intpointgrad (eid);
    break;
  }
  case linearwedget:{
    Lwt->intpointgrad (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes fluxes on integration points
   
   @param eid - element id
   
   TKr, 01/02/2010
*/
void intpointfluxes (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointflux (eid);
    break;
  }
  case barlint3d:{
    Lbt3d->intpointflux (eid);
    break;
  }
    /* case barlintax:{
    Lbat->intpointflux (eid);
    break;
    }
    case barquadt:{
    Qbt->intpointflux (eid);
    break;
    }
    case barquadtax:{
    //Qbat->intpointflux (eid);
    break;
    }
    */
  case trlint:{
    Ltt->intpointflux (eid);
    break;
  }
  case trlaxisym:{
    Ltat->intpointflux (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointflux (eid);
    break;
  }
  case quadquadt:{
    Qqt->intpointflux (eid);
    break;
  }
  case quadquadtax:{
    Qqat->intpointflux (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointflux (eid);
    break;
  }
    //     case ifacequadel:{
    //Ifcquadt->intpointflux (eid);
    // break;
    // }
    
  case lineartett:{
    Ltett->intpointflux (eid);
    break;
  }
  case linearhext:{
    Lht->intpointflux (eid);
    break;
  }
    //case quadratichext:{
    //Qht->intpointflux (eid);
    //break;
    //}
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes values of array other at integration points
   
   @param lcid - load case id, id of the medium
   @param eid - element id
   @param fl - %vector containing flux
   
   JK, 12. 8. 2014
*/
void averageflux (long lcid,long eid,vector &fl)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    //Lbt->intpointother (eid);
    break;
  }
  case barlint3d:{
    //Lbt3d->intpointother (eid);
    break;
  }
  case barlintax:{
    //Lbat->intpointother (eid);
    break;
  }
  case barquadt:{
    //Qbt->intpointother (eid);
    break;
  }
  case barquadtax:{
    //Qbat->intpointother (eid);
    break;
  }
  case trlint:{
    //Ltt->intpointother (eid);
    break;
  }
  case trlaxisym:{
    //Ltat->intpointother (eid);
    break;
  }
  case quadlint:{
    //Lqt->intpointother (eid);
    break;
  }
  case quadquadt:{
    //Qqt->intpointother (eid);
    break;
  }
  case quadquadtax:{
    //Qqat->intpointother (eid);
    break;
  }
  case quadlaxisym:{
    //Lqat->intpointother (eid);
    break;
  }
  case ifacequadel:{
    //Ifcquadt->intpointother (eid);
    break;
  }
  case lineartett:{
    //Ltett->intpointother (eid);
    break;
  }
  case linearhext:{
    Lht->average_flux (lcid,eid,fl);
    break;
  }
  case quadratichext:{
    //Qht->intpointother (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes values of array other at integration points
   
   @param eid - element id
   
   JK
*/
void intpointothers (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointother (eid);
    break;
  }
  case barlint3d:{
    Lbt3d->intpointother (eid);
    break;
  }
  case barlintax:{
    Lbat->intpointother (eid);
    break;
  }
  case barquadt:{
    Qbt->intpointother (eid);
    break;
  }
  case barquadtax:{
    Qbat->intpointother (eid);
    break;
  }
  case trlint:{
    Ltt->intpointother (eid);
    break;
  }
  case trlaxisym:{
    Ltat->intpointother (eid);
    break;
  }
  case quadlint:{
    Lqt->intpointother (eid);
    break;
  }
  case quadquadt:{
    Qqt->intpointother (eid);
    break;
  }
  case quadquadtax:{
    Qqat->intpointother (eid);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointother (eid);
    break;
  }
  case ifacequadel:{
    //Ifcquadt->intpointother (eid);
    break;
  }
  case lineartett:{
    Ltett->intpointother (eid);
    break;
  }
  case linearhext:{
    Lht->intpointother (eid);
    break;
  }
  case quadratichext:{
    Qht->intpointother (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function computes values of array eqother at integration points
   
   @param eid - element id
   
   TKr, 01/02/2010
*/
void intpointeqothers (long eid)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    //Lbt->intpointeqother (eid);
    break;
  }
  case barlint3d:{
    //Lbt3d->intpointeqother (eid);
    break;
  }
  case barlintax:{
    //Lbat->intpointeqother (eid);
    break;
  }
  case barquadt:{
    //Qbt->intpointeqother (eid);
    break;
  }
  case barquadtax:{
    //Qbat->intpointeqother (eid);
    break;
  }
  case trlint:{
    //Ltt->intpointeqother (eid);
    break;
  }
  case trlaxisym:{
    //Ltat->intpointeqother (eid);
    break;
  }
  case quadlint:{
    //Lqt->intpointeqother (eid);
    break;
  }
  case quadquadt:{
    //Qqt->intpointeqother (eid);
    break;
  }
  case quadquadtax:{
    //Qqat->intpointeqother (eid);
    break;
  }
  case quadlaxisym:{
    //Lqat->intpointeqother (eid);
    break;
  }
  case ifacequadel:{
    //Ifcquadt->intpointeqother (eid);
    break;
  }
  case lineartett:{
    //Ltett->intpointeqother (eid);
    break;
  }
  case linearhext:{
    //Lht->intpointeqother (eid);
    break;
  }
  case quadratichext:{
    //Qht->intpointeqother (eid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes contributions from source of quantity
   
   @param lcid - load case id
   @param eid - element id
   @param nodval - array of source nodal values
   @param lv - local %vector
   
   JK, 25.6.2005
*/
void source_vector (long lcid,long eid,vector &nodval,vector &lv)
{
  elemtypet te;
  
  //  type of element
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case barlint3d:{
    Lbt3d->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case barlintax:{
    Lbat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case trlint:{
    Ltt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case trlaxisym:{
    Ltat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadlint:{
    Lqt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadlaxisym:{
    Lqat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadquadt:{
    Qqt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadquadtax:{
    Qqat->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case ifacequadel:{
    break;
  }
  case lineartett:{
    Ltett->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case linearhext:{
    Lht->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case quadratichext:{
    Qht->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  case linearwedget:{
    Lwt->res_quantity_source_vector (lv,nodval,lcid,eid);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
}



/**
   The function computes fluxes at integration points.   
   All components at all integration points are evaluated.
   
   @param lcid - load case id
   
   @return The function does not return anything.
   
   Created by TKr, 01/02/2010
*/
void compute_ipfluxes ()
{
  //  number of elements
  long ne=Tt->ne;
  
#ifdef INC_OPENMP
#pragma omp parallel num_threads(4)
 {
  #pragma omp for
#endif
  for (long i=0;i<ne;i++){
    if (Gtt->leso[i]==1){
      //  only elements switched on are taken into account
      intpointfluxes (i);
    }
  }
#ifdef INC_OPENMP
 }
#endif
}

/**
   The function computes average flux
   All components in all integration points are taken into account.
   The function is used in homogenization/tiling.

   @param lcid - load case id
   
   @return The function does not return anything.
   
   JK, 12. 8. 2014
*/
void compute_average_fluxes (long lcid, vector &fl)
{
  long i,ne;
  
  //  number of elements
  ne=Tt->ne;
  
  nullv (fl.a,fl.n);

  for (i=0;i<ne;i++){
    if (Gtt->leso[i]==1){
      //  only elements switched on are taken into account
      averageflux (lcid,i,fl);
    }
  }
}


/**
   The function approximates quantity given by nodal values nv to the point given by natural coordinates
   in vector ncoord on the given element eid.

   @param eid - element id [in]
   @param xi, eta, zeta  - natural coordinates of the point for which the value is computed
   @param nv - vector of nodal values of the given quantity [in]

  @return The function returns approximated quantity value.
 
  Created by Tomas Koudelka, 19.12.2017
*/
double approx(long eid, double xi, double eta, double zeta, vector &nv)
{
  elemtypet te;
  double ret = 0.0;

  te = Tt->give_elem_type(eid);

  switch (te){
    case barlint:{
      ret = Lbt->approx(xi, nv);
      break;
    }
    case barlint3d:{
      ret = Lbt3d->approx(xi, nv);
      break;
    }
    case barlintax:{
      ret = Lbat->approx(xi, nv);
      break;
     }
    case barquadt:{
      ret = Qbt->approx(xi, nv);
      break;
    }
    case barquadtax:{
      ret = Qbat->approx(xi, nv);
      break;
    }
    case trlint:{
      ret = Ltt->approx_nat(xi, eta, nv);
      break;
    }
    case trlaxisym:{
      ret = Ltat->approx_nat(xi, eta, nv);
      break;
    }
    case quadlint:{
      ret = Lqt->approx(xi, eta, nv);
      break;
    }
    case quadquadt:{
      ret = Qqt->approx(xi, eta, nv);
      break;
    }
    case quadquadtax:{
      ret = Qqat->approx(xi, eta, nv);
     break;
    }
    case quadlaxisym:{
      ret = Lqat->approx(xi, eta, nv);
      break;
    }
    case lineartett:{
      ret = Ltett->approx_nat(xi, eta, zeta, nv);
      break;
    }
    case linearhext:{
      ret = Lht->approx(xi, eta, zeta, nv);
      break;
    }
    case quadratichext:{
      ret = Qht->approx(xi, eta, zeta, nv);
      break;
    }
    case linearwedget:{
      ret = Lwt->approx(xi, eta, zeta, nv);
      break;
    }
    default:{
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
    }
  }
  return ret;
}



/**
   The function computes gradients at integration points.   
   All components at all integration points are evaluated.
   
   @param lcid - load case id
   
   @return The function does not return anything.
   
   Created by LS, 31.8.2012
*/
void compute_ipgrads ()
{
  //  number of elements
  long ne = Tt->ne;

#ifdef INC_OPENMP
#pragma omp parallel num_threads(4)
 {
  #pragma omp for
#endif
  // loop over all elements
  for (long i=0; i<ne; i++)
    if (Gtt->leso[i])
      //  only elements switched on are taken into account
      intpointgradients (i);
#ifdef INC_OPENMP
 }
#endif
}


/**
  The function computes other values at integration points.   
  All components at all integration points are evaluated.
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by TKr, 01/02/2010
*/
void compute_ipotherst ()
{
  //  number of elements
  long ne=Tt->ne;

#ifdef INC_OPENMP
#pragma omp parallel num_threads(4)
 {
  #pragma omp for
#endif
  for (long i=0;i<ne;i++){
    if (Gtt->leso[i]==1){
      //  only elements switched on are taken into account
      intpointothers (i);
    }
  }
#ifdef INC_OPENMP
 }
#endif
}


/**
  The function computes eqother values at integration points.   
  All components at all integration points are evaluated.
   
  @param lcid - load case id
   
  @return The function does not return anything.

  Created by TKr, 01/02/2010
*/
void compute_ipeqotherst ()
{
  //  number of elements
  long ne=Tt->ne;
  
#ifdef INC_OPENMP
#pragma omp parallel num_threads(4)
 {
  #pragma omp for
#endif
  for (long i=0;i<ne;i++){
    if (Gtt->leso[i]==1){
      //  only elements switched on are taken into account
      intpointeqothers (i);
    }
  }
#ifdef INC_OPENMP
 }
#endif
}


/**
   Function computes gradients at element nodes.
   Nodal gradients are averaged. In case of bar elements 
   be carefull if averaging is reasonable.
   
   @param lcid - load case id
   
   @return The function does not return anything.
   
   TKr, 01/02/2010
*/
void compute_nodegrads ()
{
  long i;
  elemtypet te;
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].nullgrad();
  }
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      te = Tt->give_elem_type (i);
      
      switch (te){
      case barlint:{
	Lbt->nod_grads_ip (i);
	break;
      }
      case barlint3d:{
	Lbt3d->nod_grads_ip (i);
	break;
      }
      case barlintax:{
	//Lbat->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case barquadt:{
	//Qbt->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case barquadtax:{
	//Qbat->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case quadlint:{
	Lqt->nod_grads_ip (i);
	break;
      }
      case trlint:{
	Ltt->nod_grads_ip (i);
	break;
      }
      case trlaxisym:{
	//Ltat->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case quadquadt:{
	//Qqt->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case quadquadtax:{
	//Qqat->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case quadlaxisym:{
	Lqat->nod_grads_ip (i);
	break;
      }
      case ifacequadel:{
	//Ifcquadt->nod_grads_ip (lcid,i,0,0);
	break;
      }
      case lineartett:{
	Ltett->nod_grads_ip (i);
	break;
      }
      case linearhext:{
	Lht->nod_grads_ip (i);
	break;
      }
      case quadratichext:{
	//Qht->nod_grads_ip (lcid,i,0,0);
	break;
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
      }
      }
    }
  }
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].grad_averageval ();
  }
}


/**
   Function computes fluxes at element nodes.
   Nodal fluxes are averaged. In case of bar elements 
   be carefull if averaging is reasonable.
   
   @param lcid - load case id
   
   @return The function does not return anything.
   
   TKr, 01/02/2010
*/
void compute_nodefluxes ()
{
  long i;
  elemtypet te;
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].nullflux();
  }
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      te = Tt->give_elem_type (i);
      
      switch (te){
      case barlint:{
	Lbt->nod_fluxes_ip (i);
	break;
      }
      case barlint3d:{
	Lbt3d->nod_fluxes_ip (i);
	break;
      }
      case barlintax:{
	//Lbat->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case barquadt:{
	//Qbt->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case barquadtax:{
	//Qbat->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case quadlint:{
	Lqt->nod_fluxes_ip (i);
	break;
      }
      case trlint:{
	Ltt->nod_fluxes_ip (i);
	break;
      }
      case trlaxisym:{
	//Ltat->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case quadquadt:{
	//Qqt->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case quadquadtax:{
	//Qqat->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case quadlaxisym:{
	Lqat->nod_fluxes_ip (i);
	break;
      }
      case ifacequadel:{
	//Ifcquadt->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      case lineartett:{
	Ltett->nod_fluxes_ip (i);
	break;
      }
      case linearhext:{
	Lht->nod_fluxes_ip (i);
	break;
      }
      case quadratichext:{
	//Qht->nod_fluxes_ip (lcid,i,0,0);
	break;
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
      }
      }
    }
  }
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].flux_averageval ();
  }
}


/**
  Function computes other values at element nodes.
  Nodal other values are averaged. In case of bar elements 
  be carefull if averaging is reasonable.

  @param lcid - load case id
  
  @return The function does not return anything.
 
  TKr, 01/02/2010
*/
void compute_nodeotherst ()
{
  long i;
  elemtypet te;
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].nullother();
  }
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      te = Tt->give_elem_type (i);
      
      switch (te){
      case barlint:{
	Lbt->nod_others_comp (0,i,0,0);
	break;
      }
      case barlint3d:{
	Lbt3d->nod_others_comp (0,i,0,0);
	break;
      }
        //case barlintax:{
	//Lbat->nod_others_ip (lcid,i,0,0);
	//break;
        //}
      case barquadt:{
	Qbt->nod_others_comp (0,i,0,0);
	break;
      }
        //case barquadtax:{
	//Qbat->nod_others_ip (lcid,i,0,0);
	//break;
        //}
      case quadlint:{
	Lqt->nod_others_comp (0,i,0,0);
	break;
      }
      case trlint:{
	Ltt->nod_others_comp (0,i,0,0);
	break;
      }
      case trlaxisym:{
	Ltat->nod_others_comp (0,i,0,0);
	break;
      }
        //case quadquadt:{
	//Qqt->nod_others_ip (lcid,i,0,0);
        //	break;
        //}
        //case quadquadtax:{
	//Qqat->nod_others_ip (lcid,i,0,0);
        //	break;
        //}
      case quadlaxisym:{
        Lqat->nod_others_comp(0, i, 0, 0);
	break;
      }
        //case ifacequadel:{
	//Ifcquadt->nod_others_ip (lcid,i,0,0);
        //	break;
        //      }
      case lineartett:{
	Ltett->nod_others_comp (0,i,0,0);
	break;
      }
      case linearhext:{
	Lht->nod_others_comp (0,i,0,0);
	break;
      }
        //      case quadratichext:{
	//Qht->nod_others_ip (lcid,i,0,0);
        //	break;
        //      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
      }
      }
    }
  }
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].other_averageval ();
  }
}



/**
  Function computes other values directly at element nodes.
  Nodal other values are averaged. In case of bar elements 
  be carefull if averaging is reasonable.

  @param lcid - load case id
  
  @return The function does not return anything.
 
  TKr, 01/02/2010
*/
void compute_nodeotherst_comp ()
{
  long i;
  elemtypet te;
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].nullother();
  }
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      te = Tt->give_elem_type (i);
      
      switch (te){
      case barlint:{
	Lbt->nod_others_comp (0,i,0,0);
	break;
      }
      case barlint3d:{
	Lbt3d->nod_others_comp (0,i,0,0);
	break;
      }
      case barlintax:{
	//Lbat->nod_others_ip (lcid,i,0,0);
	break;
      }
      case barquadt:{
	Qbt->nod_others_comp (0,i,0,0);
	break;
      }
      case barquadtax:{
	//Qbat->nod_others_ip (lcid,i,0,0);
	break;
      }
      case quadlint:{
	Lqt->nod_others_comp (0,i,0,0);
	break;
      }
      case trlint:{
	Ltt->nod_others_comp (0,i,0,0);
	break;
      }
      case trlaxisym:{
	//Ltat->nod_others_ip (lcid,i,0,0);
	break;
      }
      case quadquadt:{
	Qqt->nod_others_comp (0,i,0,0);
	break;
      }
      case quadquadtax:{
	//Qqat->nod_others_ip (lcid,i,0,0);
	break;
      }
      case quadlaxisym:{
        Lqat->nod_others_comp(0, i, 0, 0);
	break;
      }
      case ifacequadel:{
	break;
      }
      case lineartett:{
	Ltett->nod_others_comp (0,i,0,0);
	break;
      }
      case linearhext:{
	Lht->nod_others_comp (0,i,0,0);
	break;
      }
      case quadratichext:{
	//Qht->nod_others_ip (lcid,i,0,0);
	break;
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
      }
      }
    }
  }
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].other_averageval ();
  }
}


/**
   Function computes eqother values at element nodes.
   Nodal eqother values are averaged. In case of bar elements 
   be carefull if averaging is reasonable.
   
   @param lcid - load case id
   
   @return The function does not return anything.
   
   TKr, 01/02/2010
*/
void compute_nodeeqotherst ()
{
  long i;
  elemtypet te;
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].nulleqother();
  }
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      te = Tt->give_elem_type (i);
      
      switch (te){
      case barlint:{
	//Lbt->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case barlint3d:{
	//Lbt3d->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case barlintax:{
	//Lbat->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case barquadt:{
	//Qbt->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case barquadtax:{
	//Qbat->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case quadlint:{
	Lqt->nod_eqother_ip (i);
	break;
      }
      case trlint:{
	Ltt->nod_eqother_ip (i);
	break;
      }
      case trlaxisym:{
	//Ltat->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case quadquadt:{
	//Qqt->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case quadquadtax:{
	//Qqat->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case quadlaxisym:{
	//Lqat->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      case ifacequadel:{
	break;
      }
      case lineartett:{
	Ltett->nod_eqother_ip (i);
	break;
      }
      case linearhext:{
	Lht->nod_eqother_ip (i);
	break;
      }
      case quadratichext:{
	//Qht->nod_eqothers_ip (lcid,i,0,0);
	break;
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
      }
      }
    }
  }
  
  for (i=0;i<Tt->nn;i++){
    Tt->nodes[i].eqother_averageval ();
  }
}

/**
   function assembles %vector of nodal internal fluxes

   @param lcid - load case id
   @param intflux - array containing internal nodal fluxes
   
*/
void internal_fluxes (double *intflux,long n)
{
  long i, ne, ndofe, ndofemn;
  elemtypet te;
  ivector cn;
  vector lif, alif;
  
  ne = Tt->ne;
  nullv (intflux,n);

  for (i=0;i<ne;i++){
    if (Gtt->leso[i]==1){

      te = Tt->give_elem_type (i);
      ndofe = Tt->give_ndofe (i);
      reallocv (RSTCKVEC(ndofe,lif));
      nullv (lif);
            
      switch (te){
      case barlint:{
	Lbt->res_internal_fluxes (i,lif);
	break;
      }
      case barlint3d:{
	Lbt3d->res_internal_fluxes (i,lif);
	break;
      }
      case barlintax:{
	Lbat->res_internal_fluxes (i,lif);
	break;
      }
      case barquadt:{
	Qbt->res_internal_fluxes (i,lif);
	break;
      }
      case barquadtax:{
	Qbat->res_internal_fluxes (i,lif);
	break;
      }
      case trlint:{
	Ltt->res_internal_fluxes (i,lif);
	break;
      }
      case trlaxisym:{
	Ltat->res_internal_fluxes (i,lif);
	break;
      }
      case quadlint:{
	Lqt->res_internal_fluxes (i,lif);
	break;
      }
      case quadquadt:{
	Qqt->res_internal_fluxes (i,lif);
	break;
      }
      case quadquadtax:{
	Qqat->res_internal_fluxes (i,lif);
	break;
      }
      case quadlaxisym:{
	Lqat->res_internal_fluxes (i,lif);
	break;
      }
      case ifacequadel:{
	//Ifcquadt->res_internal_fluxes (i,lif);
	break;
      }
      case lineartett:{
	Ltett->res_internal_fluxes (i,lif);
	break;
      }
      case linearhext:{
	Lht->res_internal_fluxes (i,lif);
	break;
      }
      case quadratichext:{
	Qht->res_internal_fluxes (i,lif);
	break;
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
	abort ();
      }
      }
      
      ndofemn = Gtt->give_ndofe (i);
      reallocv(RSTCKIVEC(ndofemn, cn));
      Gtt->give_code_numbers (i,cn.a);

      if (ndofe == ndofemn)
        //  localization of element values to the global vector for elements with regular nodes
        locglob (intflux,lif.a,cn.a,ndofe);
      else
      {
        //  this case means hanging nodes
        //  the element contains hanging nodes
        //  the element value vector has to be transformed
        reallocv (RSTCKVEC(ndofemn, alif));
        mtxv (*Tt->elements[i].tmat, lif, alif);
        locglob (intflux, alif.a, cn.a, ndofemn);
      }
    }
  }
  trfel_bound_flux (0,intflux,n);
}


/**
   function assembles %vector of nodal energy
   integral of B^T q is an energy in a node of a finite element

   @param nodener - node energy
   @param n - the size of the %vector nodener
   
   JK, 8.4.2019
*/
void nodal_energy (double *nodener,long n,double dt)
{
  long i, ne, ndofe, ndofemn;
  elemtypet te;
  ivector cn;
  vector lif, alif;
  
  //  the number of elements
  ne = Tt->ne;
  nullv (nodener,n);

  for (i=0;i<ne;i++){
    if (Gtt->leso[i]==1){
      //  only elements switched on are taken into account
      
      //  element type
      te = Tt->give_elem_type (i);
      //  the number of DOFs on a element
      ndofe = Tt->give_ndofe (i);
      reallocv (RSTCKVEC(ndofe,lif));
      nullv (lif);
            
      switch (te){
      case barlint:{
	Lbt->res_nodal_energy (i,lif,dt);
	break;
      }
      case barlint3d:{
	Lbt3d->res_nodal_energy (i,lif,dt);
	break;
      }
      default:{
	print_err("unknown element type is required",__FILE__,__LINE__,__func__);
	abort ();
      }
      }
      
      ndofemn = Gtt->give_ndofe (i);
      reallocv(RSTCKIVEC(ndofemn, cn));
      Gtt->give_code_numbers (i,cn.a);
      
      if (ndofe == ndofemn){
	//  localization of element values to the global vector for elements with regular nodes
	locglob (nodener,lif.a,cn.a,ndofe);
      }
      else{
	//  this case means hanging nodes
	//  the element contains hanging nodes
	//  the element value vector has to be transformed
	reallocv (RSTCKVEC(ndofemn, alif));
	mtxv (*Tt->elements[i].tmat, lif, alif);
	locglob (nodener, alif.a, cn.a, ndofemn);
      }
    }//  end of the statement if (Gtt->leso[i]==1){
  }//  end of for statement
  
}




/**
   function computes nodal values on element defined by
   prescribed fluxes (Neumann boundary condition)
   
   @param lv - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   @param i - index in the list of loaded elements
   
   JK, 22.11.2008
*/
void elem_neumann_vector (vector &lv,long lcid,long eid,long i)
{
  elemtypet te;
  //  element type
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case barlint3d:{
    Lbt3d->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case barlintax:{
    Lbat->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case barquadt:{
    Qbt->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case barquadtax:{
    Qbat->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case trlint:{
    Ltt->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case trlaxisym:{
    Ltat->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case quadlint:{
    Lqt->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case quadlaxisym:{
    Lqat->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case quadquadt:{
    Qqt->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case quadquadtax:{
    Qqat->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case lineartett:{
    Ltett->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case linearhext:{
    Lht->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case quadratichext:{
    Qht->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  case linearwedget:{
    Lwt->res_convection_vector (lv,lcid,eid,i);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes nodal values on element defined by
   prescribed transmission (Newton boundary condition)
   
   this function computes \int \kappa T_{ext}
   
   @param lv - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   @param i - index in the list of loaded elements
   
   JK, 22.11.2008
*/
void elem_newton_vector (vector &lv,long lcid,long eid,long i)
{
  elemtypet te;
  //  element type
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case barlint3d:{
    Lbt3d->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case barlintax:{
    Lbat->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case barquadt:{
    Qbt->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case barquadtax:{
    Qbat->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case trlint:{
    Ltt->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case trlaxisym:{
    Ltat->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case quadlint:{
    Lqt->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case quadlaxisym:{
    Lqat->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case quadquadt:{
    Qqt->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case quadquadtax:{
    Qqat->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case lineartett:{
    Ltett->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case linearhext:{
    Lht->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case quadratichext:{
    Qht->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  case linearwedget:{
    Lwt->res_transmission_vector (lv,lcid,eid,i);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function computes nodal values of boundary fluxes on element
   caused by transmission

   this function differs from the function elem_newton_vector
   this function computes fluxes on boundaries while the
   function elem_newton_vector computes terms for the right hand
   side %vector
   
   this function computes \int \kappa (T-T_{ext})

   @param lv - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   @param i - index in the list of loaded elements
   
   JK, 22.11.2008
*/
void elem_transmission_flux (vector &lv,long lcid,long eid,long i)
{
  elemtypet te;
  
  //  element type
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case barlint3d:{
    Lbt3d->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case barlintax:{
    Lbat->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case barquadt:{
    Qbt->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case barquadtax:{
    Qbat->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case trlint:{
    Ltt->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case trlaxisym:{
    Ltat->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case quadlint:{
    Lqt->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case quadlaxisym:{
    Lqat->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case quadquadt:{
    Qqt->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case quadquadtax:{
    Qqat->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case lineartett:{
    Ltett->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case linearhext:{
    Lht->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  case quadratichext:{
    Qht->res_boundary_flux (lv,lcid,eid,i);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function computes volume integral
   
   @param lv - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   
   JK, 22.11.2008, corrected by TKr 10/12/2013
*/
void volume_rhs_vector (vector &lv,long lcid,long eid)
{
  elemtypet te;
  //  element type
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case barlint3d:{
    Lbt3d->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case barquadt:{
    Qbt->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case barlintax:{
    Lbat->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
    /*    case barquadtax:{
	  break;
	  }
	  case trlaxisym:{
	  break;
	  }
    */
  case trlint:{
    Ltt->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case trlaxisym:{
    Ltat->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
    
  case quadlint:{
    Lqt->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case quadlaxisym:{
    Lqat->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case quadquadtax:{
    Qqat->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case ifacequadel:{
       break;
  }
  case linearhext:{
    Lht->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  case lineartett:{
    break;
  }
    
    /*     
	   
    case quadratichext:{
    break;
    }
    */
  case linearwedget:{
    Lwt->res_volume_rhs_vector (lv,eid,lcid);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes volume integral of the second type
   
   @param lv - %vector of nodal values
   @param lcid - load case id
   @param eid - element id
   
   TKr 16/05/2018
*/
void volume_rhs_vector2 (vector &lv,long lcid,long eid)
{
  elemtypet te;
  //  element type
  te=Tt->elements[eid].te;
  
  switch (te){
  case barlint:{
    Lbt->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case barlintax:{
    Lbat->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case barlint3d:{
    Lbt3d->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case trlint:{
    Ltt->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case trlaxisym:{
    Ltat->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case quadlint:{
    Lqt->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case quadlaxisym:{
    Lqat->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case quadquadtax:{
    Qqat->res_volume_rhs_vector2 (lv,eid,lcid);
    break;
  }
  case linearhext:{
    break;
  }
  case lineartett:{
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes l_t %matrix of one element
   
   @param eid - element id
   @param lcid - load case id
   @param lm - l %matrix of one element

*/
void lmat (long eid,long lcid,matrix &lm)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case trlint:{
    Ltt->res_l_matrix (eid,lcid,lm);
    break;
  }
  case quadlint:{
    Lqt->res_l_matrix (eid,lcid,lm);
    break;
  }
  case lineartett:{
    Ltett->res_l_matrix (eid,lcid,lm);
    break;
  }
  case linearhext:{
    Lht->res_l_matrix (eid,lcid,lm);
    break;
  }

  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
}


/**
   function computes l_t %matrix of one element
   
   @param eid - element id
   @param lcid - load case id
   @param lm - l_t %matrix of one element

*/
void ltmat (long eid,long lcid,matrix &lm)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case trlint:{
    Ltt->res_l_t_matrix (eid,lcid,lm);
    break;
  }
  case quadlint:{
    Lqt->res_l_t_matrix (eid,lcid,lm);
    break;
  }
  case lineartett:{
    Ltett->res_l_t_matrix (eid,lcid,lm);
    break;
  }
  case linearhext:{
    Lht->res_l_t_matrix (eid,lcid,lm);
    break;
  }

  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
}


/**
   function computes average D %matrix of one element
   
   @param eid - element id
   @param km - d %matrix of one element
   
*/
void averdmat (long eid,double &elemarea,matrix &km)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case trlint:{
    Ltt->averd_matrix (eid,km);
    elemarea = Ltt->elem_area(eid);
    break;
  }
  case quadlint:{
    Lqt->averd_matrix (eid,km);
    elemarea = Lqt->elem_area(eid);
    break;
  }
  case lineartett:{
    Ltett->averd_matrix (eid,km);
    elemarea = Ltett->elem_volume(eid);
    break;
  } 
  case linearhext:{
    Lht->averd_matrix (eid,km);
    elemarea = Lht->elem_volume(eid);
    break;
  }
    
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}


/**
   function computes average C %matrix of one element
   
   @param eid - element id
   @param cm - c %matrix of one element
   
*/
void avercmat (long eid,double &elemarea,matrix &cm)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case trlint:{
    Ltt->averc_matrix (eid,cm);
    elemarea = Ltt->elem_area(eid);
    break;
  }
  case quadlint:{
    Lqt->averc_matrix (eid,cm);
    elemarea = Lqt->elem_area(eid);
    break;
  }
  case lineartett:{
    Ltett->averc_matrix (eid,cm);
    elemarea = Ltett->elem_volume(eid);
    break;
  }
  case linearhext:{
    Lht->averc_matrix (eid,cm);
    elemarea = Lht->elem_volume(eid);
    break;
  }
    
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   function returns area/volume of one element
   
   @param eid - element id
*/
double give_elemarea (long eid)
{
  elemtypet te;
  double elemarea;

  te=Tt->give_elem_type (eid);
  
  switch (te){
  case trlint:{
    elemarea = Ltt->elem_area(eid);
    break;
  }
  case quadlint:{
    elemarea = Lqt->elem_area(eid);
    break;
  }
  case lineartett:{
    elemarea = Ltett->elem_volume(eid);
    break;
  } 
  case linearhext:{
    elemarea = Lht->elem_volume(eid);
    break;
  }
    
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }

  return elemarea;
}



/**
   function selects components from the global level
   
   @param eid - element id
   @param counter - actual position in the array buff
   @param buff - array containing components
   
   JK, 25.3.2011
*/
void higher_to_lower_level_elem (long eid,long *counter,double *buff)
{
  elemtypet te;
  
  //  type of element
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case trlint:{
    Ltt->higher_to_lower_level (eid,counter,buff);
    break;
  }
  case quadlint:{
    Lqt->higher_to_lower_level (eid,counter,buff);
    break;
  }
  case lineartett:{
    Ltett->higher_to_lower_level (eid,counter,buff);
    break;
  }
  case linearhext:{
    Lht->higher_to_lower_level (eid,counter,buff);
    break;
  }

  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
  
}



/**
  The function interpolates nodal values into integration points.

  @param gv - array containing values at all nodes of the mesh
  @param ntq - type of non-transport quantity
  @param scale - scale coefficient
   
  @return The function does not return anything.
 
  12/06/2012 TKr according to JK
*/
void intpointvalt (double *gv, nontransquant ntq, double scale)
{
  long i,j,nne,nip,ipid;
  ivector nodes;
  vector nodval,ipval;
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      nne = Tt->give_nne (i);
      nip = Tt->give_tnip (i);
      
      reallocv (nne,nodes);
      reallocv (nne,nodval);
      reallocv (nip,ipval);
      
      //  nodes on element
      Tt->give_elemnodes (i,nodes);
      
      //  nodal values on element
      for (j=0;j<nne;j++){
	nodval[j]=gv[nodes[j]]*scale;
      }
      
      //  interpolation of nodal values to integration points
      elem_intpointvalt (i,nodval,ipval);
      
      //  number of the first integration point
      ipid=Tt->elements[i].ipp[0][0];
      
      for (j=0;j<nip;j++){
	Tm->storenontransq(ntq, ipid, ipval[j]);
	ipid++;
      }
    }
  }
}

/**
   function approximates nodal values to integration points

   @param eid - element id
   @param nodval - %vector of nodal values
   @param ipval - %vector of values at integration points (output)
   
   12/06/2012 TKr according to JK
*/
void elem_intpointvalt (long eid,vector &nodval,vector &ipval)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->intpointval (eid,nodval,ipval);
    break;
  }
  case barlint3d:{
    Lbt3d->intpointval (eid,nodval,ipval);
    break;
  }
  case quadlint:{
    Lqt->intpointval (eid,nodval,ipval);
    break;
  }
  case quadlaxisym:{
    Lqat->intpointval (eid,nodval,ipval);
    break;
  }
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function approximates nontransport quantity from integration points to nodes
   
   @param gv - vector of nodal nontransport quantity
   @param nmq - type of transport (non-mechanical) quantity
 
  13/09/2012 TKr
  Modified 7.10.2013 TKo
*/
void give_transq_nodval (double *gv, nonmechquant nmq)
{
  long i,j,nne;
  ivector nodes;
  vector nodval;
  
  for (i=0;i<Tt->ne;i++){
    if (Gtt->leso[i]==1){
      
      nne = Tt->give_nne (i);
      
      reallocv (RSTCKIVEC(nne,nodes));
      reallocv (RSTCKVEC(nne,nodval));
      
      //  nodes on element
      Tt->give_elemnodes (i,nodes);
      
      //  nodal values on element
      elem_transq_nodval (i, nodval, nmq);
            
      for (j=0;j<nne;j++){
	gv[nodes[j]] = nodval[j];
      }
    }
  }
}



/**
  The function returns transport quantities from element nodes. The quantities 
  are copied from the closest integration points to element nodes.

  @param eid - element id
  @param nodval - %vector of nodal values
  @param nmq - type of transport (non-mechanical) quantity
   
  @return The function does not return anything.
 
  13/08/2012 TKr
  Modified 7.10.2013 TKo
*/
void elem_transq_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  elemtypet te;

  te = Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->transq_nodval(eid, nodval, nmq);
    break;
  }
  case barlint3d:{
    Lbt3d->transq_nodval(eid, nodval, nmq);
    break;
  }
  case quadlint:{
    Lqt->transq_nodval(eid, nodval, nmq);
    break;
  }
  case quadlaxisym:{
    Lqat->transq_nodval(eid, nodval, nmq);
    break;
  }
  case ifacequadel:
    Ifcquadt->transq_nodval(eid, nodval, nmq);
    break;
  case lineartett:{
    Ltett->transq_nodval(eid, nodval, nmq);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  The function returns initial values of transport quantities from element nodes. The quantities 
  are copied from the closest integration points to element nodes.

  @param eid - element id
  @param nodval - %vector of nodal values
  @param nmq - type of transport (non-mechanical) quantity
   
  @return The function does not return anything.
 
  13/08/2012 TKr
  Modified 7.10.2013 TKo
*/
void elem_transq_init_nodval (long eid,vector &nodval,nonmechquant nmq)
{
  elemtypet te;

  te = Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->transq_init_nodval(eid, nodval, nmq);
    break;
  }
  case barlint3d:{
    Lbt3d->transq_init_nodval(eid, nodval, nmq);
    break;
  }
  case quadlint:{
    Lqt->transq_init_nodval(eid, nodval, nmq);
    break;
  }
  case quadlaxisym:{
    Lqat->transq_init_nodval(eid, nodval, nmq);
    break;
  }
  case lineartett:{
    Ltett->transq_init_nodval(eid, nodval, nmq);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function computes transport (non-mechanical) quantities in nodes of element.

  @param eid - element id
  @param nodval - %vector of nodal values of all required quantities, i.e., 
                  nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                  is the number of calculated nodes on eid-th element.
  @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
  @param nq - number of required transport quantities
  @param qt - array of types of required transport quantities
   
  @return The function does not return anything.
 
  Created by Tomas Koudelka, 9.12.2013
*/
void elem_transq_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  elemtypet te;

  te = Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->transq_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case barlint3d:{
    Lbt3d->transq_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case quadlint:{
    Lqt->transq_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case quadlaxisym:{
    Lqat->transq_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case lineartett:{
    Ltett->transq_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
  Function computes initial transport (non-mechanical) quantities in nodes of element.

  @param eid - element id
  @param nodval - %vector of initial nodal values of all required quantities, i.e., 
                  nodal value of i-th quantity in j-th node is given by nodval[i*ncnv+j] where ncnv 
                  is the number of calculated nodes on eid-th element.
  @param ncne - number of computed nodes on element (only first ncne of nodes is calculated)
  @param nq - number of required transport quantities
  @param qt - array of types of required transport quantities
   
  @return The function does not return anything.
 
  Created by Tomas Koudelka, 5.6.2018
*/
void elem_transq_init_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt)
{
  elemtypet te;

  te = Tt->give_elem_type (eid);

  switch (te){
  case barlint:{
    Lbt->transq_init_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case barlint3d:{
    Lbt3d->transq_init_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case quadlint:{
    Lqt->transq_init_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case quadlaxisym:{
    Lqat->transq_init_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  case lineartett:{
    Ltett->transq_init_nodval_comp (eid, nodval, ncne, nq, qt);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }
}



/**
   function returns integral over a single finite element
   variable is stored in nodes
   
   @param eid - number of element
   @param nodval - array of nodal values
   
   14.3.2013, JK
*/
double elem_total_integral (long eid,vector &nodval)
{
  double value;
  elemtypet te;
  
  //  type of element
  te = Tt->give_elem_type (eid);
  
  switch (te){
  case barlint:{      value = Lbt->total_integral (eid,nodval);   break;}
  case barlint3d:{    value = Lbt3d->total_integral (eid,nodval);   break;}
  case barlintax:{    value = Lbat->total_integral (eid,nodval);  break;}
  case barquadt:{     value = Qbt->total_integral (eid,nodval);   break;}
  case barquadtax:{   value = Qbat->total_integral (eid,nodval);  break;}
  case quadlint:{     value = Lqt->total_integral (eid,nodval);   break;}
  case quadquadt:{    value = Qqt->total_integral (eid,nodval);   break;}
  case quadquadtax:{  value = Qqat->total_integral (eid,nodval);  break;}
  case lineartett:{   value = Ltett->total_integral (eid,nodval); break;}
  case linearhext:{   value = Lht->total_integral (eid,nodval);   break;}
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }
  return value;
}

/**
   function returns integral over a single finite element
   variable is stored in integration points
   
   @param eid - number of element
   @param varid - id of variable required
   
   3. 10. 2013, JK
*/
double elem_total_integral_ip (long eid,long varid)
{
  double value;
  elemtypet te;
  
  //  type of element
  te = Tt->give_elem_type (eid);
  
  switch (te){
  case barlint:{      value = Lbt->total_integral_ip (eid,varid);   break;}
  case barlint3d:{    value = Lbt3d->total_integral_ip (eid,varid);   break;}
    //case barlintax:{    value = Lbat->total_integral (eid,nodval);  break;}
    //case barquadt:{     value = Qbt->total_integral (eid,nodval);   break;}
    //case barquadtax:{   value = Qbat->total_integral (eid,nodval);  break;}
    //case quadlint:{     value = Lqt->total_integral (eid,nodval);   break;}
    //case quadquadt:{    value = Qqt->total_integral (eid,nodval);   break;}
    //case quadquadtax:{  value = Qqat->total_integral (eid,nodval);  break;}
    //case lineartett:{   value = Ltett->total_integral (eid,nodval); break;}
    //case linearhext:{   value = Lht->total_integral (eid,nodval);   break;}
  default:{
    print_err("unknown element type is required",__FILE__,__LINE__,__func__);
  }
  }//  end of the switch (te)
  
  return value;
}



/**
   The function computes global coordinates of the given integration point ipp of element eid.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ri, ci - row and column index (input)
   @param coord - %vector of global coordinates of integration point (output)

  Created by Tomas Koudelka, 1.12.2016
*/
void ipcoordt (long eid, long ipp, long ri, long ci, vector &coord)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
    case barlint:
      Lbt->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case barlint3d:
      Lbt3d->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case barlintax:
      Lbat->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case barquadt:
      Qbt->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case barquadtax:
      Qbat->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case trlint:
      Ltt->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case trlaxisym:
      Ltat->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case quadlint:
      Lqt->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case quadquadt:
      Qqt->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case quadquadtax:
      Qqat->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case quadlaxisym:
      Lqat->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case lineartett:
      Ltett->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case linearhext:
      Lht->ipcoord(eid,ipp,ri,ci,coord);
      break;
    case quadratichext:
      Qht->ipcoord(eid,ipp,ri,ci,coord);
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}



/**
   The function returns natural coordinates of the given integration point ipp of element eid.
   
   @param eid - element id (input)
   @param ipp - integration point pointer (input)
   @param ncoord - %vector of global coordinates of integration point (output)

  Created by Tomas Koudelka, 1.12.2016
*/
void ipncoordt (long eid, long ipp, vector &ncoord)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
    case barlint:
      Lbt->ipncoord(eid,ipp,ncoord);
      break;
    case barlint3d:
      Lbt3d->ipncoord(eid,ipp,ncoord);
      break;
    case barlintax:
      Lbat->ipncoord(eid,ipp,ncoord);
      break;
    case barquadt:
      Qbt->ipncoord(eid,ipp,ncoord);
      break;
    case barquadtax:
      Qbat->ipncoord(eid,ipp,ncoord);
      break;
    case trlint:
      Ltt->ipncoord(eid,ipp,ncoord);
      break;
    case trlaxisym:
      Ltat->ipncoord(eid,ipp,ncoord);
      break;
    case quadlint:
      Lqt->ipncoord(eid,ipp,ncoord);
      break;
    case quadquadt:
      Qqt->ipncoord(eid,ipp,ncoord);
      break;
    case quadquadtax:
      Qqat->ipncoord(eid,ipp,ncoord);
      break;
    case quadlaxisym:
      Lqat->ipncoord(eid,ipp,ncoord);
      break;
    case lineartett:
      Ltett->ipncoord(eid,ipp,ncoord);
      break;
    case linearhext:
      Lht->ipncoord(eid,ipp,ncoord);
      break;
    case quadratichext:
      Qht->ipncoord(eid,ipp,ncoord);
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
}


/**
  The function computes global coordinates of center of gravity for the given element eid.
 
  @param eid - element id (input)
  @param coord - %vector containing resulting coordinates (output)
   
  @return The function returns required coordinates of centroid.
 
  Created by  TKo, 1.12.2016
*/
void centroidt(long eid, vector &coord)
{
  elemtypet te = Tt->give_elem_type (eid);
  long nne = Tt->give_nne (eid);
  vector x(ASTCKVEC(nne)), y(ASTCKVEC(nne)), z(ASTCKVEC(nne));
  vector areacoord(3),volcoord(4);

  Tt->give_node_coord3d(x, y, z, eid);

  switch (te)
  {
    case barlint:
      coord(0)=Lbt->approx(0.0, x);
      coord(1)=Lbt->approx(0.0, y);
      coord(2)=Lbt->approx(0.0, z);
      break;
    case barlint3d:
      coord(0)=Lbt3d->approx(0.0, x);
      coord(1)=Lbt3d->approx(0.0, y);
      coord(2)=Lbt3d->approx(0.0, z);
      break;
    case barlintax:
      coord(0)=Lbat->approx(0.0, x);
      coord(1)=Lbat->approx(0.0, y);
      coord(2)=Lbat->approx(0.0, z);
      break;
    case barquadt:
      coord(0)=Qbt->approx(0.0, x);
      coord(1)=Qbt->approx(0.0, y);
      coord(2)=Qbt->approx(0.0, z);
      break;
    case barquadtax:
      coord(0)=Qbat->approx(0.0, x);
      coord(1)=Qbat->approx(0.0, y);
      coord(2)=Qbat->approx(0.0, z);
      break;
    case trlint:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Ltt->approx(areacoord, x);
      coord(1)=Ltt->approx(areacoord, y);
      coord(2)=Ltt->approx(areacoord, z);
      break;
    case trlaxisym:
      areacoord[0]=0.333333333333333;  areacoord[1]=0.333333333333333;  areacoord[2]=1.0-areacoord[0]-areacoord[1];
      coord(0)=Ltat->approx(areacoord, x);
      coord(1)=Ltat->approx(areacoord, y);
      coord(2)=Ltat->approx(areacoord, z);
      break;
    case quadlint:
      coord(0)=Lqt->approx(0.0, 0.0, x);
      coord(1)=Lqt->approx(0.0, 0.0, y);
      coord(2)=Lqt->approx(0.0, 0.0, z);
      break;
    case quadquadt:
      coord(0)=Qqt->approx(0.0, 0.0, x);
      coord(1)=Qqt->approx(0.0, 0.0, y);
      coord(2)=Qqt->approx(0.0, 0.0, z);
      break;
    case quadquadtax:
      coord(0)=Qqat->approx(0.0, 0.0, x);
      coord(1)=Qqat->approx(0.0, 0.0, y);
      coord(2)=Qqat->approx(0.0, 0.0, z);
      break;
    case quadlaxisym:
      coord(0)=Lqat->approx(0.0, 0.0, x);
      coord(1)=Lqat->approx(0.0, 0.0, y);
      coord(2)=Lqat->approx(0.0, 0.0, z);
      break;
    case lineartett:
      volcoord[0]=0.25;  volcoord[1]=0.25;  volcoord[2]=0.25;  volcoord[3]=1.0-volcoord[0]-volcoord[1]-volcoord[2];
      coord(0)=Ltett->approx(volcoord, x);
      coord(1)=Ltett->approx(volcoord, y);
      coord(2)=Ltett->approx(volcoord, z);
      break;
    case linearhext:
      coord(0)=Lht->approx(0.0, 0.0, 0.0, x);
      coord(1)=Lht->approx(0.0, 0.0, 0.0, y);
      coord(2)=Lht->approx(0.0, 0.0, 0.0, z);
      break;
    case quadratichext:
      coord(0)=Qht->approx(0.0, 0.0, 0.0, x);
      coord(1)=Qht->approx(0.0, 0.0, 0.0, y);
      coord(2)=Qht->approx(0.0, 0.0, 0.0, z);
      break;
    default:
      print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  return;
}



/**
  The function computes quantity values in auxiliary integration points
  values in auxiliary integration points are obtained from nodal values
   
  @param n - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]

  @return The function does not return anything but stores computed variable values 
          in particular auxiliary integration points of the transmat aip array.

  Created by Tomas Koudelka, 8.12.2017
*/
void compute_aipvals (long n, ipmap *ipm)
{

  long app, eid, i, j, k;
  long nne, ndofe;
  long **dofe, **ordering;
  double val;
  vector r, t;

  intpointst *tmp_ip;
  long       *tmp_elip;
  double    **tmp_ntq;
  double     *tmp_iv;


  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_ntq = Tm->nontransq;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->nontransq = Tm->aip_nontransq;
  Tm->initval = Tm->aip_initval;

  //  loop over elements
  for (i=0; i<n; i++)
  {
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point => aip values and their gardients have already been computed
      continue;

    eid = Tm->elip[app];
    //  only elements switched on are processed
    if (Gtt->leso[eid]==1){
      nne = Tt->give_nne(eid);
      ndofe = Tt->give_ndofe(eid);
      reallocv(RSTCKVEC(ndofe, r));
      reallocv(RSTCKVEC(nne, t));
      // nodal values of all unknowns on element
      elemvalues (eid, r);
      // array of number of element DOFs for particular transported media
      dofe = Tt->give_dofe(eid);
      // ordering of variables in the vector of unknowns
      ordering = Tt->give_ordering(eid);

      for (j=0; j<Tp->ntm; j++) // loop over all variables of transported media
      {
        //  collect nodal values of j-th single variable in the vector t
        for (k=0; k<dofe[j][j]; k++)
          t[k]=r[ordering[j][k]-1];

        // approximation of j-th variable nodal values to the auxiliary int. point
        val = approx (eid, ipm[i].xi, ipm[i].eta, ipm[i].zeta, t);
        Tm->ip[app].av[j]=val;
      }
    }
  }

  if (Tp->nvs==1){
    //  arrays of actual values are moved to arrays of previous values
    //  implementation due to saltmat4
    actual_previous_change ();
  }

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->nontransq = tmp_ntq;
  Tm->initval = tmp_iv;
}



/**
  The function computes gradients at auxiliary integration points.   
  All components at all integration points are evaluated.
   
  @param n - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]

  @return The function does not return anything but stores computed values 
          in particular auxiliary integration points of the transmat aip array.
   
  Created by Tomas Koudelka, 8.12.2017
*/
void compute_aipgrads(long n, ipmap *ipm)
{
  intpointst *tmp_ip;
  long *tmp_elip;
  double **tmp_ntq;
  double *tmp_iv;

  long app, eid, i, j, k;
  long nne, ndofe, ncomp;
  long **dofe, **ordering;
  vector x, y, z, r, t, grad;
  ivector cn;
  matrix gm;

  // swap regular integration and point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_ntq = Tm->nontransq;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->nontransq = Tm->aip_nontransq;
  Tm->initval = Tm->aip_initval;

  for (i=0; i<n; i++) // loop over all auxiliary integration points in the mapping
  {
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point => gardients have already been computed
      continue;
    // id of element connected with the given auxiliary integration point
    eid = Tm->elip[app];
    //  only elements switched on are taken into account
    if (Gtt->leso[eid] == 1)
    {
      // nodal coordinates of the given element
      nne = Tt->give_nne(eid);
      reallocv(RSTCKVEC(nne, x));
      reallocv(RSTCKVEC(nne, y));
      reallocv(RSTCKVEC(nne, z));
      Tt->give_node_coord3d(x, y, z, eid);
      // code numbers of DOFs on element
      ndofe = Tt->give_ndofe(eid);
      reallocv(RSTCKIVEC(ndofe, cn));
      Tt->give_code_numbers(eid, cn.a);
      // assemble vector of nodal values of all variables on element
      reallocv(RSTCKVEC(ndofe, r));
      reallocv(RSTCKVEC(ndofe, t));
      elemvalues(eid, r);
      // array of number of element DOFs for particular transported media
      dofe = Tt->give_dofe(eid);
      // ordering of variables in the vector of unknowns
      ordering = Tt->give_ordering(eid);
      
      for (j=0; j<Tp->ntm; j++) // loop over all variables of transported media
      {
        //  collect nodal values of j-th single variable in the vector t
        for (k=0; k<dofe[j][j]; k++)
          t[k]=r[ordering[j][k]-1];
        // number of gradient components
        ncomp = Tm->ip[app].ncompgrad;
        // allocate gradient matrix
        reallocm(RSTCKMAT(ncomp, nne, gm));
        // allocate resulting gradient vector
        reallocv(RSTCKVEC(ncomp, grad));
        // get the gradient matrix at the given auxiliary integration point
        grad_matrix(eid, x, y, z, ipm[i].xi, ipm[i].eta, ipm[i].zeta, gm);
        // grad t = gm.t
        mxv(gm, t, grad);
        // store gradients at auxiliary integration point app
        Tm->storegrad(j, app, grad);
      }
    }
  }

  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->nontransq = tmp_ntq;
  Tm->initval = tmp_iv;
}



/**
  The function computes fluxes at auxiliary integration points.   
  All components at all integration points are evaluated.
   
  @return The function does not return anything but stores computed values 
          in particular auxiliary integration points of the transmat aip array.
   
  Created by Tomas Koudelka, 8.12.2017
*/
void compute_aipfluxes()
{
  intpointst *tmp_ip;
  long *tmp_elip;
  double **tmp_ntq;
  double *tmp_iv;
  long app, i, eid;

  // swap regular integration point and auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_ntq = Tm->nontransq;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->nontransq = Tm->aip_nontransq;
  Tm->initval = Tm->aip_initval;
  
  for (app=0; app<Tm->tnaip; app++)
  {
    eid = Tm->elip[app];
    //  only elements switched on are taken into account
    if (Gtt->leso[eid] == 1)
    {
      for (i=0; i<Tp->ntm; i++)  
        Tm->computenlfluxes(i, app);
    }
  }

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->nontransq = tmp_ntq;
  Tm->initval = tmp_iv;
}



/**
  The function computes other values at auxiliary integration points.   
  All components at all integration points are evaluated.
   
  @param n - the number of components in the ipm array [in]
  @param ipm - the array of integration point mapping objects [in]

  @return The function does not return anything but stores computed values 
          in particular auxiliary integration points of the transmat aip array.

  Created by Tomas Koudelka, 8.12.2017
*/
void compute_aipotherst(long n, ipmap *ipm)
{
  intpointst *tmp_ip;
  long *tmp_elip;
  double **tmp_ntq;
  double *tmp_iv;

  long app, eid, i, j, k;
  long nne, ncompo;
  vector r, t;
  ivector nodes, cn;
  double val; 

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with AUXILIARY integration points 
  tmp_ip = Tm->ip;
  tmp_elip = Tm->elip;
  tmp_ntq = Tm->nontransq;
  tmp_iv = Tm->initval;
  Tm->ip = Tm->aip;
  Tm->elip = Tm->elaip;
  Tm->nontransq = Tm->aip_nontransq;
  Tm->initval = Tm->aip_initval;

  for (i=0; i<n; i++) // loop over all auxiliary integration points in the mapping
  {
    app = ipm[i].app;
    if (app < 0) // direct mapping to the regular integration point => gardients have already been computed
      continue;
    // id of element connected with the given auxiliary integration point
    eid = Tm->elip[app];
    //  only elements switched on are taken into account
    if (Gtt->leso[eid] == 1)
    {  
      //  nodes of required element
      nne = Tt->give_nne(eid);
      reallocv(RSTCKIVEC(nne, nodes));
      Tt->give_elemnodes (eid,nodes);
      // allocate vector of nodal values of quantities
      ncompo = Tt->nodes[nodes[0]].ncompother;
      reallocv(RSTCKVEC(ncompo*nne, r));
      //  nodal values of array other
      nodalotherval(nodes, r);
      // allocate vector of nodal values for a single other value
      reallocv(RSTCKVEC(nne, t));
      for (j=0; j<ncompo; j++)
      {
        //  nodal values of a single other value
        for (k=0; k<nne; k++)
          t[k] = r[k*ncompo+j];
        // approximation of the other value
        val = approx(eid, ipm[i].xi, ipm[i].eta, ipm[i].zeta, t);
        // storage of the j-th componet in the other array
        Tm->storeother(app, j, val);
      }
    }
  }

  // swap regular integration point auxiliary integration point arrays
  // from now, all functions will work with REGULAR integration points 
  Tm->ip = tmp_ip;
  Tm->elip = tmp_elip;
  Tm->nontransq = tmp_ntq;
  Tm->initval = tmp_iv;
}



/**
   function computes surface flux
   
   @param eid - element id (in the list of all elements)
   @param lcid - load case id
   @param beid - id of element in the list of selected elements
   @param fluxes - array of fluxes
   
   JK, 5. 3. 2018
*/
void elem_surface_flux (long eid,long lcid,long beid,double *fluxes)
{
  elemtypet te;
  
  te=Tt->give_elem_type (eid);
  
  switch (te){
  case barlint:{
    Lbt->surface_flux (lcid,eid,beid,fluxes);
    break;
  }
  case barlint3d:{
    Lbt3d->surface_flux (lcid,eid,beid,fluxes);
    break;
  }
  case quadlint:{
    Lqt->surface_flux (lcid,eid,beid,fluxes);
    break;
  }
  case quadlaxisym:{
    Lqat->surface_flux (lcid,eid,beid,fluxes);
    break;
  }
  case lineartett:{
    Ltett->surface_flux (lcid,eid,beid,fluxes);
    break;
  }
  case linearhext:{
    Lht->surface_flux (lcid,eid,beid,fluxes);
    break;
  }
  default:{
    print_err("unknown element type is required", __FILE__, __LINE__, __func__);
  }
  }

}

