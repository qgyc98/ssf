#ifndef LINWEDGET_H
#define LINWEDGET_H

#include "genfile.h"

class linwedget
{
 public:
  linwedget (void);
  ~linwedget (void);
  void codnum (long *cn,long ri);
  double element_volume (long eid);
  double approx (double xi,double eta,double zeta, vector &nodval);
  void intpointval (long eid);
  //void intpointval_puc (long eid);
  void initintpointval (long eid);
  void intpointgrad (long eid);
  //void intpointother (long eid);
  //void intpointflux (long eid);
  //void average_flux (long lcid,long eid,vector &avfl);
  
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void grad_matrix (matrix &gm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac);
  void conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  //void l_matrix (long lcid,long eid,long ri,long ci,matrix &lm);
  //void l_t_matrix (long lcid,long eid,long ri,long ci,matrix &lm);
  void capacity_matrix (long eid,long ri,long ci,matrix &cm);
  void quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci);
  //void internal_fluxes (long lcid,long eid,vector &ifl);

  void res_conductivity_matrix (long eid,long lcid,matrix &km);
  void volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs);
  void res_volume_rhs_vector (vector &f,long eid,long lcid);
  //void res_l_matrix (long eid,long lcid,matrix &lm);
  //void res_l_t_matrix (long eid,long lcid,matrix &lm);
  //void averd_matrix (long eid,matrix &lm);
  //void averc_matrix (long eid,matrix &lm);
  //double elem_volume (long eid);

  void res_capacity_matrix (long eid,matrix &cm);
  void res_convection_vector (vector &f,long lcid,long eid,long leid);
  void res_transmission_vector (vector &f,long lcid,long eid,long leid);
  void res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid);
  //void res_internal_fluxes (long eid,vector &elemif);
  //double total_integral (long eid,vector &nodval);
  void res_boundary_flux (vector &f,long lcid,long eid,long leid);
  void surface_flux (long lcid,long eid,long beid,double *fluxes);


  //void nod_grads_ip (long eid);
  //void nod_fluxes_ip (long eid);
  //void nod_others_comp (long lcid,long eid,long ri,long ci);
  //void nod_eqother_ip (long eid);


  void convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void transmission_vector (vector &v,long lcid,long eid,long leid,long cid);
  void boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void surface_integral (long surf,vector &x,vector &y,vector &z,vector &coef,matrix &km);

  void transf_flux (long surf,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_coeff (long surf,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_val (long surf,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void surfnodeval (long surf,vector &nodval,vector &list);

  void higher_to_lower_level (long eid,long *counter,double *buff);
  /// computes global coordinates of the given integration point of element
  long ipcoord (long eid, long ipp, long ri, long ci, vector &coord);
  /// function returns natural coordinates of the given integration point
  long ipncoord (long eid, long ipp, vector &ncoord);

  ///  number of transported matter
  long ntm;
  ///  total number of DOFs on the element
  long ndofe;
  ///  numbers of DOFs for particular problems
  long **dofe;
  ///  number of nodes on one element
  long nne;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of surfaces
  long nsurf;
  ///  number of nodes on one surface
  long nnsurf;
  ///  number of approximated functions
  long napfun;
  ///  problem dimension
  long ncomp;
  ///  number of integration points
  long **nip;
  ///  unknown ordering
  long **ordering;
  ///  orders of integration on triangular surfaces in conductivity matrices
  long **intordkmt;
  ///  orders of integration in the z direction in conductivity matrices
  long **intordkmz;
  ///  orders of integration on triangular surfaces in capacity matrices
  long **intordcmt;
  ///  orders of integration in the z direction in capacity matrices
  long **intordcmz;
};

#endif
