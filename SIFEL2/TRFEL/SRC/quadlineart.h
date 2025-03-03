#ifndef QUADLINEART_H
#define QUADLINEART_H

#include "genfile.h"
#include "aliast.h"

/**
   class quadlineart defines plane quadrilateral element with
   bi-linear approximation functions for transport problems
   
*/
class quadlineart
{
 public:
  quadlineart (void);
  ~quadlineart (void);
  void codnum (long *cn,long ri);
  double element_area (long eid);
  double approx (double xi,double eta,vector &nodval);
  void intpointval (long eid);
  void intpointval (long eid,vector &nodval,vector &ipval);
  void intpointval_puc (long eid);
  void initintpointval (long eid);
  void intpointgrad (long eid);
  void intpointother (long eid);

  void btf_matrix (matrix &n,double xi,double eta,matrix &d,long eid);
  void bf_matrix (matrix &n,double xi,double eta);
  void tf_matrix (matrix &w,double xi,double eta,matrix &d,long eid);
  void grad_btf_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac,long eid,long ipp,long ri,long ci);
  void grad_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac);
  void grad_tf_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,matrix &d,long eid);
  void conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void l_matrix (long lcid,long eid,long ri,long ci,matrix &lm);
  void l_t_matrix (long lcid,long eid,long ri,long ci,matrix &lm);
  void capacity_matrix (long eid,long ri,long ci,matrix &cm);
  void quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci);
  void internal_fluxes (long lcid,long eid,vector &ifl);
  void advection_matrix (long lcid,long eid,long ri,long ci,matrix &hm);
  
  
  void res_conductivity_matrix (long eid,long lcid,matrix &km);
  void res_l_matrix (long eid,long lcid,matrix &lm);
  void res_l_t_matrix (long eid,long lcid,matrix &lm);
  void averd_matrix (long eid,matrix &lm);
  void averc_matrix (long eid,matrix &lm);
  double elem_area (long eid);
  void res_capacity_matrix (long eid,matrix &cm);
  void res_reaction_matrix (long eid, matrix &rm);
  void res_convection_vector (vector &f,long lcid,long eid,long leid);
  void res_transmission_vector (vector &f,long lcid,long eid,long leid);
  void res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid);
  void res_internal_fluxes (long eid,vector &elemif);
  double total_integral(long eid,vector &nodval);
  void res_boundary_flux (vector &f,long lcid,long eid,long leid);
  void volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs);
  void volume_rhs_vector2 (long lcid,long eid,long ri,long ci,vector &vrhs);
  void res_volume_rhs_vector (vector &f,long eid,long lcid);
  void res_volume_rhs_vector2 (vector &f,long eid,long lcid);
  void res_advection_matrix (long eid,long lcid,matrix &hm);
  
  void convection_vector (vector &v,long lcid,long eid,long leid,long cid);
  void transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void transmission_vector (vector &v,long lcid,long eid,long leid,long cid);
  void boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
		      vector &t,vector &coef,matrix &km);
  void transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc, int flag);
  void transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void edgenodeval (long edg,vector &nodval,vector &list);

  void intpointflux (long eid);
  void nod_grads_ip (long eid);
  void nod_fluxes_ip (long eid);
  void nod_eqother_ip (long eid);
  void nod_other_ip (long eid);
  void nod_others_comp (long lcid,long eid,long ri,long ci);

  void higher_to_lower_level (long eid,long *counter,double *buff);
  
  double compute_error (long eid, vector *rderfull, int mattid, double &e2, double &u2, double &sizel);
  void transq_nodval (long eid,vector &nodval,nonmechquant nmq);
  void transq_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt);
  void transq_init_nodval (long eid,vector &nodval,nonmechquant nmq);
  void transq_init_nodval_comp (long eid,vector &nodval, long ncne, long nq, nonmechquant *qt);
  /// computes global coordinates of the given integration point of element
  long ipcoord (long eid, long ipp, long ri, long ci, vector &coord);
  /// computes global coordinates of the given integration point block of element
  void ipcoordblock (long eid,long ri,long ci,double **coord);
  /// function returns natural coordinates of the given integration point
  long ipncoord (long eid, long ipp, vector &ncoord);

  void surface_flux (long lcid,long eid,long beid,double *fluxes);
  
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
  ///  number of approximated functions
  long napfun;
  ///  problem dimension
  long ncomp;
  ///  number of integration points
  long **nip;
  ///  unknown ordering
  long **ordering;
  ///  orders of integration of conductivity matrices
  long **intordkm;
  ///  orders of integration of capacity matrices
  long **intordcm;
};

#endif
