#ifndef TRLINAXISYM_H
#define TRLINAXISYM_H

#include "genfile.h"

/**
   class trlinaxisym defines triangular element with linear approximation functions
   for transport problems for axisymmetric solution
*/

class trlinaxisym
{
 public:
  trlinaxisym (void);
  ~trlinaxisym (void);
  void codnum (long *cn,long ri);
  void ipcoordblock (long eid,long ri,long ci,double **coord);
  double element_area (long eid);
  double approx (vector &areacoord,vector &nodval);
  double approx_nat (double xi,double eta,vector &nodval);
  void intpointval (long eid);
  void initintpointval (long eid);
  void intpointgrad (long eid);
  void intpointother (long eid);

  void bf_matrix (matrix &n,vector &areacoord);
  void give_approx_fun (double &f,vector &areacoord,long i);
  void grad_matrix (matrix &gm,vector &b,vector &c);
  void conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void capacity_matrix (long eid,long ri,long ci,matrix &cm);
  void quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci);
  void internal_fluxes (long lcid,long eid,vector &ifl);

  void res_conductivity_matrix (long eid,long lcid,matrix &km);
  void res_capacity_matrix (long eid,matrix &cm);
  void res_convection_vector (vector &f,long lcid,long eid,long leid);
  void res_transmission_vector (vector &f,long lcid,long eid,long leid);
  void res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid);
  void res_internal_fluxes (long eid,vector &elemif);
  double total_integral (long eid,vector &nodval);
  void res_boundary_flux (vector &f,long lcid,long eid,long leid);

  void volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs);
  void volume_rhs_vector2 (long lcid,long eid,long ri,long ci,vector &vrhs);
  void res_volume_rhs_vector (vector &f,long eid,long lcid);
  void res_volume_rhs_vector2 (vector &f,long eid,long lcid);

  void intpointflux (long eid);

  void nod_others (long lcid,long eid,long ri,long ci);
  void nod_others_comp (long lcid,long eid,long ri,long ci);
    
  void convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void transmission_vector (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
		      vector &coef,matrix &km);
  void transf_flux (long edg,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_coeff (long edg,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_val (long edg,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void edgenodeval (long edg,vector &nodval,vector &list);
  /// computes global coordinates of the given integration point of element
  long ipcoord (long eid, long ipp, long ri, long ci, vector &coord);
  /// function returns natural coordinates of the given integration point
  long ipncoord (long eid, long ipp, vector &ncoord);

  ///  number of transported matter
  long ntm;
  ///  total number of DOFs on the element
  long ndofe;
  ///  number of DOFs on the element
  long **dofe;
  ///  number of nodes on one element
  long nne;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of components of stress and strain tensors
  long ncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  order of integration of conductivity matrix
  long **intordkm;
  ///  order of integration on boundaries
  long intordb;
  ///  order of integration of capacity matrix
  long **intordcm;
  ///  number of integration points for conductivity matrix
  long **nip;
  long **ordering;
};

#endif
