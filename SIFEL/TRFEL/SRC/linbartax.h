#ifndef LINBARTAX_H
#define LINBARTAX_H

#include "genfile.h"

/**
   class linbartax defines onedimensional element with linear approximation functions for axisymmetric problems
*/
class linbartax
{
 public:
  linbartax (void);
  ~linbartax (void);

  void codnum (long *cn,long ri);
  double approx (double xi,vector &nodval);
  void intpointval (long eid);
  void intpointgrad (long eid);
  void intpointother (long eid);
  
  void bf_matrix (matrix &n,double xi);
  void grad_matrix (matrix &gm,vector &x,double xi,double &jac);
  void conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void capacity_matrix (long eid,long ri,long ci,matrix &cm);
  void quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci);
  void transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void transmission_vector (vector &tmv,long lcid,long eid,long leid,long ri,long ci);
  void convection_vector (vector &f,long lcid,long eid,long leid,long ri,long ci);
  void internal_fluxes (long lcid,long eid,vector &ifl);

  void res_conductivity_matrix (long eid,long lcid,matrix &km);
  void res_capacity_matrix (long eid,matrix &cm);
  void res_convection_vector (vector &f,long lcid,long eid,long leid);
  void res_transmission_vector (vector &f,long lcid,long eid,long leid);
  void res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid);
  void res_internal_fluxes (long eid,vector &elemif);
  double total_integral(long eid,vector &nodval);
  void boundary_flux (vector &tmv,long lcid,long eid,long leid,long ri,long ci);
  void res_boundary_flux (vector &f,long lcid,long eid,long leid);  
  void nod_others (long lcid,long eid,long ri,long ci);
  /// function computes global coordinates of the given integration point
  long ipcoord (long eid, long ipp, long ri, long ci, vector &coord);
  /// function returns natural coordinates of the given integration point
  long ipncoord (long eid, long ipp, vector &ncoord);
  void volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs);
  void res_volume_rhs_vector (vector &f,long eid,long lcid);
  void volume_rhs_vector2 (long lcid,long eid,long ri,long ci,vector &vrhs);
  void res_volume_rhs_vector2 (vector &f,long eid,long lcid);

 
  ///  number of transported matters
  long ntm;
  ///  total number of DOFs on the element
  long ndofe;
  ///  numbers of DOFs for particular problems
  long **dofe;
  ///  number of nodes on one element
  long nne;
  ///  number of end nodes
  long nen;

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
