#ifndef QUADHEXT_H
#define QUADHEXT_H

#include "genfile.h"

class quadhext
{
 public:
  quadhext (void);
  ~quadhext (void);
  void codnum (long *cn,long ri);
  double approx (double xi,double eta,double zeta, vector &nodval);
  void intpointval (long eid);
  void intpointgrad (long eid);
  void intpointother (long eid);
  
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void grad_matrix (matrix &gm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac);
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
  void nod_others (long lcid,long eid,long ri,long ci);
  
  
  void convection_vector (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void transmission_vector (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void boundary_flux (vector &v,long lcid,long eid,long leid,long ri,long ci);
  void surface_integral (long surf,vector &x,vector &y,vector &z,long intord,vector &gp,vector &w,
			 vector &coef,matrix &km);
  void transf_flux (long surf,vector &coeff,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_coeff (long surf,vector &coeff,vector &list,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void transf_val (long surf,vector &nodval,vector &list,vector &trc,vector &trr,long eid,long ri,long ci,long ipp,bocontypet *bc);
  void surfnodeval (long surf,vector &nodval,vector &list);
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
  ///  orders of integration of conductivity matrices
  long **intordkm;
  ///  orders of integration of capacity matrices
  long **intordcm;
};

#endif
