#ifndef LINBART_H
#define LINBART_H

#include "genfile.h"
#include "aliast.h"

/**
   class linbart defines onedimensional element with linear approximation functions
*/
class linbart
{
 public:
  linbart (void);
  ~linbart (void);

  void codnum (long *cn,long ri);
  double approx (double xi,vector &nodval);
  void intpointval (long eid);
  void intpointval (long eid,vector &nodval,vector &ipval);
  void initintpointval (long eid);
  void intpointgrad (long eid);
  void intpointother (long eid);

  void btf_matrix (matrix &n,double xi,double v);
  void bf_matrix (matrix &n,double xi);
  void tf_matrix (matrix &w,double xi,double v);
  void grad_matrix (matrix &gm,vector &x,double xi,double &jac);
  void conductivity_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void capacity_matrix (long eid,long ri,long ci,matrix &cm);
  void quantity_source_vector (vector &sv,vector &nodval,long eid,long ri,long ci);
  void transmission_matrix (long lcid,long eid,long ri,long ci,matrix &km);
  void transmission_vector (vector &tmv,long lcid,long eid,long leid,long cid);
  void convection_vector (vector &f,long lcid,long eid,long leid,long ri,long ci);
  void internal_fluxes (long lcid,long eid,vector &ifl);
  void accumul_energy (long lcid,long eid,vector &acum);
  void advection_matrix (long lcid,long eid,long ri,long ci,matrix &hm);

  void res_conductivity_matrix (long eid,long lcid,matrix &km);
  void volume_rhs_vector (long lcid,long eid,long ri,long ci,vector &vrhs);
  void volume_rhs_vector2 (long lcid,long eid,long ri,long ci,vector &vrhs);
  void res_volume_rhs_vector (vector &f,long eid,long lcid);
  void res_volume_rhs_vector2 (vector &f,long eid,long lcid);
  void res_capacity_matrix (long eid,matrix &cm);
  void res_convection_vector (vector &f,long lcid,long eid,long leid);
  void res_transmission_vector (vector &f,long lcid,long eid,long leid);
  void res_quantity_source_vector (vector &sv,vector &nodval,long lcid,long eid);
  void res_internal_fluxes (long eid,vector &elemif);
  void res_advection_matrix (long eid,long lcid,matrix &hm);
  
  void nodal_energies (long lcid,long eid,vector &ifl);
  void res_nodal_energy (long eid,vector &elemne,double dt);
  
  ///  function computes element quantity integral
  ///  the quantity is stored in nodes
  double total_integral(long eid,vector &nodval);
  ///  function computes element quantity integral
  ///  the quantity is stored in integration points
  double total_integral_ip (long eid,long varid);


  void boundary_flux (vector &tmv,long lcid,long eid,long leid,long ri,long ci);
  void res_boundary_flux (vector &f,long lcid,long eid,long leid);  
 
  void intpointflux (long eid);
  void nod_grads_ip (long eid);
  void nod_fluxes_ip (long eid);
  void nod_others_comp (long lcid,long eid,long ri,long ci);
  void transq_nodval (long eid,vector &nodval,nonmechquant nmq);
  void transq_nodval_comp (long eid, vector &nodval, long ncne, long nq, nonmechquant *qt);
  void transq_init_nodval (long eid,vector &nodval,nonmechquant nmq);
  void transq_init_nodval_comp (long eid, vector &nodval, long ncne, long nq, nonmechquant *qt);
  /// computes global coordinates of the given integration point of element
  long ipcoord (long eid, long ipp, long ri, long ci, vector &coord);
  /// function returns natural coordinates of the given integration point
  long ipncoord (long eid, long ipp, vector &ncoord);
  
  
  void surface_flux (long lcid,long eid,long beid,double *fluxes);


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
