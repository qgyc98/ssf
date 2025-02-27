#ifndef QUADTET_H
#define QUADTET_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
  Class quadtet defines tetraheadronal finite element with
  tri-quadratic approximation functions.
   
  Created by Tomas Koudelka
*/

class quadtet
{
 public:
  quadtet (void);
  ~quadtet (void);

  /// computes approximated quantity given by the nodal values at the given point
  double approx (double xi,double eta,double zeta,vector &nodval);
  /// computes base functions %matrix
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  /// computes geometrical %matrix
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
		    double xi,double eta,double zeta,double &jac);
  /// computes block of geometrical %matrix
  void geom_matrix_block (matrix &gm,vector &x,vector &y,vector &z,
			  double xi,double eta,double zeta,double &jac);
  /// computes transformation %matrix
  void transf_matrix (ivector &nodes,matrix &tmat);
  /// computes stiffness %matrix
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  /// computes resulting stiffness %matrix
  void res_stiffness_matrix (long eid,matrix &sm);
  /// computes mass %matrix
  void mass_matrix (long eid,matrix &mm);
  /// computes resulting mass %matrix
  void res_mass_matrix (long eid,matrix &mm);
  /// computes load %matrix
  void load_matrix (long eid,matrix &lm);
  /// computes resulting load %matrix
  void res_load_matrix (long eid,matrix &lm);
  /// computes volume appropriate to integration point
  double volumeip(long eid, double w);

  /// computes resulting strains in the integration points
  void res_ip_strains (long lcid,long eid);
  /// computes strains in the integration points
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r);
  /// copies strains from integration point to nodes
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  /// computes strains at nodes
  void nod_strains_comp (long lcid,long eid,double **stra);
  /// computes required strains
  void strains (long lcid,long eid,long ri,long ci);

  /// computes resulting stresses in the integration points
  void res_ip_stresses (long lcid,long eid);
  /// computes stresses in the integration points
  void ip_stresses (long lcid,long eid,long ri,long ci);
  /// computes stresses in the integration points, elastic law is used  
  void ip_elast_stresses (long lcid,long eid,long ri,long ci);
  /// copies stresses from integration point to nodes
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  /// computes stresses at nodes
  void nod_stresses_comp (long lcid,long eid,long ri,long ci,double **stra,double **stre);
  /// computes required stresses
  void stresses (long lcid,long eid,long ri,long ci);
  
  /// copies other values from integration point to nodes
  void nod_other_ip (long eid,long ri,long ci);

  /// computes internal forces of the element
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  /// computes internal forces of the element nonlocally
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  /// computes increments of internal forces of the element
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  /// computes nodal forces of the element caused by eigenstrains
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z);

  /// computes resulting internal forces of the element
  void res_internal_forces (long lcid,long eid,vector &ifor);
  /// computes resulting internal forces of the element nonlocally
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  /// computes resulting increments of internal forces of the element
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  /// computes resulting nodal forces of the element caused by eigenstrains
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  /// computes stresses in the integration points, laws of specified material models is used for stress evaluation
  void compute_nlstress (long lcid,long eid,long ri,long ci);
  /// computes increments of stresses in the integration points, specified material model is used for stress evaluation
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  /// computes local values of averaged quantities (used in nonlocal material models)
  void local_values (long lcid,long eid,long ri,long ci);
  /// computes stresses in the integration points, nonlocal material model is used for stress evaluation
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  /// computes stresses due to eigenstrains
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  /// integrates required quantity over the element
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y,vector &z);
  
  /// returns global coordinates of integration points of the element
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  /// returns natural coordinates of integration points of the element
  void ipncoord (long eid,long ipp,vector &ncoord);
  /// sets internal data of integration points to specified initial values
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  /// computes volumes appropriated to integration points
  void ipvolume (long eid,long ri,long ci);

  /// computes nodal forces caused by surface load  
  void node_forces_surf (long lcid,long eid,long *is,double *nv,vector &nf);
  /// transforms surface load given in the local coordinate system to the global coordinate system
  void locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z);

  /// approximates nodal values to the the integration point values (quadratic approximation is used)
  void intpointval (long eid,vector &nodval,vector &ipval);
  /// approximates nodal values to the the integration point values (linear approximation is used)  
  void intpointval2 (long eid,vector &nodval,vector &ipval);
  /// computes averaged strains
  void aver_strains (long lcid,long eid,long ri,long ci,vector &averstra,double &volume);
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of surfaces
  long nsurf;
  ///  number of nodes on one surface
  long nnsurf;
  ///  array of orders of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  array of numbers of integration points in sets
  long intordb;
  ///  array of numbers of integration points on surface
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
