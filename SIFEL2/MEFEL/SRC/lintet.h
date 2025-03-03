#ifndef LINTET_H
#define LINTET_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class lintet defines tetrahedral elements with linear approximation functions
   
   JK
*/

class lintet
{
 public:
  lintet (void);
  ~lintet (void);

  double approx (vector &volcoord,vector &nodval);
  double approx_nat (double xi,double eta,double zeta,vector &nodval);
  void bf_matrix (matrix &n,vector &volcoord);
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,double &jac);
  void geom_matrix_old (matrix &gm,vector &x,vector &y,vector &z);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void bd_matrix (long eid,long ri,long ci,matrix &bd);
  void dd_matrix (long eid,long ri,long ci,matrix &dd);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);

  void strains (long lcid,long eid,long ri,long ci);

  void res_ip_stresses (long lcid,long eid);
  void ip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void nod_stresses (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);

  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z);

  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y,vector &z);
  void elem_volintegration_quant(long eid, integratedquant iq, long lcid, vector &iv);
  
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipncoord (long eid,long ipp,vector &ncoord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  double volumeip (long eid,long ri,long ci);

  void nod_other_ip (long eid);
  void node_forces_surf (long lcid,long eid,long *is,double *nv,vector &nf);
  void node_forces_surf_old (long lcid,long eid,long *is,double *nv,vector &nf);
  void locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z);
  
  void intpointval (long eid,vector &nodval,vector &ipval);
  void res_eigstrain_forces (long eid,vector &nfor);
  void eigstrain_forces (long eid,vector &nfor);
  void define_meaning (long eid);
  
  /* termitovo */  
  void ntdbr_vector (long eid,vector &ntdbr);
  void ntn_matrix (long eid,matrix &ntn);
  double compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull);
  void elchar (long eid, matrix &spsig);
  /* termitovo */  

  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges on one element
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
  ///  order of numerical interation on surfaces
  long intordb;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  array of numbers of components of blocks
  long *ncomp;
  ///  cumulative array of numbers of components of blocks
  long *cncomp;
  ///  stress/strain state
  strastrestate ssst;

};

#endif
