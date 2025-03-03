#ifndef PLELEMROTLQ_H
#define PLELEMROTLQ_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class planeelemrotlq defines plane quadrilateral finite element with
   rotational degrees of freedom
   
   JK
*/

class planeelemrotlq
{
 public:
  planeelemrotlq (void);
  ~planeelemrotlq (void);
  void auxdata (vector &x,vector &y,vector &l,vector &nx,vector &ny);
  double approx (double xi,double eta,vector &nodval);

  void bf_matrix (matrix &n,double xi,double eta,vector &l,vector &nx,vector &ny);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,
		    vector &l,vector &nx,vector &ny,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,
			  vector &l,vector &nx,vector &ny,double &jac);
  void addgeommat (matrix &gm,vector &x,vector &y,double xi,double eta,
		   vector &l,vector &nx,vector &ny,double &jac);
  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm,vector &x,vector &y);
  void res_load_matrix (long eid,matrix &lm);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid,double **stra);
  //void strains (long lcid,long eid,long ri,long ci);

  void nodecoord (vector &xi,vector &eta);

  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);


  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y);

  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void surfload (long eid,double *nv, vector &x, vector &y,vector &nf);
  void node_forces_surf (long eid,double *nodvals,vector &nf);

  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipncoord (long eid,long ipp,vector &ncoord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  number of components of stress and strain tensors
  long tncomp;
  ///  number of components for graphic purposes
  long gncomp;
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
  ///  number of surface
  long nsurf;
  ///  number of nodes on surface
  long nnsurf;
  ///  order of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration on edges
  long intordb;
  ///  number of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
