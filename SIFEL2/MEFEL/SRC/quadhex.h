#ifndef QUADHEX_H
#define QUADHEX_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class quadhex defines hexahedral finite element with
   tri-quadratic approximation functions
   
   JK
*/

class quadhex
{
 public:
  quadhex (void);
  ~quadhex (void);
  double approx (double xi,double eta,double zeta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
		    double xi,double eta,double zeta,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,vector &z,
			  double xi,double eta,double zeta,double &jac);
  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);
  void res_load_matrix (long eid,matrix &lm);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid,double **stra);
  void strains (long lcid,long eid,long ri,long ci);

  void nodecoord (vector &xi,vector &eta,vector &zeta);


  void res_ip_stresses (long lcid,long eid);
  void ip_stresses (long lcid,long eid,long ri,long ci);
  void ip_elast_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void nod_stresses_comp (long lcid,long eid,long ri,long ci,double **stra,double **stre);
  void stresses (long lcid,long eid,long ri,long ci);
  
  void nod_other_ip (long eid,long ri,long ci);
  
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
  
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipncoord (long eid,long ipp,vector &ncoord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  
  void node_forces_surf (long lcid,long eid,long *is,double *nv,vector &nf);
  void tran_mat(matrix &tran, vector &gx, vector &gy, vector &gz, long is);

  void intpointval (long eid,vector &nodval,vector &ipval);
  void intpointval2 (long eid,vector &nodval,vector &ipval);
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
