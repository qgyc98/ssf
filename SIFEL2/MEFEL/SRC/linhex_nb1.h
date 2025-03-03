#ifndef LINHEX_H
#define LINHEX_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class linhex defines hexahedral finite element with
   tri-linear approximation functions
   
   JK
*/

class linhex
{
 public:
  linhex (void);
  ~linhex (void);
  void eleminit (long eid);
  double approx (double xi,double eta,double zeta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,
		    double xi,double eta,double zeta,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,vector &z,
			  double xi,double eta,double zeta,double &jac);
  void bvectors (vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac,
		 vector &b11,vector &b12,vector &b13,
		 vector &b21,vector &b22,vector &b23,
		 vector &b31,vector &b32,vector &b33);
  void gngeom_matrix (matrix &gm,vector &r,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac);
  void gnl_grmatrix (matrix &grm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac);
  
  
  void transf_matrix (ivector &nodes,matrix &tmat);

  void gl_stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void gnl_stiffness_matrix (long lcid,long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long lcid,long eid,matrix &sm);

  void mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void res_mainip_strains (long lcid,long eid);
  void gl_mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r);
  void gnl_mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid,double **stra);
  void res_allip_strains (long lcid,long eid);
  void allip_strains (long lcid,long eid,long ri,long ci);

  //void elem_strains (double **stra,long lcid,long eid,long ri,long ci);
  void appstrain (long lcid,long eid,double xi,double eta,double zeta,long fi,long li,vector &eps);
  void strains (long lcid,long eid,long ri,long ci);

  void nodecoord (vector &xi,vector &eta,vector &zeta);
  void nodipnum (long eid,long ri,long ci,ivector &ipnum);
  //void appval (double xi,double eta,double zeta,long fi,long nc,vector &eps,double **val);

  void res_mainip_stresses (long lcid,long eid);
  void mainip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci);
  void appstress (long lcid,long eid,double xi,double eta,double zeta,long fi,long li,vector &sig);
  void res_allip_stresses (long lcid,long eid);
  void allip_stresses (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);
  
  void nod_others (long lcid,long eid,long ri,long ci);
  void nod_eqother_ip (long lcid,long eid,long ri,long ci);
  
  void gl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void gnl_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void internal_forces2 (long lcid,long eid,long ri,long ci,vector &ifor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void local_values (long lcid,long eid,long ri,long ci);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  
  void locglob_nodeval (long is,vector &nv,double *tnv,vector &x,vector &y,vector &z);
  void node_forces_surf (long lcid,long eid,long *is,double *nv,vector &nf);
  void node_forces_surf_old (long lcid,long eid,long *is,double *nv,vector &nf);
  void tran_mat(matrix &tran, vector &gx, vector &gy, vector &gz, long is);
  
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);
  void eigstrain_forces (long lcid,long eid,vector &nfor);
  void intpointval (long eid,vector &nodval,vector &ipval);

  void aver_strains (long lcid,long eid,long ri,long ci,vector &averstra,double &volume);
  
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points
  long tnip;
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
