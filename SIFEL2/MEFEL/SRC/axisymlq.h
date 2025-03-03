#ifndef AXISYMLQ_H
#define AXISYMLQ_H

#include "alias.h"
#include "galias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class axisymlq defines quadrilateral finite element with
   bi-linear approxiamtion functions for axisymmetric problem
   
   strain components are ordered in the following way:
   epsilon_x = du/dx
   epsilon_y = dv/dy
   epsilon_fi = u/r
   epsilon_xy = du/dy + dv/dx
   
   if three blocks are used, strain components are divided
   in the following way:
   1. block
   epsilon_x = du/dx
   epsilon_y = dv/dy
   2. block
   epsilon_fi = u/r
   3. block
   epsilon_xy = du/dy + dv/dx
   
   JK
*/

class axisymlq
{
 public:
  axisymlq (void);
  ~axisymlq (void);

  double approx (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac);

  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);

  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void stiffness_matrix_blocks (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);
  void res_mass_matrix (long eid,matrix &mm);
  
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  
  void edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
		      vector &coef,matrix &km);
  void edgenodeval (long edg,vector &nodval,double *list);
  void edgeload (long eid,long *le,double *list,vector &nf);


  
  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void res_allip_strains (long lcid,long eid);
  void allip_strains (long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);
  
  void nodecoord (vector &xi,vector &eta);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);
  
  void res_mainip_stresses (long lcid,long eid);
  void mainip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses (long lcid,long eid,long ri,long ci);
  void elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci);
  void appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig);
  void res_allip_stresses (long lcid,long eid);
  void allip_stresses (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);
  
  void nodipnum (long eid,long ri,long ci,ivector &ipnum);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid);
  void nod_strains (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void nod_other_ip (long eid);

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
  
  void load_matrix (long eid,matrix &lm);
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  
  // function assembles global coordinates of the given integration point of element
  void ipcoord(long eid, long ipp, long ri, long ci, vector &coord);
  // function assembles natural coordinates of the given integration point of element
  void ipncoord(long eid, long ipp, vector &ncoord);
  void intpointval(long eid, vector &nodval, vector &ipval);
  
  void mechq_nodval (long eid,vector &nodval,nontransquant qt);
  void mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt);

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges
  long ned;
  ///  number of nodes on edge
  long nned;
  ///  array of orders of numerical integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration on edges
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


