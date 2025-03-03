#ifndef PLELEMQQ_H
#define PLELEMQQ_H

#include "alias.h"
#include "galias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class planeelemqq defines plane quadrilateral element with
   bi-quadratic approximation functions
   
   JK
*/

class planeelemqq
{
 public:
  planeelemqq (void);
  ~planeelemqq (void);

  double approx (double xi,double eta,vector &nodval);

  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac);
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
  void nod_strains_comp (long lcid,long eid);
  void nod_strains (long lcid,long eid);
  void strains (long lcid,long eid,long ri,long ci);
  void nodipnum (long eid,long ri,long ci,ivector &ipnum);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);

  void res_ip_stresses (long lcid,long eid);
  void ip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void nod_stresses_comp (long lcid,long eid,long ri,long ci,double **stra,double **stre);
  void stresses (long lcid,long eid,long ri,long ci);

  void nod_other_ip (long eid,long ri,long ci);

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
  
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipcoordblock (long eid,long ri,long ci,double **coord);
  void ipncoord (long eid,long ipp,vector &ncoord);

  void res_nodeforces (long eid,long *le,double *nv,vector &nf);
  void nodeforces (long eid,long *le,double *nv,vector &nf,vector &x,vector &y);

  void midpoints (long eid,long dd,double *valel);

  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  void intpointval (long eid,vector &nodval,vector &ipval);
  void intpointval2 (long eid,vector &nodval,vector &ipval);
  

  /* termitovo */
  void ntdbr_vector (long eid,vector &ntdbr);
  void ntn_matrix (long eid,matrix &ntn);
  double compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull);
  void elchar (long eid, matrix &spsig);
  /* termitovo */

  void mechq_nodval (long eid,vector &nodval,nontransquant qt);
  void mechq_nodval2 (long eid,vector &nodval,nontransquant qt);
  void mechq_nodval_comp (long eid, vector &nodval, long ncnv, long nq, nontransquant *qt);

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  number of components for graphic purposes
  long gncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array of numbers of components
  long *ncomp;
  ///  array of cumulative numbers of components
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  array of orders of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration on edges
  long intordb;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;

};

#endif
