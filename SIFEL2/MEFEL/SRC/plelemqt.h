#ifndef PLELEMQT_H
#define PLELEMQT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class planeelemqt defines plane triangular element with quadratic
   approximation functions
   
   JK
*/

class planeelemqt
{
 public:
  planeelemqt (void);
  ~planeelemqt (void);

  double approx (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta,double &jac);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,double xi,double eta,double &jac);
  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x, vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains (long lcid,long eid,long ri,long ci);
  void elem_strains (double **stra,long lcid,long eid,long ri,long ci);
  void appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps);
  void allip_strains (double **stra,long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);
  
  void nodecoord (vector &xi,vector &eta);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);
  
  void mainip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses (long lcid,long eid,long ri,long ci);
  void elem_stresses (double **stra,double **stre,long lcid,long eid,long ri,long ci);
  void appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig);
  void allip_stresses (double **stre,long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);
  void res_ip_stresses (long lcid,long eid);

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
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);

  
  /* termitovo */
  void ntdbr_vector (long eid,vector &ntdbr);
  void ntn_matrix (long eid,matrix &ntn);
  double compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull);
  void give_sig_roof (matrix &d,double *areacoord,double *x,double *y,long *etnodes,double **xyrr,vector &sig_roof,vector &sig_star);
  //void compute_error_fine (long eid,long *etnodes,double **xyrr,vector *rsigfull,double &e2);
  void elchar (long eid, matrix &spsig);
  /* termitovo */


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
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing cumulative numbers of components of stress and strain tensors
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
