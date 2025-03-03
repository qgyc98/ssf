#ifndef PLELEMROTLT_H
#define PLELEMROTLT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class planeelemrotlt defines plane triangular element with rotational
   degrees of freedom
   
   JK
*/

class planeelemrotlt
{
 public:
  planeelemrotlt (void);
  ~planeelemrotlt (void);


  double approx (vector &areacoord,vector &nodval);
  double approx_nat (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,vector &l,vector &b,vector &c);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &areacoord);
  void geom_matrix_block (matrix &gm,long ri,vector &x,vector &y,vector &areacoord);
  void addgeommat (matrix &gm,vector &x,vector &y,vector &areacoord,double &jac);
  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void stiffness_matrix2 (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void ip_stresses (long lcid,long eid,long ri,long ci);
  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);

  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
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
  void nod_others (long lcid,long eid,long ri,long ci);

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

  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipncoord (long eid,long ipp,vector &ncoord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  
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
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  array of orders of integration of stiffness matrix
  long **intordsm;
  ///  array of numbers of integration points
  long **nip;
  ///  order of integration of mass matrix
  long intordmm;
  ///  number of blocks
  long nb;
  ///  order of integration of edge load
  long intordb;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
