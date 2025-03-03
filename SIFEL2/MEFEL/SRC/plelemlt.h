#ifndef PLELEMLT_H
#define PLELEMLT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class planeelemlt - defines plane triangular element with linear
   approximation functions
   
   
   JK
*/

class planeelemlt
{
 public:
  planeelemlt (void);
  ~planeelemlt (void);

  double approx (vector &areacoord,vector &nodval);
  double approx_nat (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void bd_matrix (long eid,long ri,long ci,matrix &mbd);
  void dd_matrix (long eid,long ri,long ci,matrix &mdd);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm,vector &x,vector &y);
  void res_load_matrix (long eid,matrix &lm);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid);
  void strains (long lcid,long eid,long ri,long ci);
  
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);

  void res_ip_stresses (long lcid,long eid);
  void ip_elast_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);

  void nod_others (long eid,long ri,long ci);
  
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y);

  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void intpointval (long eid,vector &nodval,vector &ipval);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);
  void elem_volintegration_quant (long eid, integratedquant iq, long lcid, vector &iv);
  
  void ipcoord (long eid,long ipp,long ri,long ci,vector &coord);
  void ipcoordblock (long eid,long ri,long ci,double **coord);
  void ipncoord (long eid,long ipp,vector &ncoord);
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void locglob_nodeval (long edid,vector &nv,vector &tnv,vector &x,vector &y);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void ipvolume (long eid,long ri,long ci);
  

  /* termitovo */
  void ntdbr_vector (long eid,vector &ntdbr);
  void ntn_matrix (long eid,matrix &ntn);
  double compute_error (long eid, double &e2, double &u2, double &sizel, vector *rsigfull, long flags);
  void elchar (long eid,double *&spsig,long flags);
  double error (long eid,vector &n,double &a);
  /* termitovo */

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  number of components for graphic purposes
  long gncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components
  long *ncomp;
  ///  array of cumulative numbers of components
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  array of orders of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  order of integration of edge load
  long intordb;
  ///  stress/strain state
  strastrestate ssst;
  
};

#endif
