#ifndef PLELEMSUBQT_H
#define PLELEMSUBQT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
    class planeelemsubqt defines plane triangular subparametric element
    unknown functions are approximated by quadratic functions
    geometry of the element is approximated by linear functions
    
    JK
*/

class planeelemsubqt
{
 public:
  planeelemsubqt (void);
  ~planeelemsubqt (void);

  double approx (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta);
  void geom_matrix (matrix &gm,vector &x,vector &y,double xi,double eta);
  void geom_matrix_block (matrix &gm,double ri,vector &x,vector &y,double xi,double eta);
  void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y);
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
  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);

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
