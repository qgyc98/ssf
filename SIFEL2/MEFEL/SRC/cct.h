#ifndef CCT_H
#define CCT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class cctelem defines triangular plate element with constant curvatures
   based on the Mindlin theory
   
   JK
*/

class cctelem
{
 public:
  cctelem (void);
  ~cctelem (void);

  double approx (vector &areacoord,vector &nodval);
  double approx_nat (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,vector &x,vector &y,vector &areacoord);
  void geom_matrix_block (matrix &gm,long bi,vector &x,vector &y,vector &areacoord);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &areacoord);
  void dmatblock (long ri,long ci,matrix &d,matrix &dd,double t);
  void dmat (matrix &d,double t);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid,matrix &lm);

  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void res_ip_strains (long lcid,long eid);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);

  void strains (long lcid,long eid,long ri,long ci);
  
  void nodecoord (vector &xi,vector &eta);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);
  
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
  
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);


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
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  array of orders of integration of stiffness matrix
  long **intordsm;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  array of numbers of strain/stress components
  long *ncomp;
  ///  array of cumulative numbers of strain/stress components
  long *cncomp;
  ///  number of blocks
  long nb;
  ///  order of integration for mass matrix
  long intordmm;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
