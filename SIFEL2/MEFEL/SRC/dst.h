#ifndef DST_H
#define DST_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

class dstelem
{
 public:
  dstelem (void);
  ~dstelem (void);

  double approx (vector &areacoord,vector &nodval);

  void tran_matrix (matrix &a,matrix &ct,long eid);
  void geom_matrix_bconst (matrix &gm,vector &x,vector &y);
  void geom_matrix_bending (matrix &gm,matrix &a,matrix &ct,vector &x,vector &y,vector &l);
  void geom_matrix_shear (matrix &gs,matrix &a,matrix &ct,long eid);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void dmatblock (matrix &dd,matrix &d,long ri, long ci, double t);
  void stiffness_matrix (long eid,long ri, long ci, matrix &sm,vector &x, vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  void nodecoord (vector &xi,vector &eta);
  void appval (vector &l, long fi,long nc,vector &eps,double **val);

  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void res_ip_strains (long lcid,long eid);
  void strains (long lcid,long eid,long ri, long ci);
  void res_ip_stresses (long lcid,long eid);
  void stresses (long lcid,long eid,long ri, long ci);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);

  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);

  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);

  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void areaforces (long eid,double *nv,vector &nf);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);

  //  number of DOFs on the element
  long ndofe;
  //  number of nodes on one element
  long nne;
  //  total number of components of stress and strain tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  //  number of components of stress and strain tensors
  long *ncomp;
  //  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  //  number of approximated functions on the element
  long napfun;
  //  number of edges on one element
  long ned;
  //  number of nodes on one edge
  long nned;
  //  order of integration of stiffness matrix
  long **intordsm;
  //  order of integration for mass matrix
  long intordmm;
  //  order of integration on edges
  long intordb;
  //  number of integration points in sets
  long **nip;
  //  number of blocks
  long nb;
  //  stress/strain state
  strastrestate ssst;
};

#endif
