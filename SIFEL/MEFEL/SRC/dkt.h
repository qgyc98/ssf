#ifndef DKT_H
#define DKT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

class dktelem
{
 public:
  dktelem (void);
  ~dktelem (void);

  double approx (vector &areacoord,vector &nodval);

  void geom_matrix (matrix &gm,vector &x,vector &y,vector &l);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void dbmat (matrix &d,matrix &db,double t);
  void res_stiffness_matrix (long eid,matrix &sm);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void nodecoord (vector &xi,vector &eta);
  void appval (vector &l, long fi,long nc,vector &eps,double **val);

  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void res_ip_strains (long lcid,long eid);
  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);

  void strains (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri, long ci);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);

  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);

  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y);

  void nodeforces (long eid,long *le,double *nv,vector &nf);
  void areaforces (long eid,double *nv,vector &nf);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);


  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  long *cncomp;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  /// number of surfaces
  long nsurf;
  /// number of nodes on surface
  long nnsurf;  
  //  order of integration of stiffness matrix
  long **intordsm;
  //  array of numbers of integration points in sets
  long **nip;
  //  number of blocks
  long nb;
  //  order of integration for mass matrix
  long intordmm;
  //  stress/strain state
  strastrestate ssst;
};

#endif
