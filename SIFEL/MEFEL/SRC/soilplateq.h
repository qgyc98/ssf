#ifndef SOILPLATEQ_H
#define SOILPLATEQ_H

#include "alias.h"
#include "q4plate.h"
struct matrix;
struct vector;
struct ivector;

class soilplateq
{
 public:
  soilplateq (void);
  ~soilplateq (void);

  double approx (double xi,double eta,vector &nodval);
  void atd_matrix (matrix &atd,vector &x,vector &y);
  void dbmat (matrix &db,double c11, double c22,long ri, long ci);
  void transfmat (ivector &nodes,matrix &tmat);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);
  void nodecoord (vector &xi,vector &eta);
  void geom_matrix (matrix &gm,matrix &atd,vector &x,vector &y, vector &l);
  void geom_matrix_shear (matrix &gm,matrix &atd,vector &x,vector &y,vector &l);
  void res_stiffness_matrix (long eid,matrix &sm);
  void stiffness_matrix (long eid,long ri,long ci, matrix &sm,vector &x, vector &y);

  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains (long lcid,long eid,long ri,long ci);
  void elem_strains (double **stra,long lcid,long eid,long ri,long ci);
  void appstrain (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &eps);
  void allip_strains (double **stra,long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);

  void res_allip_stresses (long lcid,long eid);
  void allip_stresses (long lcid,long eid,long ri,long ci);

  
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y, vector &z);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y, vector &z);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y, vector &z);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y, vector &z);

  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq, long lcid, long eid, long ri, long ci, vector &nv, vector &x, vector &y, vector &z);
  
//  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn) {};



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
  //  array of orders of integration of stiffness matrix
  long **intordsm;
  //  order of integration for mass matrix
  long intordmm;
  //  order of integration on edges
  long intordb;
  //  array of numbers of integration points in sets
  long **nip;
  //  number of blocks
  long nb;
  //  stress/strain state
  strastrestate ssst;
};

#endif
