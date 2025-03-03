#ifndef quadrilatkirch_H
#define quadrilatkirch_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   this class defines rectangular plate element based on Kirchhoff theory
   approximation function is bi-cubic
   
   degrees of freedom in nodes are
   w, dw/dx, dw,dy, d2w/dxdy
   
   19. 5. 2018, JK
*/
class quadrilatkirch
{
 public:
  quadrilatkirch (void);
  ~quadrilatkirch (void);
  
  double approx (double xi,double eta,vector &nodval);
  void pol_matrix (matrix &p,double &x,double &y);
  void bf_matrix (matrix &n,matrix &c,double &x,double &y);
  void lbf_matrix (matrix &nn,matrix &cc,double &x,double &y);
  void p_matrix (matrix &p,double &x,double &y);
  void s_matrix (matrix &s,vector &x,vector &y);
  void c_matrix (matrix &c,vector &x, vector &y);
  void cc_matrix (matrix &cc,vector &x, vector &y);
  void pp_matrix (matrix &pp,double &x, double &y);
  void g_matrix (matrix &g,double &x,double &y);
  void geom_matrix (matrix &gm,matrix &c,double &x,double &y);
  void h_matrix (matrix &h,double &x,double &y);
  void geom_initstiff_matrix (matrix &gm,matrix &c,double &x,double &y);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void transfmat (ivector &nodes,matrix &tmat);
  void res_stiffness_matrix (long eid,matrix &sm);
  
  void load_matrix (long eid,matrix &lm,vector &x,vector &y);
  void surfload (long lcid,long eid,double *nv,vector &nf);
  
  
  
  
  
  
  
  
  void initstr_matrix (long eid,long ri,long ci,matrix &ism);


  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);

  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);

  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
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
  
  void areaforces (long eid,double *nv,vector &lm);
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
  //  number of surfaces
  long nsurf;
  //  number of nodes on surface
  long nnsurf;

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
