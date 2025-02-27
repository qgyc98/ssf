#ifndef AXISYMQT_H
#define AXISYMQT_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class axisymqt defines triangular element for axisymmetric problems with quadratic shape functions
   
   JK, 3. 11. 2020
*/

class axisymqt
{
 public:
  axisymqt (void);
  ~axisymqt (void);

  double approx (vector &areacoord,vector &nodval);
  double approx_nat (double xi,double eta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta);
  void bf_matrix (matrix &n,vector &areacoord);
  void geom_matrix (matrix &gm,double xi, double eta,vector &x,vector &y, double &jac);
  //void geom_matrix_block (matrix &gm,long ri,vector &areacoord,vector &x,vector &y);
  //void dmatblock (long ri,long ci,matrix &d, matrix &dd);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  //void stiffness_matrix_block (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm);
  void res_mass_matrix (long eid,matrix &mm);

  void edge_integral (long edg,vector &x,vector &y,long intord,vector &gp,vector &w,
		      vector &coef,matrix &km);
  void edgenodeval (long edg,vector &nodval,double *list);
  void edgeload (long eid,long *le,double *list,vector &nf);



  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &r);
  void nod_strains_ip (long lcid,long eid);
  void nod_strains_comp (long lcid,long eid);
  void res_allip_strains (long lcid,long eid);
  void allip_strains (long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);
  void nodecoord (vector &xi,vector &eta);
  void nodipnum (long eid,ivector &ipnum);
  void appval (double xi,double eta,long fi,long nc,vector &eps,double **val);
  void res_mainip_stresses (long lcid,long eid);
  void mainip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid);
  void nod_stresses_comp (long lcid,long eid);
  void res_allip_stresses (long lcid,long eid);
  void allip_stresses (long lcid,long eid,long ri,long ci);
  void appstress (long lcid,long eid,double xi,double eta,long fi,long ncomp,vector &sig);
  void stresses (long lcid,long eid,long ri,long ci);
  void nod_other_ip (long eid);
  void load_matrix (long eid,matrix &lm);

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
  void ipncoord (long eid,long ipp,vector &coord);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);


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
  ///  total number of integration points on element
  long tnip;
  ///  array of numbers of components of stress and strain tensors
  long *ncomp;
  ///  array of cumulative numbers of components of stress and strain tensors
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
