#ifndef BAREL3D_H
#define BAREL3D_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class barel3d defines onedimensional bar element with linear approximation functions
      
   JK, 27.9.2005
*/
class barel3d
{
 public:
  barel3d (void);
  ~barel3d (void);
  void dirvect (vector &s,vector &x,vector &y,vector &z);
  void giveloccoord(vector &x,vector &y,vector &z,vector &lx);

  void bf_matrix (matrix &n,vector &s,double xi);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,double &jac);
  //void geom_matrix (matrix &gm, double xi);
  double approx (double xi, vector &nodval);
  void transf_matrix(ivector &nodes, matrix &tmat);
  void bar_transf_mat(long eid, vector &x, vector &vec, matrix &tmat, vector &gx, vector &gy, vector &gz);
  //void tran_mat (vector &x, matrix &tmat,vector &gx,vector &gy);

  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y,vector &z);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y,vector &z);
  void res_mass_matrix (long eid,matrix &mm);
  void load_matrix (long eid, matrix &lm, vector &x);
  void res_load_matrix(long eid, matrix &lm);


  void res_ip_strains (long lcid,long eid);
  void ip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);

  void res_ip_stresses (long lcid,long eid);
  void ip_stresses (long lcid,long eid,long ri,long ci);
  void ip_elast_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);
  
  void nod_eqother_ip (long lcid,long eid,long ri,long ci);


  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void ipcoord(long eid, long ipp, long ri, long ci, vector &coord);
  void ipncoord(long eid, long ipp, vector &ncoord);
  void intpointval(long eid, vector &nodval, vector &ipval);
  void ipvolume (long eid,long ri,long ci);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv);
  
  void find_extreme_strains (vector &min,vector &max,long lcid,long eid);
  void find_extreme_stresses (vector &min,vector &max,long lcid,long eid);
  

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  order of integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  array of numbers of integration points in blocks
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
  /// total number of integration point
  long tnip;
};

#endif
