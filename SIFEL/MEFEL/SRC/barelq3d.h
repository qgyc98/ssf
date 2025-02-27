#ifndef BARELQ3D_H
#define BARELQ3D_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class barelq3d defines onedimensional bar element with quadratic approximation functions
   
   
   JK, 28.9.2005
*/
class barelq3d
{
 public:
  barelq3d (void);
  ~barelq3d (void);

  void dirvect (vector &s,vector &x,vector &y,vector &z);
  void giveloccoord(vector &x,vector &y,vector &z,vector &lx);
  double approx (double xi,vector &nodval);

  void bf_matrix (matrix &n,vector &s,double xi);
  void geom_matrix (matrix &gm, vector &x,vector &y,vector &z,double xi,double &jac);

  //void give_glob_loc_tmat(vector &x, vector &y, matrix &tmat);
  //void give_loc_glob_tmat(vector &x, vector &y, matrix &tmat);
  void transf_matrix (ivector &nodes,matrix &tmat);

  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y,vector &z);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,matrix &mm,vector &x,vector &y,vector &z);
  void res_mass_matrix (long eid,matrix &sm);


  void res_mainip_strains (long lcid,long eid);
  void mainip_strains (long lcid,long eid,long ri,long ci,vector &x,vector &y,vector &z,vector &r);
  void nod_strains_ip (long lcid,long eid,long ri,long ci);
  void nod_strains_comp (long lcid,long eid,double **stra);
  void res_allip_strains (long lcid,long eid);
  void allip_strains (long lcid,long eid,long ri,long ci);
  void strains (long lcid,long eid,long ri,long ci);

  void nodecoord (vector &xi);
  void nodipnum (long eid,long ri,long ci,ivector &ipnum);

  void res_mainip_stresses (long lcid,long eid);
  void mainip_stresses (long lcid,long eid,long ri,long ci);
  void nod_stresses_ip (long lcid,long eid,long ri,long ci);
  void nod_stresses_comp (long lcid,long eid,long ri,long ci,double **stra,double **stre);
  void res_allip_stresses (long lcid,long eid);
  void allip_stresses (long lcid,long eid,long ri,long ci);
  void stresses (long lcid,long eid,long ri,long ci);


  void nod_eqother_ip (long lcid,long eid,long ri,long ci);


  void load_matrix (long eid,matrix &lm);
  void temperaturestrains (long lcid,long eid,long ri,long ci);
  void tempstrains (long lcid,long eid,long ipp,double xi,vector &eps);

  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  void nonloc_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  void incr_internal_forces (long lcid,long eid,long ri,long ci,vector &ifor,vector &x,vector &y,vector &z);
  void eigstrain_forces (long lcid,long eid,long ri,long ci,vector &nfor,vector &x,vector &y,vector &z);

  void res_internal_forces (long lcid,long eid,vector &ifor);
  void res_nonloc_internal_forces (long lcid,long eid,vector &ifor);
  void res_incr_internal_forces (long lcid,long eid,vector &ifor);
  void res_eigstrain_forces (long lcid,long eid,vector &nfor);

  void compute_nlstress (long lcid,long eid,long ri,long ci);
  void compute_nlstressincr (long lcid,long eid,long ri,long ci);
  void local_values (long lcid,long eid,long ri,long ci);
  void compute_nonloc_nlstress (long lcid,long eid,long ri,long ci);
  void compute_eigstress (long lcid,long eid,long ri,long ci);
  void elem_integration (integratedquant iq,long lcid,long eid,long ri,long ci,vector &nv,vector &x,vector &y,vector &z);


  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  void intpointval (long eid,vector &nodval,vector &ipval);
  void intpointval2 (long eid,vector &nodval,vector &ipval);
  void ipcoord(long eid, long ipp, long ri, long ci, vector &coord);
  void ipncoord(long eid, long ipp, vector &ncoord);

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  long tnip;
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
  ///  computer zero
  double zero;

};

#endif
