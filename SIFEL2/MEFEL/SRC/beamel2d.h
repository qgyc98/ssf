#ifndef BEAMEL2D_H
#define BEAMEL2D_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class beamel2d defines beam element for 2D analysis
   
   
   JK
*/
class beamel2d
{
 public:
  beamel2d (void);
  ~beamel2d (void);
  double approx (double xi,vector &nodval);
  void bf_matrix (matrix &n,double xi,double l,double kappa);
  void dbf_matrix (matrix &n,double xi,double l,double kappa);
  void geom_matrix (matrix &gm,double xi,double l,double kappa);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void beam_transf_matrix (long eid,vector &x,vector &z,matrix &tmat);

  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void stiffness_matrix_expl_local (long eid,long ri,long ci,double l,matrix &sm);
  void stiffness_matrix_expl (long eid,long ri,long ci,matrix &sm);
  void mass_matrix (long eid,long ri,long ci,matrix &mm);
  void res_mass_matrix (long eid,matrix &mm);
  void mass_matrix_expl (long eid,long ri,long ci,matrix &mm);
  void initstr_matrix (long eid,long ri,long ci,matrix &ism);
  void initstr_matrix_expl (long lcid,long eid,long ri,long ci,matrix &ism);

  void nodal_displ (long lcid,long eid);
  void nodal_forces (long lcid,long eid);
  void nodeforces(long eid, long *le, double*nv, vector &nf);
  void beamnodeforces(long eid, elloadmeaning elm, double nnv, double *la, double *lf, double*nv, vector &nf);

  void internal_forces (long lcid,long eid,vector &ifor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void define_meaning (long eid);

  
  //  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn) {};

  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
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
  ///  order of integration of initial stress matrix
  long intordism;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
  ///  coordinate system (xz, xy)
  long coordsys;
  
  long tnip;
};

#endif