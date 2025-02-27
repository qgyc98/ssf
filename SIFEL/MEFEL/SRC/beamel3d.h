#ifndef BEAMEL3D_H
#define BEAMEL3D_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class beamel3d
   
   
   PF
*/
class beamel3d
{
 public:
  beamel3d (void);
  ~beamel3d (void);
  double approx (double xi,vector &nodval);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void beam_transf_matrix (matrix &tmat,double &dl,vector &vec,vector &x,vector &y,vector &z,long eid);
  void bf_matrix (matrix &n,double s,double dl,double gy,double gz);
  void bf_matrix_transl (matrix &n,double xi,double l,double gy,double gz);
  void geom_matrix (matrix &n,double s,double dl,double gy,double gz);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void stiffness_matrix_local (long eid,long ri,long ci,matrix &sm);
  void stiffness_matrix_local2 (long eid,long ri,long ci,matrix &sm);
  void stiffness_matrix_expl (long eid,long ri,long ci,matrix &sm);
  void nodeforces(long eid, long *le, double*nv, vector &nf);
  void res_stiffness_matrix (long eid,matrix &sm);
  void mass_matrix (long eid,long ri,long ci,matrix &mm);
  void res_mass_matrix (long eid,matrix &mm);
  void mass_matrix_expl (long eid,long ri,long ci,matrix &mm);
  void initstr_matrix (long eid,long ri,long ci,matrix &ism);
  void load_matrix (long eid,matrix &lm);

  void nodal_displ (long lcid,long eid);
  void nodal_forces (long lcid,long eid);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
//  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn) {};
  void define_meaning (long eid);


  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  total number of integration points
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
  ///  order of integration of initial stress matrix
  long intordism;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
