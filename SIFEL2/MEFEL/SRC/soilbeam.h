#ifndef SOILBEAM_H
#define SOILBEAM_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

class soilbeam
{
 public:
  soilbeam (void);
  ~soilbeam (void);

  void transf_matrix (ivector &nodes,matrix &tmat);
  void geom_matrix (matrix &n,double s,double dl,double gy,double gz);
  void beam_transf_matrix (matrix &tmat,double &dl,vector &vec,vector &x,vector &y,vector &z,long eid);
  void res_stiffness_matrix (long eid,matrix &sm);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void strains (long lcid,long eid,long ri,long ci);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void internal_forces1 (long lcid,long eid,long ri,long ci,vector &ifor);
//  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn) {};



  //  number of DOFs on the element
  long ndofe;
  //  number of nodes on one element
  long nne;
  //  total number of components of stress and strain tensors
  long tncomp;
  //  array containing numbers of components of stress and strain tensors
  long *ncomp;
  //  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  //  number of approximated functions on the element
  long napfun;
  //  order of integration of stiffness matrix
  long **intordsm;
  //  order of integration of mass matrix
  long intordmm;
  //  order of integration of initial stress matrix
  long intordism;
  //  array of numbers of integration points in sets
  long **nip;
  //  numbers of integration points in sets !is sum nip
  long tnip;
  //  number of blocks
  long nb;
  //  stress/strain state
  strastrestate ssst;

  double *c1,*c2;
  double bPod;
};

#endif
