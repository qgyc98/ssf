#ifndef BEAMGEN3D_H
#define BEAMGEN3D_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class beamgen3d describes threedimensional beam with general
   cross section where centre of shear may be different from center of gravity
   
   
   PF, 2.10.2006
*/
class beamgen3d
{
 public:
  beamgen3d (void);
  ~beamgen3d (void);

  void transf_matrix (ivector &nodes,matrix &tmat);
  void beam_transf_matrix (matrix &tmat,double &dl,vector &vec,vector &x,vector &y,vector &z,long eid);
  void ck_matrix (matrix &ck,  double s,double c,double eh,double dl);
  void geom_matrix (matrix &n,double s,double dl,double gy,double gz);
  void bf_matrix (matrix &n,double s,double dl,double gy,double gz);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);
  void stiffness_matrixtor (matrix &sm, double dl,double e,double g,double gy,double gz,double *ixyz,double *ioyz);
  void stiffness_matrixtor1 (matrix &sm, double dl,double e,double g,double gy,double gz,double *ixyz,double *ioyz);
  void stiffness_matrixtor2 (matrix &sm, double dl,double e,double g,double gy,double gz,double *ixyz,double *ioyz,double *iro);
  void load_matrix (long eid,matrix &lm);

  void nodal_displ (long eid,long lcid);
  void nodal_forces (long eid,long lcid);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void internal_forces (long lcid,long eid,long ri,long ci,vector &ifor);
  void internal_forces1 (long lcid,long eid,long ri,long ci,vector &ifor);
//  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn) {};

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
