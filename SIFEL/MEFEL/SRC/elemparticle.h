#ifndef ELEMPARTICLE_H
#define ELEMPARTICLE_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class elemparticle serves for computations with particles
   it is intended for molecular dynamics
   functions from this class define particles which influence each other
      
   JK, 19.6.2005
*/
class elemparticle
{
 public:
  elemparticle (long cnne,long cdim);
  ~elemparticle (void);
  
  void direction_vector_1d (long eid,long i,long j,vector &s,vector &x,vector &u);
  void direction_vector_2d (long eid,long i,long j,vector &s,
			    vector &x,vector &y,vector &u,vector &v);
  void direction_vector_3d (long eid,long i,long j,vector &s,
			    vector &x,vector &y,vector &z,
			    vector &u,vector &v,vector &w);
  void stiffmat_1d_kii (long ipp,matrix &k,vector &s);
  void stiffmat_1d_kij (long ipp,matrix &k,vector &s);
  void stiffness_matrix_1d (long eid,matrix &sm);

  void stiffmat_2d_kii (long ipp,matrix &k,vector &s);
  void stiffmat_2d_kij (long ipp,matrix &k,vector &s);
  void stiffness_matrix_2d (long eid,matrix &sm);

  void stiffmat_3d_kii (long ipp,matrix &k,vector &s);
  void stiffmat_3d_kij (long ipp,matrix &k,vector &s);
  void stiffness_matrix_3d (long eid,matrix &sm);

  void res_stiffness_matrix (long eid,matrix &sm);
  
  void forces_1d (long ipp,vector &fij,vector &s);
  void inter_forces_1d (long eid,vector &f);
  void forces_2d (long ipp,vector &fij,vector &s);
  void inter_forces_2d (long eid,vector &f);
  void forces_3d (long ipp,vector &fij,vector &s);
  void inter_forces_3d (long eid,vector &f);
  void res_internal_forces (long eid,vector &ifor);
  
  ///  dimension of problem
  long dim;
  ///  number of nodes on one element
  long nne;
  ///  number of DOFs
  long ndofe;
  ///  number of blocks
  long nb;
  /// total number of integration point
  long tnip;
  
  ///  array of numbers of integration points in blocks
  long **nip;

};

#endif
