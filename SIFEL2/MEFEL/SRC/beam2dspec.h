#ifndef BEAM2DSPEC_H
#define BEAM2DSPEC_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;


/**
   class beam2dspec defines beam element for 2D analysis
   
   
   JK
*/
class beam2dspec
{
 public:
  beam2dspec (void);
  ~beam2dspec (void);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void beam_transf_matrix (double c,double s,matrix &tmat);
  void res_stiffness_matrix (long eid,matrix &sm);
  void stiffness_matrix_expl (long eid,long ri,long ci,matrix &sm);
  void res_mass_matrix (long eid,matrix &mm);
  void mass_matrix_expl (long eid,long ri,long ci,matrix &mm);
  void internal_forces (long lcid,long eid,vector &ifor);
  void res_internal_forces (long lcid,long eid,vector &ifor);
  void stresses (long eid,long lcid);
  void strains (long eid,long lcid);

  void res_mainip_strains (long lcid,long eid);
  void res_mainip_stresses (long lcid,long eid);

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
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration of initial stress matrix
  long intordism;
  ///  array of numbers of integration points in sets
  long **nip;
  ///   total number of integration points;
  long tnip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;
};

#endif
