#ifndef SHELLTR_H
#define SHELLTR_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   shell triangular element
   resulting matrices are composed from plane triangular element
   with rotational degrees of freedom and from constant curve triangle
   
*/
class shelltr
{
 public:
  shelltr (void);
  ~shelltr (void);

  void codnum (long *cn,long ri);
  void coord_transf_matrix (matrix &tran, vector &gx, vector &gy, vector &gz,double zero);
  void local_coordinates (vector &gx,vector &gy,vector &gz,vector &lx,vector &ly,double zero);
  void elem_transf_matrix (matrix &tran, vector &gx, vector &gy, vector &gz,double zero);
  void node_transf_matrix (ivector &nodes,matrix &tmat);

  void res_stiffness_matrix (long eid,matrix &sm);
  
  void res_ip_strains (long lcid,long eid);
  void res_ip_stresses (long lcid,long eid);
  void compute_nlstress (long lcid,long eid);
  void compute_nlstressincr (long lcid,long eid);
  void res_internal_forces (long lcid,long eid,vector &ifor);

  void nod_strains_ip (long lcid,long eid);
  void strains (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid);
  void stresses (long lcid,long eid);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);


  ///  number of nodes on one element
  long nne;
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of surfaces
  long nsurf;
  ///  number of nodes on surface
  long nnsurf;

  ///  number of DOFs on plane element
  long ndofes;
  ///  number of DOFs on plate element
  long ndofep;
  ///  number of DOFs on the element
  long ndofe;

  ///  total number of integration points on element
  long tnip,tnips,tnipp;
  ///  order of integration of stiffness matrix
  long **intordsm;
  ///  total number of components of stress and strain tensors
  long tncomp,tncomps,tncompp;
  ///  number of components of stress and strain tensors plane
  long *ncomps;
  ///  number of components of stress and strain tensors plate
  long *ncompp;
  ///  number of approximated functions on the element
  long napfun,napfuns,napfunp;
  //  array of numbers of DOFs
  //long *dofe;
  ///  number of integration points for stiffness matrix
  long **nip;
  ///  number of blocks
  long nb,nbs,nbp;
  ///  stress/strain state
  strastrestate ssst;
  ///  ordering of unknowns
  long **ordering;

};

#endif
