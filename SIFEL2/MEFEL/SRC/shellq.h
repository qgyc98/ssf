#ifndef SHELLQ_H
#define SHELLQ_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   shell triangular element
   resulting matrices are composed from plane triangular element
   with rotational degrees of freedom and from constant curve triangle
   
*/
class shellq
{
 public:
  shellq (void);
  ~shellq (void);
  
  ///  function assembles element code numbers
  void codnum (long *cn,long ri);

  ///  function assembles transformation %matrix from element coordinate system to the global system
  void coord_transf_matrix (matrix &tran, vector &gx, vector &gy, vector &gz, double zero);

  ///  function transforms node coordinates from the global system to the element coordinate system
  void local_coordinates (vector &gx,vector &gy,vector &gz,vector &lx,vector &ly, double zero);
  
  ///  function assembles transformation %matrix from local node coordinate system to the global system
  void node_transf_matrix (ivector &nodes,matrix &tmat);
  
  ///  function assembles transformation %matrix from element coordinate system to the global system
  void elem_transf_matrix (matrix &tran, vector &gx, vector &gy, vector &gz, double zero);
  
  ///  function assembles stiffness %matrix
  void res_stiffness_matrix (long eid,matrix &sm);
  
  ///  function computes strains in integration points
  void res_ip_strains (long lcid,long eid);
  
  void nod_strains_ip (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid);
  void res_ip_stresses (long lcid,long eid);

  void node_forces_surf (long eid,double *nodvals,vector &nf);
  
  void res_internal_forces (long lcid,long eid,vector &ifor);
  
  /*
  void nod_strains_ip (long lcid,long eid);
  void strains (long lcid,long eid);
  void res_ip_stresses (long lcid,long eid);
  void nod_stresses_ip (long lcid,long eid);
  void stresses (long lcid,long eid);
  void inicipval(long eid, long ri, long ci, matrix &nodval, inictype *ictn);
  */

  
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
