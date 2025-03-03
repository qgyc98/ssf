#ifndef QUADWEDGE_H
#define QUADWEDGE_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;

/**
   class quadwedge defines wedge elements with quadratic approximation functions
   
   JK
*/

class quadwedge
{
 public:
  quadwedge (void);
  ~quadwedge (void);

  double approx (double xi,double eta,double zeta,vector &nodval);
  void bf_matrix (matrix &n,double xi,double eta,double zeta);
  void geom_matrix (matrix &gm,vector &x,vector &y,vector &z,double xi,double eta,double zeta,double &jac);
  void transf_matrix (ivector &nodes,matrix &tmat);
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm);
  void res_stiffness_matrix (long eid,matrix &sm);

  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of the strain and stress tensors
  long tncomp;
  ///  total number of integration points on element
  long tnip;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges on one element
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  number of surfaces
  long nsurf;
  ///  number of nodes on one surface
  long nnsurf;
  ///  array of orders of integration of stiffness matrix
  long **intordsmt;
  long **intordsmz;
  ///  order of integration of mass matrix
  long intordmm;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  array of numbers of components of blocks
  long *ncomp;
  ///  cumulative array of numbers of components of blocks
  long *cncomp;
  ///  stress/strain state
  strastrestate ssst;

};

#endif
