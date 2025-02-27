#ifndef ARGYRISPLATE_H
#define ARGYRISPLATE_H

#include "alias.h"
struct matrix;
struct vector;
struct ivector;
struct atsel;



/**
   class argyrisplate defines triangular element with
   21 approximation functions (DOFs)
   
   basic data
   nne = 6 - number nodes on element
   ndofe = 21 - number of DOFs on element
   ncomp = 3 - number of strain (stress) tensor components
   napfun = 1 - number of approximated functions
   intordmm = 2 - order of numerical integration of the mass %matrix (2x2)
   
   7.7.2012, JK
*/
class argyrisplate
{
 public:
  argyrisplate (void);
  ~argyrisplate (void);

  
  void fx (double x,double y,vector &shapef);
  void dfdx (double x,double y,vector &shapef);
  void dfdy (double x,double y,vector &shapef);
  void dfdxdx (double x,double y,vector &shapef);
  void dfdxdy (double x,double y,vector &shapef);
  void dfdydy (double x,double y,vector &shapef);
  void shapefunctions (long eid);
  void geom_matrix (matrix &gm,double x,double y,long eid);
  
  void stiffness_matrix (long eid,long ri,long ci,matrix &sm,vector &x,vector &y);
  void res_stiffness_matrix (long eid,matrix &sm);
  
  ///  number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  total number of components of stress and strain tensors
  long tncomp;
  ///  number of components for graphic purposes
  long gncomp;
  ///  total number of integration points on element
  long tnip;
  ///  array containing numbers of components of stress and strain tensors
  long *ncomp;
  ///  array containing cumulative numbers of components of stress and strain tensors
  long *cncomp;
  ///  number of approximated functions on the element
  long napfun;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  array containing orders of numerical integration of stiffness matrix
  long **intordsm;
  ///  order of integration of mass matrix
  long intordmm;
  ///  order of integration on edges
  long intordb;
  ///  array of numbers of integration points in sets
  long **nip;
  ///  number of blocks
  long nb;
  ///  stress/strain state
  strastrestate ssst;

};

#endif
