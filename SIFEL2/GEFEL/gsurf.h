#ifndef GSURF_H
#define GSURF_H

#include <iostream>
#include <stdio.h>
#include "gnode.h"
#include "vector.h"

/**
   class gsurf
   
   this class defines general surface, which can be
   used in connection with general finite elements
   only topological informations are collected
   
   the class is motivated by solution of hemivariational
   inequalities,
   in 3D problems, there is a contact between two parts of a body
   matching mesh is assumed and therefore two elements share
   
   JK, 16.3.2009
*/

class gsurf
{
 public:
  
  gsurf ();
  ~gsurf ();
  
  void normal_vector (gnode *gnodes);
  
  void print (FILE *out);

  void give_norvect (double *v);

  void check_normal (vector &x,vector &y,ivector &nod);

  void alloc_cn (long nccnfn,long nccnln);
  
  void give_first_node_numbers (long *fnn);
  void give_last_node_numbers (long *lnn);

  void give_mult_code_numbers (long fln,long *mcn);

   ///  vypis 
    //friend std::ostream& operator<<(std::ostream &os, gsurface &gs);	

  

  ///  number of nodes on surface
  long nn;

  ///  number of DOFs defined on nodes
  ///  the same number of DOFs has to be defined on both side of the surface
  long *ndofn;

  ///  number of adjacent elements
  long nae;

  ///  number of reference element
  long re;
  ///  numbers of adjacent elements
  long *adjel;

  ///  node multiplicity
  long nm;
  ///  list of nodes on surface (on the first side and second side)
  long *nlistfn;
  long *nlistln;
  ///  array of code numbers of Lagrange Multipliers
  ///  k=cnmult[i][j] - the j-th multiplier on the i-th node has number k
  long **cnmult;

  ///  normal %vector
  double *nv;
    
  ///  threshold
  double threshold;

};

#endif
