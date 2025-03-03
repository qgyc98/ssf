#ifndef ENDNODE_H
#define ENDNODE_H

#include <stdio.h>
#include "gnode.h"
#include "vector.h"

/**
   class endnode
   
   this class defines end nodes, which can be
   used in connection with general finite elements
   only topological informations are collected
   
   the class is motivated by solution of hemivariational
   inequalities,
   in 1D problems, there is a contact between two parts of a body
   matching mesh is assumed and therefore two elements share
   one end node, the end node has to know two nodes because it will work
   with differences of displacements
   
   example:
   
   BODY 1  ---3----4 7----8----12---- BODY 2
                    3
   end node 4 from the body 1 and end node 7 from body 2 are identical,
   they are stored in one end node which has to know nodes 4, 7
   node in body 1 is ordered first, node in body 2 is ordered after
   ordering of body 1
   number of nodes collected in end node is 2

   JK, 12.10.2008
*/

class endnode
{
 public:
  
  endnode ();
  ~endnode ();
  
  void print (FILE *out);

  void alloc_cnm (long nccn);
  
  void give_node_numbers (long *nid);
  void give_mult_code_numbers (long *mcn);

  
  ///  number of nodes in end point
  long nn;
  ///  first node (number of first node)
  long fn;
  ///  last node (number of last node)
  long ln;
  ///  number of DOFs defined on first nodes
  long ndofn;

  ///  number of reference element
  long re;
  ///  numbers of adjacent elements
  long *adjel;

  ///  node multiplicity
  long nm;
  ///  array of code numbers for Lagrange multipliers
  long *cnm;

  ///  threshold
  double threshold;

};

#endif
