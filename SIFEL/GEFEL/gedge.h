#ifndef GEDGE_H
#define GEDGE_H

#include <stdio.h>
#include "gnode.h"
#include "vector.h"

/**
   class gedge
   
   this class defines general edge, which can be
   used in connection with general finite elements
   only topological informations are collected
   
   the class is motivated by solution of hemivariational
   inequalities,
   in 2D problems, there is a contact between two parts of a body
   matching mesh is assumed and therefore two elements share
   one edge, the edge has to know four nodes because it will work
   with differences of displacements
   
   example:

              BODY 1
          3        4        5
   ---2-------3--------4--------5----
   ---6-------7--------8--------9----
          2        3        4
	      BODY 2

   edge 4 from the body 1 and edge 3 from body 2 are identical,
   they are stored in one general edge which has to know nodes 3,7,4,8
   nodes in body 1 are ordered first, nodes in body 2 are ordered after
   ordering of body 1
   number of nodes on edge nn=4
   node multiplicity nm=2
   list of nodes on edge  nlist[0]=3, nlist[1]=7, nlist[2]=4, nlist[3]=8


   JK, 9.7.2007
*/

class gedge
{
 public:
  
  gedge ();
  ~gedge ();
  
  void direction_vector (gnode *gnodes);
  void normal_vector (gnode *gnodes);
  
  void print (FILE *out);

  void give_dirvect (double *v);
  void give_norvect (double *v);

  void check_normal (vector &x,vector &y,ivector &nod);

  void alloc_cn (long nccnfn,long nccnln);
  
  void give_first_node_numbers (long *fnn);
  void give_last_node_numbers (long *lnn);

  void give_mult_code_numbers (long fln,long *mcn);

  
  ///  the number of nodes on edge
  long nn;
  ///  first node (number of first node)
  long fn;
  ///  last node (number of last node)
  long ln;
  ///  the number of DOFs defined on first nodes
  long ndofnf;
  ///  the number of DOFs defined on last nodes
  long ndofnl;

  ///  previous edge
  long prev;
  ///  next edge
  long next;

  ///  number of reference element
  long re;
  ///  numbers of adjacent elements
  long *adjel;

  ///  node multiplicity
  long nm;
  ///  list of nodes on edge
  long *nlist;
  ///  array of code numbers for the first nodes (code numbers of Lagrange multipliers)
  long *cnfn;
  ///  array of code numbers for the last nodes (code numbers of Lagrange multipliers)
  long *cnln;

  ///  length of the edge
  double l;
  ///  direction %vector
  double *dv;
  ///  normal %vector
  double *nv;
    
  ///  threshold
  double threshold;

};

#endif
