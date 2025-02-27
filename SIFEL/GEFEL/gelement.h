#ifndef GELEMENT_H
#define GELEMENT_H

#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "iotools.h"
#include "galias.h"
class gnode;


/**
  Class gelement:
   
  This class defines general finite element,
  only topological informations are collected.
   
  Created by JK,
*/
class gelement
{
 public:
  gelement (void);
  ~gelement (void);
  void read (XFILE *in, long m, long n, gelemtype et, long maxnn);
  void read_gf (XFILE *in, long m, long n, gelemtype et, long maxnn);
  void print (FILE *out, long m, long n, gelemtype et);
  void print_gf (FILE *out, long m, long n);
  void initiate (long *icn, long m);
  long give_nne ();
  long give_nmne ();
  long give_ndofe ();
  long give_nmult ();
  void give_nodes (ivector &nod);
  void give_master_nodes (ivector &nod);
  long give_cne ();
  void centroid (long dim, gnode *gnodes, double *coord);
  
  ///  number of nodes on element
  long nne;
  ///  number of master nodes on element
  long nmne;
  ///  number of degrees of freedom on element
  long ndofe;
  ///  number of additional degrees of freedom on element
  ///  especially number of Lagrange multipliers in layered problems
  long nmult;
  ///  array containing node numbers
  long *nodes;
  ///  array containing master node numbers
  ///  this array is used only when the hanging nodes are used
  long *master_nodes;
  ///  code numbers on element (cne=1)
  long cne;
  ///  code numbers
  long *cn;
  ///  auxiliary information
  long auxinf;
  /// global number of edges
    //long *edgn;
  ///  general function
  ///  function is used for problems with growing number of elements in problem
  long tgf;
  ///  type of finite element with respect to geometry
  ///  get is assigned in the subroutine element::read() or elementt::read ()
  gelemtype get;

};

#endif
