#ifndef EDGET_H
#define EDGET_H

#include <stdio.h>
#include "aliast.h"
#include "iotools.h"

/**
   class edget defines general edge for transport problems

*/
class edget
{
 public:
  edget (void);
  ~edget (void);
  void init (long edid);
  void compute_jump (long edid);
  void compute_node_jump (long edid, long nodid,long lcid,double ncf);

  void init_edval (void);
  void store_edval (double *v);
  void give_edval (double *v);

  ///  material types
  mattypet *tm1,*tm2;
  ///  material id
  long *idm1,*idm2;
  
  ///  jumps in nodal values at first nodes on edge
  double *jumpfn;
  ///  jumps in nodal values at last nodes on edge
  double *jumpln;
  
  ///  the number of components in the array edval
  long ncedval;
  ///  edge values
  double *edval;

};

#endif
