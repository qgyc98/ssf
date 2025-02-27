#ifndef ENDNODET_H
#define ENDNODET_H

#include <stdio.h>
#include "aliast.h"
#include "iotools.h"

/**
   class endnodet defines general end node for transport problems
*/

class endnodet
{
 public:
  endnodet (void);
  ~endnodet (void);
  void init (long edid);
  void compute_jump (long edid);
  
  ///  material types
  mattypet *tm1,*tm2;
  ///  material id
  long *idm1,*idm2;
  
  ///  jumps in nodal values
  double *jump;
};

#endif
