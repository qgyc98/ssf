#ifndef DAMPLAST_H
#define DAMPLAST_H

#include <stdio.h>
#include "alias.h"
struct matrix;
struct vector;

/**
  This class defines artificial material model which combines
  a damage model with a plastic model
  
  
  15.1.2004
*/
class damplast
{
 public:
  damplast (void);
  ~damplast (void);
  void matstiff (matrix &d, long ipp, long ido);
  void nlstresses (long ipp);
  void updateval (long ipp);
  double give_actual_ft(long ipp, long im, long ido);
};

#endif
