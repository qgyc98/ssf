#ifndef CREEPDAM_H
#define CREEPDAM_H

#include <stdio.h>
#include "alias.h"
struct matrix;
struct vectror;

/**
  This class defines artificial material model which combines
  a creep model with damage model
  
  Order of internal variables in the other array :
  -----------------------------------------------
  The first components of eqother array belong to creep model,
  followed by the damage eqother components and
  optionally followed by the thermisomatttime components

  Exact order of the components depends on the used material 
  models.  
  
  15.1.2004
*/
class creepdam
{
 public:
  creepdam (void);
  ~creepdam (void);
  void matstiff (matrix &d,long ipp,long im,long ido);
  void nlstressesincr (long ipp, long im, long ido);
  void nlstresses (long ipp, long im, long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  void initvalues (long lcid, long ipp, long im, long ido, bool rinit);
  double give_actual_ft (long ipp, long im, long ido);
  double give_actual_fc (long ipp, long im, long ido);  
};

#endif
