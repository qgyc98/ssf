#ifndef VISPLAST_H
#define VISPLAST_H

#include <stdio.h>
#include "alias.h"
struct matrix;
struct vector;

/**
   artificial material model for combination of viscous and plastic material models
   viscous material is the first one while plastic material is the second one
   
   structure of eqother array of a visco-plastic material:
   increments of stress components, previous total strains, increments of irreversible strains,
   consistency parameter, hardening parameters
   
   direct access from viscous material: increments of stress components, previous total strains, previous total stress
   direct access from plastic material: total irreversible strains, consistency parameter,
                                        hardening parameters
   
   25.6.2004,
*/
class visplast
{
 public:
  visplast (void);
  ~visplast (void);
  void matstiff (matrix &d,long ipp,long im,long ido);
  
  void givestressincr (long ipp,long ido,long fi,vector &sig);
  void storestressincr (long ipp,long ido,vector &sig);

  void nlstressesincr (long ipp,long im,long ido);
  void nlstresses (long ipp,long ido);

  void updateval (long ipp, long im, long ido);

};

#endif
