#ifndef SIMVISCOUS_H
#define SIMVISCOUS_H

//#include <stdio.h>
#include "iotools.h"


/**
   structure of the array eqother:
   total stresses, viscoplastic (irreversible) strain increments,
   previous total strains,

*/
class simviscous
{
 public:
  
  simviscous (void);
  ~simviscous (void);
  void read(XFILE *in);
  void print(FILE *out);
  double gfun (double f);
  double dergfun (double f);
  
  ///  viscosity
  double eta;
};

#endif
