#ifndef RELAXEUROC_H
#define RELAXEUROC_H

#include "alias.h"
#include "iotools.h"
#include "vector.h"

/**
   The class relaxeuroc contains material model of stress relaxation defined in Eurocode 2
   
   Created by JK, 11. 6. 2013
*/
class relaxeuroc
{
 public:
  relaxeuroc (void);
  ~relaxeuroc (void);
  void read (XFILE *in);
  void print (FILE *out);

  double stress_decrement (void);
  void stress (vector &sig,vector &eps,strastrestate ssst);
  
  ///  initial prestress \sigma_{pm0}
  double siginit;
  ///  characteristic strength f_{pk}
  double fpk;
  ///  model coefficient A (determined by regression)
  double a;
  ///  model coefficient B (determined by regression)
  double b;
  ///  Young modulus of elasticity
  double e;
    
};

#endif
