#ifndef LEMAITRE_H
#define LEMAITRE_H

#include "iotools.h"
struct matrix;


class lemaitre
{
 public:
  lemaitre (void);
  ~lemaitre (void);

  void read (XFILE *in);
  void print (FILE *out);
  double gfun (double f,double cs);
  //void nlstresses (long ipp,long ido);
  //void matstiff (matrix &d,long ipp,long ido);
  
  //  viscosity coefficient
  double eta;
  //  hardening parameter
  double m;
  //  stress exponent
  double n;

};

#endif
