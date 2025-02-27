#ifndef NONLOCMICROM4_H
#define NONLOCMICROM4_H

#include "iotools.h"
#include "matrix.h"
#include "vector.h"


class nonlocmicroM4
{
 public:
  nonlocmicroM4 (void);
  ~nonlocmicroM4 (void);
  void read (XFILE *in);
  long give_num_averq (long ipp);
  void average (long ipp);

  //material parameters
  double r;
  
 protected:
};

#endif
