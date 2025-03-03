#ifndef ITABLEFUNCT_H
#define ITABLEFUNCT_H

#include <stdio.h>
#include "iotools.h"


/**
 This file declares the class of tablefunct, which implements various types of interpolation.
 The values for interpolation are given by the table.
*/
class itablefunct
{
 public:
  itablefunct();
  ~itablefunct();
  
  void read(XFILE *in);
  void print(FILE *out);
  void read_prop(FILE *in);
  long getval(double temp);

  void copy(itablefunct &tf);
  long compare(itablefunct &tf);


  ///number of function values
  long asize;
  ///array of time values
  double *x;
  ///array of function values
  long *y;
  
};
#endif
