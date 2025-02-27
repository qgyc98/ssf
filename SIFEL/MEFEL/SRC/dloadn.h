#ifndef DLOADN_H
#define DLOADN_H

#include "xfile.h"
#include <stdio.h>

class gfunct;

/**
  This class defines time dependent nodal load or the load which is the function of time or displacement.
  The load may be described with constant values, expression defined by the string or set of strings and
  finally it may be described with table.
*/
class dloadn
{
 public:
  dloadn();
  ~dloadn();
  long read(XFILE *in);
  long print(FILE *out);
  double getval(double t, long id);

  ///  number of loaded node
  long idn;
  /// array with time functions for particular load components
  gfunct *gf;
};

#endif
