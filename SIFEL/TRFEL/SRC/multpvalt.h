#ifndef MULTPVALT_H
#define MULTPVALT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"

class multpvalt
{
 public:
  multpvalt ();
  ~multpvalt ();
  long read(FILE *in);
  double getval(void);
  long read_prop(FILE *in, long lc);
  
  // prescribed value
  double v;

};

#endif
