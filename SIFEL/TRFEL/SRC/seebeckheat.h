#ifndef SEEBECKHEAT_H
#define SEEBECKHEAT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"

/**
   class seebeckheat sets the heat conduction effect from electrical filed temporarily as a source

   TKr+JM, 15.6.2016

*/
class seebeckheat
{
 public:
  seebeckheat (void);
  ~seebeckheat (void);
  void read (XFILE *in);
  void print (FILE *out);
  double give_value ();
  
  double as,sigma;
};

#endif
