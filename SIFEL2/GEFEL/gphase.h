#ifndef GPHASE_H
#define GPHASE_H

#include <stdio.h>

/**
   class gphase
   
   this class contains general topological informations about phase
   
   TKr
*/

class gphase
{
 public:
  gphase (void);
  ~gphase (void);
  void read (FILE *in);

  ///  phase fraction
  double f;
};

#endif
