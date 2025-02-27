#ifndef CRSECNOD_H
#define CRSECNOD_H

#include "iotools.h"

/**
   class crsecnod defines cross section for layered nodes used in layered problems
   
   JK
*/

class crsecnod
{
 public:
  crsecnod (void);
  ~crsecnod (void);
  void read (XFILE *in);
  void print (FILE *out);
  
  ///  height of the cross section
  ///  it is used in analysis of layered beams
  double t;
  
  ///  concentrated weight at node
  double m;

};

#endif
