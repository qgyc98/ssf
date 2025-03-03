#ifndef SOURCET_H
#define SOURCET_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "hydrationheat.h"
#include "seebeckheat.h"
//#include "cemhyd.h"

/**
   class contains model of quantity source for transport problems
   the class can describe, e.g., source of heat
   
   JK, revised 25.6.2005
*/
class sourcet
{
 public:
  sourcet ();
  ~sourcet ();
  long read(XFILE *in);
  long print(FILE *out);
  double giveval (long eid);
  long read_prop (FILE *in, long lc);
  long compare(sourcet &src);
  
  ///  type of quantity source
  sourcetype sourtype;

  ///  model of source of general quantity described by mathematical function
  gfunct *gf;
  ///  model of heat source caused by cement hydration, amount of heat is given by a special function
  hydrationheat *hydrh;
  
  ///  model of heat source caused by cement hydration
  //cemhyd *cemh;
  
  ///  heat source for seebeck effect
  seebeckheat *seebh;
  

  ///  constant value of source of general quantity for stationary problems only
  double v;
  
};

#endif
