#ifndef BNODVALT_H
#define BNODVALT_H

#include <stdio.h>
#include "iotools.h"
#include "vector.h"
#include "gfunct.h"

/**
   this class serves for data storage
   the class stores nodal values defined on boundary
   stored values may represent prescribed values, transmission
   coefficients, radiation coeffcients, prescribed fluxes, etc.
   
   JK, 24.11.2008
*/
class bnodvalt
{
 public:
  bnodvalt ();
  ~bnodvalt ();
  void read (XFILE *in);
  void print (FILE *out);
  void give_val (double t,vector &nv);

  ///  number of stored components
  long nsc;
  ///  nodal values defined on boundaries
  gfunct *nodval;
  
};

#endif
