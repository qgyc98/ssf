#ifndef GEN2DELEM_H
#define GEN2DELEM_H

#include "genfile.h"

/**
   The class gen2delement defines a fictious general 2D element
   for transport problems which is used as an interface for material
   models in HERMES adaptivity.
   
   Created by Tomas Koudelka 03.2011
*/
class gen2delem
{
 public:
  gen2delem (void);
  ~gen2delem (void);

  ///  number of transported matter
  long ntm;
  ///  total number of DOFs on the element
  long ndofe;
  ///  number of nodes on one element
  long nne;
  ///  number of edges
  long ned;
  ///  number of nodes on one edge
  long nned;
  ///  problem dimension
  long ncomp;
  ///  number of integration points
  long **nip;
};

#endif
