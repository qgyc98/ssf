#ifndef DAMPERMEABILITY_H
#define DAMPERMEABILITY_H

#include <stdio.h>
#include "genfile.h"
struct vector;

class dampermeability
{
 public:
  dampermeability (void);
  ~dampermeability (void);
 
  void read (XFILE *in);
  void print (FILE *out);

  void matcond (matrix &d,long ipp);
  void give_reqntq(long *antq);
  
  ///  kinematic viscosity of water (nu)
  ///  nu = eta/rho
  ///  eta - dynamic viscosity
  ///  rho - density of water
  gfunct *kinvis;
};

#endif
