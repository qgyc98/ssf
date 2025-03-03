#ifndef ISOVISCOUS_H
#define ISOVISCOUS_H

#include "alias.h"
#include "iotools.h"
struct matrix;
struct vector;


/**
   material model for isotropic viscous material

*/
class isoviscous
{
 public:
  
  isoviscous (void);
  ~isoviscous (void);
  void read(XFILE *in);
  void print(FILE *out);
  
  void matdamp (matrix &d,strastrestate ssst);
  
  void matdamp_plstress (matrix &d);
  
  void matdamp_spacestr (matrix &d);

  ///  coefficient of volumetric viscosity
  double xi;
  ///  coefficient of dynamic viscosity
  double eta;
};

#endif
