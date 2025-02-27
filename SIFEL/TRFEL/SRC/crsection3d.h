#ifndef CRSECTION3D_H
#define CRSECTION3D_H

#include <stdio.h>
#include "iotools.h"
struct vector;
struct atsel;

/**
   class crsection3d defines cross section for threedimensional transport problems
*/

class crsection3d
{
 public:
  crsection3d (void);
  ~crsection3d (void);
  void read (XFILE *in);
  void print (FILE *out);
  void changeparam (atsel &atcs,vector &val);
  
  ///  density of material (necessary for nonstationary problems)
  double rho;

};

#endif
