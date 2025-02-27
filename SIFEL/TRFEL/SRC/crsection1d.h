#ifndef CRSECTION1D_H
#define CRSECTION1D_H

#include <stdio.h>
#include "iotools.h"
struct vector;
struct atsel;

/**
   class crsection1d for onedimensional transport problems
*/

class crsection1d
{
 public:
  crsection1d (void);
  ~crsection1d (void);
  void read (XFILE *in);
  void print (FILE *out);
  void changeparam (atsel &atcs,vector &val);
  
  //  necessary components for every problem type
  ///  area of cross section
  double a;
  
  ///  density of material (necessary for nonstationary problems)
  double rho;

};

#endif
