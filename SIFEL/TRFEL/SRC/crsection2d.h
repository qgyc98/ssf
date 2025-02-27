#ifndef CRSECTION2D_H
#define CRSECTION2D_H

#include <stdio.h>
#include "iotools.h"
struct vector;
struct atsel;

/**
   class crsection2d defines cross section for twodimensional transport problems
*/

class crsection2d
{
 public:
  crsection2d (void);
  ~crsection2d (void);
  void read (XFILE *in);
  void print (FILE *out);
  void changeparam (atsel &atcs,vector &val);
  
  //  necessary components for every problem type
  ///  thickness of the problem
  double t;
  
  ///  density of material (necessary for nonstationary problems)
  double rho;

};

#endif
