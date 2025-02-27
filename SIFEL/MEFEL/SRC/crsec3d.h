#ifndef CRSEC3D_H
#define CRSEC3D_H

#include "xfile.h"
#include <stdio.h>
struct vector;
struct atsel;



/**
  Class crsec3d defines cross section for 3D problems.
   
  Creatde by JK,
*/
class crsec3d
{
 public:
  crsec3d (void);
  ~crsec3d (void);
  void read (XFILE *in);
  void print (FILE *out);
  void changeparam (atsel &atcs,vector &val);
  
  //  optional components
  ///  density of material (necessary for dynamics)
  double rho;
  
};

#endif
