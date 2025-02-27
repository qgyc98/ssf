#ifndef ELASTGMAT3D_H
#define ELASTGMAT3D_H

#include<stdio.h>
#include "iotools.h"

struct matrix;

/**
  Class elastgmat3d defines elastic fully anisotropic material model.
   
  Created by JK,
*/

class elastgmat3d
{
 public:
  elastgmat3d (void);
  ~elastgmat3d (void);
  void read (XFILE *in);
  void matstiff (matrix &d);
  void elmatstiff (matrix &d);
  void nlstresses(long ipp);

  ///  stiffness coefficients
  double d11,d12,d13,d14,d15,d16;
  double d22,d23,d24,d25,d26;
  double d33,d34,d35,d36;
  double d44,d45,d46;
  double d55,d56;
  double d66;
};

#endif
