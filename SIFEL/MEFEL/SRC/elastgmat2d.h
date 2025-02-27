#ifndef ELASTGMAT2D_H
#define ELASTGMAT2D_H

#include<stdio.h>
#include "alias.h"
#include "iotools.h"

struct matrix;

/**
  Class elastgmat2d defines elastic fully anisotropic material model.
   
  Created by JK,
*/

class elastgmat2d
{
 public:
  elastgmat2d (void);
  ~elastgmat2d (void);
  void read (XFILE *in);
  void matstiff (matrix &d,strastrestate ssst);
  void elmatstiff (matrix &d,strastrestate ssst);
  void nlstresses(long ipp);
  void matcompl (matrix &c,strastrestate ssst);

  ///  stiffness coefficients
  double d11,d12,d13;
  double d22,d23;
  double d33;
};

#endif
