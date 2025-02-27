#ifndef ELASTORTOMAT_H
#define ELASTORTOMAT_H

#include <stdio.h>
#include "iotools.h"
#include "alias.h"
#include "gfmatrix.h"
#include "matrix.h"



/**
  The class defines elastic orthotropic material given by the 
  stiffness %matrix de in the local coordinate system of the 
  material. There must be also defined transformation %matrix t
  whose components represent general functions (gfunct). The 
  components of the %matrix t must be evaluated with respect 
  to space coordinates of the given integration point.

  Created by Tomas Koudelka, 12.2014
*/
class elastortomat
{
 public:
  elastortomat (void);
  ~elastortomat (void);
  /// reads material from the opened text file
  void read (XFILE *in);
  /// printf material to the opened text file
  void print (FILE *out);

  /// computes the stiffness %matrix of the elastic orthotropic material in the global coordinate system
  void matstiff (matrix &d, long ipp);
  /// computes the stiffness %matrix of the elastic orthotropic material in the global coordinate system
  void matstiff (matrix &d, long ipp, strastrestate ssst);
  /// returns the stiffness %matrix of the elastic orthotropic material in the local coordinate system
  void loc_matstiff (matrix &d, strastrestate ssst);
  /// returns transformation %matrix from the local material coordinate system to the global one
  void give_transf_mat (matrix &tmat, long ipp);
  /// calculates stresses
  void nlstresses (long ipp);

  /// stiffness %matrix in the local coordinate system of the material
  matrix de;

  /** Transformation %matrix whose components are evaluated 
      with respect to coordinates of the given integration point,
      It must be given in the form x_g = T x_l where 
      x_g and x_l are %vectors in the global and local coordinate system. */
  gfmatrix t;
};

#endif
