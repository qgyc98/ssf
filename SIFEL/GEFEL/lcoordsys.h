#ifndef LCOORDSYS_H
#define LCOORDSYS_H

#include "galias.h"
#include <stdio.h>

class gfmatrix;
class matrix;
class gnode;
class vector;
struct XFILE;

class lcoordsys
{
 public:
  /// default constructor
  lcoordsys (void);

  /// constructor with problem dimension specified
  lcoordsys (long pdim);

  /// destructor
  ~lcoordsys (void);

  /// function reads definition of a local coordinate system from the opened text file
  void read (XFILE *in);

  /// function prints definition of a local coordinate system from the opened text file
  void print (FILE *out);

  /// function assembles transformation matrix
  void give_transfmat(matrix &t, vector &pc, double time);

  /// function assembles transformation matrix for a cylindrical domain
  void give_cyl_transfmat(matrix &t, vector &pc);

  static const double zero;
  static const char *namepar[];

  ///  dimensions of the transformation matrix
  long dim;

  /// transformation matrix type
  tmat_type tt;

  /// base vector in the direction of cylinder axis (for a cylindrical domains, tm = gen_tr)
  vector *caxis;
  /// a point on the cylinder axis
  vector *paxis;

  /// transformation matrix from the local to global coordinate system for general case (tm = gen_tr)
  matrix *tmat; 

  /// transformation matrix from the local to global coordinate system for general case (tm = gfgen_tr)
  gfmatrix *gftmat; 
};

#endif
