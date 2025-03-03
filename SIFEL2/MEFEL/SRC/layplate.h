#ifndef LAYPLATE_H
#define LAYPLATE_H

#include "alias.h"
#include "iotools.h"
#include "strretalg.h"
#include "intpoints.h"
#include "matrix.h"
struct matrix;
struct vector;
struct atsel;

/**
   Class layplate defines layered model for plates.
   
   Josef Fiedler, 01/2015
*/

class layplate
{
 public:
  layplate (void);
  ~layplate (void);
  void read (XFILE *in);
  
  void   matstiff       (matrix &d, long ipp,long ido);
  void   nlstresses     (long ipp, long im, long ido);
  long   compeqother    (long ipp);
  long   compother      (long ipp);
  void   backup         (long ipp, double *&k, long j);
  void   restore_values (long ipp, double *k, long j);
  void   stress_calc    (vector &df, vector &eps, double *k, vector &intfor, double *layz, double *layth, long ipp, long ido);
  void   updateval      (long ipp, long im, long ido);
  double dstep_red      (long ipp, long im, long ido);


  /// number of layers
  long nl;
  /// number of materials on current layer
  long *nm;
  /// array of material types
  mattype **tm;
  /// array of number of appropriate material type
  long **idm;
  /// number of iteration steps for cross-section eqilibrium
  long nli;
  /// error for iteration of cross-section eqilibrium
  double err;
  /// backup of informations stored on integration point
  intpoints bcup;
  /// compliance matrix for determination of normal forces eqilibrium
  matrix c;
  /// a matrix providing relation between eps and curvatures;
  /// it is used in case of uneven distribution of stiffness along height of cross-section
  matrix cn;

};

#endif
