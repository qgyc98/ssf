#ifndef FIXORTODAM_H
#define FIXORTODAM_H

#include <stdio.h>
#include "gfunct.h"
#include "iotools.h"
#include "vector.h"

struct vector;
struct matrix;

class fixortodam
{
 public:
  fixortodam (void);
  ~fixortodam (void);
  void read (XFILE *in);
  void print (FILE *out);
  void initval(long ipp, long ido);
  void matstiff (long ipp,long ido,matrix &d);
  void tmatstiff (long ipp,long ido,matrix &d);

  void compute_eqdispl(long ipp, matrix &epst, vector &xeq);
  void secstiffmat(long ipp, vector &omega, matrix &d);
  void compute_dam(long ipp, vector &xeq, vector &xeq0, vector &omegao, vector &omega);

  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp,long im,long ido);


  /// tensile stregth in longitudinal direction
  double x1t;
  
  /// tensile stregth in the first transversal direction
  double x2t;
  
  /// tensile stregth in the second transversal direction
  double x3t;
	
  double x1c;
  double x2c;
  double x3c;

  /** equivalent displacement in fully damaged state in longitudinal, the first transversal and the second 
      transversal directions (tension)*/
  vector xeqtf;
  vector xeqcf;

  static double L_MAX;
  static double L_MIN;
};

#endif
