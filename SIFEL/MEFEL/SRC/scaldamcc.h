#ifndef SCALDAMCC_H
#define SCALDAMCC_H

#include "iotools.h"
#include "alias.h"
struct matrix;
struct vector;


class scaldamcc
/**
  This class defines scalar damage with crack closure material model.
  The different type of norms for the computing parameters of
  the damage function can be used.
*/
/**
*/
{
 public:
  scaldamcc (void);
  ~scaldamcc (void);
  void   read (XFILE *in);
  void   print (FILE *out);
  void   elmatstiff (matrix &d,long ipp);
  double damfuncpar(long ipp, vector &eps);
  double damfunction(long ipp, double kappa, vector &eps, double omegao, double &dt, double &dc, double &alphat, double &alphac);
  void   matstiff(matrix &d, long ipp, long ido);
  void   nlstresses(long ipp, long im, long ido);
  void   updateval(long ipp, long im, long ido);
  double epsefunction(long ipp);
  

  // Parameters for the damage function defined by the simple exponential function (Jirasek)
  /// tensile strength
  double ft;
  /// softening slope
  double uf;

  // Parameters for the damage function defined by the exponential function (Mazars)
  /// damage evolution treshold
  double k0;
  ///
  double at;
  /// tensile strength
  double bt;
  ///
  double ac;
  /// strength at compression
  double bc;
  /// 
  double beta;
  ///
  double k;
  /// function type for evaluation damage funtion parameter
  paramf_type eqepsnt;
  damfunc_type damfunct;
};

#endif
