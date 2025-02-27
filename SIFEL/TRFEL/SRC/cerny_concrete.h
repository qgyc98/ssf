#ifndef CERNYMAT_H
#define CERNYMAT_H

#include <stdio.h>
#include "genfile.h"

class cernymat
{
 public:
  cernymat (void);
  ~cernymat (void);

 
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  

  void read (XFILE *in);
  void print (FILE *out);
  double get_k(double t);
  double get_dk_dt(double t);
  double get_c(double t);
  double get_dc_dt(double t);

  //density
  double rho;
  //  thermal conductivity
  double k;
  //  cpecific heat
  double c;

};

#endif
