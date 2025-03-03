#ifndef DISCISOTRMAT_H
#define DISCISOTRMAT_H

#include <stdio.h>
#include "genfile.h"
struct vector;
struct atsel;

class discisotrmat
{
 public:
  discisotrmat (void);
  ~discisotrmat (void);
 
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);
  void print (FILE *out);
  double get_k();
  double get_c();
  void changeparam (atsel &atm,vector &val);

  double correct_val (double *in,double *iin);

  double compute_rel (double in);
  double compute_abs (double in);


  ///  coefficient of conductivity
  double k;
  ///  coefficient of capacity
  double c;
  ///  jump in quantity at discontinuity
  double ju;

};

#endif
