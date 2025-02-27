#ifndef GENERAL3MAT_H
#define GENERAL3MAT_H

#include <stdio.h>
#include "genfile.h"

class general3mat
{
 public:
  general3mat (void);    //constructor
  ~general3mat (void);   //destructor
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void read (FILE *in);

  double k11 (double x1,double x2,double x3);
  double k12 (double x1,double x2,double x3);
  double k13 (double x1,double x2,double x3);

  double k21 (double x1,double x2,double x3);
  double k22 (double x1,double x2,double x3);
  double k23 (double x1,double x2,double x3);

  double k31 (double x1,double x2,double x3);
  double k32 (double x1,double x2,double x3);
  double k33 (double x1,double x2,double x3);

  double c11 (double x1,double x2,double x3);
  double c12 (double x1,double x2,double x3);
  double c13 (double x1,double x2,double x3);

  double c21 (double x1,double x2,double x3);
  double c22 (double x1,double x2,double x3);
  double c23 (double x1,double x2,double x3);

  double c31 (double x1,double x2,double x3);
  double c32 (double x1,double x2,double x3);
  double c33 (double x1,double x2,double x3);

  void auxiliarydata (double x1,double x2,double x3);
  
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_nodval_11(double bv,double x1,double x2,double x3,long bc,long ipp);
  double get_transmission_transcoeff_11(double x1,double x2,double x3,long bc,long ipp);
  double get_transmission_flux_11(double bv,double x1,double x2,double x3,long bc,long ipp);

  
  double get_othervalue(long compother,long ipp, double x1,double x2,double x3);
  void print_othervalue_name(FILE *out,long compother);

 private:

};

#endif

