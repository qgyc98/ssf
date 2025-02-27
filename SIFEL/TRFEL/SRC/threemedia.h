#ifndef THREEMEDIA_H
#define THREEMEDIA_H

#include "genfile.h"

class med3
{
 public:
  med3();    //constructor
  ~med3();   //destructor

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  void rhs_volume (matrix &d,long ri, long ci,long ipp);
  void rhs_volume2 (double &c,long ri, long ci,long ipp);
  
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  double compute_othervalues (long compother,long ipp,double *r);
  void print_othervaluesnames (FILE *out,long ipp,long compother);

  
};  

#endif
