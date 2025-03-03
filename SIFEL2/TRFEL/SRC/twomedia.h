#ifndef TWOMEDIA_H
#define TWOMEDIA_H

#include "genfile.h"

class med2
{
 public:
  med2();    //constructor
  ~med2();   //destructor

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);

  void rhs_volume (matrix &d,long ri, long ci,long ipp);
  void rhs_volume2 (double &c,long ri, long ci,long ipp);
  
  double transmission_transcoeff (double trc,long ri,long ci,long nid,long bc,long ipp);
  double transmission_nodval (double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux (double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  void print_othervalues (FILE *out);

  double compute_othervalues (long compother,long ipp,double *r);
  void print_othervaluesnames (FILE *out,long ipp,long compother);

 private: 

  double scale_w,scale_t, scale_cf;

};  

#endif
