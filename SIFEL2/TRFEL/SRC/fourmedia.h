#ifndef FOURMEDIA_H
#define FOURMEDIA_H

#include "genfile.h"

class med4
{
 public:
  med4();    //constructor
  ~med4();   //destructor
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  double transmission_transcoeff (double trc,long ri,long ci,long nid,long bc,long ipp);
  double transmission_nodval (double nodval,double trc2,long ri,long ci,long nid,long bc,long ipp);
  double transmission_flux (double nodval,double trc2,long ri,long ci,long nid,long bc,long ipp);
  
  double compute_othervalues (long compother,long ipp,double *r);
  void print_othervaluesnames (FILE *out,long ipp,long compother);
  
 private: 
  double scale_w,scale_t, scale_cf;
  
};  

#endif
