#ifndef ONEMEDIUM_H
#define ONEMEDIUM_H

#include "genfile.h"

class med1
{
 public:
  med1();    //constructor
  ~med1();   //destructor

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  void matreact (double &r,long ri,long ci,long ipp);

  double cond_k (long ipp);
  double cap_c (long ipp);

  void rhs_volume (matrix &d,long ri, long ci,long ipp);
  void rhs_volume2 (double &c,long ri, long ci,long ipp);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  double compute_othervalues (long compother,long ipp,double *r);
  void print_othervaluesnames (FILE *out,long ipp,long compother);
  void eigstrains (long ipp);

 private: 

  double scale;
};  

#endif
