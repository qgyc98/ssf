#ifndef SEJTKRMATC_H
#define SEJTKRMATC_H

#include "genfile.h"

class sejtkrmatc
{
 public:
  sejtkrmatc();    //constructor
  ~sejtkrmatc();   //destructor

  void read(XFILE *in);

  double get_kuw(double pw);
  double get_kwu(double pw);
  double get_kww(double pw);

  double get_capuw(double pw);
  double get_capwu(double pw);
  double get_capww(double pw);

  double get_fw1(double pw);
  double get_fu1(double pw);

 private:
  
  double emod,nu,alpha,ks,kt,kw,phi0, kintr,rhos,rhow,muw;
};  

#endif
