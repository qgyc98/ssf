#ifndef SEJTKRMATC_H
#define SEJTKRMATC_H

#include "genfile.h"

class sejtkrmatc
{
 public:
  sejtkrmatc();    //constructor
  ~sejtkrmatc();   //destructor

  void read(XFILE *in);

  void rhs_volume (matrix &d,long ri,long ci,long ipp);
  void rhs1d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d1 (matrix &d,long ri,long ci,long ipp);

  double get_kuw(double pw);
  double get_kwu(double pw);
  double get_kww(double pw);

  double get_capuw(double pw);
  double get_capwu(double pw);
  double get_capww(double pw);

  double get_fw1(double pw);
  double get_fu1(double pw);

  double get_othervalue(long compother,double pw, long ipp);
  void print_othervalue_name(FILE *out,long compother);
  void values_correction (vector &nv, long ipp);
  void water_pressure_check(double &pw,long ipp);

 private:
  
  double emod,nu,alpha,kz,ks,kk,phi0,k,rhok,g,lambda,c;
};  

#endif
