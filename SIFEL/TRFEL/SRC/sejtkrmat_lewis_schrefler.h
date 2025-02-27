#ifndef SEJTKRMAT_H
#define SEJTKRMAT_H

#include "genfile.h"

class sejtkrmat
{
 public:
  sejtkrmat();    //constructor
  ~sejtkrmat();   //destructor

  void read(XFILE *in);

  double get_sw(double pw);
  double get_dsw_dpw(double pw);
  double get_krw(double pw);
  double get_phi();
  double get_kintr();
  double get_rhos();
  double get_alpha();
  double get_emod();
  double get_nu();

  double get_ks();
  double get_kw();
  double get_muw();

  double get_kuw(double pw);
  double get_kwu(double pw);
  double get_kww(double pw);

  double get_capuw(double pw);
  double get_capwu(double pw);
  double get_capww(double pw);

  double get_fw1(double pw);
  double get_fu1(double pw);

 private:
  
  double ab,bb;
  double emod,nu,alpha,ks,phi0,rhos;
  double kw,rhow,muw;
  double kintr;
};  

#endif
