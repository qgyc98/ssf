#ifndef VAN_GENUCHTEN_RETEN_H
#define VAN_GENUCHTEN_RETEN_H

#include "genfile.h"

/**
   This class defines retention curve of Van_Genuchten and Schrefler's model
*/

class van_genuchten_reten
{
 public:
  van_genuchten_reten();    //constructor
  ~van_genuchten_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double sw(double pw, double t);
  double dsw_dpw(double pw, double t);
  double dsw_dt(double pw, double t);
  double get_krw(double pw, double t);
  double get_k(double pw);

 private:
  long vg_ret_type;
  double ssat,sirr,ksat,gamaw,delta,expn,expm;
  double p0,pd,lambda0,lambdap;
  double sig0,t0;

};  

#endif
