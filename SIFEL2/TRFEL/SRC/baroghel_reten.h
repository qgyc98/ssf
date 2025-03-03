#ifndef BAROGHEL_RETEN_H
#define BAROGHEL_RETEN_H

#include "genfile.h"

class baroghel_reten
{
 public:
  baroghel_reten();    //constructor
  ~baroghel_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);
  double sw(double pc,double t);
  double dsw_dpc(double pc,double t);
  double dsw_dt(double pc,double t);

 private:
  
  double t0,tcr,b,q2,q3,n,z;
};  

#endif
