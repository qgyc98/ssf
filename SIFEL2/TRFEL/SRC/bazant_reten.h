#ifndef BAZANT_RETEN_H
#define BAZANT_RETEN_H

#include "genfile.h"

class bazant_reten
{
 public:
  bazant_reten();    //constructor
  ~bazant_reten();   //destructor

  double sat(double pc,double t);
  double dsat_dpc(double pc,double t);
  double dsat_dt(double pc,double t);

 private:
  double t0;
  double p0;
  double tcr;
};  

#endif
