#ifndef FEBEX_GRANIT_RETEN_H
#define FEBEX_GRANIT_RETEN_H

#include "genfile.h"

/**
   This class defines retention curve of Van_Genuchten and Schrefler's model
*/

class febex_granit_reten
{
 public:
  febex_granit_reten();    //constructor
  ~febex_granit_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double sw(double pw);
  double dsw_ds(double pw);
  double dsw_dt(double pw);
  double get_krw(double sw);

 private:
  double consta,n1,n2;
  //FEBEX granit krw:
  double m1,m2,m3,m4;
};  

#endif
