#ifndef LEWIS_RETEN_H
#define LEWIS_RETEN_H

#include "genfile.h"

/**
   This class defines retention curve of Lewis and Schrefler's model
*/

class lewis_reten
{
 public:
  lewis_reten();    //constructor
  ~lewis_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double sw(double pw);
  double dsw_dpw(double pw);
  double dsw_dt(double pw);
  double krw(double sw);
  double get_se(double sw);
  double get_krg(double pw);

 private:
  double sr;
};  

#endif
