#ifndef GARDNER_RETEN_H
#define GARDNER_RETEN_H

#include "genfile.h"

/**
   This class defines retention curve of Gardner and Schrefler's model
*/

class gardner_reten
{
 public:
  gardner_reten();    //constructor
  ~gardner_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double sw(double pw);
  double dsw_dpw(double pw);
  double dsw_dt(double pw);
  double krw(double pw);
  double get_k(double pw);

 private:
  double ssat,sirr,alfa,ksat,gamaw;
};  

#endif
