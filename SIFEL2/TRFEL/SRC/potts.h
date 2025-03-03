#ifndef POTTS_RETEN_H
#define POTTS_RETEN_H

#include "genfile.h"

/**
   This class defines retention curve of Potts and Schrefler's model
*/

class potts_reten
{
 public:
  potts_reten();    //constructor
  ~potts_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double sw(double pw);
  double dsw_dpw(double pw);
  double dsw_dt(double pw);
  double krw(double pw);
  double get_k(double pw);

 private:
  double ssat,sirr,hp_min,htz,r,gammaw,ksat;

};  

#endif
