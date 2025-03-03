#ifndef MASIN_RETEN_H
#define MASIN_RETEN_H

#include "genfile.h"

/**
   This class defines retention curve of Masin and Schrefler's model
*/

class masin_reten
{
 public:
  masin_reten();    //constructor
  ~masin_reten();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double psi(double suction, double dsuc, double e, double temp);
  double get_lambdap(double suction, double dsuc, double e, double temp);
  double get_se(double suction, double dsuc, double e, double temp);
  double get_dse_de(double suction, double dsuc, double e, double temp);
  double sw(double suction, double dsuc, double e, double temp);
  double dsw_dpw(double suction, double dsuc, double e, double temp);
  double dsw_dt(double suction, double dsuc, double e, double temp);
  double dsw_depsv(double suction, double dsuc, double e, double temp);

 private:
  long sr_type; /// type of retention curve and effective stress factor
  double se0; /// the air entry value of suction for the reference macrostructural void ratio e0
  double e0;  /// reference void ratio
  double ae,lambdap0,at,bt,tref,gamma0; ///parameters
};  

#endif
