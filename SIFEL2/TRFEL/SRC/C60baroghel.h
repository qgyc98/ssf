#ifndef C60_BAROGHEL_H
#define C60_BAROGHEL_H

#include "genfile.h"

class C60barmat
{
 public:
  C60barmat();    //constructor
  ~C60barmat();   //destructor

  double sat(double pc,double t);
  double dsat_dpc(double pc,double t);
  double dsat_dt(double pc,double t);
  double ssp();
  double C60bar_phi(double t);
  double C60bar_kintr(double pc,double pg,double t,double dam);
  double C60bar_krg(double pc,double t);
  double C60bar_krw(double pc,double t,double rh);
  double C60bar_dd(double pc,double t);
  double C60bar_deff(double pc,double pg,double t);
  double C60bar_cps(double t);
  double C60bar_rhocp(double pc,double pg,double t);
  double C60bar_tau(double pc,double t);
  double C60bar_lambdaeff(double pc,double pg,double t);
  double C60bar_lambdas(double t);
  double C60bar_rhos();
  double C60bar_betas();

  double C60bar_dehydw_dt(double pc,double pg,double t);
  double C60bar_hydw(double pc,double pg,double t);
  double C60bar_hydren(double pc, double pg, double t);
  double C60bar_fste(double pc,double pg,double t);
  double C60bar_ddbw(double pc,double pg,double t);

  double C60bar_emod(double pc,double pg,double t);
  double C60bar_fct(double pc,double pg,double t);
  double C60bar_xk0(double pc,double pg,double t);
  double C60bar_bcc(double pc,double pg,double t);
  double C60bar_alpha();
  double C60bar_nu();
  void read(XFILE *in);
  void print(FILE *out);
  void give_reqntq(long *antq);

 private:
  double mw;
  double ma;
  double gasr;
  
  double t0;
  double p0;
  double tcr;
  //free water content at 20°C for Bazant isotherms
  double w1;
  //cement content
  double c1;
  
  // porosity
  double phi0;
  double aphi;
  //intrinsic permeability
  double k0;
  double ak;
  //relative permeability
  double scr,ag,sir,aw,bw;
  //skeleton density
  double rhos;
  //structure coefficient
  double fs;
  //thermal conductivity of solid skeleton
  double lambdas0;
  double alam;
  //thermal capacity of solid skeleton
  double ac;
  double cps0;
  //Hydration energy
  double hydren;
  //finv= aging factor
  double finv;
  //fste= Water/Cement ratio
  double fste;
  //diffusion of bound water at refference temperature (295K)
  double ddbw0;
  //chracteristic length
  double dld;
  //Young's modulus
  double emod0;
  //Poisson's ceofficient
  double vcoeff;
  //cubic thermal expansion coeffcient
  double betas;
  //Biot's constant
  double alpha;
  //MAZAR'S COEFFICIENTS
  //tensile coefficients
  double at;
  double bt;
  //compressive coefficients
  double acc;
};  

#endif
