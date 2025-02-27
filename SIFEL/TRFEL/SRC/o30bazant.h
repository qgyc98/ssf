#ifndef o30_BAZANT_H
#define o30_BAZANT_H

#include "genfile.h"

class o30bazmat
{
 public:
  o30bazmat();    //constructor
  ~o30bazmat();   //destructor

  double sat(double pc,double t);
  double dsat_dpc(double pc,double t);
  double dsat_dt(double pc,double t);
  double ssp();
  double o30baz_phi(double t);
  double o30baz_kintr(double pc,double pg,double t,double dam);
  double o30baz_krg(double pc,double t);
  double o30baz_krw(double pc,double t,double rh);
  double o30baz_dd();
  double o30baz_deff(double pc,double pg,double t);
  double o30baz_cps(double t);
  double o30baz_rhocp(double pc,double pg,double t,long ipp);
  double o30baz_fs(double pc,double pg,double t);
  double o30baz_tau();
  double o30baz_lambdaeff(double pc,double pg,double t);
  double o30baz_lambdas(double t);
  double o30baz_rhos(double t);
  double o30baz_betas();

  double o30baz_dehydw_dt(double pc,double pg,double t);
  double o30baz_hydw(double pc,double pg,double t);
  double o30baz_hydren(double pc, double pg, double t);
  double o30baz_fste(double pc,double pg,double t);
  double o30baz_ddbw(double pc,double pg,double t);

  double o30baz_emod();
  double o30baz_fct(double pc,double pg,double t);
  double o30baz_xk0(double pc,double pg,double t);
  double o30baz_bcc();
  double o30baz_alpha();
  double o30baz_nu();
  void read(XFILE *in);
  void print(FILE *out);
  void give_reqntq(long *antq);

 private:
  double mw;
  double ma;
  double gasr;
  
  double t0;
  double p0;
  double t00;
  double tcr;
  //free water content at 20°C for Bazant isotherms
  double w1;
  //cement conctent
  double c1;
  
  // porosity
  double phi0;
  double aphi;
  //intrinsic permeability
  double k0;
  double ak;
  //relative permeability
  double scr,ag,sir,aw,bw;
  //structure coefficient
  //double fs;
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
  //Biot's constant
  double alpha;
  //Poisson's ceofficient
  double vcoeff;
  //cubic thermal expansion coeffcient
  double betas;
  //MAZAR'S COEFFICIENTS
  //tensile coefficients
  double at;
  double bt;
  //compressive coefficients
  double acc;
};  

#endif
