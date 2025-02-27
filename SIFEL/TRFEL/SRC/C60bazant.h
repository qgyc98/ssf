#ifndef C60_BAZANT_H
#define C60_BAZANT_H

#include "aliast.h"
#include "genfile.h"

class C60bazmat
{
 public:
  C60bazmat();    //constructor
  ~C60bazmat();   //destructor

  double sat(double pc,double t);
  double dsat_dpc(double pc,double t);
  double dsat_dt(double pc,double t);
  double ssp();
  double C60baz_krg(double pc,double t);
  double C60baz_krw(double pc, double t);
  double C60baz_phi();
  double C60baz_kintr();
  double C60baz_dd(double pc,double t);
  double C60baz_deff(double pc,double pg,double t);
  double C60baz_cps();
  double C60baz_rhocp(double pc,double pg,double t);
  double C60baz_tau(double pc,double t);
  double C60baz_lambdaeff(double pc,double t);
  double C60baz_lambdas();
  double C60baz_rhos();
  double C60baz_betas();
  double C60baz_dmdh_dt(double pc,double pg,double t);

  double C60baz_dehydw_dt(double pc,double pg,double t);
  double C60baz_hydw(double pc,double pg,double t);
  double C60baz_hydren(double pc, double pg, double t);
  double C60baz_fste(double pc,double pg,double t);
  double C60baz_ddbw(double pc,double pg,double t);

  double C60baz_emod(double pc,double pg,double t);
  double C60baz_fct(double pc,double pg,double t);
  double C60baz_xk0(double pc,double pg,double t);
  double C60baz_bcc(double pc,double pg,double t);
  double C60baz_alpha();
  double C60baz_nu();
  void read(XFILE *in);
  void print(FILE *out);

  void give_dof_names(namevart *dofname, long ntm);

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
  //chracteristic length
  double dld;

  //diffusion of bound water at refference temperature (295K)
  double ddbw0;

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
