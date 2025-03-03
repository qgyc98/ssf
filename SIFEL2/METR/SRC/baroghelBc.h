#ifndef BAROGHELBC_H
#define BAROGHELBC_H

#include "genfile.h"

class baroghelmatc
{
 public:
  baroghelmatc();    //constructor
  ~baroghelmatc();   //destructor

  double baroghel_sw(double pc,double t);
  double baroghel_dsw_dpc(double pc,double t);
  double baroghel_dsw_dt(double pc,double t);
  double baroghel_ssp();
  double baroghel_krg(double s);
  double baroghel_krw(double pc, double t);
  double baroghel_phi();
  double baroghel_kintr();
  double baroghel_cps();
  double baroghel_rhocp(double pc,double pg,double t);
  double baroghel_cp(double pc,double pg,double t,long ipp);
  double baroghel_tau(double pc,double t);
  double baroghel_dd(double pc,double t);
  double baroghel_deff(double pc,double pg,double t);
  double baroghel_lambdaeff(double pc,double pg,double t);
  double baroghel_rhos();
  double baroghel_betas();
  double baroghel_alpha();
  double baroghel_emod();
  double baroghel_nu();

  double baroghel_dehydw_dt(double pc,double pg,double t);
  double baroghel_hydw(double pc,double pg,double t);
  double baroghel_hydren(double pc, double pg, double t);
  double baroghel_fste(double pc,double pg,double t);
  double baroghel_ddbw(double pc,double pg,double t);

  void read(XFILE *in);

 private:
  
  double gasr;
  double ma;
  double mw;
  double t0,p0,tcr;
  double betas;
  double scr,ag;
  double av,fs;
  double ab, bb, c, lambdab, cpb, rhosb, phib;
  double alpha;
  double emod,nu;

  //Hydration energy
  double hydren;
  //finv= aging factor
  double finv;
  //fste= Water/Cement ratio
  double fste;
  //diffusion of bound water at refference temperature (295K)
  double ddbw0;
  //cement content
  double c1;
  //free water content
  double w1;
};  

#endif
