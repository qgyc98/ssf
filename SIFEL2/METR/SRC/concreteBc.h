#ifndef CONCRETEBC_H
#define CONCRETEBC_H

#include "genfile.h"

class concreteBmatc
{
 public:
  concreteBmatc();    //constructor
  ~concreteBmatc();   //destructor

  double concreteB_sw(double pc,double t);
  double concreteB_dsw_dpc(double pc,double t);
  double concreteB_dsw_dt(double pc,double t);
  double concreteB_ssp();
  double concreteB_krg(double s);
  double concreteB_krw(double s,double rh);
  double concreteB_phi(double t);
  double concreteB_kintr(double pg,double t,double dam);
  double concreteB_cps(double t);
  double concreteB_rhocp(double pc,double pg,double t,long ipp);
  double concreteB_cp(double pc,double pg,double t,long ipp);
  double concreteB_fs(double pc,double t);
  double concreteB_tau(double pc,double t);
  double concreteB_deff(double pc,double pg,double t);
  double concreteB_lambdaeff(double pc,double pg,double t);
  double concreteB_kt(double pc,double pg,double t);
  double concreteB_ks(double pc,double pg,double t);
  double concreteB_rhos(double t);
  double concreteB_betas();
  double concreteB_emod();
  double concreteB_nu();
  double concreteB_dmdh_dt(double pc,double pg,double t);
  double concreteB_dhdehydr(double pc,double pg,double t);
  double concreteB_drhos_dgammadh(double pc,double pg,double t);
  double concreteB_dgammadh_dt(double pc,double pg,double t);
  void read(XFILE *in);

 private:
  
  double gasr;//universal gas constant
  double ma;//molar mass of dry air
  double mw;//molar mass of water
  
  //cement content
  double c1;
  //finv= aging factor
  double finv;
  //fste= Water/Cement ratio
  double fste;
  //Hydration energy
  double hydren;

  double emod,nu;
  double t0,p0,tcr;
  double rhos_th0,betas;//temp.
  double scr,ag,sir,aw,bw,phi0,aphi,k0,ak,bk;
  double cps0,ac,tref;
  double ads,bds,nds;
  double av,fs;
};  

#endif
