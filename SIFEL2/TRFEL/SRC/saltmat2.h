#ifndef SALTMAT2_H
#define SALTMAT2_H

#include <stdio.h>
#include "genfile.h"

class saltmat2
{
 public:
  saltmat2 (void);    //constructor
  ~saltmat2 (void);   //destructor

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  void matcond2d2 (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);

  double k11 (double x1,double x2,double x3);
  double k12 (double x1,double x2,double x3);
  double k13 (double x1,double x2,double x3);

  double k21 (double x1,double x2,double x3);
  double k22 (double x1,double x2,double x3);
  double k23 (double x1,double x2,double x3);

  double k31 (double x1,double x2,double x3);
  double k32 (double x1,double x2,double x3);
  double k33 (double x1,double x2,double x3);

  double c11 (double x1,double x2,double x3);
  double c12 (double x1,double x2,double x3);
  double c13 (double x1,double x2,double x3);

  double c21 (double x1,double x2,double x3);
  double c22 (double x1,double x2,double x3);
  double c23 (double x1,double x2,double x3);

  double c31 (double x1,double x2,double x3);
  double c32 (double x1,double x2,double x3);
  double c33 (double x1,double x2,double x3);

  void auxiliarydata (double x1,double x2,double x3);
  
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_nodval_11(double bv,double x1,double x2,double x3,long bc,long ipp);
  double get_transmission_transcoeff_11(double x1,double x2,double x3,long bc,long ipp);
  double get_transmission_flux_11(double bv,double x1,double x2,double x3,long bc,long ipp);

  
  double get_othervalue(long compother,long ipp, double x1,double x2,double x3);
  void print_othervalue_name(FILE *out,long compother);
  double pgws(double t);
  double permeabilitavodnipary(double t, double  p);
  void inverze_sorption_izoterm_data(double x1, double &fi, double &dfi);
  void sisotherm(int kod, double x1,double &fiw, double &dfdw);
  double linear_data(int kod, double x1,double x2, double x3);
  double kapa(int kod, double x1,double x2, double x3);
  double diffcoefiont(int kod, double x1,double x2, double x3);
  double mw, ma, gasr;
     void values_correction (vector & nv);
     void Cf_check(double cf, double w, long ipp);
     void binding_izoterm_derivation(double x2, double & derbi);
 private:
  int kd;
  double a1,a2,a3;
  int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  double MatConst [20];     // constanti hodnoty mat. vlastnosti
  double MatData [20][3][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych
  void CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2);


};

#endif

