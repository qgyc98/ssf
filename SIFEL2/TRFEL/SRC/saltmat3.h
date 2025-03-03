#ifndef saltmat3_H
#define saltmat3_H

#include <stdio.h>
#include "genfile.h"

/**
   class describes material model which deals with simultaneous
   transport of moisture, salt and salt crystals
   
   model uses array eqother
   eqother[0] - moisture diffusivity (kappa)
   eqother[1] - binding izoterm
   eqother[2] - salt diffusivity (D)
   eqother[3] - ???
   
   JM, 29.5.2007
*/
class saltmat3
{
 public:
  saltmat3 (void);    //constructor
  ~saltmat3 (void);   //destructor
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  void matcond2d2 (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);

  double k11 (double x1,double x2,double x3,long ipp);
  double k12 (double x1,double x2,double x3,long ipp);
  double k13 (double x1,double x2,double x3,long ipp);

  double k21 (double x1,double x2,double x3,long ipp);
  double k22 (double x1,double x2,double x3,long ipp);
  double k23 (double x1,double x2,double x3,long ipp);

  double k31 (double x1,double x2,double x3,long ipp);
  double k32 (double x1,double x2,double x3,long ipp);
  double k33 (double x1,double x2,double x3,long ipp);

  double c11 (double x1,double x2,double x3,long ipp);
  double c12 (double x1,double x2,double x3,long ipp);
  double c13 (double x1,double x2,double x3,long ipp);

  double c21 (double x1,double x2,double x3,long ipp);
  double c22 (double x1,double x2,double x3,long ipp);
  double c23 (double x1,double x2,double x3,long ipp);

  double c31 (double x1,double x2,double x3,long ipp);
  double c32 (double x1,double x2,double x3,long ipp);
  double c33 (double x1,double x2,double x3,long ipp);

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
  void values_correction (vector & nv,long ipp);
  void Cf_check(double cf, double w, long ipp);
  void binding_izoterm_derivation(double x2, double & derbi);
  void kapa_values (int kod, long ipp, double &kapa, double av);
  void salt_diffusivity_values (int kod, long ipp, double &diff, double x2);
  void binding_izoterm_derivation_values (int kod, long ipp, double &dcbdcf, double x2);
  void initvalues (long ipp,long ido);
  void aux_values (long ipp);
  void hystereze (int matchar, int matchar2, double & outvalue, long ipp);
  void der_value_hyst (int matchar,int kod, double pv, double & outvalue,double & outvalue2, long ipp);
  double get_moisture(double rh);
  double get_rel_hum(double w);
  double derivation_dy_dx (int matchar, double prom, int pomk1, int pomk2);


 private:
  int kd;
  double a1,a2,a3;
  int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  double MatConst [20];     // constanti hodnoty mat. vlastnosti
  double MatData [20][4][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakteristika zadana nejakou funkci - 5 promennych
  void CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2);
  double Init[20];
  double k1[20],k2[20],k3[20],k4[20],k5[20],k6[20];
  
  
};

#endif

