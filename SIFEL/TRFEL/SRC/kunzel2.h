#ifndef KUNMAT2_H
#define KUNMAT2_H

#include "genfile.h"

class kunmat2
{
public:
  kunmat2();    //constructor
  ~kunmat2();   //destructor

  double give_hum (long nn);
  double give_temp (long nn);

  void print(FILE *out);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void values_correction (vector &nv, long ipp);
  void relhum_check(double &x1,double x2,long ipp);

  double DerivaceTlakuNasycenychParNaTeplote(double x1, double x2);
  double PermeabilitaVodniPary( double x1, double x2, long ipp);
  double tokJ1(double x1, double x2, long ipp);
  double tokJ2(double x1, double x2, long ipp);
  double tokJ4(double x1, double x2, long ipp);
  double tokJ3(double x1, double x2, long ipp);
  double TlakNasycenychVodnichParNaTeplote(double x1, double x2);
  double LatentHeatofEvaporationOfWater(double x1, double x2);
  double DerivaceHustotyEntalpiePodleTeploty(double x1, double x2, long ipp);
  double HygroscopicMoisture(double x1, double x2);
  double DerivativeOfTheSorptionIsotherm(double x1, double x2);
  double DerivativeOfTheRetentionCurve(double x1, double x2);
  double DerivativeOfTheMoistureRetentionCharacteristik(double x1, double x2, long ipp);
  double soptionizothermDerivation(double x1, double x2);
  void read(XFILE *in);
 
  void sorption_izoterms_values(int kod, long ipp, double x1,double xpv, double ineq1, double & w, double & dwdf);
  void kapa_values(int kod, long ipp,double x1, double xpv, double ineq1, double &kapa);
  void hystereze2 (int matchar, double x, double xpv, double ineq1,double & outvalue,double & outvalue2, long ipp);


  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_transcoeff_hh(double x1,double x2,long bc,long ipp);
  double get_transmission_nodval_hh(double bv,double x1,double x2,long bc,long ipp);
  double get_transmission_flux_hh(double bv,double x1,double x2,long bc,long ipp);

  double get_transmission_transcoeff_tt(double x1,double x2,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double x1,double x2,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double x1,double x2,long bc,long ipp);

  double get_othervalue(long compother,double x1,double x2, long ipp);
  void print_othervalue_name(FILE *out,long compother);

  void give_data(double rh,double Mhmc, double Smc, double Mhrh,double & moistakt);

  double kapa_exp(double a, double b,double x1w, double x2, long ipp);
  void CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2);
 // void sorption_izoterms_giva_data(int kod,double x1, double x2, double & moist, double & dmoistdrh, long ipp);
  double si_kk_hansen(double x1, double x2,double u, double a, double n);

  double derivation_sorption_izoterm_data(double x1, double x2, double x3, long ipp);
  double derivation_dy_dx (int matchar, double prom, int pomk1, int pomk2);
  void der_value_hyst (int matchar,int kod, double pv, double & outvalue,double & outvalue2, long ipp);
  
  void give_values (long ipp,double *av,double *inp, double *ineq);
  void aux_values (long ipp,double *in,double *inp, double *ineq,double *out);
  void save_values (long ipp,double *out);
  
  void initvalues (long ipp,long ido);

 private:
  double J1; // tok vlhkosti podle fi
  double J2; // tok vlhkosti podle T
  double J3; // tok energie podle Fi
  double J4; // tok energie podle T

  int madripom;
  int kd;
  double a1,a2,a3;
  double k1[20],k2[20],k3[20],k4[20],k5[20],k6[20];


  // nove udelane materialove charakterisitky

  int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  double MatConst [20];     // constanti hodnoty mat. vlastnosti
  double MatData [20][4][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych
  double Init[20];
  
  
  //void aux_values (double x1, double x2, long ipp);
  
  double rho_m,rho_w,moist,dmoistdrh;
};

#endif
