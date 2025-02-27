#ifndef KUNMAT_H
#define KUNMAT_H

#include "genfile.h"

class kunmat
{
public:
  kunmat();    //constructor
  ~kunmat();   //destructor

  void print(FILE *out);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void values_correction (vector &nv, long ipp);
  void relhum_check(double &rh,double &t,long ipp);

  double DerivaceTlakuNasycenychParNaTeplote(double Fi, double TK);
  double PermeabilitaVodniPary( double Fi, double TK, long ipp);
  double tokJ2(double Fi, double TK, long ipp);
  double tokJ4(double Fi, double TK, long ipp);
  double tokJ3(double Fi, double TK, long ipp);
  double TlakNasycenychVodnichParNaTeplote(double Fi, double TK);
  double LatentHeatofEvaporationOfWater(double Fi, double TK);
  double DerivaceHustotyEntalpiePodleTeploty(double Fi, double TK, long ipp);
  double HygroscopicMoisture(double Fi, double TK);
  double DerivativeOfTheSorptionIsotherm(double Fi, double TK);
  double DerivativeOfTheRetentionCurve(double rh, double tk, long ipp);
  double DerivativeOfTheMoistureRetentionCharacteristik(double Fi, double TK, long ipp);
  double soptionizothermDerivation(double Fi, double TK);
  double tokJ1(double Fi, double TK, long ipp);
  void read(XFILE *in);


  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_transcoeff_hh(double rh,double t,long bc,long ipp);
  double get_transmission_nodval_hh(double bv,double rh,double t,long bc,long ipp);
  double get_transmission_flux_hh(double bv,double rh,double t,long bc,long ipp);

  double get_transmission_transcoeff_tt(double h,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double h,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double h,double t,long bc,long ipp);

  double get_othervalue(long compother,double rh,double t, long ipp);
  void print_othervalue_name(FILE *out,long compother);

  void give_data(double rh,double Mhmc, double Smc, double Mhrh,double & moistakt);

  double kapa_exp(double a, double b,double rh, double tk, long ipp);
  void CorD(int cislochar, int &kvyhl,double in, double x, double & y, double & z, double &z2);
  void sorption_izoterms_giva_data(int kod,double rh, double tk, double & moist, double & dmoistdrh, long ipp);
  double si_kk_hansen(double rh, double tk,double u, double a, double n);

  double derivation_sorption_izoterm_data(double x1, double x2, double x3);
  void sorption_izoterms_values(int kod, long ipp, double x1,double xpv, double ineq1, double & w, double & dwdf);
  void aux_values (long ipp,double *in,double *inp, double *ineq,double *out);
  void save_values (long ipp,double *out);
  void give_values (long ipp,double *av,double *inp, double *ineq);
  void initvalues (long ipp,long ido);

   void kapa_values (int kod, long ipp,double x1, double xpv, double ineq1, double &kapa);

 private:
  double J1; // tok vlhkosti podle fi
  double J2; // tok vlhkosti podle T
  double J3; // tok energie podle Fi
  double J4; // tok energie podle T

  int madripom;
  int kd;
  double a1,a2,a3;


  // nove udelane materialove charakterisitky

  int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  double MatConst [20];     // constanti hodnoty mat. vlastnosti
  double MatData [20][3][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych

  
  
  void aux_values (long ipp);
  
  double rho_m,rho_w,moist,dmoistdrh;
};

#endif
