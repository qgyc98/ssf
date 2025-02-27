#ifndef KUNMAT_H
#define KUNMAT_H

#include "genfile.h"

/**
   class contains Kunzel material model for coupled heat and moisture transfer
   
   list of material parameters used in the model:
   2 - density
   3 - porosity
   4 - faktor difusniho odporu
   5 - kapa
   6 - sorption izoterm
   7 - saturated moisture
   8 - hydraulic conductivity
   9 - Cecko
   10 - Lambda
   11 - water retention curve
   12 - none
   13 - none
   14 - Dcoef
   15 - binding isotherm
   16 - cfmax
   17 - ws
   18 - none
   19 - none    
   
   values stored in eqother array
   Tm->ip[ipp].eqother[0] - moisture content
   Tm->ip[ipp].eqother[1] - derivative of the sorption isotherm
   Tm->ip[ipp].eqother[2] - maximalni saturace
   Tm->ip[ipp].eqother[3] - vlhkostni vodivost
   Tm->ip[ipp].eqother[4] - 
   Tm->ip[ipp].eqother[5] - 
   
   
   JM
*/
class kunmat
{
public:
  kunmat();    //constructor
  ~kunmat();   //destructor

  void read(XFILE *in);
  void print(FILE *out);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void values_correction (vector &nv);

  double DerivaceTlakuNasycenychParNaTeplote(double rh, double tk);
  double PermeabilitaVodniPary( double rh, double tk, long ipp);
  double tokJ1(double rh, double tk, long ipp);
  double tokJ2(double rh, double tk, long ipp);
  double tokJ4(double rh, double tk, long ipp);
  double tokJ3(double rh, double tk, long ipp);
  double TlakNasycenychVodnichParNaTeplote(double rh, double tk);
  double LatentHeatofEvaporationOfWater(double tk);
  double DerivaceHustotyEntalpiePodleTeploty(double rh, double tk, long ipp);
  double HygroscopicMoisture(double rh, double tk);
  double DerivativeOfTheSorptionIsotherm(double rh, double tk);
  double DerivativeOfTheRetentionCurve(double rh, double tk, long ipp);
  double sorptionizothermDerivation(double rh, double tk);


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

  double kapa_exp(double a, double b,double rh);
  void CorD(long charid,long &kvyhl,double x, double &y, double &z,double &z2);
  void sorption_izoterms_giva_data(long kod,double rh, double tk, double & moist, double & dmoistdrh, long ipp);
  double si_kk_hansen(double rh, double tk,double u, double a, double n);

  double derivative_sorption_izoterm_data(double x1);
  void aux_values (long ipp,double *in,double *inp, double *ineq,double *out);
  void save_values (long ipp,double *out);
  void give_values (long ipp,double *av,double *eq);
  void initvalues (long ipp,long ido);

   void kapa_values (long kod, long ipp,double x1, double xpv, double ineq1, double &kapa);

 private:
  double J1; // tok vlhkosti podle rh
  double J2; // tok vlhkosti podle T
  double J3; // tok energie podle rh
  double J4; // tok energie podle T

  long kd;
  double a1,a2,a3;


  // nove udelane materialove charakterisitky

  long MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  //double MatConst [20];     // constanti hodnoty mat. vlastnosti
  //double MatData [20][3][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych

  gfunct *data[20];

  
  void aux_values (long ipp);
  
  double rho_m,rho_w,moist,dmoistdrh;
};

#endif
