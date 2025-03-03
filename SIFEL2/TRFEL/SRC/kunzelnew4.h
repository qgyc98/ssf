#ifndef KUNMAT_H
#define KUNMAT_H

#include "genfile.h"

/**
   class contains Kunzel material model for coupled heat and moisture transfer
   
   list of material parameters used in the model:
    position CORD:  2 - density
		    3 - porosity
		    4 - water vapour diffusion resistance factor
		    5 - moisture diffusivity
		    6 - sorption isoterm
		    7 - saturated moisture
		    8 - none
		    9 - specific heat capacity
		    10 - thermal conductivity
		    11 - 13 - none
		    14 Dcoef
		    15 - binding isotherm
		    16 - cfmax
		    17 ws
		    18 - none
		    19 - kunzeltype
   
   values stored in eqother array
   Tm->ip[ipp].eqother[0] - volumetric moisture content
   Tm->ip[ipp].eqother[1] - derivative of the sorption isotherm
   Tm->ip[ipp].eqother[2] - saturated volumetric moisture content
   Tm->ip[ipp].eqother[3] - moisture diffusivity
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

  double derivative_of_saturated_water_vapor_pressure(double tk);
  double water_vapour_permeability(double tk, long ipp);
  double kmm(long ipp);
  double kmt(long ipp);
  double khm(long ipp);
  double kht(long ipp);
  double cmm(long ipp);
  double cmt(long ipp);
  double chm(long ipp);
  double cht(long ipp);
  double saturated_water_vapor_pressure(double tk);
  double latent_heat_of_evaporation_of_water(double tk);
  double derivative_of_the_enthalpy_density(long ipp);
  double hygroscopic_moisture(double rh, double tk);
  double derivative_of_the_sorption_isotherm_root(double rh, double tk);
  double derivative_of_retention_curve_root(double rh, double tk, long ipp);
  double sorption_izotherm_derivative_kk_hansen(double rh, double tk);


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

   void kapa_values (long kod, long ipp,double x1, double &kapa);

 private:

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
  double M, R;
  
  //  type of Kunzel model
  //  kunzeltype=1 - relative humidity and temperature
  //  kunzeltype=2 - partial pressure and temperature
  //  kunzeltype=3 - partial pressure and temperature, without liquid moisture
  long kunzeltype;
};

#endif
