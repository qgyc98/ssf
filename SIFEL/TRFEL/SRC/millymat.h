#ifndef MILLYMAT_H
#define MILLYMAT_H

#include <stdio.h>
#include "genfile.h"

class millymat
{
 public:
  millymat (void);    //constructor
  ~millymat (void);   //destructor

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);

  double k11 (double x1,double x2,long ipp);
  double k12 (double x1,double x2,long ipp);

  double k21 (double x1,double x2,long ipp);
  double k22 (double x1,double x2,long ipp);

  double c11 (double x1,double x2,long ipp);
  double c12 (double x1,double x2,long ipp);

  double c21 (double x1,double x2,long ipp);
  double c22 (double x1,double x2,long ipp);

  void auxiliarydata (double x1,double x2);
  
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_nodval_11(double bv,double x1,double x2,long bc,long ipp);
  double get_transmission_transcoeff_11(double x1,double x2,long bc,long ipp);
  double get_transmission_flux_11(double bv,double x1,double x2,long bc,long ipp);

  void values_correction (vector &nv);
  double get_othervalue(long compother,long ipp, double x1,double x2);
  void print_othervalue_name(FILE *out,long compother);
  double saturation_water_vapour_pressure(double x1, double x2);
  double derivation_saturation_water_vapour_pressure_temperature(double x1, double x2);
  double relative_volume_ration_a(double x1, double x2, double ul);
  double relative_humidity_pc(double x1, double x2);
  double relative_humidity_psi(double x1, double x2);
  double coeef_zaporneA (double x1, double x2, long ipp);
  double coeef_B (double x1, double x2, long ipp);
  double coeef_C (double x1, double x2, long ipp);
  double diffusion_coefficient_of_water_vapor_in_air(double x1, double x2);
  double partial_water_vapour_pressure_function(double x1, double x2, long ipp);
  double coeff_Kv (double x1, double x2, long ipp);
  double coeff_Dtv (double x1, double x2, long ipp);
  double gaseous_moisture_content_by_mass (double x1, double x2, long ipp);
  double specific_enthalpy_of_porous_matrix (double x1, double x2);
  double specific_enthalpy_of_liquid_phase (double x1, double x2);
  double specific_enthalpy_of_gaseous_phase (double x1, double x2);
  double latent_heat_of_evaporation_of_water(double x1, double x2);
  double derivation_water_retention_curve_by_pressure_head (double x1, double x2);
  double derivation_water_retention_curve_by_temperature (double x1, double x2);
  double give_ul (double x1, double x2, double &fi,double &dmoistdrh);
  double give_ul_sorption_izothemrs (double x1, double x2,double rh, double &dmoistdrh );
  double give_ul_retention_curve (double x1, double x2, double rh, double &dmoistdrh );
  double capilar_pressure (double x1, double x2);
  double ul_to_w (double ul);
  double derivation_sorption_izoterm_data(double x1, double x2, int kod);
  double DerivativeOfTheSorptionIsotherm(double x1, double x2,double rh,double &moistakt1);
  void CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2);
  void initvalues (long ipp,long ido);
  void save_values (long ipp,double *out);
  void give_values (long ipp,double *av, double *pv, double *eq);
  void aux_values (long ipp,double *in,double *inp, double *ineq,double *out);


 private:

     int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
     double MatConst [20];     // constanti hodnoty mat. vlastnosti
     double MatData [20][3][150]; // maximalni pocet 150 radku na jednu charakteristiku
     double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych

     int madripom;
     int kd;
     double a1,a2,a3;


};

#endif

