#ifndef DEVRIESMAT_H
#define DEVRIESMAT_H

#include <stdio.h>
#include "genfile.h"

class devriesmat
{
 public:
  devriesmat (void);    //constructor
  ~devriesmat (void);   //destructor
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);
  void print(FILE *out);

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
  void CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2);
  double saturation_water_vapour_pressure(double x1, double x2, double x3);
  double derivation_saturation_water_vapour_pressure_temperature(double x1, double x2, double x3);
  double partial_water_vapour_pressure_function(double x1, double x2, double x3);
  void sorption_izotherm_derivation(double x1, double x2, double x3, double & derfi);
  double derivation_specific_internal_energy_of_water_vapour(double x1, double x2, double x3);
  double get_rel_hum(double w);
  double sortpion_isotherm_root_shifted(double x1, double w_hyg, double rh_hyg);
  void give_data_si_root_dfidw(double x1, double x2, double x3, double rh_hyg, double w_sat, double w_hyg, double shift_w, double & dfdw);
  double get_moisture(double rh);
  void sorption_izoterms_giva_data(int kod,double x1, double x2, double x3, double & fi, double & dfdw);
  void give_data_si_root_fi(double x1, double x2, double x3, double w_hyg, double w_sat, double rh_hyg, double shift_w,  double & relh);
  double pressure_head(double x1, double x2, double x3);
     double derivation_pressure_head(double x1, double x2, double x3);
     double relative_volume_ration_a(double x1, double x2, double x3);
     double surface_tension(double x1, double x2, double x3);
     double derivation_surface_tension_on_temperature(double x1, double x2, double x3);
     double coeff_of_water_diffusion_by_moisture_grad(double x1, double x2, double x3);
     double coeff_of_water_diffusion_by_temperature_grad(double x1, double x2, double x3);
     double coeff_dwv(double x1, double x2, double x3);
     double coeff_dtv(double x1, double x2, double x3);
     double specific_heat_capacities_star(double x1, double x2, double x3);
     double partial_density_of_water_vapor(double x1, double x2, double x3);
     double latent_heat_of_evaporation_of_water(double x1, double x2, double x3);
     double derivation_drov_dw(double x1, double x2, double x3);
     double diffusion_coefficient_of_water_vapor_in_air(double x1, double x2, double x3);



 private:
        // nove udelane materialove charakterisitky

     int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
     double MatConst [20];     // constanti hodnoty mat. vlastnosti
     double MatData [20][3][150]; // maximalni pocet 150 radku na jednu charakteristiku
     double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych

     int madripom;
     int kd;
     double a1,a2,a3;

};

#endif


