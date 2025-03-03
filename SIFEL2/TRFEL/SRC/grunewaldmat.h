#ifndef GRUNEWALDMAT_H
#define GRUNEWALDMAT_H

#include <stdio.h>
#include "genfile.h"

class grunewaldmat
{
 public:
  grunewaldmat (void);    //constructor
  ~grunewaldmat (void);   //destructor
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void read (XFILE *in);
  void print(FILE *out);

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

  void auxiliarydata (double x1,double x2,double x3,long ipp);
  
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_nodval_11(double bv,double x1,double x2,double x3,long bc,long ipp);
  double get_transmission_transcoeff_11(double x1,double x2,double x3,long bc,long ipp);
  double get_transmission_flux_11(double bv,double x1,double x2,double x3,long bc,long ipp);

  
  double get_othervalue(long compother,long ipp, double x1,double x2,double x3);
  void print_othervalue_name(FILE *out,long compother);

  void give_data_si_root_fi(double x1, double x2, double x3, double w_hyg, double w_sat, double rh_hyg, double shift_w,  double & relh);

  double kapa_exp(double a, double b,double x1, double x2, double x3);
  void CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2);
  void sorption_izoterms_giva_data(int kod,double x1, double x2, double x3, double & moist, double & dmoistdrh);
  double si_kk_hansen(double x1, double x2, double x3,double u, double a, double n);

  void aux_values (long ipp,double *in,double *inp, double *ineq,double *out);
  void save_values (long ipp,double *out);
  void give_values (long ipp,double *av,double *inp, double *ineq);
  void initvalues (long ipp,long ido);

  double density_lql(double x1,double x2, double x3);
     double derivation_density_droldt(double x1, double x2, double x3);
     double density_lqw(double x1, double x2, double x3);
     double derivation_density_drowdt(double x1, double x2, double x3);
     double saturation_water_vapour_pressure(double x1, double x2, double x3);
     double derivation_saturation_water_vapour_pressure_temperature(double x1, double x2, double x3);
        double specific_internal_energy_of_the_liquid_phase(double x1, double x2, double x3);
        double specific_internal_energy_of_the_liquid_water(double x1, double x2, double x3);
        double specific_internal_energy_of_water_vapour(double x1, double x2, double x3);
        double derivation_of_specific_internal_energy_of_the_solid_material_dependence_temperature(double x1, double x2, double x3);
        double derivation_of_specific_internal_energy_of_the_liquid_phase_dependence_temperature(double x1, double x2, double x3);
        double partial_water_vapour_pressure_function(double x1, double x2, double x3);
        double diffusion_number_function(double x1, double x2, double x3);
     double specific_enthalpy_of_water_vapour_hv(double x1, double x2, double x3);
        void sorption_izotherm_derivation(double x1, double x2, double x3, double & derfi);
     double derivation_specific_internal_energy_of_water_vapour(double x1, double x2, double x3);
     void get_rel_hum(double w, double &fi, double &dfdw);
     double sortpion_isotherm_root_shifted(double x1, double w_hyg, double rh_hyg);
     void give_data_si_root_dfidw(double x1, double x2, double x3, double rh_hyg, double w_sat, double w_hyg, double shift_w, double & dfdw);
     double get_moisture(double rh);
	 double give_kapa(double x1, double x2);

 private:
     double rom; // mass density of the solid material
     double row; //density of the liquid water
     double rov; // water vapour density
     double pvs; // saturation water vapour pressure
     double dpvsdt;
     double drowdt; // droldt
     double droldt; // derivation of the density of the liquid phase dependence of temperature
     double dfidt; // rel. vlhkost na teplote
     double pv; //partial water vapour pressure
     double rol; // density of the liquid phase
     double ul; //specific internal energy of the liqiud phase
     double uv; //specific internal energy of water vapour
     double dumdt; //derivation of specific internal energy of the solid material dependence of temperature
     double duvdt; //derivation of specific internal energy of water vapour dependence of temperature
     double duldt; //derivation of specific internal energy of the liquid phase dependence of temperature
     double hv;
     double pir; // porovitost
     double relhum;
     double dp;
     double dvn;
     double mi;
     double dfidw;
     double kapa;
     double lambda;
     double cmat;
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


