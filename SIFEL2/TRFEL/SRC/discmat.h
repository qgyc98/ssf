#ifndef DISCMAT_H
#define DISCMAT_H

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
class discmat
{
 public:
  discmat (void);    //constructor
  ~discmat (void);   //destructor
  
  double give_x1 (long nn);
  double give_x2 (long nn);

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  void matcond4d (matrix &d,long ri,long ci,long ipp);
  void matcond2d2 (matrix &d,long ri,long ci,long ipp);
  
  void read (XFILE *in);
  void print(FILE *out);
  
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

  
  double get_othervalue(long compother,long ipp, double x1,double x2);
  void print_othervalue_name(FILE *out,long compother);
  double pgws(double x1,double x2);
  double permeabilitavodnipary(double x1,double x2);
  void inverze_sorption_izoterm_data(double x1, double &fi, double &dfi);
  void sisotherm(int kod, double x1,double &fiw, double &dfdw);
  double linear_data(int kod, double x1,double x2);
  double kapa(int kod, double x1,double x2);
  double mw, ma, gasr;
  void values_correction (vector & nv,long ipp);
  
  void kapa_values (int kod, double x1, double xpv, double ineq1, double &kapa);
  
  
  void water_content_relhum (long nid,double *in,double *inp,double *ineq, double *out);
  void get_moisture (long nid,double in,double *inp, double &out);

  void initvalues (long ipp,long ido);

  void give_values (long ipp,double *av, double *pv, double *eq);
  void aux_values(double * in, double * inp, double * ineq, double * out);
  void save_values (long ipp,double *out);
  void aux_values_elements (double *in,double *inp, double *ineq,double *out);

  double derivation_saturation_water_vapour_pressure_temperature(double x1, double x2);
  double latent_heat_of_evaporation_of_water(double x1, double x2);
  void sorption_izotherm_derivation(double x1, double x2, double & derfi);
  void sorption_izoterms_values(int kod, double x1, double xpv, double ineq1, double & fi, double & dfdw);
  void give_data_si_root_dfidw(double x1, double x2,double rh_hyg, double w_sat, double w_hyg, double shift_w, double & dfdw);
  double sortpion_isotherm_root_shifted(double x1, double w_hyg, double rh_hyg);
  double get_rel_hum(double w);
  void give_data_si_root_fi(double x1, double x2, double w_hyg, double w_sat, double rh_hyg, double shift_w,  double & relh);
  double derivation_dy_dx (int matchar, double prom, int pomk1, int pomk2);
  
  double get_moisture2 ( double rh);
  void get_rel_hum2 (double w,double &fi, double &dfdw);
  
 private:
  int kd;
  double a1,a2,a3;
  int MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  double MatConst [20];     // constanti hodnoty mat. vlastnosti
  double MatData [20][4][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakteristika zadana nejakou funkci - 5 promennych
  void CorD(int cislochar, int &kvyhl,double in,int rhw, double x, double & y, double & z, double &z2);
  
};

#endif

