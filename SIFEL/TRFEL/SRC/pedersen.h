#ifndef PED_H
#define PED_H

#include "genfile.h"

class pedmat
{
 public:
  pedmat();    //constructor
  ~pedmat();   //destructor
  
  double give_hum (long nn);
  double give_temp (long nn);

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void matcond2 (matrix &d,long ri,long ci,long ipp);
  
  void matcond1d_2 (matrix &d,long ri,long ci,long ipp);
  void matcond2d_2 (matrix &d,long ri,long ci,long ipp);
  void matcond3d_2 (matrix &d,long ri,long ci,long ipp);

  double perm_ww(double w,double t);
  double perm_wt(double w,double t);
  double perm_tw(double w,double t);  
  double perm_tt(double w,double t);
  double c_ww (double w,double t);
  double c_wt (double w,double t);
  double c_tw (double w,double t);
  double c_tt (double w,double t);

  void values_correction (vector &nv);
  void moisture_check(double &w,double t,long ipp);

  double get_delta_gw(double phi,double w);
  double get_kw_g(double w);
  double sorption_isotherm(double phi);
  double inverse_sorption_isotherm(double w);
  double suction_curve(double s);
  double inverse_suction_curve(double w,double t);
  double get_p_gws(double t);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  
  double get_transmission_transcoeff_ww(double w,double t,long bc,long ipp);
  double get_transmission_nodval_ww(double bv,double w,double t,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double w,double t,long bc,long ipp);
  
  double get_transmission_transcoeff_tt(double w,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double w,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double w,double t,long bc,long ipp);

  double get_othervalue(long compother,double w,double t);
  void print_othervalue_name(FILE *out,long compother);

  void read(XFILE *in);
  void print(FILE *out);


 private: 
  double por;         //porosity
  double a_0;         //constant (obtained from experiments) a_0 [Bazant and Najjar, 1972]
  double nn;          //constant-exponent (obtained from experiments) n [Bazant and Najjar, 1972]
  double phi_c;       //constant-relative humidity  (obtained from experiments) phi_c [Bazant and Najjar, 1972]
  double delta_wet;   //constant-water vapor permeability (obtained from experiments) delta_wet [Bazant and Najjar, 1972]
  double delta_dry;   //constant-water vapor permeability (obtained from experiments) delta_dry
  double rho;         //volume density
  
  double ceff;       //effective thermal capacity
  double chieff;     //effective thermal conductivity
  
  double w_h_sorp;         //constant water content (obtained from experiments) w_h [Pedersen, 1990]
  double n_sorp;           //constant-exponent (obtained from experiments) n [Pedersen, 1990]
  double a_sorp;           //constant (obtained from experiments) A [Pedersen, 1990]
  double mw;          //molar mass of water kg.mol-1
  double gasr;        //universal gas constant J.mol-1.K-1
  double rhow;        //kg/m^3 = water density
  double awet,bwet;   //2 coefficients for suction curve
  double k_wg;         //hydraulic conductivity (maximum)
  double ak,bk;       //2 coefficients for hydraulic conductivity
  double nk;          //coefficient for hydraulic conductivity  
  double w_cr,w_cap,w_vac,w_98; //critical water content, capillary water content,vacuum water content,water content accord. 98% rel.hum
  double dhvap;       //enthalpy of evaporation (latent heat of vaporization)

  double get_dphi_dw(double w);
  double get_dphi_dw_hansen(double w);
  double sorption_isotherm_hansen(double w_h_s, double a_s, double n_s, double phi);
  double inverse_sorption_isotherm_hansen(double w);
  double inverse_sorption_isotherm_hansen(double w_h_s, double a_s, double n_s, double w);
  double sorption_isotherm_data_get_w(double phi);
  double inverse_sorption_isotherm_data_get_phi(double w);
  double inverse_sorption_isotherm_data_get_derivative(double w);
  void sorption_isotherm_data_read(XFILE *in);
  void sorption_isotherm_data_print(FILE *out);
  double get_dpgw_dt(double t, double phi);
  double get_dpc_dw(double w, double t);
  double get_dpc_dt(double w, double t);


  ///  model of sorption isotherm from data
  tablefunct *gf;
  int s_type;
};  

#endif
