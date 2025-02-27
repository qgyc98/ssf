#ifndef BAZPED_H
#define BAZPED_H

#include "genfile.h"

class bazpedmat
{
 public:
  bazpedmat();    //constructor
  ~bazpedmat();   //destructor
  
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
  double sorption_isotherm(double w_h_s, double a_s, double n_s, double phi);
  double inverse_sorption_isotherm(double w);
  double inverse_sorption_isotherm(double w_h_s, double a_s, double n_s, double w);
  double suction_curve(double s);
  double inverse_suction_curve(double w,double t);
  double get_p_gws(double t);

  double transmission_nodval (double nodval,double trc2,long ri,long ci,long nid,long bc);
  double get_transmission_nodval_ww (double bv,double w,double t,long bc);
  double get_transmission_nodval_tt (double bv,double t,long bc);

  double transmission_flux (double nodval,double trc2,long ri,long ci,long nid,long bc);
  double get_transmission_flux_ww (double bv,double w,double t,long bc);
  double get_transmission_flux_tt (double bv,double w,double t,long bc);
  
  double transmission_transcoeff (double trc,long ri,long ci,long nid,long bc);
  double get_transmission_transcoeff_ww (double w,double t,long bc);
  double get_transmission_transcoeff_tt (double w,double t,long bc);

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
  
  double w_h_sorp;         //water content (obtained from experiments) w_h [Pedersen, 1990] = Hansen's sorption isotherm
  double n_sorp;           //exponent (obtained from experiments) n [Pedersen, 1990] = Hansen's sorption isotherm
  double a_sorp;           //parameter (obtained from experiments) A [Pedersen, 1990] = Hansen's sorption isotherm
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
  double get_dpgw_dt(double t, double phi);
  double get_dpc_dw(double w, double t);
  double get_dpc_dt(double w, double t);
};  

#endif
