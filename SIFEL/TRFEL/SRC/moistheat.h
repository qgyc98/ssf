#ifndef MOISTHEATMAT_H
#define MOISTHEATMAT_H

#include "genfile.h"
#include "isotherm.h"
#include "dampermeability.h"

/**
   class contains Kunzel material model for coupled heat and moisture transfer
   
   ip[ipp].av[0] - actual partial pressure of water vapour
   ip[ipp].av[1] - actual temperature
   
   ip[ipp].pv[0] - the partial pressure of water vapour from the previous time step
   ip[ipp].pv[1] - the temperature from the previous time step
   
   values stored in eqother array
   Tm->ip[ipp].eqother[0] - volumetric moisture content
   Tm->ip[ipp].eqother[1] - derivative of the sorption isotherm
   Tm->ip[ipp].eqother[2] - saturated volumetric moisture content
   Tm->ip[ipp].eqother[3] - temperature of transport parameter
   Tm->ip[ipp].eqother[4] - volume increase of ice in pores
   Tm->ip[ipp].eqother[5] - lower radius of ice filling a pore
   Tm->ip[ipp].eqother[6] - upper radius of ice filling a pore
   Tm->ip[ipp].eqother[7] - volumetric moisture content with unfrozen water
   
   ice volume is defined by lower and upper radius, 
   the upper radius is obtained from the volumetric moisture
   the lower radius is obtained from the temperature
   both are calculated for each integration point
   
   JM, revised by JK 10. 10. 2013
*/
class moistheatmat
{
 public:
  moistheatmat();    //constructor
  ~moistheatmat();   //destructor
  
  void read(XFILE *in);
  void print(FILE *out);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  
  void values_correction (vector &nv);
  void values_correction_ipp (long ipp);
  
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
  double derivative_of_the_enthalpy_density_h(long ipp);
  void gibbs_thomson(long ipp);
  double freezing_volume_increment (double w, double t, double pt, double &rl, double &ru, double &icevol);

  double transmission_nodval (double nodval,long ri,long ci,long nid,long bc);
  double get_transmission_nodval_hh (double bv,double rh,double t,long bc);
  double get_transmission_nodval_th (double bv,long bc);
  double get_transmission_nodval_tt (double bv,long bc);


  double transmission_transcoeff (double trc,long ri,long ci,long nid,long bc);
  double get_transmission_transcoeff_hh (double t,long bc);
  double get_transmission_transcoeff_tt (double trcp,long bc);
  double get_transmission_transcoeff_th (double t,long bc);

  double transmission_flux (double nodval,long ri,long ci,long nid,long bc);
  double get_transmission_flux_hh (double bv,double rh,double t,long bc);
  double get_transmission_flux_tt (double bv,double t,long bc);
  double get_transmission_flux_th (double bv,double pv,double t,long bc);

  double get_othervalue(long compother,double rh,double t, long ipp);
  void print_othervalue_name(FILE *out,long compother);

  void aux_values (long ipp,double *inv,double *inp,double *ine,double *out);
  void save_values (long ipp,double *out);
  void give_values (long ipp,double *av,double *pv,double *eq);
  void initvalues (long ipp,long ido);
  
  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);

  ///  function returns temperature in integration point
  double give_temperature (long ipp);

  ///  function returns initial temperature in integration point
  double give_inittemperature (long ipp);

  ///  function returns relative humidity in integration point
  double give_rel_hum (long ipp);

  ///  function returns partial pressure of water vapour in integration point
  double give_press_water_vapor (long ipp);
  
  ///  function returns volumetric moisture content
  double give_vol_moist (long ipp);

  ///  function returns volume change
  double give_volume_change (long ipp);

  /// function returns required non-transport quantities
  void give_reqntq(long *antq);
  
  
  long kd;
  double a1,a2,a3;
  
  
  // nove udelane materialove charakterisitky

  gfunct *data[20];
  
  
  double rho_m,rho_w,moist,dmoistdrh;
  double M, R;
  double cs, cl, c_ice, ro_ice, beta, t_1, t_2, del_l;

  ///  sorption isotherm
  isotherm isoth;
  
  ///  density
  gfunct rho;
  ///  sorption isotherm
  isotherm sorpiso;
  ///  saturated volumetric moisture content
  gfunct sm;
  ///  specific heat capacity
  gfunct c;
  ///  thermal conductivity
  gfunct lambda;

   ///  permeabilita vodni pary
  gfunct dpp;

  ///  transportni clen
  gfunct transpar;

  ///  temperature of transport parametr
  gfunct teptr;

  /// porozitemtrie
  gfunct poroz;

  ///  influence of damage on permeability
  static dampermeability damper;
  ///  flag for influence of damage on permeability
  flagsw daminfl;
};

#endif
