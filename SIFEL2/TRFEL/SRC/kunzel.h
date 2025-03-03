#ifndef KUNMAT_H
#define KUNMAT_H

#include "genfile.h"
#include "isotherm.h"
#include "dampermeability.h"

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
		    14 - Dcoef
		    15 - binding isotherm
		    16 - cfmax
		    17 - ws
		    18 - kmm type
		    19 - kunzeltype

   ip[ipp].av[0] - actual relative humidity
   ip[ipp].av[1] - actual temperature
   
   ip[ipp].pv[0] - the relative humidity from the previous time step
   ip[ipp].pv[1] - the temperature from the previous time step
   
   values stored in eqother array
   Tm->ip[ipp].eqother[0] - volumetric moisture content
   Tm->ip[ipp].eqother[1] - derivative of the sorption isotherm
   Tm->ip[ipp].eqother[2] - saturated volumetric moisture content
   Tm->ip[ipp].eqother[3] - moisture diffusivity
   Tm->ip[ipp].eqother[4] - type of Kunzel model

   
   JM, revised by JK 10. 10. 2013
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

  double get_othervalue(long compother,double rh,double t, long ipp);
  void print_othervalue_name(FILE *out,long compother);

  double kapa_exp(double a, double b,double rh);
  void CorD(long charid,long &kvyhl,double x, double &y, double &z,double &z2);

  //void sorption_izoterms_giva_data(long kod,double rh, double tk, double & moist, double & dmoistdrh, long ipp);

  void aux_values (long ipp,double *inv,double *inp,double *ine,double *out);
  void save_values (long ipp,double *out);
  void give_values (long ipp,double *av,double *pv,double *eq);
  void initvalues (long ipp,long ido);
  
  ///  function evaluates the sorption isotherm
    //double sorption_isotherm_value (double in);
  
  void hystereze (long ipp,double *inv,double *inp,double *ine,double *out);
  
  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);

  ///  function returns temperature in integration point
  double give_temperature (long ipp);

  ///  function returns initial temperature in integration point
  double give_inittemperature (long ipp);

  ///  function returns relative humidity in integration point
  double give_rel_hum (long ipp);
  
  ///  function returns volumetric moisture content
  double give_vol_moist (long ipp);

  /// function returns required non-transport quantities
  void give_reqntq(long *antq);
  
  
  long kd;
  double a1,a2,a3;
  
  
  // nove udelane materialove charakterisitky
  
  long MatChar [20];		// popisuje pro jednotlivy material typ modelu a typ jednotlivych vlastnosti
  //double MatConst [20];     // constanti hodnoty mat. vlastnosti
  //double MatData [20][3][150]; // maximalni pocet 150 radku na jednu charakteristiku
  double MatFunce [20][5];   // charakterisitka zadana nejakou funkci - 5 promennych
  
  gfunct *data[20];
  
  
  double rho_m,rho_w,moist,dmoistdrh;
  double M, R;
  
  //  type of Kunzel model
  //  kunzeltype=1 - relative humidity and temperature
  //  kunzeltype=2 - partial pressure and temperature
  //  kunzeltype=3 - partial pressure and temperature, without liquid moisture
  //  kunzeltype=11 - relative humidity and temperature, hysteresis
  long kunzeltype;
  
  //  type of 
  long kmmtype;  
  
  ///  sorption isotherm
  isotherm isoth;
  ///  influence of damage on permeability
  static dampermeability damper;
  
  ///  density
  gfunct rho;
  ///  porosity
  gfunct por;
  ///  water vapour diffusion resistance factor
  gfunct mu;
  ///  moisture diffusivity
  gfunct kappa;
  ///  moisture diffusivity
  gfunct kappadry;
  ///  sorption isotherm
  isotherm sorpiso;
  ///  desorption isotherm
  isotherm desorpiso;
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
  
  ///  correction coefficient for the sorption isotherm
  double asi;
  ///  correction coefficient for the desorption isotherm
  double adsi;
  ///  correction coefficient for the moisture diffusivity
  double akappa;
  ///  correction coefficient for the moisture diffusivity of drying
  double adkappa;
  
  ///  minimum time before hysteresis
  double time_hyst;
  
  ///  flag for influence of damage on permeability
  flagsw daminfl;
};

#endif
