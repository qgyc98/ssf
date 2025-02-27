#ifndef saltmat4_H
#define saltmat4_H

#include <stdio.h>
#include "genfile.h"
#include "dampermeability.h"

/**
   class describes material model which deals with simultaneous
   transport of heat, moisture, salt and salt crystals
   
   
   ordering of unknowns:
   w - the volumetric moisture content (m^3/m^3)
   C_f - the concentration of free salts in water (kg/m^3 of solution)
   C_c - the amount of crystallized salt (kg/m^3 of sample)
   T - temperature (K)
   
   ordering in the array eqother
   eqother[0] - the water vapour diffusion permeability
   eqother[1] - relative humidity
   eqother[2] - derivative of the relative humidity with repsect to the moisture content
   eqother[3] - saturated volumetric moisture content
   eqother[4] - maximum concentration
   eqother[5] - total salt concentration



   components in CORD:
   2 - density
   3 - porosity
   4 - faktor difusniho odporu
   5 - kappa
   6 - sorption izoterm
   7 - saturated moisture
   8 - none
   9 - Cecko
   10 - Lambda - the thermal conductivity (W/m/K)
   11 - not used
   12 - not used
   13 - not used
   14 - Dcoef - the salt diffusion coefficient (m$^2$/s)
   15 - binding isotherm
   16 - cfmax - the saturated free salt concentration (kg/m$^3$ of solution)
   17 - ws
   18 - not used
   19 - not used
   

   JM, 29.5.2007, revised 2. 10. 2013
*/
class saltmat4
{
 public:
  saltmat4 (void);    //constructor
  ~saltmat4 (void);   //destructor
  
  void read (XFILE *in);
  void print(FILE *out);


  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcond2 (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  void matcond4d (matrix &d,long ri,long ci,long ipp);
  void matcond2d2 (matrix &d,long ri,long ci,long ipp);


  double k11 (long ipp);
  double k12 ();
  double k13 ();
  double k14 (long ipp);

  double k21 (long ipp);
  double k22 (long ipp);
  double k23 ();
  double k24 ();

  double k31 ();
  double k32 ();
  double k33 ();
  double k34 ();

  double k41 (long ipp);
  double k42 ();
  double k43 ();
  double k44 (long ipp);
  
  
  double c11 (long ipp);
  double c12 ();
  double c13 ();
  double c14 ();

  double c21 (long ipp);
  double c22 (long ipp);
  double c23 ();
  double c24 ();

  double c31 (long ipp);
  double c32 (long ipp);
  double c33 ();
  double c34 ();

  double c41 ();
  double c42 ();
  double c43 ();
  double c44 (long ipp);

  
  double transmission_transcoeff(double trc,long ri,long ci,long nid,long bc);
  //double get_transmission_transcoeff_11(double w,double cf,double cc,double t,long bc,long ipp);

  double transmission_nodval(double nodval,long ri,long ci,long nid,long bc);
  //double get_transmission_nodval_11(double bv,double w,double cf,double cc,double t,long bc,long ipp);

  double transmission_flux(double nodval,long ri,long ci,long nid,long bc);
  double get_transmission_flux_ww (double bv,double w,long bc);

  
  double get_othervalue(long compother,long ipp, double w,double cf,double cc,double t);
  void print_othervalue_name(FILE *out,long compother);
  
  
  
  ///  pressure of saturated water vapour
  double pgws (double t);
  ///  water vapour diffusion permeability
  double permeabilitavodnipary(double w,double t);
  ///  water density
  double water_density ();
  ///  derivative of saturation water vapour pressure with respect to temperature
  double derivative_saturation_water_vapour_pressure_temperature(double t);
  ///  latent heat of evaporation of water (J/kg)
  double latent_heat_of_evaporation_of_water(double t);


  ///  derivative of C_b with respect to C_f
    //double binding_isotherm_derivative (double cf);
  
  ///  relative humidity computed from the volumetric moisture content
    //double relative_humidity (double w);
  







  void inverze_sorption_isotherm_data(double w, double &fi, double &dfi);
  void sisotherm(int kod, double w,double &fiw, double &dfdw);
  double linear_data(int kod, double w,double cf, double cc,double t);
  double diffcoefiont(int kod, double w,double cf, double cc,double t);
  double mw, ma, gasr;
  void values_correction (vector & nv,long ipp);
  void Cf_check(double cf, double w, long ipp);
  

  void salt_diffusivity_values (int kod, long ipp, double cf, double xpv,double ineq1, double &diff);
  
  
  void water_content_relhum (long nid,double *in,double *inp,double *ineq, double *out);
  double inverse_hystereze_sorption_isotherms (long ipp, double in, double inp, double ineq);
  void get_moisture (long nid,double in,double *inp, double &out);


  void der_value_hyst (int matchar,int kod, double pv, double & outvalue,double & outvalue2, long ipp);
  



  long cycle_detection (double *r,double *pr, double *ppr);
  
  void hystereze2 (int matchar, double x, double xpv, double ineq1,double & outvalue,double & outvalue2, long ipp, long timeH);
  double get_moisture2 ( double rh);
  void get_rel_hum2 (double w, double &fi, double &dfdw);
  void read_Sourcet(XFILE *in);
  double getval_source (double t);

  void give_values (long ipp,double *av, double *pv, double *eq);
  void aux_values (long ipp,double *in,double *inp,double *ine,double *out);
  void save_values (long ipp,double *out);

  void give_reqntq(long *antq);

  ///  function returns temperature in integration point
  double give_temperature (long ipp);
  
  ///  function returns initial temperature in integration point
  double give_inittemperature (long ipp);

  ///  function returns relative humidity in integration point
  double give_rel_hum (long ipp);

  ///  function returns pore pressure in integration point
  double give_pore_pressure (long ipp);


  double **source;
  long pocet_radku;

  //  general function for material parameters previously stored in MatData
  gfunct *data[20];

  
  //void CorD (int cislochar, double in,int rhw, double x, double &y, double &z, double &z2);

  ///  sorption isotherm
  isotherm isoth;
  ///  influence of damage on permeability
  dampermeability damper;
  
  
  ///  density
  gfunct rho;
  ///  porosity
  gfunct por;
  ///  water vapour diffusion resistance factor
  gfunct mu;
  ///  moisture diffusivity
  gfunct kappa;
  ///  sorption isotherm
  isotherm sorpiso;
  ///  saturated volumetric moisture content
  gfunct sm;
  ///  specific heat capacity
  gfunct c;
  ///  thermal conductivity
  gfunct lambda;
  ///  Dcoef
  gfunct dcoef;
  ///  binding isotherm
  isotherm bindiso;
  ///  cfmax
  gfunct cfmax;
  ///  ws
  gfunct ws;
  
  ///  parameter for the generalized Heaviside function
  double eps;

  ///  influence of damage on permeability
  flagsw daminfl;
};

#endif

