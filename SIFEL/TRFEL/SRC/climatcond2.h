#ifndef CLIMATCOND2_H
#define CLIMATCOND2_H

#include <stdio.h>
#include "iotools.h"
#include "gfunct.h"
#include "isotherm.h"

class climatcond2
{
 public:
  climatcond2 (void);
  
  ~climatcond2 (void);
  
  ///  function reads data from input files
  void read (XFILE *in, long readclimfile=1);

  ///  function prints data to output file
  void print(FILE *out);

  ///  function reads files with climatic data
  void read_climatic_files (char *Cfile, long cod);

  ///  function prepares the actual values
  void prepare_values (double time);

  ///  function computes values of all external variables
  void compute_values (double time);

  ///  function computes flux densities
  double give_flux (long flid, double *nval,long ipp,long nid);

  ///  function computes the heat flux
  double give_heat_flux (double nval_0,double temp,long ipp, long nid);
  
  ///  function computes short wave radiation flux
  double condit_short_wave_radiat_flux (long swkind, double time);

  ///  function computes long wave radiation flux
  double condit_long_wave_radiat_flux (long lwkind);

   ///  function computes internal energy of convective liquid water flux (J/m2/s)
  double condit_rain_heat_flux (long rainkind,long nid);

  ///  function returns transmission coefficient for heat flux
  double transmission_coefficient_tt (long VDKind);
  
  ///  function computes temperature of exterior
  double transmission_exterior_value_tt (long VDKind);
  
  ///  function computes transmission coefficient tw
  double transmission_coefficient_tw (long VDKind,double temp);


  ///  function computes the moisture flux density
  double give_moisture_flux (double nval_0, double nval_1,long ipp, long nid);
  
  ///  function computes mass of convective liqiud water flux
  double condit_rain_moisture_flux (long rainkind,long nid);
  
  ///  function returns transmission coefficient for moisture flux defined by moisture difference
  double transmission_coefficient_ww (long VDKind);
  
  ///  function computes water vapour pressure in exterior
  double transmission_exterior_value_ww (long VDKind);
  
  ///  function returns external value of the cid-th variable
  double external_nodval (long cid);
  
  ///  function returns transmission coefficient for the lcid-th flux caused by the cid-th gradient
  double transmission_coeff (long lcid,long cid,double temp);


 private:
  
  ///  type of heat conduction
  ///  KIND1 =  0 - no heat conduction
  ///  KIND1 = 11 - prescribed heat flux (on boundary)
  ///  KIND1 = 12 - heat transmission
  ///  KIND1 = 13 - heat transmission extended by the wind factor
  ///  KIND1 = 16 - constant heat transmission
  long KIND1;

  ///  type of short wave radiation heat flux
  ///  KIND2 = 0 - no short wave radiation heat flux
  ///  KIND2 = 21 - short wave radiation heat flux
  ///  KIND2 = 22 - short wave radiation heat flux extended by the wind
  long KIND2;

  ///  type of long wave radiation heat flux
  ///  KIND3 = 0 - no long wave radiation heat flux
  ///  KIND3 = 31 - long wave radiation heat flux
  ///  KIND3 = 32 - long wave radiation heat flux extended by sky radiation
  long KIND3;


  ///  type of water contact
  ///  KIND4 = 0 - no water contact
  long KIND4;

  ///  type of rain - mass of convective liquid water flux
  ///  KIND5 =  0 - no rain
  ///  KIND5 = 51 - prescribed water flux caused by rain
  ///  KIND5 = 52 - prescribed water flux caused by rain extended by wind effect
  long KIND5;


  ///  type of mass of diffusive water vapour flux
  ///  KIND6 = 0 - no mass of diffusive water vapour flux
  ///  KIND6 = 61 - mass of diffusive water vapour flux
  ///  KIND6 = 62 - mass of diffusive water vapour flux - extended by wind effect
  ///  KIND6 = 63 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  ///  KIND6 = 64 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  ///  KIND6 = 65 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by function
  ///  KIND6 = 66 - mass of diffusive water vapour flux - water vapor exchange coefficient is given by a table
  ///  KIND6 = 67 - mass of diffusive water vapour flux - constant temperature, relative humidity, water vapor exchange coefficient
  long KIND6;

  
  ///  input files
  char *files0,*files1,*files2,*files3,*files4,*files5,*files6,*files7,*files8;

  ///  last time when values were computed
  double last_time;
  ///  computer zero
  double zero;

  ///  the number of lines in input file
  long d;
  ///  the number of rows in input files with climatic data
  long tn;
  ///  cc = 0 - the classical referential year
  ///  cc = 1 - measured data with variable length
  long cc;

  ///  type of material model (see alias.h)
  long mattyp;
  ///  wall orientation in radians (rotation around vertical axis)
  double basicdirect;
  ///  wall inclination in radians (rotation around horizontal axis)
  double basicinclin;
  
  ///  the Ludolf number
    //double M_PI;
  ///  water molecular weight
  double M;
  /// universal gas constant
  double R;
  ///  Specific heat capacity of pure liquid water
  double Cwat;
  ///  Specific heat capacity of water vapour
  double cvap;
  ///  Specific enthalpy of water vapour at 0 C
  double hvap;
  ///  Reference temperature for enthalpy (0 C)
  double T0;

  double rhint, rhalfa, rhtemp;
  double tint;

  ///  heat exchange coefficent (W/m^2/K)
  double coefAlf;
  ///  wind factor (-) (greater or equal to zero)
  double WindFact;

  ///  absorption coefficient for short wave radiation (-) (it is between 0 and 1)
  double AbsCoef;
  ///  ground reflextion coefficient
  double Albedo;
  ///  geographical latitude of location in radians
  double Latitude;

  ///  emission coefficient for long wave radiation (-) (it is between 0 and 1)
  double EmisCoeff;

  ///  rain exposure coefficient kg/m^2/s
  double ExchCoeff;
  ///  rain coefficient (it is between 0 and 1)
  double RainCoeff;
  ///  minimum rain temperature from -40 C to 0 C (K)
  double MinRainTemp;
  ///  minimum normal rain intensity (l/m^2/s)
  double MinRainFlux;

  ///  wind factor (-) (greater or equal to zero)
  double WindFactVD;
  /// water vapor exchange coefficient
  double ExchCoeffVD;
  /// coeffcient for function of water vapor exchange coefficient
  double alfaK, Ax, Bx, Cx;

  ///  temperature
  double TepVal;
  ///  relative humidity
  double RHVal;
  ///  diffuse radiation flux on a horizontal plane
  double DifSVal;
  ///  direct radiation flux on a horizontal plane
  double DirSVal;
  ///  wind direction
  double WDVal;
  ///  wind velocity
  double WSVal;
  ///  vertical rain flow density, normal to the ground
  double VRVal;
  ///  long wave radiation emitted by a horizontal surface
  double LWVal;
  ///long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
  double SLWVal;
  /// rain flow density normal to the surface area of the wall
  double NorRainVal;
  /// imposed gaseous mass flux (dry air+ water vapor)
  double GasFluxVal;
  /// imposed heat flux
  double HeatFluxVal;
  /// imposed short wave radiation flux normal to the wall surface
  double ShWRadVal;
  /// imposed long wave radiation flux normal to the wall surface
  double LoWRadVal;
 
  ///  temperature
  gfunct TepValue;
  ///  relative humidity
  gfunct RHValue;
  ///  diffuse radiation flux on a horizontal plane
  gfunct DifSValue;
  ///  direct radiation flux on a horizontal plane
  gfunct DirSValue;
  ///  wind direction
  gfunct WDValue;
  ///  wind velocity
  gfunct WSValue;
  ///  vertical rain flow density, normal to the ground
  gfunct VRValue;
  ///  long wave radiation emitted by a horizontal surface
  gfunct LWValue;
  ///long wave radiation reflected on a horizontal surface by the clouds and particles in the atmosfere
  gfunct SLWValue;
  /// rain flow density normal to the surface area of the wall
  gfunct NorRainValue;
  /// imposed gaseous mass flux (dry air+ water vapor)
  gfunct GasFluxValue;
  /// imposed heat flux
  gfunct HeatFluxValue;
  /// imposed short wave radiation flux normal to the wall surface
  gfunct ShWRadValue;
  /// imposed long wave radiation flux normal to the wall surface
  gfunct LoWRadValue;
  /// water vapor exchange coefficient
  gfunct exchangecoeff;
  

};

#endif
