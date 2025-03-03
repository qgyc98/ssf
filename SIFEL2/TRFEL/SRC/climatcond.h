#ifndef CLIMATCOND_H
#define CLIMATCOND_H

#include <stdio.h>
#include "iotools.h"
#include "gfunct.h"
#include "isotherm.h"

# ifndef M_PI
#   define M_PI 3.14159265358979323846  
# endif

class climatcond
{
 public:
  climatcond(void);
  ~climatcond(void);
  void read(XFILE *in, long readclimfile=1);
  void print(FILE *out);
  void read_climatic_files (char * Cfile, long cod);
  
  double give_flux(long icid, double * nval,long ipp,long nid);

  double give_heat_flux(double nval_0,double nval_1,long ipp, long nid);
  double condit_heat_flux(long tkind, double temperstr);
  double condit_short_wave_radiat_flux(long swkind,double time);
  double condit_long_wave_radiat_flux (long lwkind);
  double condit_rain_heat_flux(long rainkind, double nval_0,long nid);
  double condit_vapour_diffusion_heat_flux(long vdkind,double nval_0 , double nval_1,long ipp, long nid);


  double give_moisture_flux(double nval_0, double nval_1,long ipp, long nid);
  double condit_water_contact ();
  double condit_rain_moisture_flux(long rainkind, double nval_0,long nid);
  double condit_vapour_diffusion_moisture_flux(long VDKind,double nval_0, double nval_1,long ipp, long nid);



  double give_salt_flux(double w, double cf, double t, long nid);


  void compute_values (double time);
  void prepare_values (double time);
  
 
  double condit_vapour_diffusion_moisture_flux_salt(long VDKind,double w, double cf, double temp, long nid);



  long MatCharC[3];
  double MatConstC[3];
  long numd;
  

 private:
  
  ///  type of heat conduction
  ///  KIND1 =  0 - no heat conduction
  ///  KIND1 = 11 - prescribed heat flux (on boundary)
  ///  KIND1 = 12 - heat transmission
  ///  KIND1 = 13 - heat transmission extended by the wind factor
  ///  KIND1 = 14 - JK+Tomaschko
  ///  KIND1 = 15 - JK+Klier
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


  ///  last time when values were computed
  double last_time;
  ///  computer zero
  double zero;
  ///  the Ludolf number
    //double M_PI;


  ///  heat exchange coefficent (W/m^2/K)
  double coefAlf;
  ///  wind factor (-) (greater or equal to zero)
  double WindFact;
  ///  coefficients for computation of transmission coefficent
  double aa,b;
  ///  coefficients for Stefan-Boltzmann law
  double eps,sig;
  ///  coefficients for radiation
  double alpha,m;

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

  double rhint, rhalfa, rhtemp;
  double tint;
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
  ///  array containing radiation (JK+Tomaschko)
  gfunct RadValue;
  ///  array containing temperature Tb (JK+Klier)
  gfunct TempbValue;
  
  /// water vapor exchange coefficient
  gfunct exchangecoeff;
  
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
  ///  array containing radiation (JK+Tomaschko)
  double RadVal;
  ///  array containing temperature Tb (JK+Klier)
  double TempbVal;


  ///  type of material model (see alias.h)
  long mattyp;
  ///  wall orientation in radians (rotation around vertical axis)
  double basicdirect;
  ///  wall inclination in radians (rotation around horizontal axis)
  double basicinclin;




	double Jkq;
	double Cwat;
	double T0;
	double hvap;
	double cvap;

	//.... toto musim nacist ze vstup souboru

	long WCorPH;//  pouziva se to, ale necte se to!
	long RHorVP;

	double cfmax;

	double u, a, n;
	//double tint; // pro interier HEMOT
	//double rhint,rhtemp; // pro interier HEMOT

  //    double rhint, rhalfa, rhtemp;

    double **SIdata;

    long SI, numberSI;
    long d;
    long tn,cc; // tn je pocet radku climatickych hodnot, cc je 0 pro TRY  a 1 pro namerena klim data ruzne delky

    char *files0,*files1,*files2,*files3,*files4,*files5,*files6,*files7,*files8,*files14,*files15;
    
    
};


#endif
