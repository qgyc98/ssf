#ifndef SALTMAT1_H
#define SALTMAT1_H

#include <stdio.h>
#include "genfile.h"
#include "dampermeability.h"

/**
   class describes material model which deals with simultaneous
   transport of moisture and salt
   
   
   ordering of unknowns:
   w - the volumetric moisture content (m^3/m^3)
   C_f - the concentration of free salts in water (kg/m^3 of solution)
   
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
   

   JM, 29.5.2007, revised 18. 11. 2013
*/
class saltmat1
{
 public:
  saltmat1 (void);
  ~saltmat1 (void);
  void read (XFILE *in);
  void print (FILE *out);

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  
  double k11 (double w);
  double k12 ();
  
  double k21 (double w,double cf);
  double k22 (double w,double cf);

  double c11 ();
  double c12 ();

  double c21 (double cf);
  double c22 (double w,double cf);
  
  
  double transmission_transcoeff (double trc,long ri,long ci,long nid,long bc);
  double transmission_nodval (double nodval,long ri,long ci,long nid,long bc);
  double transmission_flux (double nodval,long ri,long ci,long nid,long bc);
  double get_transmission_flux_ww (double bv,double w,long bc);


  double get_othervalue(long compother,long ipp, double x1,double x2,double x3);
  void print_othervalue_name(FILE *out,long compother);

  void give_dof_names (namevart *dofname, long ntm);

  ///  moisture diffusivity
  gfunct kappa;
  ///  saturated volumetric moisture content
  gfunct sm;
  ///  Dcoef
  gfunct dcoef;
  ///  binding isotherm
  isotherm bindiso;

  ///  influence of damage on permeability
  static dampermeability damper;

  ///  flag for influence of damage on permeability
  flagsw daminfl;
};

#endif

