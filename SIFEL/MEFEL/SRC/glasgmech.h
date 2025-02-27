#ifndef GLASGMECH_H
#define GLASGMECH_H

#include "iotools.h"
#include "alias.h"
struct vector;
struct matrix;

/**
   class contains mechanical part of the Glasgow material model of
   hygro-thermo-mechanical problems
   
   strain rate is the sum of rates of elastic, plastic, free thermal, creep
   and load induced thermal strains
   
   stored values in the array eqother:
   previous total stresses, previous total strains, total irreversible strains,
   previous temperature, equivalent strain, mechanical damage parameter,
   maximum reached normalized temperature, thermal damage parameter
   
*/
class glasgmech
{
 public:
  glasgmech ();
  ~glasgmech ();
  void read (XFILE *in);
  
  double betacoeff (long ipp,double t);
  void damfuncpar(long ipp, vector &eps, vector &kappa);
  double damfunction(long ipp,double tempr,double chi,vector &kappa);
  double domegadkmd (long ipp,double tempr,double chi, vector &kappa);
  double domegadt (long ipp,double tempr,double chi, vector &kappa);
  void matstiff (matrix &d,long ipp,long ido);
  void compute_thermdilat (double t,double dt,matrix &eps);
  double thermdamfunction (long ipp,double tempr,vector &kappa);
  double dtdfdt (long ipp,double tempr,vector &kappa);
  void compute_lits (long ipp,vector &epstm,matrix &sigt,double told,double tnew);
  void depseqdepsel(long ipp, vector &eps, vector &deeqdeel);
  void nlstresses (long ipp,long im,long ido);


  ///  type of norm of equivalent strains
  paramf_type ft;

  ///  thermal damage threshold
  double tkappa0;
  ///  tensile strength
  double st;
  ///  fracture energy at room temperature
  double gf0;
  ///  ratio between compressive and tensile strengths
  double k;
  ///  characteristic length
  double lc;
  ///  parameters of LITS (Load Induced Thermal Strains)
  double a,b,c;
    
};

#endif
