#ifndef GLASGOWDAM_H
#define GLASGOWDAM_H

#include "iotools.h"
#include "alias.h"
struct matrix;
struct vector;

/**
   class contains mechanical part of the Glasgow material model of
   hygro-thermo-mechanical problems
   
   strain rate is the sum of rates of elastic, plastic, free thermal, creep
   and load induced thermal strains
   
   stored values in the array eqother:
   previous total stresses, previous total strains, previous total irreversible strains,
   increment of irreversible strains,
   previous temperature, equivalent strain, mechanical damage parameter,
   maximum reached normalized temperature, thermal damage parameter, mechanical damage increment
   
*/
class glasgowdam
{
 public:
  glasgowdam();
  ~glasgowdam();
  void   read(XFILE *in);
  void   damfuncpar(long ipp, vector &eps, vector &kappa);
  double damfunction(long ipp,double tempr,double chi,vector &kappa);
  void   compute_thermdilat(double t, double dt, vector &eps);
  double thermdamfunction(long ipp, double tempr, vector &kappa);
  void   matstiff(matrix &d, long ipp, long ido);
  void   nlstresses(long ipp, long im, long ido);
  void   updateval(long ipp, long im, long ido);

  ///  type of norm of equivalent strains
  paramf_type ft;

  ///  tensile strength
  double st;
  ///  fracture energy at room temperature
  double gf0;
  ///  ratio between compressive and tensile strengths
  double k;
  ///  characteristic length of localization zone
  double lc;
  /// reference temperature
  double reftemp;
};

#endif
