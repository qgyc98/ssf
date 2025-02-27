#ifndef ACI78_H
#define ACI78_H

#include "xfile.h"
struct matrix;

class aci78
/**< 
   This file declares the class for the calculation of the compliance function according to the model ACI 1978
 */
{
 public:
  aci78 ();  
  ~aci78 (); 
  void read (XFILE *in);
  void compliance (double t_current,double &fi_t_t_dash,double &fcyl_t_dash,double &eps_shr_t);
  void matstiff (matrix &d, long ipp, long ido);
  
 private:
  enum curing {
    WATER_CURING=1,
    AIR_CURING=2,
    STEAM_CURING=3
  };
  double t_end_curing;     ///<   time and curing [days]   
  double t_loading;        ///<   age at loading [days]
  double slump;            ///<   slump of fresh concrete [m]
  double density;          ///<   volume density [kg/m3]
  double ratio_ac;         ///<   aggregate/cement ratio
  double ratio_wc;         ///<   water/cement ratio
  double ratio_as;         ///<   aggregate/send ratio
  double humidity;         ///<   relative humidity <0,1> 
  double cs_thickness;     ///<   cross-sectional thickness [m]
  double air_content;      ///<   air content   
  double p6;               ///<   coefficient of correction which is dependent on lab measurement, normally 1.0 
  double fcyl28;           ///<   28-day cylinder strength of concrete [kPa]
  long   curing;           ///<   type of curing
  long   concrete_type;    ///<   type of concrete
  
};

#endif
