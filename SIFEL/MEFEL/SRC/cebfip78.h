#ifndef CEBFIP78_H
#define CEBFIP78_H

#include "xfile.h"
struct matrix;


class cebfip78
/**< 
   This file declares the class for the calculation of the compliance function according to the model CEB FIP 1978
 */ 
{
 public:
  cebfip78();  
  ~cebfip78(); 
  void read (XFILE *in);
  void compliance (double t_current, double &fi_t_t_dash, double &fcyl_t_dash, double &eps_shr_t);
  void matstiff (matrix &d, long ipp, long ido);
  
 private:
  double t_end_curing; ///<time and curing [days]
  double t_loading;    ///<age at loading [days]
  double humidity;     ///<relative humidity <0,1> 
  double cs_thickness; ///<cross-sectional thickness [m]
  double fcyl28;       ///<28-day cylinder strength of concrete [kPa]
  double p6;           ///<coefficient of correction which is dependent on lab measurement, normally 1.0
  
};

#endif



