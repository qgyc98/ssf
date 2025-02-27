#ifndef SORPISOROORGRUN_H
#define SORPISOROORGRUN_H

#include <stdio.h>
#include "genfile.h"

/**
   class contains Root (Grunewald) sorption isotherm
   
   Grunewald, J., 2000, DELPHIN 4.1 - Documentation, Theoretical Fundamentals.
   TU Dresden, Dresden
   
   19. 11. 2012, JK
*/
class sorpisorootgrun
{
 public:
  sorpisorootgrun (void);
  ~sorpisorootgrun (void);
 
  void read (XFILE *in);
  void print (FILE *out);

  ///  sorption isotherm value
  double sorption_isotherm (double rh);
  ///  inverse sorption isotherm value
  double inverse_sorption_isotherm (double w);
  ///  derivative of the sorption isotherm with respect to relative humidity
  double derivative_sorption_isotherm (double rh);
  ///  derivative of the inverse sorption isotherm with respect to volumetric moisture content
  double derivative_inverse_sorption_isotherm (double w);
  
  ///  the hygroscopic moisture content
  double whyg;
  ///  the maximum hygroscopic relative humidity
  double rhhyg;
  ///  the saturated volumetric moisture content
  double wsat;
};

#endif
