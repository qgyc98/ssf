#ifndef SORPISOHANSEN_H
#define SORPISOHANSEN_H

#include <stdio.h>
#include "genfile.h"

/**
   class contains Hansen sorption isotherm
   
   K.K. Hansen: Sorption Isotherms: A Catalog and a Data Base.
   Water Vapor Transmisson Through Building Materials and Systems:
   Mechanisms and Measurement. ASTM STP 1039
   
   19. 11. 2012, JK
*/
class sorpisohansen
{
 public:
  sorpisohansen (void);
  ~sorpisohansen (void);
 
  void read (XFILE *in);
  void print (FILE *out);

  double derivative_relhum (double rh);
  double hansen_sorption_isotherm (double rh);
  double hansen_inverse_sorption_isotherm (double w);
  
  ///  maximum hygroscopically bound water by adsorption
  double uh;
  ///  empirical fixed exponent
  double n;
  ///  coefficient A = u_n/u_h (u_n - non-evaporable water content)
  double a;
};

#endif
