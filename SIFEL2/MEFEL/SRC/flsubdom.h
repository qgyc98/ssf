#ifndef FLSUBDOM_H
#define FLSUBDOM_H

#include <stdio.h>
#include "alias.h"
#include "gmatrix.h"


/**
   class flsubdom serves for problems with floating subdomains
   solved by the FETI method
   
   JK, 20.6.2006, modified 3. 8. 2012
*/

class flsubdom
{
 public:
  flsubdom ();
  ~flsubdom ();
  
  void initiation (gmatrix *gm,FILE *out);
  void solve_lin_alg_system (double *lhs,double *rhs,FILE *out);

  ///  the %vector of Lagrange multipliers
  ///  it contains increments of multipliers in nonlinear problems
  double *lambda;
  ///  the %vector of total Lagrange multipliers
  ///  it is used only in nonlinear problems
  double *totlambda;
  
  ///  the number of changes in the %vector of Lagrange multipliers
  long nch;
  ///  compliances
  double *compli;

};

#endif
