#ifndef GLASGOWCOUP_H
#define GLASGOWCOUP_H

#include <stdio.h>
#include "alias.h"
#include "genfile.h"

/**
   ordering of variables in eqother array
   
   temperature [K],

*/
class glasgowcoup
{
 public:
  
  glasgowcoup ();
  ~glasgowcoup ();
  void read (XFILE *in);
  
  void material_matrix_fts (long ipp,matrix &d,strastrestate ssst);
  void material_matrix_td (long ipp,matrix &d,strastrestate ssst);
  
  ///  normalized temperature
  double norm_tempr;
  ///  coefficient of free thermal strain
  double alpha;
  ///  the highest reached normalized temperature
  double t_hat;
  
};

#endif
