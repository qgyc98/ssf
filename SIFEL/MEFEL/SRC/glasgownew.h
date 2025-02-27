#ifndef GLASGOWNEW_H
#define GLASGOWNEW_H

#include "iotools.h"
#include "alias.h"

struct matrix;
struct vector;

/**
   ordering of variables in eqother array
   
   temperature [K],

*/
class glasgownew
{
 public:
  
  glasgownew ();
  ~glasgownew ();
  void read (XFILE *in);
  
  void material_matrix_fts (long ipp,matrix &d,strastrestate ssst);
  void material_matrix_td (long ipp,matrix &d,strastrestate ssst);

  void free_thermal_strains (long ipp, vector &epsft);
  double thermdamfunction (long ipp,double tempr,vector &kappa);
  void nlstresses(long ipp, long ido);
  void updateval(long ipp, long ido);
  
  ///  normalized temperature
  double norm_tempr;
  ///  coefficient of free thermal strain
  double alpha;
  ///  the highest reached normalized temperature
  double t_hat;
  
};

#endif
