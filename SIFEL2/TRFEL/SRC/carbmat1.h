#ifndef CARBMAT1_H
#define CARBMAT1_H

#include "genfile.h"

class carbmat1
{
  public:
  carbmat1();    //constructor
  ~carbmat1();   //destructor
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  
  void read(XFILE *in);
  double Equivalent_diffusion(double A);
  double get_conductivity(void);
  
  int global_variable;
  
 private:
  double param;
};

#endif
