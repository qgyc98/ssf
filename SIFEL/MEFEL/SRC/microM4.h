#ifndef MICROM4_H
#define MICROM4_H

#include "iotools.h"
#include "matrix.h"
#include "vector.h"


class microM4
{
 public:
  microM4 (void);
  ~microM4 (void);
  void read (XFILE *in);
  void matstiff (matrix &d,long ipp,long ido);
  void nlstresses (long ipp,long ido);
  void nonloc_nlstresses (long ipp,long ido);
  void updateval(long ipp,long im,long ido);
  void changeparam (atsel &atm,vector &val);

  //material parameters
  long numberOfMicroplanes;
  double e,nu;

 protected:
  vector  microplaneWeights; //(numberOfMicroplanes);
  matrix projN,projM,projL;
  vector kronecker; //(6);
  double k1,k2,k3,k4,c3,c20,c1,c2,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,mu,ev,ed,et;
  void initializeData (long numberOfMicroplanes );
  inline double macbra(double x);
  inline double FVplus(double epsV);
  inline double FVminus (double epsV);
  inline double FDminus(double epsD);
  inline double FDplus(double epsD);
  inline double FN(double epsN,double sigmaV);
  inline double FT (double sigmaN,double epsV);
  inline double maxim (double a,double b);
  inline double minim (double a,double b);
};

#endif
