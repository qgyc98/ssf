#ifndef BOERMAT_H
#define BOERMAT_H

#include "xfile.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
  This class implements Boer plasticity model, which is similar to Mohr-Coulomb but
  the singularities on the cone edges are removed by the smoothing. Smoothing is
  controled by the attribute n which has default valu 0.229.
*/
class boermat
{
 public:
  boermat (void);
  ~boermat (void);
  void read (XFILE *in);
  double yieldfunction (vector &sig);
  void deryieldfsigma (vector &sig, vector &dfds);
  void derpotsigma (vector &sig, vector &dgds);
  void matstiff (matrix &d, long ipp,long ido);
  void nlstresses (long ipp,long im,long ido);
  void nonloc_nlstresses (long ipp,long im,long ido);
  void updateval (long ipp,long im,long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  double give_consparam (long ipp,long ido);
  void changeparam (atsel &atm, vector &val);

  ///  friction angle
  double phi;
  ///  cohesion
  double c;
  ///  dilation
  double psi;
  ///  exponent for smoothing
  double n;
  
  
  //  1 - hat
  //  2 - prime
  double alpha,alpha1,alpha2,beta,beta1,delta,a;

  ///  stress return algorithm
  strretalg sra;
};

#endif
