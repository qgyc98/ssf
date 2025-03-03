#ifndef MOHRCPARAB_H
#define MOHRCPARAB_H

#include "iotools.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
  This class defines plastic material model with parabolic Mohr-Coulomb yield criterion.

  Created by Tomas Koudelka,
*/
class mohrcoulombpar
{
 public:
  mohrcoulombpar (void);
  ~mohrcoulombpar (void);
  void read (XFILE *in);
  double yieldfunction (vector &sig);
  void deryieldfsigma (vector &sig,vector &dfds);
  void derplaspotsigma (vector &sig,vector &dfds);
  void matstiff (matrix &d,long ipp,long ido);
  void nlstresses (long ipp, long im, long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  void updateval (long ipp, long im, long ido);
  double give_consparam (long ipp, long ido);
  void changeparam (atsel &atm,vector &val);

  ///  friction angle
  double phi;
  ///  cohesion
  double c;
  ///  dilation
  double psi;
  // Internal parameters
  /// parameter alpha for phi
  double alphaphi;
  /// parameter beta for phi
  double betaphi;
  /// parameter sig cap for phi
  double sigcphi;
  /// parameter alpha for psi
  double alphapsi;
  /// parameter beta for psi
  double betapsi;
  /// parameter sig cap for psi
  double sigcpsi;

  ///  stress return algorithm
  strretalg sra;
};

#endif

