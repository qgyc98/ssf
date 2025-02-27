#ifndef HISSPLAS_H
#define HISSPLAS_H

#include "iotools.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;


class hissplas
{
 public:
  hissplas(void);
  ~hissplas(void);
  void read(XFILE *in);
  /// yield function (computed from stress invariants)
  double yieldfunction (double i1, double j2, double cos3t, vector &q);
  /// yield function (computed from stress tensor)
  double yieldfunction (vector &sigt, vector &q);
  /// condition of modified flow rule (g1 function)
  double g1function(double depsvp, double depsdp, double i1, double j2, double alpha, double gamma);
  /// derivatives of yield function with respect of stresses
  void deryieldfsigma (long ipp, vector &dfds);
  /// derivatives of yield function with respect to the increment of plastic volumetric strain
  double dfdepsvp (long ipp, double i1, double alpha, double gamma, double ksi, vector &depsp);
  /// derivatives of yield function with respect to the increment of plastic deviatoric strain
  double dfdepspd (long ipp, double i1, double j2, vector &dev, double ksi, vector &depsp);
  /// derivatives of g1 = (df/d[sqrt(j2)])*deps_vp - df/dI1*deps_dp with respect to increment of plastic volumetric strain
  double dg1depsvp (long ipp, double i1, double j2, double alpha, double gamma, vector &dev, double depspd, 
                    double ksi, vector &depsp);
  /// derivatives of g1 = (df/d[sqrt(j2)])*deps_vp - df/dI1*deps_dp with respect to increment of plastic deviatoric strain
  double dg1depsdp (long ipp, double i1, double j2, double alpha, double gamma, vector &dev, double depsvp, 
                    double depsdp, double ksi, vector &depsp);
  /// stiffness matrix for given integration point
  void matstiff (matrix &d,long ipp,long ido);
  /// computes stresses at the given integration point from the nonlocal strain values
  void nonloc_nlstresses (long ipp,long ido);
  /// computes stresses at the given integration point
  void nlstresses (long ipp,long ido);
  /// updates hardening parameters
  void updateq (long ipp, vector &depspt, vector &sigt, vector &depst, vector &q);
  /// updates values of the integration point eqother array with the new reached
  void updateval (long ipp, long ido);
  /// returns plastic strains
  void giveirrstrains (long ipp, long ido, vector &epsp);
  /// returns consistency parameter
  double give_consparam (long ipp,long ido);
  /// function changes material parameters (used in stochastic calculations)
  void changeparam (atsel &atm,vector &val);

  /// shape parameter (beta=0 circular corss-section, beta > 0.0 cross-section -> triangular)
  double beta;
  // Dilation parameter determines apex of the yield surface on the I1-sqrt(J2) plane.
  double n;
  /// triaxial tensile strength
  double r;
  /// atmospheric pressure
  double pa;
  /// coefficient influencing limit value of hardening parameter alpha
  double alpha1;
  /// parameter controlling hardening
  double eta1;
  /// limit value of norm of plastic strains used in hardening
  double ksi_lim;
  /// maximum value of isotropic measure of response surface degradation (gamma)
  double gamma_max;
  /// softening (degradation) parameter
  double k;
  ///  stress return algorithm
  strretalg sra;
};


#endif
