#ifndef MOHRC_H
#define MOHRC_H

#include "iotools.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;



/**
  The class defines nonassociated plastic material model with Mohr-Coulomb yield criterion.
  Single-surface stress return cutting plane algorithm is used for detection of active yield surfaces and 
  if it is necessary, the multi-surface cutting plane algorithm is used.

  Created by Tomas Koudelka,
*/
class mohrcoulomb
{
 public:
  mohrcoulomb (void);
  ~mohrcoulomb (void);
  /// reads material parameters from the file
  void read (XFILE *in);
  /// prints material parameters to the file
  void print (FILE *out);
  /** yield function from principal stresses for standard stress ordering 
      psig[0] < psig[1] < psig[2] */
  double yieldfunction (vector &psig);
  /** derivatives of yield function with respect of  principal stresses for standard stress ordering 
      sig1 < sig2 < sig3 */
  void deryieldfsigma (vector &dfds);
  /** derivatives of plastic potential function with respect of  principal stresses for standard stress ordering 
      sig1 < sig2 < sig3 */
  void derplaspotsigma (vector &dgds);
  /// stiffness matrix for given intgration point
  void matstiff (matrix &d,long ipp,long ido);
  /// tangent stiffness matrix for given integration point
  void tangentstiff (matrix &d,matrix &td,long ipp,long ido);
  /// elastic stiffness matrix for elastic isotropic material computed at principal stress directions
  void pelmatstiff (long ipp, matrix &d, double e, double nu);
  /// computes stresses at the given integration point 
  void nlstresses (long ipp, long ido);
  /// computes stresses at the given integration point from the nonlocal strain values
  void nonloc_nlstresses (long ipp,long ido);
  /// cutting plane algorithm at principal stresses
  long cutting_plane(long ipp, double &gamma, vector &epsn, vector &epsp, vector &q, long ni, double err);
  /// returns hardening contribution for cutting_plane method
  double plasmodscalar(long ipp, vector &q);
  /// updates hardening parameters
  void updateq(long ipp, vector &epsp, vector &q);
  /// checks type of singularity region and returns singularity indicator
  long checkpsig(vector &psig);
  /// checks for index of out of plane principal stress in case plane stress problem
  long checkzeropsig(vector &psig);
  /// multisurface cutting plane algorithm at principal stresses
  void mc_msurf_cp(long ipp,double &gamma,vector &epsn,vector &epsp,vector &q, long mu,long ni,double err);
  /// returns hardening contribution for multisurface cutting plane method
  void plasmodscalar(long ipp, vector &q, long mu, matrix &hcpm);
  /** yield function from principal stresses for standard stress ordering and
      ordering given by mu and stat parameter */
  void yieldfunction(vector &psig, long mu, long *stat, vector &f);
  /** yield function from principal stresses for standard stress ordering and
      ordering given by mu parameter */
  void yieldfunction(vector &psig, long mu, vector &f);
  /** derivatives of yield function with respect of  principal stresses for standard stress ordering and
      ordering given by mu parameter */
  void dfdsigma(long *stat, long mu, matrix &dfds);
  /** derivatives of plastic potential function with respect of  principal stresses for standard stress ordering and
      ordering given by mu parameter */
  void dgdsigma(long *stat, long mu, matrix &dgds);
  /// updates values of the integration point eqother array with the new reached
  void updateval (long ipp,long im, long ido);
  /// returns plastic strains
  void giveirrstrains (long ipp, long ido, vector &epsp);
  /// returns consistency parameter
  double give_consparam (long ipp,long ido);
  /// function changes material parameters (used in stochastic calculations)
  void changeparam (atsel &atm,vector &val);
  
  ///  friction angle
  double phi;
  ///  cohesion
  double c;
  ///  dilatation
  double psi;

  ///  stress return algorithm
  strretalg sra;
};

#endif
