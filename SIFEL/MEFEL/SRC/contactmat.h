#ifndef CONTACTMAT_H
#define CONTACTMAT_H

#include "alias.h"
#include "iotools.h"
struct matrix;
struct vector;
struct atsel;

/**
   class contactmat defines material for contact problems
   
   JK, 11.6.2006
*/

class contactmat
{
 public:
  contactmat (void);
  ~contactmat (void);

  /// read material parameters
  void read (XFILE *in);
  /// assemble stiffness %matrix of the interface material
  void matstiff (matrix &d,long ipp);
  /// computes actual stress values
  void nlstresses (long ipp, long im, long ido);
  /// changes material parameters according to stochastic driver setup
  void changeparam (atsel &atm,vector &val);
  /// updates values of internal variables in eqother array
  void updateval (long ipp, long im, long ido);
  
  ///  stiffness along slip
  double ks;
  ///  stiffness orthogonal to slip
  double kn;
};

#endif
