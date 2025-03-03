#ifndef PLAST_IFACE_MAT
#define PLAST_IFACE_MAT
#include "strretalg.h"
#include "gfunct.h"
#include "galias.h"
#include <stdio.h>


struct XFILE;
struct matrix;
struct vector;
struct atsel;


/**
  The class defines material model for 2D interface elements which 
  involves the damage model applied to the normal stress component and
  Mohr-Coulomb plasticity model influenced by the damage for the shear 
  stress component.
  The order of array eqother:
  eqother[0] - irreversible part of the normal relative displacement (rdnp)
  eqother[1] - cumulative value of onsistency parameter (gamma)
  eqother[2] - first component of irreversible part of relative displacements in tangential direction (rdsp)
  .
  .
  eqother[1+ncompstr] - last component of the irreversible part of relative displacements in tangential direction (rdsp)

  Created by Tomas Koudelka, koudelka@fsv.cvut.cz, 4.4.2024
*/
class plastifacemat
{
 public:
  plastifacemat();
  ~plastifacemat();
  void read(XFILE *in);
  void print(FILE *out);
  void elmatstiff(matrix &d, long ipp);
  void matstiff(matrix &d, long ipp, long ido);
  void nlstresses(long ipp, long im, long ido);
  void updateval(long ipp, long im, long ido);
  double compute_normal_stress(long ipp, double rdn, double &rdnp);
  double compute_shear_stress(long ipp, double rds, double sigma, double &rdsp, double &gamma);
  double give_consparam(long ipp, long ido);
  void giveirrstrains(long ipp, long ido, vector &rdp);
  void changeparam(atsel &atm, vector &val);
  
  double ks;  /// shear stiffnes coefficient of the material
  double kn;  /// normal stiffness coefficient of the material
  double sig0t;  /// tensile yield stress in the normal direction
  double sig0c;  /// compressive yield stress in the normal direction
  double rdn_lt; /// limit relative displacement in tension in the normal direction
  double rdn_lc; /// limit relative displacement in compression in the normal direction
  double phideg; /// friction angle in degrees
  double c;      /// cohesion in shear
  double tanphi; /// tangent of friction angle
  strretalg sra; /// stress algorithm settings for damage model
};

#endif
