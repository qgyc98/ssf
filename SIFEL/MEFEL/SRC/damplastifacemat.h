#ifndef DAMPLAST_IFACE_MAT
#define DAMPLAST_IFACE_MAT
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
  eqother[0] - maximum attained relative normal displacement in history (kappa)
  eqother[1] - damage parameter (omega)
  eqother[2] - actual crack width (w)
  eqother[3] - cumulative value of onsistency parameter (gamma)
  eqother[4] - cumulative irreversible part of relative displacements in tangential direction (dup)

  Created by Tomas Koudelka, koudelka@cml.fsv.cvut.cz, 17.12.2015
*/
class damplastifacemat
{
 public:
  damplastifacemat();
  ~damplastifacemat();
  void read(XFILE *in);
  void print(FILE *out);
  void elmatstiff(matrix &d, long ipp);
  void matstiff(matrix &d, long ipp, long ido);
  void nlstresses(long ipp, long im, long ido);
  void updateval(long ipp, long im, long ido);
  void changeparam(atsel &atm, vector &val);
  double compute_normal_stress(long ipp, double dv, double aft, double &kappa, double &omega, double &w);
  double compute_shear_stress(long ipp, double du, double aft, double omega, double kappa, double sigma, double &dup, double &gamma);
  double give_actual_ft();
  double epsefunction(long ipp);
  double givedamage(long ipp, long ido);
  double give_crackwidth(long ipp, long ido);
  double give_consparam(long ipp, long ido);
  void giveirrstrains(long ipp, long ido, vector &epsp);
  
  double ks;  /// shear stiffnes coefficient of the material
  double kn;  /// normal stiffness coefficient of the material
  double ft;  /// tensile strength of the interface (in normal direction)
  double wf;  /// initial slope of the softening branch
  double phideg; /// friction angle in degrees
  double c;      /// cohesion in shear
  double tanphi; /// tangent of friction angle
  strretalg sra; /// stress algorithm settings for damage model

  answertype coref; /// flag for the treatment of the corrosion effect
  
  gfunct alpha_c;       /// evolution function of the expansion factor for corrosion products

  gfunct icor_func;     /// evolution function of the corrosion current density [uA/cm^2]

  gfunct flux_cor_func; /// evolution function of the corrosion product flux
};

#endif
