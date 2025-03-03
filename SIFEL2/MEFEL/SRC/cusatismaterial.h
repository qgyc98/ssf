#ifndef CUSATISMATERIAL_H
#define CUSATISMATERIAL_H

#include "alias.h"
#include "iotools.h"
struct matrix;
struct vector;
struct atsel;

/**
   The class cusatismaterial defines material model for lattice discrete particle models

   Material parameters are:
   - Young's modulus (modulus of elasticity)
   - Poisson's ratio
   
   JK, 7. 3. 2021
*/
class cusatismaterial
{
 public:
  cusatismaterial (void);
  ~cusatismaterial (void);
  void read (XFILE *in);
  void init();
  void print (FILE *out);

  void matstiff (matrix &d,strastrestate ssst);
  void elmatstiff (matrix &d,strastrestate ssst);
  void matstiff_spacestr (matrix &d);
  
  void matcompl (matrix &c,strastrestate ssst);
  void matcompl_spacestr (matrix &c);

  void nlstresses (long ipp);
  void updateval (long ipp,long ido);
  void changeparam (atsel &atm,vector &val);
  
  ///  Young's modulus
  double e;
  ///  Poisson's ratio
  double nu;

  /// Normal modulus
  double eN;

  /// Normal-Shear coupling
  double alpha;

  /// Shear modulus
  double eT;

  /// Tensile strength
  double ft;
  
  /// Tensile fracture energy 
  double gft;

  /// Tensile characteristic length
  double lt;

  /// Shear strength
  double fs;
  
  /// Shear fracture energy
  double gfs;
  
  /// Shear characteristic length
  double ls; 
  
  /// Tension-shear interaction exponent
  double nt;
  
  /// Initial hardening moduls
  double hc0;
  
  /// Transitional strain ratio
  double kc0;
  
  /// Deviatoric-volumetric strain ratio treshold
  double kc1;
  
  /// Deviatoric damage parameter - reduction of hc for increasing rdv
  double kc2; 
  
  /// Volumetric-deviatoric coupling coefficient
  double beta;
  
  /// Yielding compressive strength
  double sigC0;
  
  /// Densified normal modulus
  double ed;

  /// Initial friction coefficient 
  double m0;
  
  /// Final friction coefficient
  double mInf;
  
  /// Transitional normal stress for friction coefficient
  double sigN0; 


};

#endif
