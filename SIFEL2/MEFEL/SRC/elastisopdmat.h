#ifndef ELASTISOPDMAT_H
#define ELASTISOPDMAT_H

#include "alias.h"
#include "iotools.h"
#include "gfunct.h"
struct matrix;
struct vector;
struct atsel;

/**
  The class elastisopdmat defines elastic isotropic material model with pressure dependent 
  stiffness. It is used in geotechnical problems to simulate increasing elastic modulus
  with mean stress. 
  Material parameters are:
  - initial Young's modulus (modulus of elasticity) 
  - Poisson's ratio
  - mean stress treshold
  - parser expression for the calculation of actual value of elastic modulus
   
  Created by TKo
*/
class elastisopdmat
{
 public:
  elastisopdmat();
  ~elastisopdmat();
  void read (XFILE *in);
  void print (FILE *out);
  void initval(long ipp, long ido);

  void matstiff (matrix &d, long ipp, long ido);
  void elmatstiff (matrix &d, long ipp, long ido);
  void matstiff_bar (matrix &d, double e_c);
  void matstiff_plbeam (matrix &d, double e_c);
  void matstiff_spacebeam (matrix &d, double e_c);
  void matstiff_plstress (matrix &d, double e_c);
  void matstiff_plstrain (matrix &d, double e_c);
  void matstiff_axi (matrix &d, double e_c);
  void matstiff_platek (matrix &d, double e_c);
  void matstiff_plates (matrix &d, double e_c);
  void matstiff_spacestr (matrix &d, double e_c);
  
  void matcompl (matrix &c, long ipp, long ido);
  void matcompl_bar (matrix &c, double e_c);
  void matcompl_plbeam (matrix &c, double e_c);
  void matcompl_plstress (matrix &c, double e_c);
  void matcompl_plstrain (matrix &c, double e_c);
  void matcompl_axi (matrix &c, double e_c);
  void matcompl_spacestr (matrix &c, double e_c);

  void   nlstresses (long ipp, long ido);
  void   updateval (long ipp, long ido);
  double give_actual_ym(long ipp, long ido);
  void   changeparam (atsel &atm,vector &val);
  
  ///  Young's modulus
  double e;
  ///  Poisson's number
  double nu;
  /** mean stress threshold, if the mean stress is greater than the threshold, expression pf
      will be applied to modify elastic moduls */
  double p0;
  /// expression for the determinantion of elastic modulus with respect to actual mean stress
  gfunct pf;  
};

#endif
