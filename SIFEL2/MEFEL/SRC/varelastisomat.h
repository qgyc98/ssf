#ifndef VARELASTISOMAT_H
#define VARELASTISOMAT_H

#include "alias.h"
#include "iotools.h"
struct matrix;
struct vector;
struct atsel;

/**
  The class varelastisomat defines elastic isotropic material model with variable 
  Young modulus and Poisson's ratio. These parameters are assumed to be calculated outside of 
  this model, e.g. in the Cam-Clay model. Only initial values of Young modulus and Poisson's ratio 
  are required to be specified in the input file. The current values are calculated by the calling 
  of the mechmat give_actual_ym(ipp) and give_actual_nu(ipp) which are called with the default values 
  of arguments and thus the specific calculation procedure implemented in the control advanced material 
  model (e.g. plasticity) can be invoked.

  Material parameters are:
  - initial Young's modulus (modulus of elasticity)
  - initial Poisson's ratio
   
  Created by TKo, 02.2013,
*/
class varelastisomat
{
 public:
  varelastisomat (void);
  ~varelastisomat (void);
  void read (XFILE *in);
  void print(FILE *in);

  void matstiff (matrix &d, long ipp);
  void elmatstiff (matrix &d, long ipp);
  void elmatstiff (matrix &d, long ipp, strastrestate ssst);
  void matstiff_bar (matrix &d, double e_c);
  void matstiff_plbeam (matrix &d, double e_c, double nu_c);
  void matstiff_spacebeam (matrix &d, double e_c, double nu_c);
  void matstiff_plstress (matrix &d, double e_c, double nu_c);
  void matstiff_plstrain (matrix &d, double e_c, double nu_c);
  void matstiff_axi (matrix &d, double e_c, double nu_c);
  void matstiff_platek (matrix &d, double e_c, double nu_c);
  void matstiff_plates (matrix &d, double e_c, double nu_c);
  void matstiff_spacestr (matrix &d, double e_c, double nu_c);
  
  void matcompl (matrix &c, long ipp);
  void matcompl_bar (matrix &c, double e_c);
  void matcompl_plbeam (matrix &c, double e_c, double nu_c);
  void matcompl_plstress (matrix &c, double e_c, double nu_c);
  void matcompl_plstrain (matrix &c, double e_c, double nu_c);
  void matcompl_axi (matrix &c, double e_c, double nu_c);
  void matcompl_spacestr (matrix &c, double e_c, double nu_c);

  void nlstresses (long ipp);
  void changeparam (atsel &atm,vector &val);
  
  ///  initial Young's modulus
  double e;
  ///  initial Poisson's number
  double nu;
};

#endif
