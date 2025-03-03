#ifndef ELASTISOMAT_H
#define ELASTISOMAT_H

#include "alias.h"
#include "iotools.h"
struct matrix;
struct vector;
struct atsel;

/**
  The class elastisomat defines elastic isotropic material model.
  Material parameters are:
  - Young's modulus (modulus of elasticity)
  - Poisson's ratio
   
  Created by JK,
*/
class elastisomat
{
 public:
  elastisomat (void);
  ~elastisomat (void);
  void read (XFILE *in);
  void print (FILE *out);

  void matstiff (matrix &d,strastrestate ssst);
  void elmatstiff (matrix &d,strastrestate ssst);
  void matstiff_bar (matrix &d);
  void matstiff_plbeam (matrix &d);
  void matstiff_spacebeam (matrix &d);
  void matstiff_plstress (matrix &d);
  void matstiff_plstrain (matrix &d);
  void matstiff_axi (matrix &d);
  void matstiff_platek (matrix &d);
  void matstiff_plates (matrix &d);
  void matstiff_spacestr (matrix &d);
  
  void matcompl (matrix &c,strastrestate ssst);
  void matcompl_bar (matrix &c);
  void matcompl_plbeam (matrix &c);
  void matcompl_plstress (matrix &c);
  void matcompl_plstrain (matrix &c);
  void matcompl_axi (matrix &c);
  void matcompl_spacestr (matrix &c);

  void initval(long ipp, long im, long ido);
  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp,long ido);
  double give_strain_vol(long ipp, long ido);
  void changeparam (atsel &atm,vector &val);
  
  ///  Young's modulus
  double e;
  ///  Poisson's number
  double nu;
};

#endif
