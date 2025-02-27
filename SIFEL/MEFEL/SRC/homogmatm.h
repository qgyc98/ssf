#ifndef HOMOGMATM_H
#define HOMOGMATM_H

#include "alias.h"
#include "genfile.h"
#include "iotools.h"
struct matrix;
struct vector;
struct atsel;

/**
  The class homogmatm defines elastic isotropic material model.
  Material parameters are:
  - Young's modulus (modulus of elasticity)
  - Poisson's ratio
   
  Created by TKr,
*/
class homogmatm
{
 public:
  homogmatm (void);
  ~homogmatm (void);
  void read (XFILE *in);
  void print (FILE *out);

  void matstiff (matrix &d,strastrestate ssst);
  void matstiff_plstress (matrix &d);
  void matstiff_plstrain (matrix &d);
  void matstiff_axi (matrix &d);
  void matstiff_spacestr (matrix &d);
  
  void matcompl (matrix &c,strastrestate ssst);

  void assemble_matrices (double *d,long ncomp, long dim);

  void initval(long ipp, long im, long ido);
  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp,long ido);
  double give_strain_vol(long ipp, long ido);
  
  matrix dd;

  long hom_mattypem; //material type for mogenization
  long hom_mattypem_number; //number of material for homogenization

};

#endif
