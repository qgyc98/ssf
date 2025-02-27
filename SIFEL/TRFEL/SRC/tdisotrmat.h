#ifndef TDISOTRMAT_H
#define TDISOTRMAT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "aliast.h"
#include "gfunct.h"
//#include "dampermeability.h"

struct vector;
struct atsel;

class tdisotrmat
{
 public:
  tdisotrmat (void);
  ~tdisotrmat (void);
 
  /// computes conductivity %matrix 
  void matcond (matrix &d,long ri,long ci,long ipp);
  /// computes capacity %matrix 
  void matcap (double &cc,long ri,long ci,long ipp);
  
  /// computes conductivity %matrix for the 1D problems  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  /// computes conductivity %matrix for the 2D problems  
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  /// computes conductivity %matrix for the 3D problems  
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  /// reads material parameters from the file
  void read (XFILE *in);
  /// prints material parameters to the file
  void print (FILE *out);

  /// returns the conductivity coefficient
  double get_k();
  /// returns the capacity coefficient
  double get_c();
  /// fills dofnames used in the material model
  void give_dof_names(namevart *dofname, long ntm);
  void changeparam (atsel &atm,vector &val);
  /// returns volumetric moisture content in the material point
  double give_vol_moist(long ipp);
  /// fills array with non-transport quantities required by the model
  void give_reqntq(long *antq);

  //  coefficient of conductivity
  gfunct k;
  //  coefficient of capacity
  double c;
};

#endif
