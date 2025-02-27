#ifndef ISOTRMAT_H
#define ISOTRMAT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "aliast.h"
#include "dampermeability.h"

struct vector;
struct atsel;

class isotrmat
{
 public:
  isotrmat (void);
  ~isotrmat (void);
 
  /// computes conductivity %matrix 
  void matcond (matrix &d,long ri,long ci,long ipp);
  /// computes capacity %matrix 
  void matcap (double &cc,long ri,long ci,long ipp);
  /// computes reaction %matrix 
  void matreact (double &r,long ri,long ci,long ipp);
  
  /// computes conductivity %matrix for the 1D problems  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  /// computes conductivity %matrix for the 2D problems  
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  /// computes conductivity %matrix for the 3D problems  
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  /// reads material parameters from the file
  void read (XFILE *in);
  /// prints material parameters to the file
  void print (FILE *out);

  /// returns the conductivity coefficient
  double get_k();
  /// returns the capacity coefficient
  double get_c();
   /// returns the reaction coefficient
  double get_r();
 /// fills dofnames used in the material model
  void give_dof_names(namevart *dofname, long ntm);
  double give_temperature(long ipp);
  void changeparam (atsel &atm,vector &val);
  /// returns volumetric moisture content in the material point
  double give_vol_moist(long ipp);
  /// fills array with non-transport quantities required by the model
  void give_reqntq(long *antq);
  
  double get_othervalue(long compother,double t, long ipp);
  void print_othervalue_name(FILE *out,long compother);

  //  coefficient of conductivity
  double k;
  //  coefficient of capacity
  double c;
  //  coefficient of reaction
  double r;
  
};

#endif
