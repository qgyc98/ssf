#ifndef INTERFACEMATT_H
#define INTERFACEMATT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "aliast.h"

struct vector;
struct atsel;

class interfacematt
{
 public:
  interfacematt (void);
  ~interfacematt (void);
 
  /// computes conductivity %matrix 
  void matcond (matrix &d,long ri,long ci,long ipp);
  /// computes conductivity %matrix for the 2D problems  
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  /// computes capacity %matrix 
  void matcap (double &cc,long ri,long ci,long ipp);

  /// reads material parameters from the file
  void read (XFILE *in);
  /// prints material parameters to the file
  void print (FILE *out);

  /// returns conductivity along the element
  double get_ks ();
  /// returns conductivity in the normal direction
  double get_kn ();

  /// fills dofnames used in the material model
  void give_dof_names (namevart *dofname, long ntm);
  /// returns volumetric moisture content in the material point
  double give_vol_moist(long ipp);
  /// fills array with non-transport quantities required by the model
  void give_reqntq(long *antq);

  ///  coefficient of conductivity in the direction of element
  double ks;
  ///  coefficient of conductivity in the normal direction
  double kn;
  /**  
    h - fictitious width of the element
    h and kn are correlated and they should be selected properly
  */
  double h;
};

#endif
