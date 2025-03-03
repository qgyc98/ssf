#ifndef CRSEC3DBEAM_H
#define CRSEC3DBEAM_H

#include "xfile.h"
#include "vector.h"
#include <stdio.h>


/**
   class crsec3dbeam defines cross section for 3D beam elements
   
   JK
*/
class crsec3dbeam
{
 public:
  crsec3dbeam (void);
  ~crsec3dbeam (void);
  void read (XFILE *in);
  void print (FILE *out);

  double give_area ();
  double give_ix ();
  double give_iy ();
  double give_iz ();
  void give_moments (double *i);
  void give_shearcoeff (double *sc);
  double give_density ();  

  void changeparam (atsel &atcs,vector &val);
  
  //  necessary components for every problem type
  ///  area of cross section of the beam
  double a;
  ///  moment of inertia of the cross section of the beam
  double ix,iy,iz;
  ///  shear coefficients
  double shearcoeffy,shearcoeffz;
  /// vector which determines local coordinate system on 3d beam element
  vector lcs;
  //  optional components
  ///  density of material (necessary for dynamics)
  double rho;
};

#endif
