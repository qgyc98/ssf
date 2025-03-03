#ifndef CRSEC2DBEAM_H
#define CRSEC2DBEAM_H

#include "iotools.h"
struct vector;
struct atsel;

/**
  Class crsec2dbeam defines cross section for 2D beam elements.
   
  Created by JK,
*/

class crsec2dbeam
{
 public:
  crsec2dbeam (void);
  ~crsec2dbeam (void);
  void read (XFILE *in);
  void print (FILE *out);
  
  double give_area ();
  void give_moments (double *i);
  void give_shearcoeff (double *sc);
  double give_density ();
  
  void changeparam (atsel &atcs,vector &val);
  
  //  necessary components for every problem type
  ///  area of cross section of the beam
  double a;
  ///  moment of inertia of the cross section of the beam
  double iy;
  ///  shear coefficient
  double shearcoeff;
  ///  height of the cross section
  ///  it is used in analysis of layered beams
  double t;
  
  //  optional components
  ///  density of material (necessary for dynamics)
  double rho;

};

#endif
