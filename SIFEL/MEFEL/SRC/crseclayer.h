#ifndef CRSECLAYER_H
#define CRSECLAYER_H

#include "iotools.h"
struct vector;
struct atsel;




/**
  Class crseclayer defines cross section for layered model of planar structure.
   
  Created by JF 10/2014,
*/
class crseclayer
{
 public:
  crseclayer (void);
  ~crseclayer (void);
  
  void read (XFILE *in);
  void zcoordinates ();
  void changeparam (atsel &atcs,vector &val);
  
  ///  number of layers
  long nl;
  
  ///  array with thicknesses of layers
  double *layth;

  /// array with z/coordinates of layers 
  double *layz;
  
  /// thickness of plate
  double th;
  
  ///  optional components
  ///  density of material (necessary for dynamics)
  double rho;
  
  ///  optional parameter
  ///  concentrated mass (only for dynamic problems)
  ///  it should be used only for cross sections defined at nodes
  double m;

};

#endif
