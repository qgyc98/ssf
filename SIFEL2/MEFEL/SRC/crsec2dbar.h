#ifndef CRSEC2DBAR_H
#define CRSEC2DBAR_H

#include "iotools.h"
struct vector;
struct atsel;

/**
  Class crsec2dbar defines cross section for bar elements in 2D.
   
  Created by JK
*/

class crsec2dbar
{
 public:
  crsec2dbar (void);
  ~crsec2dbar (void);
  void read (XFILE *in);
  void print (FILE *out);
  void changeparam (atsel &atcs,vector &val);
 
  //  necessary components for every problem type
  ///  area of cross section of the bar
  double a;

  //  optional components
  ///  density of material (necessary for dynamics)
  double rho;

  ///  optional parameter
  ///  concentrated mass (only for dynamic problems)
  ///  it should be used only for cross sections defined at nodes
  double m;

};

#endif
