#ifndef CRSECPLSTR_H
#define CRSECPLSTR_H

#include "iotools.h"
struct vector;
struct atsel;




/**
  Class crsecplstr defines cross section for 2D plane problems.
   
  Creatde by JK,
*/
class crsecplstr
{
 public:
  crsecplstr (void);
  ~crsecplstr (void);
  void read (XFILE *in);
  void print (FILE *out);
  void changeparam (atsel &atcs,vector &val);
  
  ///  necessary components for every problem type
  ///  thickness of the plate or plane stress problem
  double t;
  
  ///  optional components
  ///  density of material (necessary for dynamics)
  double rho;
  
  ///  optional parameter
  ///  concentrated mass (only for dynamic problems)
  ///  it should be used only for cross sections defined at nodes
  double m;

};

#endif
