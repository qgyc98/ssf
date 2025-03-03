#ifndef GNODVALVM_H
#define GNODVALVM_H



/**
  The class stores pointers to the global vectors with nodal values.
  The dimensions of all %vector stored must be the toatl number of DOFs in the problem, i.e Ndofm.
  The instances of the class are used in functions for the quantity retrieving.
  
  Created by Tomas Koudelka, 09.2023
*/
class gnodvalvm
{
 public:
  gnodvalvm();  
  gnodvalvm(double *d, double *fl, double *fi, double *fr);
  /// %vector of nodal displacements
  double *displv; 
  /// load %vector
  double *loadv; 
  /// internal force %vector
  double *iforv; 
  /// residual %vector
  double *residv; 
};

#endif
