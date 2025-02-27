#ifndef EFFSTRESS_H
#define EFFSTRESS_H

#include <stdio.h>
#include "alias.h"
struct XFILE;
class  gfunct;

/**
  This class defines artificial material model which combines
  a damage/(visco)plastictiy/elasticity models with the concept 
  of effective stresses. It can be used also for other types of
  'pore' pressures like a crystallisation pressure of water, salt, etc.
  
  Order of internal variables in the other array :
  -----------------------------------------------
  eqother[0] = other[0] = u_w,0 - initial effective pore pressure (e.g. Sr,w*u_w + Sr,a*u_a)
  eqother[1,...,n], other[1,...,n]  the order of remaining eqother array components
         is given by the used damage/plasticty models
  
  19.9.2012
    04.2022
*/
class effstress
{
 public:
  effstress (void);
  ~effstress (void);

  /// reading of parameters from file
  void read(XFILE *in); 

  /// computation of nonlinera stresses
  void nlstresses (long ipp, long im, long ido); 

  /// initializes mechmat::nonmechq array by constant value of pore pressure
    void initvalues(long lcid, long ipp, long im, long ido, bool rinit);

  /// marks required non-mechanical quantities
  void give_reqnmq(long *anmq);

  pore_press_comp ppct; /// type of pore pressure computation, there can be specified to take values from TRFEL in the coupled problems
  double ppv;           /// prescribed constant pore pressure value assumed for the given material
  double srv;           /// prescribed constant value of degree of saturation assumed for the given material
  gfunct *varppv;       /// prescribed variable pore pressure value assumed for the given material
  gfunct *varsrv;       /// prescribed variable value of degree of saturation assumed for the given material
};

#endif
