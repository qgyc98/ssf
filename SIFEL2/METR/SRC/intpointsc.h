#ifndef INTPOINTSTC_H
#define INTPOINTSTC_H

#include <stdio.h>
#include "aliasc.h"
#include "alias.h"
#include "genfile.h"

/**
   class intpointsc
*/

class intpointsc
{
 public:
  intpointsc (void);
  ~intpointsc (void);
  void alloc (long ipp);
  
  ///  material type
  mattypec tm;
  ///  material id
  long idm;
  
  ///  the number of transported media
  long ntm;
  ///  the number of components of gradient/flux
  long ncompgrad;
  ///  array of actual values
  double *av;
  ///  array of previous values
  double *pv;
  ///  array of gradients
  double **grad;
  ///  array of fluxes
  double **fluxes;

  ///  stress-strain state
  strastrestate ssst;
  ///  number of components of displacement
  long ncompdispl;
  ///  number of components of stress/strain array
  long ncompstr;
  ///  number of components of eqother array
  long ncompeqother;
  ///  number of components of other array
  long ncompother;
  ///  array of displacements
  double *displ;
  ///  array of stresses
  double *stress;
  ///  array of strains
  double *strain;
  ///  other components
  double *other;
  ///  equilibriated components of other array
  double *eqother;
  
};

#endif
