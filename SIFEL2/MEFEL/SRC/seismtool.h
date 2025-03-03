#ifndef SEISMTOOL_H
#define SEISMTOOL_H

#include "iotools.h"
#include "alias.h"
#include "genfile.h"

/**
   class seismtool contains tools used in seismic analysis
   
   JK, 21.8.2005
*/

class seismtool
{
 public:
  seismtool (void);
  ~seismtool (void);

  void read (XFILE *in);
  void print (FILE *out);

  void seisminit (double *seism);
  void assemble (double *rhs,double time);
  void assemble (double *v,double *w);
  
  
  ///  number of seismic acceleration components
  long nsac;

  ///  directions of seismic loads
  dirdynload *direction;

  ///  array containing components of right hand side
  double *seism;
  
  ///  amplitude function or response spectrum
  ///  gf describes amplitude functions depending on time in case of step-by-step integration methods
  ///  gf describes response spectrum depending on periods of structure in case of response spectrum analysis
  gfunct *gf;

  // name of file with prescribed set of amplitude functions or response spectra
  char *name;

};

#endif
