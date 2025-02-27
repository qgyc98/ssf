#ifndef TRANSBCLC_H
#define TRANSBCLC_H

#include <stdio.h>
#include "iotools.h"
#include "loadcaset.h"
#include "boundfluxes.h"
#include "genfile.h"

class transbclc
{
 public:
  transbclc (void);
  ~transbclc (void);
  void read (XFILE *in);
  void print (FILE *out);
  void elemsource ();
  
  /// number of load cases
  long nlc;
  
  ///  load cases
  loadcaset *lc;
  
  ///  the number of boundaries with output fluxes
  long nbf;

  ///  output boundary fluxes
  boundfluxes *bf;
  ///  array of fluxes
  double *fluxes;
};

#endif
