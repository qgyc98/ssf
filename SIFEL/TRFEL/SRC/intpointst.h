#ifndef INTPOINTST_H
#define INTPOINTST_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"
#include "selection.h"

/**
   class intpoints
*/
class intpointst
{
 public:
  intpointst (void);
  ~intpointst (void);
  void read (FILE *in);
  void alloc ();
  void clean ();
  void storegrad (long lcid, double *gradv);
  void copy(intpointst &ip, long ntm, long realloc);
  void save_data_txt    (FILE *aux,sel &selother);
  void restore_data_txt (FILE *aux,long ncompo, sel &selother, long *selid);
  void actual_previous_change ();

  void save_data_bin    (FILE *aux,sel &selother);
  void restore_data_bin (FILE *aux,long ncompo, sel &selother, long *selid);

  ///  material type
  mattypet tm;
  ///  material id
  long idm;

  ///  number of component of flux/grad array
  long ncompgrad;
  ///  number of component of other array
  long ncompother;
  ///  number of component of eqother array
  long ncompeqother;
  
  ///  array of influence of particular gradients
  long **infl;

  ///  array of actual values
  double *av;
  ///  array of previous values
  double *pv;
  ///  array of gradients
  double **grad;
  ///  array of fluxes
  double **fluxes;
  double *cap;
  ///  other components
  double *other;
  ///  eqother components
  double *eqother;
};

#endif
