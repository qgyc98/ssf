#ifndef OUTDIAGT_H
#define OUTDIAGT_H
#include <stdio.h>
#include "galias.h"
#include "xfile.h"
#include "aliast.h"
#include "selection.h"


class outdiagt
{
  public :
  /// constructor
  outdiagt();
  /// destructor
  ~outdiagt();
  /// reads data from the input file
  long read(XFILE *in);
  /// prints data to the input file
  long print(FILE *out);
  /// prints header to the output diagram file
  long print_header(FILE *out);
  /// prints diagram to output diagram file
  long printval(FILE *out, long lcid, double lambda, long istep, double *fi);
  /// the forced print of diagram to output diagram file
  long printval_forced(FILE *out, long lcid, double lambda, long istep, double *fi);
  /// prints unknowns
  long print_unknowns(FILE *out, long lcid, long idp);
  /// prints gradients
  long print_gradients(FILE *out, long lcid, long idp);
  /// prints fluxes
  long print_fluxes(FILE *out, long lcid, long idp);
  /// prints surface fluxes
  long print_surffluxes(FILE *out, long lcid, long idp);
  /// prints others
  long print_others(FILE *out, long lcid, long idp);
  /// prints eqothers
  long print_eqothers(FILE *out, long lcid, long idp);

  /// number of printed unknowns
  long npun;
  /// print step i.e. printout at every pstep-th step
  sel dstep;
  /// pid is node or ipp flag
  nodip *nif;
  /// point id
  long *pid;
  /// element id
  long *eid;
  /// coordinates of selected point 
  double *x, *y, *z;
  /// order of integration point on element
  long *ipeid;
  /// array of type of printed unknowns
  prunkt *pu;
  /// indeces of pu components
  long *ipu;

};

#endif
