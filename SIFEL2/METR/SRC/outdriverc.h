#ifndef OUTDRIVERC_H
#define OUTDRIVERC_H

#include <stdio.h>
#include "galias.h"
#include "aliasc.h"
#include "xfile.h"

#define FNAMELEN 1001

class outdriverc;

class outdriverc
{
  public :
    /// constructor
    outdriverc();
  /// destructor
  ~outdriverc();
  /// reads data from file
  long read(XFILE *in);
  /// prints description data to file
  void print(FILE *out);
  /// prints header for output file
  void print_header(FILE *out);
  /// prints header for new step
  void print_newstep(FILE *out, long step, double time);
  /// prints values at nodes
  void print_out(FILE *out, long lcid, long istep, double time);
  /// prints diagrams
  void print_diags(long lcid, double lambda, long istep, double *fi);
  /// prints graphics
  void print_graphics(FILE *out, long lcid, double lambda, long istep, double *fi);
  
  /// output text filename
  char outfn[FNAMELEN];
};

#endif
