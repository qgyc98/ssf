#ifndef LHSRHST_H
#define LHSRHST_H

#include <stdio.h>
#include "galias.h"
#include "iotools.h"

/**
   class lhsrhst defines left and right hand side %vectors of transport problems
*/

class lhsrhst
{
 public:
  lhsrhst (void);
  ~lhsrhst (void);
  void alloc ();
  double *give_lhs (long lcid);
  double *give_lhsi (long lcid);
  double *give_tdlhs (long lcid);
  double *give_rhs (long lcid);
  void initcond (XFILE *in);
  void initcondprint (FILE *out);
  
  ///  number of loading cases
  long nlc;
  ///  array containing left hand sides
  double *lhs;
  ///  array containing initial values
  double *lhsi;
  ///  array containing first time derivative of lhs array
  double *tdlhs;
  ///  array containing right hand sides
  double *rhs;
  /// deallocation flag - in the case of fully coupled problems arrays lhs, tdlhs, lhsi ans rhs may be just references to arrays allocated somewhere else
  answertype deallocf;
};

#endif
