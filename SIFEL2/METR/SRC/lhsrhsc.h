#ifndef LHSRHSC_H
#define LHSRHSC_H

#include <stdio.h>

class lhsrhsc
{
 public:
  lhsrhsc (void);
  ~lhsrhsc (void);
  void alloc ();

  double *give_lhs (long lcid);
  double *give_lhsi (long lcid);
  double *give_tdlhs (long lcid);
  double *give_rhs (long lcid);

  void initcond (FILE *in);
  
  ///  number of loading cases
  long nlc;
  ///  array containing vector of solution
  double *lhs;
  ///  array containing initial values
  double *lhsi;
  ///  array containing first time derivative of lhs array
  double *tdlhs;
  ///  array containing right hand sides
  double *rhs;
};

#endif
