#ifndef LHSRHS_H
#define LHSRHS_H

#include <stdio.h>
#include "galias.h"



class lhsrhs
{
 public:
  lhsrhs ();
  ~lhsrhs ();
  void alloc ();
  double *give_lhs (long i);
  double *give_tdlhs (long i);
  double *give_stdlhs (long i);
  double *give_rhs (long i);
  //void output (FILE *out,long lcid);
  void initcond ();
  void clean_lhs ();
  
  ///  dimension of LHS and RHS
  long ndof;
  ///  number of loading cases
  long nlc;
  ///  array containing left hand sides
  double *lhs;
  ///  array containing time derivative of unknowns
  double *tdlhs;
  ///  array containing second time derivative of unknowns
  double *stdlhs;
  ///  array containing initial values of left hand sides
  double *lhsi;
  ///  array containing initial values of time derivative of unknowns
  double *tdlhsi;
  ///  array containing right hand sides
  double *rhs;
  ///  array containing eigenvalues
  double *w;
  /// deallocation flag - in the case of fully coupled problems arrays lhs, tdlhs, lhsi ans rhs may be just references to arrays allocated somewhere else
  answertype deallocf;
};

#endif
