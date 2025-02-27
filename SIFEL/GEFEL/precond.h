#ifndef PRECOND_H
#define PRECOND_H

#include <stdio.h>
#include "galias.h"

class gmatrix;
class aggregator;
class gtopology;
struct XFILE;

/**
   class deals with preconditioning
   
   JK, 25.2.2007
*/
class precond
{
 public:
  precond ();
  ~precond ();
  
  void read (gtopology *gt,XFILE *in,long mespr);
  void print (FILE *out);

  void setting (precondtype p,double ift,long numaggr);
  
  void initiation (gtopology *gt,gmatrix *gm,FILE *out);

  void preconditioning (double *r,double *h,double *rhs);
  
  
  ///  type of preconditioning
  precondtype pt;
  
  ///  number of unknowns/equations
  long n;
    
  ///  treshold for incomplete factorization
  double incompltresh;

  ///  class for aggregation method
  aggregator *agg;
  
  
  gmatrix *agm;
  
  double ssoromega;
  double indegamma;
};


#endif
