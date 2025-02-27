#ifndef LINCOUPMAT_H
#define LINCOUPMAT_H

#include <stdio.h>
#include "genfile.h"
struct vector;
struct atsel;

/**
   class %lincoupmat defines an linear material model for coupled
   transport processes
   
   JK, 15.11.2008
*/
class lincoupmat
{
 public:
  lincoupmat (void);
  ~lincoupmat (void);
  
  void read (XFILE *in);
  void print (FILE *out);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);
  
  void matcap (double &cc,long ri,long ci,long ipp);
  
  double give_k (long bri,long bci, long ri,long ci);
  double give_c (long bri,long bci, long ri,long ci);
  
  
  ///  number of transported matters
  long ntm;
  ///  geometrical dimension of the problem
  long dim;

  ///  array of conductivity coefficients
  ///  conductivity %matrix has a block form, there are ntm rows of blocks
  ///  and ntm columns of blocks, each block is a (dim x dim) %matrix
  ///  example in 2D and two matters
  ///   k1  k2     k5  k6
  ///   k3  k4     k7  k8
  ///  
  ///   k9 k10    k13 k14
  ///  k11 k12    k15 k16
  double *k;
  ///  array of capacity coefficients
  ///  storage of coefficient is the same as for the conductivity coeffcients
  double *c;

};

#endif
