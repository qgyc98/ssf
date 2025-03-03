#ifndef PARDISO_H
#define PARDISO_H

/**
   class defining Intel sparse direct solver
   
   IS, 7.12.2005
*/
class pardiso  
{
 public:
  
  pardiso ();
  ~pardiso ();

  void symbfact (double *a,long *ci,long *adr,long ndof);
  void numfact (double *a,long *ci,long *adr,long ndof);
  void backsubst (double *a,long *ci,long *adr,long ndof,double *x,double *y);
  
  pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*,int*, int*, int*, double*, double*, int*);

  //void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*,int*, int*, int*, double*, double*, int*);
  
  ///  number of rows
  int nRows;
  ///  number of columns
  int nCols;
  ///  number of nonzeros
  int nNonZeros;
  ///  number of right hand sides
  int nrhs;

  //double *rhs,*solValues;

  ///  type of solved matrix
  ///  mtype=11 - real unsymmetric matrix
  int mtype;
  
  ///  internal solver memory pointer pt
  ///  32-bit: int pt[64]; 64-bit: long int pt[64]
  ///  or void *pt[64] should be OK on both architectures
  void *pt[64];

  /// pardiso control parameters
  int iparm[64];
  int maxfct, mnum, phase, error, msglvl;
  ///  auxiliary variables. */
  double ddum; /* Double dummy */
  int idum; /* Integer dummy. */
  

};

#endif
