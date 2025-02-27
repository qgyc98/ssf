#ifndef DIAGMAT_H
#define DIAGMAT_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"

class precond;

/**
   class diagmat
   
   this class serves for storage of diagonal matrices

   the main puropose of this class is to unify matrix-vector multiplication
   in the class gmatrix

   JK, 25.5.2019
*/

class diagmat
{
 public:
  diagmat (void);
  ~diagmat (void);

  void alloc (long m);
  void dealloc(void);
  void copy (diagmat *dm);
  double* status ();
  long decomp ();
  void changedecomp ();
  void setfact ();
  void setnotfact ();
  void nullmat ();
  void localize (matrix &b,long *cn);
  void localized (double *b,long *cn,long m);
  void initiate (long ndof,long mespr);
  
  void mxv_diag (double *b,double *c);
  void addmat_diag (double c,diagmat &dm);
  void scalmat_diag (double c);
  
  void printmat (FILE *out);
  void printdiag (FILE *out);
  long give_negm ();
  
  double give_entry (long ri);
  
  double estim_spect_radius ();

  void gemp (double *x,double *y,double zero);
  void diagonal_solver (double *x,double *y,double zero);
  
  ///  number of rows (columns)
  long n;
  ///  number of entries in the global %matrix
  long negm;
  ///  array containing global %matrix
  double *a;
  ///  decomposition (factorization) indicator
  long decompid;
  
};

#endif
