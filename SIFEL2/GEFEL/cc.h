#ifndef COMPCOL_H
#define COMPCOL_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"
#include "csv.h"
#include "densemat.h"

#ifdef INC_CHMD
 #include "cholmod.h"
#endif

class precond;
class diagmat;

/**
   class compcol
   
   serves for %matrix storage called compressed columns
   compressed columns storage is available for real matrices
   only non-zero %matrix entries are stored column by column
   
   basic data
   double array a where %matrix entries are stored
   long array ri where row indices of stored non-zero %matrix entries are stored
   long array adr where the first column entries are stored
   long number n stands for number of columns (rows) of %matrix
   long number negm stands for number of entries of %matrix
   
   basic relationships
   array a has negm components
   array ri has negm components
   array adr has n+1 components
   adr[n] is equal to negm
   
   JK, 24. 3. 2023
*/
class compcol
{
 public:
  compcol ();
  ~compcol ();
  void allocadr (long m);
  double* status ();
  long decomp ();
  void changedecomp ();
  void nullmat ();
  void numcontr_elem (long *cn,long ndofe);
  void numcontr (gtopology *top);
  void fillarray_elem (long *cn,long ndofe);
  void fillarray (gtopology *top);
  void addresses (void);
  void sort_and_colindex (void);
  void localize (matrix &b,long *cn);
  void localized (double *b,long *cn,long nc);
  long minimize (double limit);
  void initiate (gtopology *top,long ndof,long mespr);
  void printmat (FILE *out);
  void printdiag (FILE *out);
  double give_entry (long rri,long ci);
  void add_entry (double e,long rri,long ci);
  
  void mxv_cc (double *b,double *c);
  void addmat_cc (double c, compcol &cc);
  void addmat_diag (double c, diagmat &d);
  void scalmat_cc (double c);
  void solve(double *x, double *y);
 
  void copy_cr (compcol &cc);
  
  
  void fill_data (long nsub,long *smadr,long *smri,double *sma);
  //void select_submatrix (long *li,long nsub,comprow *smcr);
  //void select_submatrix (long *li,long nsub,densemat *smdm);
  void crxcv_cv (compvect *cvi,compvect *cvo);

  long give_negm ();

  double estim_spect_radius ();
  
  
  ///  number of rows of the matrix
  long n;
  ///  number of entries in the global matrix
  long negm;
  ///  array containing addresses of the beginnings of rows
  long *adr;
  ///  auxiliary array containing addresses
  long *adra;
  ///  auxiliary array containing column indices
  long *aux;
  ///  array containing row indices
  long *ri;
  ///  array containing nonzero entries of the matrix
  double *a;
  ///  threshold for entries rejection
  double limit;
  ///  decomposition indicator
  long decompid;

#ifdef INC_CHMD  
  // the following data members relate to the CHOLMOD structures

  /// pointer to the structure with the symbolic decomposition of the given %matrix
  void *symbolic;
  /// pointer to the structure with the decomposed %matrix
  void *numeric;

  /*
  /// structure with sparse matrix stored by compressed columns
  cholmod_sparse  chmd_spr;
  // structure with common data required by CHOLMOD (parameters for symbolic/numeric factorization and update)
  cholmod_common  chmd_com;
  /// pointer to CHOLMOD structure with the sparse Cholesky factorization data (simplicial and supernodal)
  cholmod_factor* pchmd_fact;
  /// dense matrix with left hand sides (vectors of unknowns)
  cholmod_dense   chmd_lhs;
  /// dense matrix with righ hand sides (rhs vectors) 
  cholmod_dense   chmd_rhs;*/
#endif  
};

#endif
