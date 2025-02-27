#ifndef SYMCOMPCOL_H
#define SYMCOMPCOL_H

#include "galias.h"
#include "matrix.h"
#include "vector.h"
#include "gtopology.h"

#ifdef INC_CHMD
 #include "cholmod.h"
#endif

class diagmat;

/**
   class symcompcol
   
   serves for %matrix storage called symmetric compressed rows
   compressed rows storage is available for symmetric real matrices
   only non-zero matrix entries are stored row by row from the first entry in the row to the diagonal one
   
   basic data
   double array a where matrix entries are stored
   long array ri where row indices of stored non-zero matrix entries are stored
   long array adr where the first column entries are stored
   long number n stands for number of rows (columns) of matrix
   long number negm stands for number of entries of matrix
   
   basic relations
   array a has negm components
   array ci has negm components
   array adr has n+1 components
   adr[n] is equal to negm
   
   JK
*/
class symcompcol
{
 public:
  symcompcol(void);

  ~symcompcol(void);

  /// function returns allocation status of array s
  const double* status() {return a;};

  ///  function returns indicator of decomposition (factorization)
  long decomp() {return decompid;};

  /// function changes indicator of decomposition (factorization)
  void changedecomp(){decompid == 0 ? decompid=1 : decompid = 0;};

  /// function fills array a by zeros
  void nullmat() { decompid = 0; memset(a, 0, negm*sizeof(*a));};

  void allocadr(long m);

  void numcontr_elem(const ivector &cn);

  void numcontr(gtopology *top);

  void addresses(void);

  void fillarray_elem(const ivector &cn);

  void fillarray(gtopology *top);

  void sort_and_colindex(void);

  void sort_and_colindex_tko (void);

  long initiate (gtopology *top, long ndof, long mespr);
  
  void localize (matrix &b, long *cn);

  void localized(double *b, long *cn, long nc);

  long minimize(double limit);

  const char*chmd_ord_method(int i);
  
  double give_entry(long ir, long ic);

  void add_entry(double e,long ir,long ic);
  
  void factorize(void);

  void solve(double *x, double *y);

  void mxv_scc(double *b, double *c);

  void addmat_scc(double c, symcompcol &scc);

  void addmat_diag (double c, diagmat &d);

  void scalmat_scc(double c);

  void printmat(FILE *out);

  void printdiag(FILE *out);

  void copy_scc(symcompcol &scc);

  long give_negm();


  ///  number of rows of the %matrix
  long n;
  ///  number of entries in the global %matrix
  long negm;
  ///  array containing addresses of the beginnings of rows
  long *adr;
  ///  auxiliary array containing addresses
  long *adra;
  ///  auxiliary array containing column indices
  long *aux;
  ///  array containing column indices
  long *ri;
  ///  array containing nonzero entries of the %matrix
  double *a;
  ///  threshold for entries rejection
  double limit;
  ///  decomposition indicator
  long decompid;

#ifdef INC_CHMD  
  // the following data members relate to the CHOLMOD structures

  /// structure with sparse matrix stored by compressed columns
  cholmod_sparse* chmd_spr;
  // structure with common data required by CHOLMOD (parameters for symbolic/numeric factorization and update)
  cholmod_common  chmd_com;
  /// pointer to CHOLMOD structure with the sparse Cholesky factorization data (simplicial and supernodal)
  cholmod_factor* chmd_fact;
  /// dense matrix with left hand sides (vectors of unknowns)
  cholmod_dense*  chmd_lhs;
  /// dense matrix with righ hand sides (rhs vectors) 
  cholmod_dense   chmd_rhs;
#endif  
};

#endif
