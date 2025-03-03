#ifndef SYMCOMPROW_H
#define SYMCOMPROW_H

#include "galias.h"
#include "matrix.h"
#include "vector.h"
#include "gtopology.h"
#include "SPARSE/DSSolver.h"
class gmatrix;

/**
   class symcomprow
   
   serves for %matrix storage called symmetric compressed rows
   compressed rows storage is available for symmetric real matrices
   only non-zero matrix entries are stored row by row from the first entry in the row to the diagonal one, i.e. 
   lower triangular %matrix is stored
   
   basic data
   double array a where matrix entries are stored
   long array ci where column indices of stored non-zero matrix entries are stored
   long array adr where the first row entries are stored
   long number n stands for number of rows (columns) of matrix
   long number negm stands for number of entries of matrix
   
   basic relations
   array a has negm components
   array ci has negm components
   array adr has n+1 components
   adr[n] is equal to negm
   
   JK
*/
class symcomprow
{
 public:
  symcomprow (void);
  ~symcomprow (void);
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
  void sort_and_colindex_tko (void);
  void localize (matrix &b,long *cn);
  void localized (double *b,long *cn,long nc);
  long minimize (double limit);
  void initiate (gtopology *top,long ndof,long mespr);
  
  double give_entry (long ri,long rci);
  void add_entry (double e,long ri,long rci);
  
  void mxv_scr (double *b,double *c);
  void addmat_scr (double c,symcomprow &scr);
  void addmat_gm (double c,gmatrix &gm);
  void scalmat_scr (double c);
  void cg (double *x,double *y,long ni,double res,long &ani,double &ares,double limit,long iv);

  void prec_diag_scr (double *x,double *y);
  void prec_ssor_scr (double *x,double *y,double omega);
  double rows_mult_ldl_scr (long m,long o);
  void incomplete_ldl (double gamma);
  void prec_ildl_scr (double *x,double *y);

  void cg_prec (double *x,double *y,long ni,double res,long &ani,double &ares,
		double zero,long iv,long tprec,double par,ISolver *sdirect);
  
  void printmat (FILE *out);
  void printdiag (FILE *out);
  void copy_scr (symcomprow &scr);

  long give_negm ();
  
  void select_submatrix(symcomprow *smscr,long nsdof,long *sdof);


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
  long *ci;
  ///  array containing nonzero entries of the %matrix
  double *a;
  ///  array containing incomplete decomposition
  double *incdec;
  ///  threshold for entries rejection
  double limit;
  ///  decomposition indicator
  long decompid;
};

#endif
