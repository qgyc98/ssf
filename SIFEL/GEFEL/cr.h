#ifndef COMPROW_H
#define COMPROW_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"
#include "csv.h"
#include "densemat.h"

class precond;

/**
   class comprow
   
   serves for %matrix storage called compressed rows
   compressed rows storage is available for real matrices
   only non-zero %matrix entries are stored row by row
   
   basic data
   double array a where %matrix entries are stored
   long array ci where column indices of stored non-zero %matrix entries are stored
   long array adr where the first row entries are stored
   long number n stands for number of rows (columns) of %matrix
   long number negm stands for number of entries of %matrix
   
   basic relationships
   array a has negm components
   array ci has negm components
   array adr has n+1 components
   adr[n] is equal to negm
   
   JK
*/
class comprow
{
 public:
  comprow ();
  ~comprow ();
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
  double give_entry (long ri,long rci);
  void add_entry (double e,long ri,long rci);
  
  void mxv_cr (double *b,double *c);
  void mxv_cr_new (double *b,double *c);
  double mxv_cr_new2 (double *b,double *c);
  double mxv_cr_pom (double *b1,double *c1,double *b2,double *c2);
  
  void mv_cr15 (double *b,double *c);

  void mtxv_cr (double *b,double *c);
  void addmat_cr (double c,comprow &cr);
  void scalmat_cr (double c);
  void cg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv);
  void cg_new (double *x,double *y,
	       long ni,double res,long &ani,double &ares,double zero,long iv);

  void cg_prec (precond &pr,double *x,double *y,
		long ni,double res,long &ani,double &ares,double zero,long iv);

  void bicg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv);
  void bicg_new (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv);
  
  double mv_cr15_rev (double *b,double *c);
  void cg_cr_rev (double *x,double *y,
		  long ni,double res,long &ani,double &ares,double limit,long iv);
  
  void copy_cr (comprow &cr);
  
  
  void fill_data (long nsub,long *smadr,long *smci,double *sma);
  void select_submatrix (long *li,long nsub,comprow *smcr);
  void select_submatrix (long *li,long nsub,densemat *smdm);
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
  ///  array containing column indices
  long *ci;
  ///  array containing nonzero entries of the matrix
  double *a;
  ///  threshold for entries rejection
  double limit;
  ///  decomposition indicator
  long decompid;
};

#endif
