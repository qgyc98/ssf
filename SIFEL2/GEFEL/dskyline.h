#ifndef DSKYLINE_H
#define DSKYLINE_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"

/**
   class dskyline
   
   serves for %matrix storage called double skyline or double profile
   double skyline storage is available generaly for non-symmetric real matrices
   
   basic data
   double array a where %matrix entries are stored
   long array adr where addresses of diagonal entries are stored
   long number n stands for number of rows (columns) of the matrix
   long number negm stands for number of entries of matrix
   
   let skynegm be the number of stored %matrix entries in the %skyline
   storage scheme
   if the %matrix is factorized into form LU, the %matrix L is stored
   in the array a from 0 to skynegm and the %matrix U is stored
   in the array a from skynegm to 2 skynegm
   in other words, the matrix entries with ri>ci are between 0 and skynegm
   while the matrix entries with ri<ci are between skynegm and 2 skynegm
   
   basic relations
   array a has negm components
   array adr has n+1 components
   adr[n] is equal to negm
   
   matrices stored in dskyline format can be eliminated
   
   JK
*/
class dskyline
{
 public:
  dskyline (void);
  ~dskyline (void);
  void allocadr (long m);
  double* status ();
  long decomp ();
  void changedecomp ();
  void setfact ();
  void setnotfact ();

  void column_lengths_elem (long *cn,long ndofe);
  void column_lengths_mult (long *ncn1,long *ncn2,long *mcn,long nm);
  void column_lengths(gtopology *top);
  void addresses ();
  void neglobmat ();
  void allocglomat ();
  void nullsky ();
  void localize (matrix &b,long *cn);
  void localized (double *b,long *cn,long k);
  void glocalize (matrix &b,long *rcn,long *ccn);
  void mult_localize (long nm,long *ncn1,long *ncn2,long *mcn);
  void initiate (gtopology *top,long ndof,long mespr);

  void mxv_dsky (double *x,double *y);
  void mtxv_dsky (double *x,double *y);
  void addmat_dsky (double c,dskyline &dsky);
  void scalmat_dsky (double c);
  void copy_dsky (dskyline &dsky);
  void lu_dsky (double *x,double *y,double zero,long tc);
  void lukon_dsky (double *b,double *c,double *x,double *y,double zero,long nr,long tc);
  void bicg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv);
  void printmat (FILE *out);
  void printdiag (FILE *out);
  
  long give_negm ();
  void diag_check (double thr,double *rhs);

  ///  number of rows of the %matrix
  long n;
  ///  number of entries in the dskyline, it is equal to the size of the array a
  long negm;
  ///  addresses of diagonal entries
  long *adr;
  ///  global %matrix
  double *a;
  ///  decomposition (factorization) indicator
  long decompid;
};

#endif
