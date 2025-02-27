#ifndef GFMATRIX_H
#define GFMATRIX_H

#include "gfunct.h"
#include "xfile.h"

struct vector;
struct matrix;

struct gfmatrix
/** This file declares struct gfmatrix, which implements
    %matrix with %elements type of gfunct.
    There are also declarations of the functions for the %matrix and %vector computing */
{
  long m;    ///< number of rows
  long n;    ///< number of columns
  gfunct *a; ///< pointer to onedimensional array with matrix elements stored in the rows
  long size; ///< real length of array a (due to reallocation of memory)
  
  inline gfmatrix() {m = n = size = 0L; a = NULL;}; ///< default constructor
  gfmatrix(long m, long n);          ///< allocating constructor
  gfmatrix(long n);                  ///< allocating constructor for the square matrix
  
  /**
   The operator enables access to the elements of the member array a.

    @param i is the index of row
    @param j is the index of column

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline gfunct& operator () (long i, long j) const
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("matrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);
     if (j >= n)
       print_err("matrix column index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, j, n-1);

     assert((i < m) && (j < n));
    #endif
    return (a[i*n+j]);
  };

  /// evaluates individual %matrix components and store them into the real %matrix res
  void evaluate(vector &p, const char **namevar, matrix &res);
  /// evaluates individual %matrix components and store them in the transposed form to the real %matrix res
  void evaluate_t(vector &p, const char **namevar, matrix &res);

  inline ~gfmatrix() ///< destructor
  {
    m = 0;
    n = 0;
    delete [] a;
    a = NULL;
    size = 0;
  };
};



/// allocates gfmatrix to the dimensions m x n
long allocm(long m, long n, gfmatrix &mat);

/// allocates and initializes identity gfmatrix of dimension n x n
long allocim(long n, gfmatrix &mat);

/// reads gfmatrix from the file in
long readm(XFILE *in, gfmatrix &a);

/// prints the gfmatrix to the file
long printm(FILE *out, const gfmatrix &a);

#endif
