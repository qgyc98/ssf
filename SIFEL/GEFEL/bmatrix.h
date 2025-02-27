#ifndef DEF_BMATRIX
#define DEF_BMATRIX
#include <assert.h>
#include "iotools.h"
#include "matrix.h"
struct vector;
struct ivector;
//#define DEBUG_BMATRIX


#ifdef DEBUG_BMATRIX
 static unsigned long Acbm; ///< alocation counter
#endif 

struct bmatrix
/** This file declares struct bmatrix, which implements 
    %matrix with block structure and with %elements type of double.
    It can be used as a fourth order tensor or default matrix.
    There are also declarations of the functions for the basic %matrix operations
*/
{
  long m;      ///< number of rows of block matrices
  long n;      ///< number of columns of block matrices
  long tm, tn; ///< total number of double elements in rows and columns
  matrix *a;   ///< pointer to onedimensional array with block matrices elements stored in rows
  ivector row; ///< %vector of row indices
  ivector col; ///< %vector of column indices
  
  bmatrix() {m = n = 0; a = NULL;};  ///< default constructor
  bmatrix(long m, long n);           ///< allocating constructor
  bmatrix(const bmatrix &mat);       ///< copy constructor
  
  /**
   The operator enables access to the rows of elements of the member array a.

    @param i is the number of row

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  matrix* operator [] (long i) const
  {
    #ifdef DEBUG_BMATRIX
     if (i >= m)
       print_err("block row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);

     assert(i < m);
    #endif
    return (&a[i * n]);
  };


  /**
   The operator enables access to the elements of the member array a.

    @param i is the index of row
    @param j is the index of column

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  matrix& operator () (long i, long j) const
  {
    #ifdef DEBUG_BMATRIX
     if (i >= m)
       print_err("block row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);

     if (j >= n)
       print_err("block column index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert((i < m) && (j < n));
    #endif
    return (a[i * n + j]);
  };

  /**
   The operator enables access to the elements of the member array a.

    @param i is the index of row
    @param j is the index of column

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  double& operator () (long i, long j, long k, long l) const
  {
    #ifdef DEBUG_BMATRIX
     if (i >= m)
       print_err("block row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);

     if (j >= n)
       print_err("block column index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert((i < m) && (j < n));
    #endif
    return (a[i * n + j][k][l]);
  };

  /// generates  global indices vectors
  long gen_indices(); 

  /// returns total number of rows
  long give_totm() const;   
 
  /// returns total number of columns
  long give_totn() const;
  
  /// returns total indices from the given local indices
  void give_totid(long i, long j, long k, long l, long &ti, long &tj) const;

  /// returns local indices from the given total indices
  void give_locid(long ti, long tj, long &i, long &j, long &k, long &l) const;
 
  /// destructor
  ~bmatrix();
};




/// allocates matrix
long allocm(long m, long n, bmatrix &mat);

/// copies contents of matrix
long copym(const bmatrix &src, bmatrix &dest);

/// fills contents of matrix with given value
long fillm(double c, bmatrix &mat);
long fillrow(double c, long i, matrix &bmat);
long fillcol(double c, long i, matrix &bmat);

/// deallocates matrix
long destrm(bmatrix &mat);

/// adds 2 matrices
long addm(const bmatrix &a, const bmatrix &b, bmatrix &c);

/// substracts 2 matrices
long subm(const bmatrix &a, const bmatrix &b, bmatrix &c);

/// performs tensor product of two vectors
void tensprd (vector &a,vector &b,matrix &c);

/// multiplies two matrices A.B=C
long mxm(const bmatrix &a, const bmatrix &b, bmatrix &c);

/// multiplies matrix A with real constant c 
long cmulm(double c, bmatrix &a);

#endif
