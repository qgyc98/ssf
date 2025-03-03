#include <stdio.h>
#include <math.h>
#include <string.h>
#include "bmatrix.h"
#include "matrix.h"
#include "vector.h"


/**
   The constructor allocates memory form the heap to the member array a.
   
   @param m is the number of rows
   @param n is the number of columns
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
bmatrix::bmatrix(long m, long n)
{
  bmatrix::m = m;
  bmatrix::n = n;
  a = new matrix[m*n];
  if (a == NULL)
    print_err("allocating memory for bmatrix", __FILE__, __LINE__, __func__);

  #ifdef DEBUG_BMATRIX
   Acbm++;
  #endif
}


/**
   The copy constructor creates copy of the object given by the constant reference parameter.
   
   @param v is reference to the object of the vector which should be copied
   
   created  16.10.2003 by Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
bmatrix::bmatrix (const bmatrix &mat)
{
  long i, j;
  if ((mat.m != m) || (mat.n != n) || (mat.a == NULL))
  {
    m = mat.m;
    n = mat.n;
    a = new matrix[m*n];
    if (a == NULL)
      print_err("allocating memory for bmatrix", __FILE__, __LINE__, __func__);

  }
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      a[i*n+j] = mat.a[i*n+j];
  if (mat.row.n)
    row = mat.row;
  if (mat.col.n)
    col = mat.col;
  #ifdef DEBUG_BMATRIX
   print_warning("copy constructor of bmatrix is called.\n"
                 "Please check your code and make sure you want use copy constructor.",__FILE__, __LINE__, __func__);
   Acbm++;
  #endif
}


/**
   The destructor deallocates the memory occupied by the member array a.
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
bmatrix::~bmatrix()
{
  m = 0;
  n = 0;
  #ifdef DEBUG_BMATRIX
  if (a != NULL)
    Acbm--;
  #endif
  delete [] a;
  a = NULL;
}



/**
  Function checks all block matrices in rows and columns of bmatrix 
  for dimension consistency and generates vector of global indices in rows and columns.
 
  @retval 0 - on success
  @retval 1 - error in rows
  @retval 2 - error in columns
*/
long bmatrix::gen_indices()
{
  long i, j, tid, d;
  allocv(m, row);
  allocv(n, col);
  for (i=0; i < m; i++)
  {
    for(j=0; j < n; j++)
    {
      if (j == 0)
      {
        row[i] = a[i*n+j].m;
        continue;
      }
      if (row[i] != a[i*n+j].m)
      {
        print_err("wrong number of rows in matrices in %ld. row of block matrix\n", __FILE__, __LINE__, __func__, i+1);
        return 1;
      }
    }
  }
  for (i=0; i < m; i++)
  {
    for(j=0; j < n; j++)
    {
      if (i == 0)
      {
        col[j] = a[i*n+j].n;
        continue;
      }
      if (col[j] != a[i*n+j].n)
      {
        print_err("wrong number of columns in matrices in %ld.column of block matrix", __FILE__, __LINE__, __func__, j+1);
        return 2;
      }
    }
  }
  tid = 0;
  for(i=0; i < m; i++)
  {
    d = row[i];
    row[i] = tid;
    tid += d;;
  }
  tm = tid;
  tid = 0;
  for(i=0; i < n; i++)
  {
    d = col[i];
    col[i] = tid;
    tid += d;;
  }
  tn = tid;
  return 0;
}



/**
  Function computes the total number of rows in block matrix.

  @retval The function returns total number of rows.
*/
long bmatrix::give_totm() const
{
  long i, totm = 0;
  for (i=0; i < m; i++)
    totm += a[i*n].m;
  return totm;
}



/**
  Function computes the total number of columns in block matrix.

  @retval The function returns total number of columns.
*/
long bmatrix::give_totn() const
{
  long i, totn = 0;
  for (i=0; i < n; i++)
    totn += a[i].n;
  return totn;
}

/**
  Function computes the total index of double element from indices of block (i,j)
  and element indices (k,l).

  @param i - row index of block matrix
  @param j - column index of block matrix
  @param k - row index of element in block(i,j)
  @param k - column index of element in block(i,j)
  @param ti - total row index of element - output parameter
  @param tj - total column index of element - output parameter

  @retval The total indices are returned via parameters ti and tj
*/
void bmatrix::give_totid(long i, long j, long k, long l, long &ti, long &tj) const
{
  #ifdef DEBUG_BMATRIX
   if ((i >= m) || (j >=n) || (a[i][j].m >= k) || (a[i][j]>=l) ||
       (i < 0) || (j < 0) || (k < 0) || (l < 0))
   {
     print_err("required indices are out of range", __FILE__, __LINE__, __func__);
     abort();
   }
  #endif
  if (row.n && col.n)
  {
    ti = row[i]+k;
    tj = col[j]+l;
    return;
  }

  long ii;
  ti = 0;
  for(ii=0; ii<i; ii++)
    ti += a[ii*n].m;
  ti += k;
  tj = 0;
  for(ii=0; ii<j; ii++)
    tj += a[ii].n;
  tj += l;
}

/**
  Function computes the total index of double element from indices of block (i,j)
  and element indices (k,l).

  @param i - row index of block matrix
  @param j - column index of block matrix
  @param k - row index of element in block(i,j)
  @param k - column index of element in block(i,j)
  @param ti - total row index of element - output parameter
  @param tj - total column index of element - output parameter

  @retval The total indices are returned via parameters ti and tj
  
*/
void bmatrix::give_locid(long ti, long tj, long &i, long &j, long &k, long &l) const
{
  #ifdef DEBUG_BMATRIX
   if (tm == 0)
     tm = give_totm();
   if (tn == 0)
     tn = give_totn();
   if ((ti >= tm) || (tj >= tn) || (ti < 0) || (tj < 0))
   {
     print_err("required indices are out of range", __FILE__, __LINE__, __func__);
     abort();
   }
  #endif
  if (row.n && col.n)
  {
    for(i = 0; i < m; i++)
      if (ti < row[i])  break;
    for(j = 0; j < n; j++)
      if (tj < col[j])  break;
    if (i > 0)
      k = ti - row[i-1];
    else
      k = ti;
    if (j > 0)
      l = tj - col[j-1];
    else
      l = tj;
    return;
  }

  long  d;
  d = 0;
  for(i=0; i < m; i++)
  {
    if (ti < d)  break;
    d += a[i*n].m;
  }
  if (i > 0)
    k = ti - d;
  else
    k = ti;
  
  d = 0;
  for(j=0; j < n; j++)
  {
    if (tj < d)  break;
    d += a[j].n;
  }
  if (j > 0)
    l = tj - d;
  else
    l = tj;
}


/**
   The function allocates memory from heap to the mat's member array a.
   
   @param m is the number of rows
   @param n is the number of columns
   @param mat is the structure for allocated matrix
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long allocm(long m, long n, bmatrix &mat)
{
  mat.m = m;
  mat.n = n;
  mat.a = new matrix[m*n];
  if (mat.a == NULL)
  {
    print_err("allocating memory for bmatrix", __FILE__, __LINE__, __func__);
    return (1);
  }
  #ifdef DEBUG_BMATRIX
   Acbm++;
  #endif
  return (0);
}



/**
   The function copies matrix given by src to dest.
   
   @param src is the structure of source bmatrix to copy
   @param dest is the structure of destination bmatrix to which the contents of src will be copied 
   
   @b Requests:
   dest has to be setuped dimensions and allocated memory array for elements which is enough
   large to hold all contents of the matrix src.
   
   @retval 0: on succes
   @retval 1: in case incompatibility sizes of the src and dest matrices
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long copym(const bmatrix &src, bmatrix &dest)
{
  long i, j;
  if ((src.m != dest.m) || (src.n != dest.n))
  {
    print_err("copying bmatrix - incompatible size of matrices", __FILE__, __LINE__, __func__);
    return (1);
  }
  for (i = 0; i < src.m; i++)
    for (j = 0; j < src.n; j++)
      dest[i][j] = src[i][j];

  if (src.row.n)
    dest.row = src.row;    
  if (src.col.n)
    dest.col = src.col;    

  return (0);
}


/**
   The function fills memory of mat's member array a with constant c.
   
   @param c   is constant, which will be used for filling entire matrix
   @param mat is the structure for the allocated bmatrix
   
   @return always zero
   
   created  22.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long fillm(double c, bmatrix &mat)
{
  long i, j, k, l;
  for (i = 0; i < mat.m; i++)
    for (j = 0; j < mat.n; j++)
      for (k = 0; k < mat[i][j].m; k++)
        for (l = 0; l < mat[i][j].n; l++)
	  mat[i][j][k][l] = c;
  return (0);
}



/**
   The function deallocates memory occupied by mat's member array a.
   
   @param mat is the structure for the deallocated bmatrix
   
   @return always zero
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long destrm(bmatrix &mat)
{
  #ifdef DEBUG_BMATRIX
   if (mat.a != NULL)
     Acbm--;
  #endif
  mat.m = 0;
  mat.n = 0;
  delete [] mat.a;
  mat.a = NULL;
  return (0);
}



/**
   The function adds bmatrix given by a to bmatrix given by b, the result is stored in c.
   
   @param a is the structure of the first added bmatrix
   @param b is the structure of the second added bmatrix
   @param c is the structure of the result bmatrix
   
   @b Requests:
   a, b and c have to be same dimensions, c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0: on succes.
   @retval 1: in case incompatibility sizes of a, b and c matrices.
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long addm(const bmatrix &a, const bmatrix &b, bmatrix &c)
{
  long i,j,k,l;
  long ti, tj;
  long ai, aj, ak, al;
  long bi, bj, bk, bl;
  
  if ((a.tm != b.tm) || (a.tn != b.tn) ||
      (c.tm != b.tm) || (c.tn != b.tn))
  {
    print_err("adding bmatrix - incompatible size", __FILE__, __LINE__, __func__);
    return (1);
  }  
  for (i = 0; i < c.m; i++)
  {
    for (j = 0; j < c.n; j++)
    {
      for (k = 0; k < c[i][j].m; k++)
      {
        for (l = 0; l < c[i][j].n; l++)
	{
          c.give_totid(i, j, k, l, ti, tj);
          a.give_locid(ti, tj, ai, aj, ak, al);
          b.give_locid(ti, tj, bi, bj, bk, bl);
          c[i][j][k][l] = a[ai][aj][ak][al] + b[bi][bj][bk][bl];
	}
      }
    }
  }
  return(0);
}


/**
   The function substracts matrix given by b from matrix given by a, the result is stored in the c
   
   @param a is the structure of the bmatrix, from which is sbustracted
   @param b is the structure of substracted bmatrix
   @param c is the structure of the result bmatrix
   
   @b Requests:
   a, b and c have to be same dimensions, c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0: on succes
   @retval 1: in case incompatibility sizes of a, b and c matrices
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long subm(const bmatrix &a, const bmatrix &b, bmatrix &c)
{
  long i,j,k,l;
  long ti, tj;
  long ai, aj, ak, al;
  long bi, bj, bk, bl;
  
  if ((a.tm != b.tm) || (a.tn != b.tn) ||
      (c.tm != b.tm) || (c.tn != b.tn))
  {
    print_err("adding bmatrix - incompatible size", __FILE__, __LINE__, __func__);
    return (1);
  }  
  for (i = 0; i < c.m; i++)
  {
    for (j = 0; j < c.n; j++)
    {
      for (k = 0; k < c[i][j].m; k++)
      {
        for (l = 0; l < c[i][j].n; l++)
	{
          c.give_totid(i, j, k, l, ti, tj);
          a.give_locid(ti, tj, ai, aj, ak, al);
          b.give_locid(ti, tj, bi, bj, bk, bl);
          c[i][j][k][l] = a[ai][aj][ak][al] + b[bi][bj][bk][bl];
	}
      }
    }
  }
  return(0);
}




/**
   The function multiplies matrix given by a from right with matrix given by b,
   the result is stored in c
   
   @param a is the structure of the bmatrix, from which is multiplied
   @param b is the structure of the multiplicating bmatrix
   @param c is the structure of the result bmatrix
   
   @b Requests:
   a, b and c have to be following total dimensions
   a (tm,tn),  b (tn,tp),  c (tm,tp)
   c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0: on succes
   @retval 1: in case incompatibility sizes of a, b and c matrices
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long mxm(const bmatrix &a, const bmatrix &b, bmatrix &c)
{
  
  long i,j,k;
  long ai, aj, ak, al;
  long bi, bj, bk, bl;
  long ci, cj, ck, cl;

  if ((a.tn != b.tm) ||
      (c.tm != a.tm) || (c.tn != b.tn))
  {
    print_err("multiplying matrix - incompatible size", __FILE__, __LINE__, __func__);
    return (1);
  }
  for (i=0; i < a.tm; i++)
    for (j=0; j < b.tn; j++)
    {
      c.give_locid(i, j, ci, cj, ck, cl);
      c[ci][cj][ck][cl] = 0.0;
      for (k=0; k < a.tn; k++)
      {
        a.give_locid(i, k, ai, aj, ak, al);
        b.give_locid(k, j, bi, bj, bk, bl);
        c[ci][cj][ck][cl] += a[ai][aj][ak][al] * b[bi][bj][bk][bl];
      }
    }
  return(0);
}



/**
   The function multiplies matrix given by a with constant c,
   the result is stored in a.
   
   @param c is the real number
   @param a is the structure of the bmatrix, from which is multiplied
   
   @return always zero
   
   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
*/
long cmulm(double c, bmatrix &a)
{
  
  long i, j, k, l;
  for (i = 0; i < a.m; i++)
    for (j = 0; j < a.n; j++)
      for (k = 0; k < a[i][j].m; k++)
        for (l = 0; l < a[i][j].n; l++)
	  a[i][j][k][l] *= c;
  
  return (0);
}
