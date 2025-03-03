#include <stdlib.h>
#include <stdio.h>
#include "gfmatrix.h"
#include "gfunct.h"
#include "iotools.h"
#include "vector.h"
#include "matrix.h"

/**
   The constructor allocates memory form the heap to the member array a for the matrix mxn.
   
   @param m is the number of rows
   @param n is the number of columns
   
   created  10.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
gfmatrix :: gfmatrix(long m, long n)
{
  gfmatrix::m = m;
  gfmatrix::n = n;
  a = new gfunct[m*n];
  if (a == NULL)
  {
    print_err("cannot allocate memory for gfmatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
  }
}



/**
   The constructor allocates memory form the heap to the member array a for the quare %matrix nxn.
   
   @param n is the number of rows and columns
   
   created  5.1.2022, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
gfmatrix :: gfmatrix(long n)
{
  gfmatrix::m = n;
  gfmatrix::n = n;
  a = new gfunct[n*n];
  if (a == NULL)
  {
    print_err("cannot allocate memory for gfmatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, n, n);
  }
}



/**
  The function evaluates individual %matrix components and store them
  in the real %matrix res. Actual values of variables defined in the 
  componet general functions are given by %vector p and string with names 
  of variables are given by namevar.

  @param p - %vector with actual values of variables which the %matrix components depends on 
  @param namevar - array of pointers to strings with names of variables whose values are defined by p
  @param res - %matrix with resulting values of components

  Created  10.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void gfmatrix::evaluate(vector &p, const char **namevar, matrix &res)
{
  long i,j;

  for (i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
      res(i,j) = a[i*n+j].getval(p, namevar);
  }
}



/**
  The function evaluates individual %matrix components and store them
  in transposed form in the real %matrix res. Actual values of variables defined in the 
  componet general functions are given by %vector p and string with names 
  of variables are given by namevar.

  @param p - %vector with actual values of variables which the %matrix components depends on 
  @param namevar - array of pointers to strings with names of variables whose values are defined by p
  @param res - %matrix with resulting values of components

  Created  10.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void gfmatrix::evaluate_t(vector &p, const char **namevar, matrix &res)
{
  long i,j;

  for (i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
    {
      res(j,i) = a[i*n+j].getval(p, namevar);
    }
  }
}



/**
   The function allocates memory from heap to the mat's member array a.
   
   @param m is the number of rows
   @param n is the number of columns
   @param mat is the structure for allocated gfmatrix
   
   created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocm(long m, long n, gfmatrix &mat)
{
  mat.m = m;
  mat.n = n;
  mat.a = new gfunct[m*n];
  if (mat.a)
    return 0;
  else
    print_err("cannot allocate memory for gfmatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);

  return (1);
}



/**
   The function allocates memory from heap to the mat's member array a and
   it is initialized to be identity gfmatrix, i.e. all %matrix
   component general functions are stat type with 0 value except 
   of diagonal elements whose value is set to 1 .
   
   @param m is the number of rows
   @param n is the number of columns
   @param mat is the structure for allocated gfmatrix
   
   created  27.11.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocim(long n, gfmatrix &mat)
{
  mat.m = n;
  mat.n = n;
  mat.a = new gfunct[n*n];
  if (mat.a)
  {
    for (long i=0; i<n; i++)
      mat.a[i*n+i].f = 1.0;

    return (0);
  }
  else
    print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);

  return (1);
}



/**
   The function reads the contents of the gfmatrix a from the opened text file.

   @param in - pointer to the opened XFILE   
   @param a  - is the structure of the matrix which is read
   
   @return always zero.
   
   Created  10.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long readm(XFILE *in, gfmatrix &a)
{
  long i,j;

  for (i = 0; i < a.m; i++)
  {
    for (j = 0; j < a.n; j++)
      a(i,j).read(in);
  }
  return (0);
}



/**
   The function prints out the contents of the %matrix a to the file given by out.
   
   @param out   is the structure with opened file for %matrix output
   @param a     is the structure of the %matrix which is printed
   
   @return always zero.
   
   Created  10.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printm(FILE *out, const gfmatrix &a)
{
  long i,j;

  for (i = 0; i < a.m; i++)
  {
    for (j = 0; j < a.n; j++)
    {
      a(i,j).print(out);
      fprintf(out, "\n");
    }
  }
  return (0);
}

