#include "diagmat.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "precond.h"

diagmat::diagmat (void)
{
  //  the number of rows in the matrix (it is equal to the number of unknowns)
  n=0;
  //  array containing matrix entries
  a=NULL;
  //  indicator of matrix factorization
  decompid=0;
}

diagmat::~diagmat (void)
{
  delete [] a;
}

/**
   function allocates array containing %matrix
   
   @param m - number of rows/columns
   
   JK
*/
void diagmat::alloc (long m)
{
  n=m;
  negm=m;
  a = new double [negm];
  memset (a,0,negm*sizeof(double));
}

void diagmat::dealloc (void)
{
  delete [] a;
  a = NULL;
}

void diagmat::copy (diagmat *dm)
{
  long i;
  
  if (a != NULL)
    delete [] a;  
  
  n=dm->n;
  negm=n;
  a = new double [negm];
  
  for (i=0;i<n;i++){
    a[i]=dm->a[i];
  }
}

/**
   function returns status of array a
   
   JK
*/
double* diagmat::status ()
{
  return a;
}

/**
   function returns indicator of decomposition (factorization)
   
   JK
*/
long diagmat::decomp ()
{
  return decompid;
}

/**
   function changes indicator of decomposition (factorization)
   
   JK
*/
void diagmat::changedecomp ()
{
  if (decompid==0)  decompid=1;
  else              decompid=0;
}

void diagmat::setfact ()
{
  decompid=1;
}
void diagmat::setnotfact ()
{
  decompid=0;
}

/**
   function fills array by zeros
   
   JK
*/
void diagmat::nullmat ()
{
  long i;

  for (i=0;i<n;i++){
    a[i]=0.0;
  }
}

/**
   function localizes local %matrix b into global %matrix
   
   @param b - local %matrix
   @param cn - array containing code numbers
   
   JK
*/
void diagmat::localize (matrix &b,long *cn)
{
  long i,ii;
  
  for (i=0;i<b.m;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    a[ii]+=b[i][i];
  }
}

/**
   function localizes local %matrix stored in array into global %matrix
   
   @param b - array containing local %matrix
   @param cn - array containing code numbers
   @param m - number of rows/columns of the %matrix b
   
   JK
*/
void diagmat::localized (double *b,long *cn,long m)
{
  long i,ii;
  
  for (i=0;i<m;i++){
    ii=cn[i]-1;
    if (ii<0)  continue;
    a[ii]+=b[i*m+i];
  }
}

/**
   function initiates diagmat storage
   
   @param ndof - number of rows/columns of the %matrix
   @param mespr - message printing indicator
   
   JK
*/
void diagmat::initiate (long ndof,long mespr)
{
  if (status ()==NULL){
    alloc (ndof);
  }
  else{
    nullmat ();
  }
  
  if (mespr==1)  fprintf (stdout,"\n number of matrix entries   %ld",negm);
}

/**
   function multiplies %matrix stored as dense %matrix by %vector b,
   resulting %vector is c
   
   @param b - array containing %vector b
   @param c - array containing %vector c
   
   JK
*/
void diagmat::mxv_diag (double *b,double *c)
{
  long i;
  
  for (i=0;i<n;i++){
    c[i]=a[i]*b[i];
  }
}

/**
   function adds premultiplied %matrix stored in dm by coefficient c to actual %matrix
   
   @param c - multiplicative coefficient
   @param dm - another dense %matrix storage
   
   JK
*/
void diagmat::addmat_diag (double c,diagmat &dm)
{
  long i;

  for (i=0;i<n;i++){
    a[i]+=c*dm.a[i];
  }
}

/**
   function multiplies components of the %matrix by coefficient c
   
   @param c - multiplicative coefficient
   
   JK
*/
void diagmat::scalmat_diag (double c)
{
  long i;

  for (i=0;i<n;i++){
    a[i]*=c;
  }
}





/**
   function prints %matrix into output file
   
   @param out - output stream
   
   JK
*/
void diagmat::printmat (FILE *out)
{
  long i;
  
  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"%3ld %3ld  %20.15le\n",i,i,a[i]);
  }
  
}


/**
   function prints diagonal entries of the %matrix
   
   @param out - output stream
   
   JK
*/
void diagmat::printdiag (FILE *out)
{
  long i;
  
  fprintf (out,"\n\n");
  for (i=0;i<n;i++){
    fprintf (out,"%5ld  %.10e\n",i,a[i]);
  }
}


/**
   function returns number of stored %matrix entries
   
   JK, 25.5.2019
*/
long diagmat::give_negm ()
{
  return negm;
}

double diagmat::give_entry (long ri)
{
  return a[ri];
}



/**
   function estimates spectral radius

   the estimates are based on the Gershgorin circle theorem
   for details see G.H. Golub, C.F. Van Loan: Matrix computations.
   The Johns Hopkins University Press, 3rd edition, page 320, 1996
   
   JK, 27.8.2008
*/
double diagmat::estim_spect_radius ()
{
  long i,j;
  double s,sr;
  
  //  estimate of the spectral radius
  sr=0.0;
  for (i=0;i<n;i++){
    s=0.0;
    for (j=0;j<n;j++){
      s+=fabs(a[i*n+j]);
    }
    if (s>sr)
      sr=s;
  }
  
  return sr;
}

/**
   function solves system of linear algebraic equations with diagonal %matrix
   
   @param x - array containing vectors of unknowns
   @param y - array containing right hand sides
   @param zero - computer zero
   
   JK, 19. 7. 2019
*/
void diagmat::gemp (double *x,double *y,double zero)
{
  long i;

  
  for (i=0;i<n;i++){
    if (fabs(a[i])<zero){
      print_err("singular matrix is detected a[%ld]=%le",__FILE__, __LINE__, __func__,i,a[i]);
      abort ();
    }
    x[i]=y[i]/a[i];
  }
}

/**
   function solves system of linear algebraic equations with diagonal %matrix
   
   @param x - array containing vectors of unknowns
   @param y - array containing right hand sides
   @param zero - computer zero
   
   JK, 5. 3. 2022
*/
void diagmat::diagonal_solver (double *x,double *y,double zero)
{
  long i;

  
  for (i=0;i<n;i++){
    if (fabs(a[i])<zero){
      print_err("singular matrix is detected a[%ld]=%le",__FILE__, __LINE__, __func__,i,a[i]);
      abort ();
    }
    x[i]=y[i]/a[i];
  }
}
