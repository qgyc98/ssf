#ifndef ELEMMAT_H
#define ELEMMAT_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"

/**
   class elemmat serves for storage of matrices of finite elements
   no global matrix is used, only iterative methods are applicable,
   maximum exploatation of cache is considered
   
   JK
*/

class elemmat
{
 public:
  elemmat (void);
  ~elemmat (void);

  void alloc (gtopology *top);
  double** status ();
  long decomp ();
  void changedecomp ();
  void nullmat ();
  void localize (matrix &b,long *cn,long eid);
  void localized (double *b,long *cn,long eid,long m);
  void initiate (gtopology *top,long ndof,long mespr);
  void mxv_em (double *b,double *c);
  void addmat_em (double c,elemmat &dm);
  void scalmat_em (double c);
  void printmat (FILE *out);
  void cg (double *x,double *y,
	   long ni,double err,long &ani,double &ares,double zero,long iv);

  ///  number of unknowns (degrees of freedom)
  long n;
  ///  number of finite elements
  long ne;
  ///  total number of stored entries
  long mem;
  ///  maximum number of components in local %vector
  long maxarr;
  ///  array containing number of DOFs on elements
  long *andofe;
  ///  array containing number of matrix entries on elements
  long *neme;
  ///  array containing element matrices
  double **a;
  ///  array containing code numbers of finite elements
  long **acn;
  ///  decomposition indicator
  long decompid;

};

#endif
