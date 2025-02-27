#ifndef LOCMATRIX_H
#define LOCMATRIX_H

#include "gmatrix.h"

/**
   class is considered for localization %matrix
   localization %matrix is used for localization of local matrices to
   global matrices and vice versa
   
   localization %matrix is stored in twodimensional array with
   variable lengths, only nonzero entries are stored
   
   JK, 23.12.2006
*/
class locmatrix
{
 public:
  locmatrix ();
  ~locmatrix ();
  
  void initiate_var (long nrows,long ncolumns);
  void addresses ();
  void allocate_ci ();
  void allocate_lm ();
  void initiate_nncr (long *a);
  void initiate_ci (long *colind);
  void initiate_ci (long id,long *colind);
  void initiate_lm (double *a);
  void initiate_lm (long id,double *a);
  
  void lmxv (double *a,double *b);
  void lmtxv (double *a,double *b);
  void lm01xv (double *a,double *b);
  void lmt01xv (double *a,double *b);
  
  void lm01xm (gmatrix &a,gmatrix &b);
  //void lmt01xm (double *a,double *b);
 
  void lmxmxlmt01 (gmatrix &a,gmatrix &b);
  void lmxmxlmt (gmatrix &a,gmatrix &b);

  ///  number of rows
  long nr;
  ///  number of columns
  long nc;
  ///  number of nonzero components
  long nnc;
  ///  threshold for zero
  double threshold;

  ///  list of the numbers of nonzero entries in particular rows
  long *nncr;
  ///  addresses of the first components in particular rows
  long *adr;
    
  ///  array containing localization %matrix
  double *lm;
  
  ///  column indices
  long *ci;

};

#endif
