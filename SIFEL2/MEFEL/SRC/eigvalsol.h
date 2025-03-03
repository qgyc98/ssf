#ifndef EIGVALSOL_H
#define EIGVALSOL_H

#include "alias.h"
#include "galias.h"
#include "xfile.h"
#include <stdio.h>

/**
   class eigvalsol defines solver of eigenvalues and eigenvectors
   
   JK, 20.8.2005
*/

class eigvalsol
{
 public:
  
  eigvalsol (void);
  ~eigvalsol (void);
  void read (XFILE *in);
  void print (FILE *out);
  
  
  ///  type of solver of eigenvalues and eigenvectors
  eigensolver teigsol;

  ///  number of required eigenvectors
  long neigv;
  ///  number of vectors used in computation
  long nev;
  ///  maximum number of iterations
  long nies;
  ///  number of performed iterations
  long anies;
  ///  required error
  double erres;
  ///  attained error
  double aerres;
  ///  maximum number of iteration in Jacobi' method
  long nijmr;
  ///  number of thresholds in Jacobi' method
  long njacthr;
  ///  array containing thresholds in Jacobi' method
  double *jacthr;
  
  ///  shift
  double shift;

};


#endif
