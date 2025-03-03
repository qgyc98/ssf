#ifndef SLESOLV_H
#define SLESOLV_H

#include <stdio.h>
#include "galias.h"
#include "xfile.h"
#include "precond.h"
#include "seqfeti.h"
#include "seqschur.h"

class gmatrix;
class gtopology;
class comprow;
class fllopsif;
class saddpoint;
class comprow;
class fllopsif;


/**
   class slesolv defines type of solver of systems of linear algebraic equations
   
   JK
*/

class slesolv
{
 public:
  slesolv ();
  ~slesolv ();
  void read (gtopology *top,XFILE *in,long mespr);
  void print (FILE *out);
  void initiate (slesolv *ssle);

  void prepare_data (long *ndof,gtopology *top,FILE *out);
  
  void solve_system (gtopology *gt,gmatrix *gm,double *lhs,double *rhs,FILE *out);

  int solver_fllop (comprow *compr,double *rhs, double *lhs);
  
  ///  type of solver of system of linear algebraic equations
  linsolvertype tlinsol;

  ///  if tlinsol is conden, additional information is needed
  linsolvertype tsol;
  
  ///  block size in direct sparse solver
  unsigned char bsize;
    
  ///  computer zero
  double zero;
  
  ///  maximum number of iterations in iterative methods
  long ni;
  ///  required norm of residual %vector
  double res;
  ///  number of performed iterations (atual number of iterations)
  long ani;
  ///  attained norm of residual %vector
  double ares;
  
  ///  initial values in iterative methods
  ///  iv=0 - the input %vector is multiplied by zero
  ///       the iteration starts from zero %vector
  ///  iv=1 - the input %vector is used as a start %vector
  long iv;

  ///  number of rows which are not condensed
  long nbdof;

  ///  preconditioner
  precond prec;
  
  ///  FETI solver
  seqfeti feti;
  
  ///  solver based on the Schur complement method
  seqschur schur;

  //  saddle point solver
  saddpoint *sp;

  /// FLLOP solver
  fllopsif *fllopptr;
};

#endif
