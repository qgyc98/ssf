#ifndef GMATRIX_H
#define GMATRIX_H

//#include "LAPACK/spmatrix.h"
#include "SPARSE/DSSolver.h"

#include "galias.h"
#include "densemat.h"
#include "diagmat.h"
#include "dskyline.h"
#include "skyline.h"
#include "cr.h"
#include "scr.h"
#include "cc.h"
#include "scc.h"
//#include "spasol.h"
#include "elemmat.h"
#include "gtopology.h"

#include "slesolv.h"

class precond;

/**
   class describing a general %matrix
   
   JK
*/
class gmatrix
{
 public:
  gmatrix ();
  ~gmatrix ();
  void alloc ();
  void dealloc ();
  void initiate (gtopology *top,long ndof,storagetype ats,long mespr,FILE *out);

  //void setval (linsolvertype ls,precondtype p,long ni,double err);
  void setval (slesolv *ssle);

  void localize (matrix &lm,ivector &cn,long eid);
  void localized (double *lm,long *cn,long nc,long eid);
  void glocalize (matrix &lm,ivector &rcn,ivector &ccn);
  void mult_localize (long nm,long *ncn1,long *ncn2,long *mcn);

  void prepmat (double limit,long mespr);
  void prepmat2 (gtopology *top,FILE *out);
  void auxdatsparsesolver (gtopology *top,FILE *out);

  void solve_system (gtopology *top,precond &prec,double *lhs,double *rhs,FILE *out);
  void decompose_matrix ();
  void back_substitution (double *lhs,double *rhs);
  
  void incomplete_fact (double incompltresh);
  void back_incomplete_fact (double *x,double *y,double incompltresh);

  void condense (gtopology *top,double *condmat,double *condvect,double *lhs,double *rhs,long nrdof,long tc,FILE *out);
  void kernel (double *rbm,long &nse,long *se,long ense,double limit,long tc);
  void ldl_feti (double *lhs,double *rhs,long nse,long *se,double zero);
  void gmxv (double *a,double *b);
  void decompgmxv (double *a,double *b);
  void addgm (double a,gmatrix &gm);
  void scalgm (double a);
  void copygm (gmatrix &gm);
  void a12block (gtopology *top,double *block,long nrdof,FILE *out);
  void printmat (FILE *out);
  void printdiag (FILE *out);

  void changedecomp ();
  //void setfact ();
  //void setnotfact ();
  long decomp ();
  
  double give_entry (long ri,long ci);
  void add_entry (double e,long ri,long ci);


  long give_negm ();

  void diag_scale (double *d);
  void diag_check (double thr,double *rhs);

  double power_method (double *v,long ni,double err);
  double inverse_iteration (double *v,long ni,double err);

  double estim_spect_radius ();
  /// turns off solver messages on console output 
  void mute();

  
  ///  type of matrix storage
  storagetype ts;
  ///  type of solver of system of linear algebraic equations
  linsolvertype tlinsol;
  ///  type of condensation if tlinsol indicates condensation
  linsolvertype tsol;
  ///  type of preconditioning
  precondtype tprec;
  
  ///  block size in direct sparse solver
  unsigned char bsize;
  
  ///  number of unknowns = number of rows/columns
  long n;
  ///  computer zero
  double zero;
  ///  maximum number of iterations in an iterative method
  long ni;
  ///  required norm of the residual %vector in an iterative method
  double res;
  ///  number of performed iterations in an iterative method
  long ani;
  ///  attained norm of residual %vector in an iterative method
    double ares;
  ///  initial values in iterative methods
  ///  iv=0 - the input %vector is multiplied by zero
  ///       the iteration starts from zero %vector
  ///  iv=1 - the input %vector is used as a start %vector
  long iv;

  ///  parameters of preconditioning
  double omega;
  double indegamma;
  
  ///  number of rows which are not condensed
  long nbdof;

  ///  threshold for rejection from compressed storages
  double limit;
  
  ///  auxiliary information about block structure for sparse direct solver
  long *auxsd;
  
  ///  describes status of cr/scr - sparse exchange
  long pmat;
  
  densemat *dm;
  diagmat *diagm;
  skyline *sky;
  dskyline *dsky;
  comprow *cr;
  symcomprow *scr;
  compcol *cc;
  symcompcol *scc;
  elemmat *em;
  //MatrixPtr Matrix;
  ISolver *sdirect;
  //spasol *iss;

  //  PETSC
  //Mat petscmat;
  //Vec x;
};


#endif

