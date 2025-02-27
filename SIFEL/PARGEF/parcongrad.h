#ifndef PARCONGRAD_H
#define PARCONGRAD_H

#include <stdio.h>
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"
#include "pgalias.h"
#include "partop.h"

/**

*/
class parcongrad
{
 public:
  parcongrad(int np,int mr,long nd,long mes);
  ~parcongrad ();
  void initiate(partop *ptop,gtopology *top,FILE *out);
  void subdomain_ordering (partop *ptop,gtopology *top,FILE *out);
  void coarse_problem_ordering (long *domproc,gtopology *top,partop *ptop,FILE *out);

  void local_buff (double *lv,double *buff);
  void buff_coarse (double *cv,double *buff,long nd);
  void coarse_buff (double *cv,double *buff,long nd);
  void buff_local (double *lv,double *buff);
  
  double mod_scal_prod (double *u,double *v);
  void solve_system (partop *ptop,gtopology *top,gmatrix *gm,long *domproc,double *lhs,double *rhs,FILE *out,long iv);
  
  ///message printing
  long mespr;
  
  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  number of subdomain
  long ndom;
  
  ///  number of degrees of freedom on subdomain
  long ndof;
  ///  number of internal degrees of freedom on subdomain
  long nidof;
  ///  number of boundary degrees of freedom on subdomain
  long nbdof;
  ///  number of boundary nodes on subdomain
  long nbn;
  ///  maximum number of boundary DOFs on one subdomain
  long maxnbdof;
  
  ///  list of boundary DOFs
  long *lbdof;
  ///  list of internal DOFs
  long *lidof;
  
  ///  DATA DEFINED ONLY ON MASTER
    
  ///  number of global DOFs = total number of boundary DOFs
  ///  number of DOFs of coarse problem
  long ndofcp;

  ///  array containing numbers of boundary DOFs on subdomains
  long *nbdofdom;

  long **bcn;

  //  nova data
  long maxndof;
  long tndof;
  long *ndofdom;
  long ni;
  double err;
};


#endif
