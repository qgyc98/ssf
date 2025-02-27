#ifndef BOUNDPARCONGRAD_H
#define BOUNDPARCONGRAD_H

#include <stdio.h>
#include "gtopology.h"
#include "gmatrix.h"
#include "pgalias.h"
#include "partop.h"
#include "selnodes.h"
class PETSC_CONT;

/**

*/
class boundparcongrad
{
 public:
  boundparcongrad(int np,int mr,long nd,long mes);
  ~boundparcongrad ();
  
  void initiate (selnodes *selnodschur,partop *ptop,gtopology *top,FILE *out);
  double v_ixu_i(double *u,double *v);
  double uxv(double *u,double *v, long n);
  void select_bound(double *buff,double *v);
  void add_bound(double *buff,double *v);
  void add_master_buff (double *buff,double *cv,long nd);
  void select_master_buff(double *cv,double *buff,long nd);
  void ilu_mat_petsc(gmatrix *gm,FILE *out);
  void solve_fact_mat_petsc(double *r,double *z,long ndof,FILE *out);
  void initialize_diag_prec(gmatrix *gm,double *precvec,FILE *out);
  void jacobi_precondition(double *precvec,double *invec,double *outvec,long ndof,FILE *out);
  void solve_system (partop *ptop,gtopology *top,gmatrix *gm,long *domproc,double *lhs,double *rhs,FILE *out,long iv);
  
  /// message printing
  long mespr;
  
  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  number of subdomain
  long ndom;
  
  // number of iterations
  long ni;
  // error of computation
  double err;
  // type of precondition
  long prec;

  double *factmat;
      
  /// TOPOLOGY
  /// number of  all nodes on subdomain
  long nn; 
  ///  number of boundary nodes on subdomain
  long nbn;
  //  number of internal nodes on subdomain
  long nin;
  ///  number of degrees of freedom on subdomain
  long ndof;
  ///  number of internal degrees of freedom on subdomain
  long nidof;
  ///  number of boundary degrees of freedom on subdomain
  long nbdof;
  ///  maximum number of boundary DOFs on one subdomain
  long maxnbdof;
  ///  list of boundary DOFs
  /// lbdof[i] = j - i-th boundary DOF has local number j
  long *lbdof;
  ///  list of internal DOFs
  /// lidof[i] = j - i-th internal DOF has local number j
  long *lidof;
  
  ///  DATA DEFINED ONLY ON MASTER
  
  ///  TOPOLOGY
  
  ///  number of global DOFs = total number of boundary DOFs
  ///  number of DOFs of coarse problem
  long ndofcp;
  
  
  ///  array containing numbers of boundary DOFs on subdomains
  ///  nbdofdom[i] = j - i-th domain has j-th  numbers of boundary DOFs  
  long *nbdofdom;

  /// array containing boundary global code numbers
  /// bcn[i][j] = k - on i-th domain has j-th boundary DOF global code numbers k 
  long **bcn;
  PETSC_CONT *petscpoint;  


  int Argc;
  char **Argv;
  
};


#endif
