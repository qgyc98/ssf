#ifndef DPFETI_NEW_H
#define DPFETI_NEW_H

#include <stdio.h>
#include "mpi.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"
#include "pgalias.h"

/**
   class dpfeti contains information about local-global ordering
   
   totmaxndofn, maxncn, maxncdof, maxnbn
   maximum number of degrees of freedom on node, maximum number of corner
   noeds on subdomain, maximum number of corner unknowns on subdomain and
   maximum number of bounadary nodes on subdomain are established after
   after execution of function globcnnum_dpfeti
   
   
   
*/
class dpfeti
{
 public:
  dpfeti (int np,int mr,int nd);
  ~dpfeti();
  void dpfetiordering (gtopology *top,long *ltg);
  void globcnnum_dpfeti (gtopology *top,long *ltg,long *domproc,FILE *out);
  void extract_from_local_vector (double *ev,double *lv);
  void put_into_local_vector (double *ev,double *lv);
  void extract_from_global_vector (double *ev,double *gv,long nsub);
  void put_into_global_vector (double *ev,double *gv,long nsub);

  void arrmatrix (double *condmat);
  void vectors_br_bm (gtopology *top,gmatrix *gm,
		      double *condmat,double *condvect,double *rhs);
  
  void rhs_dpfeti (gtopology *top,long *domproc,gmatrix *gm);
  
  void matxvect (gtopology *top,long *domproc,
		 gmatrix *gm,double *input,double *output);
  
  void cg (gtopology *top,long *domproc,gmatrix *gm,long iv,FILE *out);
  
  void corner_displ (gtopology *top,long *domproc,gmatrix *gm);

  void compute_displ (gtopology *top,gmatrix *gm,long *domproc,
		      double *subdispl,double *rhs);

  void solve_system (gtopology *top,gmatrix *gm,
		     long *domproc,double *lhs,double *rhs,FILE *out,long mespr);
  
  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  number of subdomain
  int ndom;
  
  ///  number of degrees of freedom on subdomain
  long ndof;
  ///  number of internal degrees of freedom on subdomain
  long nidof;
  ///  number of corner degrees of freedom on subdomain
  long ncdof;
  ///  number of all corner unknowns
  long tncdof;
  ///  number of all multipliers
  long tnmdof;


  //
  long maxncdof;
  long maxnbdof;



  
  


  //  maximum number of iterations in conjugate gradient method
  long nicg;
  //  number of performed iterations in conjugate gradient method
  long anicg;
  //  required error
  double errcg;
  //  attained error
  double aerrcg;
  //  computer zero
  double zero;
  //
  double limit,omega,indegamma;
  
  
  
  
  //  matrix K_rc
  double *krc;
  
  //
  //  VARIABLES DEFINED ONLY ON MASTER PROCESSOR
  //
  
  //
  long **loccn;
  //
  long **globcn;
  //
  long *nlbdof;
  //
  long *ncdofdom;
  
  

  //  type of reduced system solver
  redsystsolver trssol;
  
  //  storage of reduced system matrix
  storagetype rsmstor;
  
  //  type of linear system solver
  linsolvertype tlinsol;
  
  //  type of preconditioner
  precondtype tprec;

  gtopology *ptop;

  //  vector b_r
  double *br;
  //  vector b_m
  double *bm;
  
  //  matrix A_rr
  gmatrix *arr;
  
  double *cgrhs;
  double *cglhs;
  double *displ;
};


#endif
