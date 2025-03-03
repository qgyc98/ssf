#ifndef DPFETI_H
#define DPFETI_H

#include <stdio.h>
#include "mpi.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"
#include "pgalias.h"
#include "selnodes.h"

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
  dpfeti (int np,int mr,long nd);
  ~dpfeti();
  void initiate (selnodes *selnodschur,selnodes *selnodfeti);
  void cornodedetection (gtopology *top,long *ltg,long *domproc,FILE *out);
  void dpfetiordering (gtopology *top,long *ltg);
  void globcnnum_dpfeti (gtopology *top,long *ltg,long *domproc,FILE *out);
  void extract_from_local_vector (double *ev,double *lv);
  void put_into_local_vector (double *ev,double *lv);
  void extract_from_global_vector (double *ev,double *gv,long nsub);
  void put_into_global_vector (double *ev,double *gv,long nsub);

  void arrmatrix (double *condmat,long *domproc,FILE *out);
  void vectors_br_bm (gtopology *top,gmatrix *gm,
		      double *condmat,double *condvect,double *rhs,long *domproc,FILE *out);
  
  void rhs_dpfeti (gtopology *top,long *domproc,gmatrix *gm,FILE *out);
  
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
  long ndom;
  
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
  long nlbdof;
  //

  long *localcn;
  long *globalcn;

  
  


  ///  maximum number of iterations in conjugate gradient method
  long nicgdpfeti;
  ///  number of performed iterations in conjugate gradient method
  long anicgdpfeti;
  ///  required error
  double errcgdpfeti;
  ///  attained error
  double aerrcgdpfeti;
  ///  computer zero
  double zero;
  
  
  ///  matrix K_rc
  double *krc;
  
  //
  //  VARIABLES DEFINED ONLY ON MASTER PROCESSOR
  //
  
  //
  long **loccn;
  //
  long **globcn;
  //
  long *ncdofdom;
  //
  long *nlbdofdom;
  

  ///  type of reduced system solver
  redsystsolver trssol;
  
  ///  storage of reduced system matrix
  storagetype rsmstor;
  
  ///  type of linear system solver
  slesolv ssle;

  
  gtopology *ptop;

  ///  vector b_r
  double *br;
  ///  vector b_m
  double *bm;
  
  ///  matrix A_rr
  gmatrix *arr;
  
  double *cgrhs;
  double *cglhs;
  double *displ;
};


#endif
