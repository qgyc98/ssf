#ifndef FETI1_H
#define FETI1_H

#include <stdio.h>
#include "mpi.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"

/**
   
   *************************************************
   FINITE ELEMENT TEARING AND INTERCONNECTING METHOD
   *************************************************
   nproc, myrank, ndom
   number of processors, myrank and number of domain are established
   after constructor execution
   
   ndof
   total number of local unknowns is established after execution
   of function schurordering
   
   nbndom (M)
   array of numbers of boundary nodes on subdomains is established after
   execution of function globcnnum_feti
   
   maxnbn, totmaxndofn
   maximum number of boundary nodes and maximum number of degrees of freedom
   on node are established after execution of globcnnum_feti
   
   maxinc
   maximum number of incidencies is established after execution of globcnnum_feti

   ngdof (M)
   total number of global unknowns is established after execution of
   function globcnnum_feti
   
   nodnum
   array of legth nbn containing global node numbers; components correspond
   only to boundary nodes, internal nodes are excluded
   array is established after execution of function globcnnum_feti
   
   inc, lcngcn
   array of numbers of incidencies of nodes to subdomains and
   global code numbers in FETI method are established after execution
   of function globcnnum_feti
   
   
*/
class feti1
{
 public:
  feti1 (int np,int mr,int nd);
  ~feti1();
  void globcnnum_feti (gtopology *top,long *ltg,long *domproc,FILE *out);
  void locbuff (double *buff,double *lv);
  void buffloc (double *buff,double *lv);
  void globbuff (double *buff,double *gv,long nd);
  void buffglob (double *buff,double *gv,long nd);

  void hmatrixsize (double *rbm,long *domproc,FILE *out);
  void hmatrix (double *h,double *rbm);
  void qvector (double *q,double *rbm,double *f);

  void feti_projection (double *v,double *h,double *h1);
  void mpcg (gtopology *top,gmatrix *gm,
	     double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out);
  void lagrmultdispl (gtopology *top,gmatrix *gm,long *domproc,
		      double *w,double *d,double *f,double *rbm,long *rbmi,
		      double *h,double *h1);
  void lumpedprec (gtopology *top,double *v,double *pv);

  
  void solve_system (gtopology *top,gmatrix *gm,
		     long *domproc,double *lhs,double *rhs,FILE *out);


  //  number of processors
  int nproc;
  //  rank of processor
  int myrank;
  //  number of subdomain
  int ndom;
  
  //  number of degrees of freedom on subdomain
  long ndof;
  //  number of boundary nodes
  long nbn;
  //  total number of boundary nodes in problem
  long tnbn;
  //  number of boundary unknowns on subdomain
  long nbdof;
  
  //  maximum number of components in local-global arrays
  long maxlggl;
  //  number of components in local-global arrays
  long nclggl;
  //  maximum number of RBM
  long maxnrbm;
  
  
  //  array containing local code numbers of boundary unknonws
  long *lcn;
  
  //  array containing global code numbers in FETI method
  //long **lcngcn;
  //  array containing numbers of incidencies of nodes to subdomains
  //long *inc;
  //  array containing global node numbers of boundary nodes
  //long *nodnum;

  //  estimate of number of rigid body motions
  long enrbm;
  //  threshold for linear dependence
  double lithr;
  //  number of ridig body motions
  long nrbm;
  //  size of matrix H
  long hsize;

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

  //********************************************
  //  VARIABLES DEFINED ONLY ON MASTER PROCESSOR
  //
  //  number of global DOFs = total number of boundary DOFs
  long ngdof;
  
  //  array containing numbers of boundary nodes on subdomains
  long *nbndom;
  //  array containing numbers of RBM on subdomains
  long *nrbmdom;
  //  array containing addresses of first RBM in coarse matrix
  long *rbmadr;

};


#endif
