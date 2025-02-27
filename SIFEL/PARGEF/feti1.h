#ifndef FETI1_H
#define FETI1_H

#include <stdio.h>
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"
#include "partop.h"
#include "selnodes.h"

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
  feti1 (int np,int mr,long nd);
  ~feti1();
  
  void initiate (selnodes *selnodfeti,FILE *out);

  void coarse_dofs (partop *ptop,FILE *out);
  void number_contributing_nodes_dofs (gtopology *top,partop *ptop,long *domproc,FILE *out);
  void contributing_nodes_dofs (gtopology *top,partop *ptop,long *domproc,FILE *out);
  void coarse_problem_ordering (gtopology *top,partop *ptop,long *domproc,FILE *out);
  void subdomain_ordering (gtopology *top,long *ltg,FILE *out);
  
  void local_buff (double *lv,double *buff);
  void buff_coarse (double *cv,double *buff,long nd);
  void coarse_buff (double *cv,double *buff,long nd);
  void buff_local (double *lv,double *buff);


  //void globcnnum_feti (gtopology *top,long *ltg,long *domproc,FILE *out);
  //void locbuff (double *buff,double *lv);
  //void buffloc (double *buff,double *lv);
  //void globbuff (double *buff,double *gv,long nd);
  //void buffglob (double *buff,double *gv,long nd);

  void hmatrixsize (double *rbm,long *domproc,FILE *out);
  void hmatrix (double *h,double *rbm,long *domproc);
  void qvector (double *q,double *rbm,double *f,long *domproc);

  void feti_projection (double *v,double *h,double *h1);
  void mpcg (gtopology *top,gmatrix *gm,long *domproc,
	     double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out);
  void lagrmultdispl (gtopology *top,gmatrix *gm,long *domproc,
		      double *w,double *d,double *f,double *rbm,long *rbmi,
		      double *h,double *h1,FILE *out); 
  void subdomain_matrix (gmatrix *gm,long *domproc,FILE *out);
  void scaling(double *invect,double *outvect,long n,FILE *out);
  void locscaling (double *invect,double *outvect,FILE *out);
  void lumpedprec (gmatrix *gm,double *dd,double *pp,FILE *out);
  void dirichletprec (double *dd,double *pp,FILE *out);

  
  void solve_system (gtopology *top,gmatrix *gm,
		     long *domproc,double *lhs,double *rhs,FILE *out);


  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  number of subdomain
  long ndom;
  
  
  

  //  maximum number of components in local-global arrays
  long maxlggl;
  //  number of components in local-global arrays
  long nclggl;
  
  
  //  array containing local code numbers of boundary unknonws
  long *lcn;
  
  
  //  estimate of number of rigid body motions
  long enrbm;
  //  threshold for linear dependence
  double lithr;

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
  // type of preconditioner of FETI method
  long fetiprecond;
  
  // ********************************************
  //  VARIABLES DEFINED ONLY ON MASTER PROCESSOR
  //
  



  
  //
  long *lggl;
  //
  long **glcor;

  //  array containing numbers of RBM on subdomains
  long *nrbmdom;
  //  array containing addresses of first RBM in coarse matrix
  long *rbmadr;
  
  
  
  
  // ***************************
  //  sobota 30.7.2005
  // *************************
  
    
  
  ///  M
  
    
  
  
  
  // **********************************
  // **********************************
  // **********************************
  // **********************************
  // **********************************
  //  nedele 31.7.2005
  
  ///  number of degrees of freedom on subdomain
  long ndof;

  ///  M
  ///  number of DOFs (unknowns) in coarse problem = total number of boundary DOFs
  long ndofcp;

  ///  maximum number of nodes contributing to the coarse problem
  long maxncnd;
  
  ///  maximum number of DOFs contributing to the coarse problem
  long maxncdofd;

  ///  number of nodes contributing to the coarse problem
  long ncn;

  ///  number of contributions (unknowns) from particular subdomains to the coarse problem
  long ncdof;
  
  ///  number of rigid body motions (RBM)
  long nrbm;
 
  ///  maximum number of rigid body motions (RBM)
  long maxnrbm;

  ///  M
  ///  size of matrix H
  long hsize;

  ///  M
  ///  array containing numbers of DOFs at corner nodes
  ///  it contains tnbn rows and bnmultip[i] columns
  ///  ndofncn[i][j]=k - the j-th node shared by the i-th coarse node contains k DOFs
  ///  array is allocated in function coarse_dofs
  long **ndofncn;
  
  ///  M
  ///  array of code numbers at coarse nodes
  ///  it contains tnbn components, each contains bnmultip[i] subcomponents and each subcomponents contains ndofncn[i][j] subsubcomponents
  ///  cncn[i][j][k]=m - the k-th DOF at the j-th node shared by the i-th coarse node has number m
  ///  array is allocated in function coarse_dofs
  long ***cncn;

  ///  M
  ///  array of numbers of nodes contributing to coarse problem
  ///  ncnd contains nproc components
  ///  ncnd[i]=j - the i-th subdomains contributes to coarse problem by j nodes
  ///  defined on the master processor, array is allocated in function number_contributing_nodes_dofs
  long *ncnd;

  ///  M
  ///  array of numbers of unknowns (DOFs) contributing to the coarse problem
  ///  ncdofd contains nproc components
  ///  ncdofd[i]=j - the i-th subdomains contributes to coarse problem by j contributions
  ///  defined on the master processor, array is allocated in function number_contributing_nodes_dofs
  long *ncdofd;
  
  ///  array containing code numbers contributing to the coarse problem
  ///  extracted values from subdomains to the coarse problem
  ///  edofs[i]=j - the i-th components contributing to the coarse problem has number j
  ///  array is allocated in function contributing_nodes_dofs
  long *edofs;

  ///  M
  ///  array of coarse code numbers
  ///  it contains tnbn rows and ncdofd[i] columns
  ///  ccn[i][j]=k - the j-th contribution from the i-th subdomain goes to the k-th coarse unknown
  ///  array is allocated in function contributing_nodes_dofs
  long **ccn;

  /// M
  double *wscalmat;
  
  // each processor
  double *lwscalmat;
  
  /// each processor
  ///  subdomain matrix for preconditioning
  skyline *smsky;
  ///  %matrix in the dense storage scheme (it contains Schur complements which are dense)
  densemat *smdm;
  //
  long  ndofprec;
  
  //
  long *cnprec;
  
  //
  long *cpreccn;





  
  // **********************
  //  REVISION 16.10.2011
  // **********************

  ///  type of FETI implementation
  ///  fetiimpl=no_impl=0 - no implementation is defined
  ///  fetiimpl=boolean_matrices=1 - Boolean %matrix is assembled, it is used for tests only, it is not efficient implementation
  ///  fetiimpl=nonredundant=2 - nonredundant constraints are defined, %matrix B has linearly independent rows
  ///  fetiimpl=redundant=3 - redundant constraints are defined, %matrix B has linearly dependent rows
  fetiimplem fetiimpl;

};


#endif
