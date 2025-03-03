#ifndef SCHURCOMPL_H
#define SCHURCOMPL_H

#include "pgalias.h"
#include "partop.h"
#include "partopjk.h"
#include "selnodes.h"
#include "slesolv.h"
#include "gtopology.h"
#include "gmatrix.h"
#include <stdio.h>

//class selnodes;

/**
   class paral contains information about local-global ordering
   
   **********************************
   PRIMAL DOMAIN DECOMPOSITION METHOD
   **********************************
   nproc, myrank, ndom
   number of processors, myrank and number of domain are established
   after constructor execution
   
   ndof, nidof
   total number of local unknowns and number of internal unknowns on subdomains
   are established after execution of function schurordering
   
   ngdof (M)
   total number of global unknowns is established after execution of
   function globcnnum_pdd
   
   maxnbn, maxnrdof, totmaxndofn
   maximum number of boundary nodes, maximum number of boundary unknowns on
   subdomain and maximum number of degrees of freedom on node are established
   after execution of globcnnum_pdd
   
   gcn
   array of local to global code numbers correspondence is established
   after execution of function globcnnum_pdd
   
   nbndom (M)
   array of numbers of boundary nodes on subdomains is established after
   execution of function globcnnum_pdd
   
   nrdofdom (M)
   array of numbers of reduced unknowns on subdomains is established after
   execution of function globcnnum_pdd
   
   masgcn (M)
   array of local to global code numbers of all subdomains is established after
   execution of function globcnnum_pdd
*/
class schurcompl
{
 public:
  schurcompl(int np,int mr,long nd);
  ~schurcompl ();
  //void search_comcn(gtopology *top, long *ltg, long *idgnn, long nbn, long **gnnc, long id, long ccn);
  //void subdomain_ordering (gtopology *top,partop *ptop);
  //void subdomain_ordering (gtopology *top,long *ltg,FILE *out);
  void initiate (partop *ptop,selnodes *selnodschur);
  void initiate (partopjk *ptop,FILE *out);
  //void coarse_problem_ordering (partop *ptop,FILE *out);
  void schurordering_new (gtopology *top, long *nodeidentif);
  void globcnnum_pdd (gtopology *top,long *ltg,long *domproc,FILE *out);

  void solve_red_sys_iter (long *domproc,double *condmat,double *condvect,FILE *out);
  void solve_red_sys_fin (long *domproc,double *condmat,double *condvect,long decomp,FILE *out);
  void solve_red_sys (long *domproc,double *condmat,double *condvect,long decomp,FILE *out);

  void solve_system (gtopology *top,gmatrix *gm,long *domproc,double *lhs,double *rhs,FILE *out);
  
  void gather_bound_vect (double *lv,double *gv,long *domproc);
  double pss_gather_bound_vect (double *lv1,double *gv1,double *lv2,double *gv2,long *domproc);
  long selected_norm_calculation (long cid,double err,double thresh,long sc,partopjk *ptop,gtopology *top,double *lv1,double *gv1,double *rhs,double *grhs,long *domproc,long n,FILE *out, double &norfv);
  double pss_gather_bound_vect_old (double *lv,double *gv,long *domproc);
  void scatter_bound_vect (double *lv,double *gv,long *domproc);
  
  double unbalanced_values (double *lv,long *domproc,FILE *out);
  
  ///  number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  ///  number of subdomain
  long ndom;
  
  ///  number of degrees of freedom on subdomain
  long ndof;
  ///  number of internal degrees of freedom (unknowns) on subdomain
  long nidof;
  ///  number of boundary degrees of freedom (unknowns) on subdomain
  long nbdof;

  //  number of edge degrees of freedom
  //long nedof;
  //  number of vertex degrees of freedom
  //long nvdof;
  

  //  maximum number of reduced DOFs on subdomain
  long maxnbdof;

  //  maximum number of iterations in conjugate gradient method
  long nicgsch;
  //  number of performed iterations in conjugate gradient method
  long anicgsch;
  //  required error
  double errcgsch;
  //  attained error
  double aerrcgsch;
  //  computer zero
  double zero;
  
  //long *dofinc;

  
  // ********************************************
  //  VARIABLES DEFINED ONLY ON MASTER PROCESSOR
  //

  ///  number of global DOFs = total number of boundary DOFs
  ///  number of DOFs of coarse problem
  long ndofcp;
  
  ///  type of reduced system solver
  redsystsolver trssol;
  
  ///  storage of reduced system matrix
  storagetype rsmstor;
  
  ///  type of linear system solver
  slesolv *ssle;

  ///  type of storage of domain matrix
  storagetype dmstor;

  //  type of preconditioner
  //precondtype tprec;
  
  //  array containing numbers of boundary DOFs on subdomains
  long *nbdofmas;
  // flag for deallocation of nbdofmas array in the destructor (1=delete array nbdofmas in the destructor)
  long destr_nbdofmas;
  
  gtopology *gtop;
  gmatrix *arr;
  
};


#endif
