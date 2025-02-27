#ifndef SEQSCHUR_H
#define SEQSCHUR_H

#include "skyline.h"
#include "seqselnodes.h"
#include "scr.h"
#include "densemat.h"

class gtopology;
class gmatrix;
class slesolv;

/**
   class seqschur deals with the Schur complement method
   this is single-processor implementation, it is different
   from schur in PARGEF which is multi-processor implementation
   
   JK, 13.5.2009
*/
class seqschur
{
 public:
  seqschur ();
  ~seqschur (); 
  
  void read (gtopology *top,long mespr,XFILE *in);
  void print (FILE *out);

  //void search_comcn(gtopology *top, long *ltg, long *idgnn, long nbn, long **gnnc, long id, long ccn);
  //void subdomain_ordering (gtopology *top,partop *ptop);
  //void subdomain_ordering (gtopology *top,long *ltg,FILE *out);
  void initiate (seqselnodes *selnodschur);


  //void coarse_problem_ordering (FILE *out);
  void schurordering_new (gtopology *top, long *nodeidentif);
  void globcnnum_pdd (gtopology *top,long *ltg,long *domproc,FILE *out);

  void assemble_subdom_unknowns (gtopology *top,FILE *out);
  void subdomain_matrices (gmatrix *gm,FILE *out);

  void solve_red_sys_iter (double **condmat,double **condvect,FILE *out);
  void solve_red_sys_fin (double **condmat,double **condvect,FILE *out);
  void solve_red_sys (double **condmat,double **condvect,FILE *out);

  void solve_system (gtopology *top,gmatrix *gm,double *lhs,double *rhs,FILE *out);
  
  //void gather_bound_vect (double *lv,double *gv,long *domproc);
  //double pss_gather_bound_vect (double *lv,double *gv,long *domproc);
  //void scatter_bound_vect (double *lv,double *gv,long *domproc);
  
  //double unbalanced_values (double *lv,long *domproc,FILE *out);
  
  ///  number of subdomains
  long ns;

  ///  type of reduced system solver
  redsystsolver trssol;
  
  ///  type of linear system solver
  slesolv *ssle;
  
  ///  storage of reduced system %matrix
  storagetype rsmstor;

  ///  maximum number of iterations in conjugate gradient method
  long nicg;
  ///  number of performed iterations in conjugate gradient method
  long anicg;
  ///  required error
  double errcg;
  ///  attained error
  double aerrcg;
  ///  computer zero
  double zero;
  

  ///  number of degrees of freedom on subdomain
  long ndof;
  ///  number of internal degrees of freedom (unknowns) on subdomain
  long nidof;
  ///  number of boundary degrees of freedom (unknowns) on subdomain
  long nbdof;

  //  maximum number of reduced DOFs on subdomain
  long maxnbdof;

  ///  number of coarse DOFs = total number of boundary/interface DOFs
  ///  number of DOFs of the coarse problem
  long ndofcp;

  //  array containing numbers of boundary/interface DOFs on subdomains
  long *nbdofmas;
  
  ///  generalized topology, each subdomain is assumed as generalized element (superelement)
  gtopology *gtop;

  ///  subdomain matrices stored in the %skyline storage
  skyline *smsky;

  // *******************






  



  ///  node-subdomain correspondence
  ///  nsid[i]=j - the i-th node belongs to the j-th subdomain
  long *nsid;
  
  
  ///  numbers of DOFs on subdomains
  ///  ndofdom[i]=j - the i-th subdomain contains j DOFs
  long *ndofdom;
  
  ///  list of DOFs on subdomains
  ///  cndom[i][j]=k - the j-th DOF on the i-th subdomain has number k
  long **cndom;


  
  
  

  ///  type of storage of domain matrix
  storagetype dmstor;

  //  type of preconditioner
  //precondtype tprec;
  
  
  

  gmatrix *arr;

};

#endif
