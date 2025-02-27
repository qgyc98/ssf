#ifndef PARAL_H
#define PARAL_H

#include <stdio.h>
#include "mpi.h"
#include "../GEFEL/gtopology.h"

/**
   class paral contains information about local-global ordering
   
   **********************************
   PRIMAL DOMAIN DECOMPOSITION METHOD
   **********************************
   nproc, myrank, ndom
   number of processors, myrank and number of domain are established
   after constructor execution
   
   ndof, indof
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
   
   
   ***************************************************************
   DUAL-PRIMAL FINITE ELEMENT TEARING AND INTERCONNECTING METHOD
   ***************************************************************
   
   totmaxndofn, maxncn, maxncdof, maxnbn
   maximum number of degrees of freedom on node, maximum number of corner
   noeds on subdomain, maximum number of corner unknowns on subdomain and
   maximum number of bounadary nodes on subdomain are established after
   after execution of function globcnnum_dpfeti
   
   
   
*/
class paral
{
 public:
  paral (int np,int mr,int nd);
  ~paral();
  void read (FILE *in,gtopology *top);
  void procdomcorr ();
  void schurordering (gtopology *top);
  void globcnnum_pdd (gtopology *top,FILE *out);

  void globcnnum_feti (gtopology *top,FILE *out);
  void locglobfeti (gtopology *top,double *gv,double *lv);
  void globlocfeti (gtopology *top,double *gv,double *lv);
  
  void dpfetiordering (gtopology *top);
  void globcnnum_dpfeti (gtopology *top,FILE *out);
  
  
  //  number of processors
  int nproc;
  //  rank of processor
  int myrank;
  //  number of subdomain
  int ndom;
  
  //  domain-processor correspondence
  long *domproc;
  
  //  number of degrees of freedom on subdomain
  long ndof;
  //  number of internal degrees of freedom on subdomain
  long indof;
  //  number of boundary nodes
  long nbn;
  //  maximum number of reduced DOFs on subdomain
  long maxnrdof;
  //  maximum number of degrees of freedom on node
  long totmaxndofn;

  //  local to global correspondence of node numbers
  long *ltg;

  //  PDD
  //  local to global code numbers correspondence in PDD
  long *gcn;
  //  asi to smazu

  
  //  FETI
  //  array containing global code numbers in FETI method
  long **lcngcn;
  //  array containing numbers of incidencies of nodes to subdomains
  long *inc;
  //  array containing global node numbers of boundary nodes
  long *nodnum;


  //********************************************
  //  VARIABLES DEFINED ONLY ON MASTER PROCESSOR
  //
  //  number of global DOFs = total number of boundary DOFs
  long ngdof;
  
  //  array containing numbers of boundary nodes on subdomains
  long *nbndom;
  //  array containing numbers of reduced DOFs on subdomains
  long *nrdofdom;
  //  array containing global code numbers of all subdomains PDD
  long **masgcn;

  //  array containing numbers of RBM on subdomains
  long *nrbmdom;
  //  array containing addresses of first RBM in coarse matrix
  long *rbmadr;
};


#endif
