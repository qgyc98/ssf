#ifndef FLOATSUB_H
#define FLOATSUB_H

#include "gnode.h"
#include "lgnode.h"

/**
   class floatsub
   
   class deals with floating subdomains
   floating subdomains are used in e.g. problem of reinforcement in matrix, pull-out test, etc.
   
   JK, 7.11.2005
*/
class floatsub
{
 public:
  floatsub (void);
  ~floatsub (void);
  
  void ndof_subdom (gnode *gnodes);
  long number_of_lagr_mult (long nlgn,lgnode *lgnodes,gnode *gnodes);
  void displ_extract (long nlgn,lgnode *lgnodes,gnode *gnodes);
  
  void local_coarse (double *cv,double *lv);
  void coarse_local (double *cv,double *lv);
  
  ///  the number of subdomains (for problems with floating parts)
  long nsd;
  ///  the number of nodes on particular subdomains
  long *nnsd;
  ///  first node numbers on subdomains
  long *fnnsd;
  ///  number of elements on particular subdomains
  long *nesd;
  ///  first element numbers on subdomains
  long *fensd;
  ///  array containing numbers of DOFs on particular subdomains
  long *ndofd;
  ///  array containing number of the first unknown on subdomains
  long *fdofd;
  
  ///  the number of Lagrange multipliers
  long nlm;
  ///  extraction table, contains nonzero entries of constraint %matrix
  long **extrtab;
    
    
};

#endif
