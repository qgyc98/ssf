#ifndef SADDPOINT_H
#define SADDPOINT_H

#include "densemat.h"
#include "gtopology.h"
#include "gmatrix.h"

/**
   class deals with saddle point problems
   
   JK, 28.11.2007
*/
class saddpoint
{
 public:
  saddpoint ();
  ~saddpoint ();
  
  void read (XFILE *in);
  void initiate (seqselnodes *selnodfeti,gtopology *top,FILE *out);
  void matrix_e ();
  
  void get_jumps (long nc,double *jumps);
  void solve_system (gtopology *top,gmatrix *gm,double *lhs,double *rhs,FILE *out);
  
  ///  number of subdomains
  long ns;
  
  long n;
  
  ///  number of dual unknowns (number of Lagrange multipliers)
  long nm;
  ///  number of reduced/boundary/interface DOFs
  long nrdof;
  ///  number of internal DOFs
  long nidof;
  
  ///  dense %matrix for the final system of equations
  densemat *dm;
  
  double *c;


  ///  number of DOFs (unknowns) in coarse problem
  long ndofcp;

  ///  array of numbers of unknowns (DOFs) contributing to the coarse problem
  ///  ncdofd contains nproc components
  ///  ncdofd[i]=j - the i-th subdomains contributes to coarse problem by j contributions
  long *ncdofd;

  ///  array containing code numbers contributing to the coarse problem
  ///  extracted values from subdomains to the coarse problem
  ///  edofs[i][j]=k - the j-th components contributing to the coarse problem from the i-th subdomains has number k
  long **edofs;

  ///  array of coarse code numbers
  ///  it contains tnbn rows and ncdofd[i] columns
  ///  ccn[i][j]=k - the j-th contribution from the i-th subdomain goes to the k-th coarse unknown
  long **ccn;

  ///  node-subdomain correspondence
  ///  nsid[i]=j - the i-th node belongs to the j-th subdomain
  long *nsid;
  
  
  ///  numbers of DOFs on subdomains
  ///  ndofdom[i]=j - the i-th subdomain contains j DOFs
  long *ndofdom;
  
  ///  list of DOFs on subdomains
  ///  cndom[i][j]=k - the j-th DOF on the i-th subdomain has number k
  long **cndom;


  ///  matrix E
  double *e;
};

#endif
