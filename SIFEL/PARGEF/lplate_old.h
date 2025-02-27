#ifndef LPLATE_H
#define LPLATE_H

#include "mpi.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"

class lplate
{
 public:
  lplate (int np,int mr,int nd,gtopology *top);
  ~lplate ();
  void globcnnum_lpp (gtopology *top,long *domproc,FILE *out);
  void constraintmat (double *th,long *domproc,FILE *out);
  void maxndofdom ();
  void assemble_cn (long *buff,long ndom);
  void assemble_thicknesses (double *buff,long ndom);
  void codenum_mult (long &nm);
  void assemble_constr ();
  void assemble_constr_dm ();
  void orthonormalization ();
  void orthonormalization_dm ();
  void locglob (double *gv,double *lv,long ns);
  void globloc (double *gv,double *lv,long ns);
  void printconstrmat (FILE *out);
  void printconstrmat_dm (FILE *out);

  void cg (gmatrix *gm,double *lhs,double *rhs,long *domproc,FILE *out);
  void nodaldisplacements (gmatrix *gm,double *lhs,double *rhs,double *w,long *domproc);

  void solve_system (gtopology *top,gmatrix *gm,
		     long *domproc,double *lhs,double *rhs,FILE *out);
  
  
  ///  number of processors
  int nproc;
  ///  my rank
  int myrank;
  ///  number of domain
  long ndom;
  ///  number of nodes in one layer
  long nn;
  ///  number of layers
  long nl;
  ///  number of degrees of freedom of one node
  long ndofn;
  ///  number of multipliers on one node
  long nnmult;
  ///  number of degrees of freedom of one subdomain
  long ndof;
  ///  number of Lagrange multipliers
  long nmult;
  
  ///  maximum number of DOFs on subdomain
  long maxndof;
  ///  required error in conjugate gradient method
  double errcg;
  ///  maximum naumber of iterations
  long nicg;
  ///  residuum
  double aerrcg;
  ///  number of performed iterations
  long anicg;
  ///  computer zero
  double zero;
  
  ///  array containing nodal code numbers
  ///  (only on the master processor)
  ///  cn[number of layer][number of node][number of DOF]
  long ***cn;
  ///  array containing thicknesses in nodes
  ///  (only on the master processor)
  ///  thick[number of layer][number of node]
  double **thick;
  ///  array containing constraint matrix
  ///  (only on the master processor)
  ///  constrmat[number of multiplier][first or second component]
  double **constrmat;
  
  
  lgnode *lgnodes;
  
  matrix dcm;
};

#endif
