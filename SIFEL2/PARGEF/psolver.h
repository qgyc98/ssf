#ifndef PSOLVER_H
#define PSOLVER_H

#include "pgalias.h"
#include "schurcompl.h"
#include "feti1.h"
#include "dpfeti.h"
#include "lplate.h"
#include "partop.h"
#include "partopjk.h"
#include "parcongrad.h"
#include "boundparcongrad.h"
#include "gtopology.h"
#include "gmatrix.h"
#include "xfile.h"
#include "selnodes.h"
#include "fixnodesel.h"
#include "vector.h"
#include <stdio.h>

/**
   class psolver
   
   this class controls parallel computation
   
   JK
*/

class psolver
{
 public:
  psolver (int np,int mr,char *nameproc,int nl);
  ~psolver ();
  void initiate (gtopology *top,int argc,const char**argv);
  void read (XFILE *in,gtopology *gt,long nd,storagetype tsm,long mes);
  void print (FILE *out);
  void read_ltg (XFILE *in,gtopology *top,FILE *out);
  void procdomcorr ();
  
  //void nodesplit (gtopology *top,FILE *out);
  
  long ordering (gtopology *top,FILE *out,char *proc_name);
  
  void constr_mat (double *th,FILE *out);
  
  void par_linear_solver (gtopology *top,gmatrix *gm,
			  double *lhs,double *rhs,FILE *out,long mespr);
  
  void gather_bound_vect (double *lv,double *gv);
  void scatter_bound_vect (double *lv,double *gv);
  double unbalanced_values (double *lv,FILE *out);
  
  double pss (double *lv1,double *lv2,FILE *out);
  long selected_norm_calculation (long cid,double err,double thresh,long sc,gtopology *top,double *lv1,double *gv1,double *rhs,double *grhs,long n,FILE *out, double &norfv);
  double compute_quantfluxres_norm(vector &lv, FILE *Out);

  void computation_stat (gtopology *top,gmatrix *gm);
  void computation_stat_print (FILE *out);
  long compare_on_master(double val, double limit);

 
  ///  DOMAIN DECOMPOSITION TYPE
  domdectype tdd;
  
  ///  type of reduced system solver
  redsystsolver trssol;
  
  ///  storage of redused system %matrix
  storagetype rsmstor;
  
  ///  data about linear system solver
  slesolv *ssle;
  
  ///  type of storage of domain %matrix
  storagetype dmstor;
  
  ///  description of the mesh
    ///  md=1 - nodes are denoted by the original numbers (numbers before mesh partitioning)
    ///  md=2 - boundary nodes are denoted by 
  meshdescription md;
  
  ///  topology for parallel computation
  partop *ptop;
  partopjk *ptopjk;
  
  selnodes *selnodschur;
  selnodes *selnodfeti;
  
  
  //  message printing
  long mespr;
  

  ///  local to global correspondence of node numbers
  long *ltg;
  ///  domain-processor correspondence
  ///  domproc[i]=j - the i-th processor contains the j-th subdomain
  long *domproc;


  ///  total number of processors
  int nproc;
  ///  rank of processor
  int myrank;
  /// name of processor
  char procName[10000];
  //char procName[MPI_MAX_PROCESSOR_NAME];
  ///
  int nameLength;
  
  ///  number of subdomain
  long ndom;
  
  ///  estimate of number of rigid body motions
  long enrbm;
  ///  threshold for linear dependence
  double lithr;
  ///  maximum number of iterations in conjugate gradient method
  long nicg;
  ///  number of performed iterations in conjugate gradient method
  long anicg;
  ///  required error
  double errcg;
  ///  attained error
  double aerrcg;
  /// type of preconditioner of FETI method
  long fetiprecond;

  long prec;
  
  /// condensation of corner node in FETI-DP method
  long condcorner;
  /// type automatic and user defined condensation of corner node in FETI-DP method
  long typecondcorner;

  long nmember;
   
  ///  computer zero
  double zero;
  
  
  ///  minimum number of nodes on subdomain
  long minnn;
  ///  maximum number of nodes on subdomain
  long maxnn;
  ///  minimum number of elements on subdomain
  long minne;
  ///  maximum number of elements on subdomain
  long maxne;
  ///  minimum number of stored matrix entries on subdomain
  long minnegm;
  ///  maximum number of stored matrix entries on subdomain
  long maxnegm;
  
  ///  list of numbers of nodes on subdomain
  long *nndom;
  ///  list of numbers of elements on subdomain
  long *nedom;
  ///  list of numbers of stored matrix entries on subdomain
  long *negmdom;
  
  
  ///  schur complement method
  schurcompl *schcom;
  ///  one-level FETI method
  feti1 *f1;
  ///  dual-primal FETI method
  dpfeti *dpf;
  ///  layered plate method
  lplate *lpp;
  ///  parallel conjugate gradient method
  parcongrad *parcg;
  ///  boundary version of parallel conjugate gradient method
  boundparcongrad *boundparcg;

  fixnodesel *fixnodes;
    /// condensation of fixing node in FETI-DP method
  long condfixing;
  //// method of condensation
  long methodcondcor;
  /// type automatic and user defined condensation of fixing node in FETI-DP method
  long typecondcur;
  long typecondsurf;
  

  long nmembercur;
  long nmembersurf;
  // order of numerical intergration
  long intpointocur;
  long nuserdefnod;
  long *userdefnod;
  long nring;
  long *ring;

  double **gquantfluxres;


};

#endif
