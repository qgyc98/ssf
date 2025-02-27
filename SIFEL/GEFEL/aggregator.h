#ifndef AGGREGATOR_H
#define AGGREGATOR_H

#include "gtopology.h"
#include "gmatrix.h"
#include "csv.h"
#include "densemat.h"
#include "skyline.h"
#include "cr.h"
#include "slesolv.h"
#include "SPARSE/DSSolver.h"

/**
   class aggregator deals with preconditioning of solvers of systems of linear algebraic equations
   the preconditioner is based on the BOSS method (Black-box Overlapping Schwarz with Smoothed Coarse Space)
   the BOSS method was developed ba M. Brezina in his Ph.D. thesis
   
   JK
*/
class aggregator
{
 public:
  //konstrukce a destrukce
  aggregator(void);
  ~aggregator(void);
  
  void read (gtopology *gt,XFILE *in,long mespr);
  void print (FILE *out);
  

  //   JK -> PM

  ///  defines number of aggregates
  void define_na (long i);

  /// vytvori graf matice
  /// nn        - pocet uzlu
  /// nneigh    - pocty sousednich uzlu
  /// listneigh - seznamy sousednich uzlu
  void define_graph(long nnod, long *nneigh, long **listneigh);

  
  //  PM -> JK

  ///  assembles list of numbers of nodes in aggregates without overlap
  ///  pocet uzlu v agregatech bez prekryvu
  ///  assembles node numbers in aggregates without overlap
  ///  seznamy uzlu v jednotlivych agregatech bez prekryvu
  ///  lnntagr[i] = k -  v i-tem agregatu je k uzlu
  ///  lntagr[i][j] = k - j-ty uzel v i-tem agregatu ma cislo k
  void prepare_tlnagr(void);

  ///  assembles list of numbers of nodes in aggregates with overlap
  ///  pocet uzlu v agregatech s prekryvem
  ///  assembles node numbers in aggregates with overlap
  ///  seznamy uzlu v jednotlivych agregatech s prekryvem
  ///  lnnagr[i] = k -  v i-tem agregatu je k uzlu
  ///  lnagr[i][j] = k - j-ty uzel v i-tem agregatu ma cislo k
  void prepare_lnagr ();
  

  ///  function sets up degree of recursion
  ///  real degree of polynomial is equal to (3^dg-1)/2
  void define_degree ();
  

  //  sestavuje tentative prolongator
  //  agrid - cislo pozadovaneho agregatu
  //  cid - cislo pozadovaneho sloupce v agregatu
  //  i - pole indexu
  //  a - pole hodnot (v pripade rigid body motions v mechanice)
  //void assemble_tentative_prol (long agrid,long cid,long *i,double *a);
  
  //  vraci skutecny stupen prolongatoru
  long give_deg ();


  //  sestavuje smoothed prolongator
  //  agrid - cislo pozadovaneho agregatu
  //  cid - cislo pozadovaneho sloupce v agregatu
  //  i - pole indexu
  //  a - pole hodnot (v pripade rigid body motions v mechanice)
  void assemble_smoothed_prol (gtopology *gt,gmatrix *mtx);
  void assemble_smoothed_prol2 (long agrid,long cid,long *i,compvect *cvi,compvect *cvo,gmatrix *mtx);
  
  
  
  ///  function copies arrays lnntagr and lntagr from the
  ///  object stop (class seqtop) allocated in class gtopology
  ///  it is an alternative to the function void prepare_tlnagr(void);
  void metis_aggr (gtopology *gt);
  
  void gener_rbm_1 (gtopology *gt,long agid,double *v);
  void gener_rbm_2 (gtopology *gt,long agid,long cnid,double *v);
  void coarse_matrix (gmatrix *gm);
  
  
  void assemble_aggr_unknowns (gtopology *gt);

  void local_matrices (long bsize,gmatrix *gm);

  void clean_memory ();
  
  void prepare_boss (gtopology *gt,gmatrix *gm,FILE *out);
  
  //void boss (gmatrix *gm,double *u,double *v,double *rhs);
  void boss (gmatrix *gm,double *u,double *v);



 /// lokalni metody - casem asi private
  /**
     provedena aplikaci s(st) na pozadovany vekor
     st    .... stupen S, t.j. dolni index
     x     .... vektor
     mtx   .... matice
  **/
  void mulS(long st,double *x, gmatrix *mtx);
  /** 
     st    .... stupen S, t.j. dolni index
     x     .... vektor
     mtx   .... matice
  **/
  void mulA(long st, double *x, gmatrix *mtx);

  
  ///  number of nodes in the finite element mesh
  long nn;
  ///  number of aggregates
  long na;
  ///  maximum number of unknowns in aggregate
  long maxnu;
  ///  number of unknowns in the whole problem
  long n;
  ///  size of the coarse matrix
  long cms;

  ///  type of BOSS algorithm
  ///  impl=1 - own implementation of the BOSS algorithm
  ///  impl=2 - implementation based on METIS
  long impl;
  
  ///  type of solver - exact or inexact
  ///  exinex=1 - exact solver is used
  ///  exinex=2 - inexact solver is used
  long exinex;
  
  
  ///  type of factorization of system of linear algebraic equations
    //linsolvertype tlinsol;
  
  ///  type of rigid body modes (kernels)
  long trbm;
  ///  number of rigid body modes of one aggregate
  long nrbm;
    
  ///  data about solver of system of linear equations
  slesolv *ssle;

  ///  list of numbers of nodes in aggregates with overlap
  ///  lnnagr[i]=j - the i-th aggregate contains j nodes
  long *lnnagr;
  ///  list of node numbers in aggregates with overlap
  ///  lnagr[i][j]=k - the j-th node in the i-th aggregate has the number k
  long **lnagr;

  ///  list of numbers of nodes in tentative aggregates (without overlap)
  long *lnntagr;
  ///  list of node numbers in tentative aggregates (without overlap)
  long **lntagr;
   
  ///  numbers of adjacent nodes to nodes
  ///  nadjnodnod[i]=j - the i-th node has j adjacent nodes
  long *nadjnodnod;
  ///  array of adjacent nodes to nodes
  ///  adjnodnod[i][j]=k - the j-th adjacent node to the i-th node has number k
  long **adjnodnod;

  ///  list of numbers of unknowns on aggregates with overlap
  ///  lnuaggr[i]=j - the i-th aggregate contains j unknowns
  long *lnuaggr;
  ///  list of unknown numbers on aggregates with overlap
  ///  luaggr[i][j]=k - the j-th unknown on the i-th aggregate has number k
  long **luaggr;
  

  // odhad spektralniho polomeru matice 
  double specrad;
  
  // degree of smoothing aggregator - real
  long degree_r;
  long degree_k;
  
  //  local matrices
  comprow *lmcr;
  
  skyline *lmsky;
  
  ///  smoothed prolongator
  ///  it is a %matrix with ndof rows and na*nrbm columns
  ///  the matrix is stored column after column in onedimensional array
  double *p;

  ///  coarse %matrix
  densemat *cm;

  //  provizorne pro ladeni
  densemat *lmdm;

  ISolver *sdirect;

};

#endif
