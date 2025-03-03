#ifndef SEQFETI_H
#define SEQFETI_H

#include "skyline.h"
#include "seqselnodes.h"
#include "scr.h"
#include "densemat.h"

class gtopology;
class gmatrix;

/**
   class seqfeti deals with the FETI method
   (Finite Element Tearing and Interconnecting method)
   this is single-processor implementation, it is different
   from feti1 in PARGEF which is multi-processor implementation
   
   this class contains three different implementations of the FETI method
   1. implementation based on the Boolean matrices, the matrices
   define ordering of unknowns and multipliers
   2. implementation without the Boolean matrices with nonredundant conditions
   correspondence among nodes is used for deinition of unknowns
   3. implementation without the Boolean matrices with redundant conditions
   correspondence among nodes is used for deinition of unknowns
   
   
   JK, 21.8.2007
*/
class seqfeti
{
 public:
  seqfeti ();
  ~seqfeti (); 
  
  void read (gtopology *top,XFILE *in);
  void print (FILE *out);

  void read_booldata (XFILE *in);
  void initiate (seqselnodes *selnodfeti,gtopology *top,FILE *out);

  void det_ndofmax ();
 
  void assemble_subdom_unknowns (gtopology *top,FILE *out);
  void subdomain_matrices (gmatrix *gm,FILE *out);
  

  void boolean_matrix (long sdid,double *a);
  void local_coarse (long nd,double *lv,double *cv);
  void coarse_local (long nd,double *lv,double *cv);

  void kernel (FILE *out);


  void g_matrixsize (FILE *out);
  void g_matrix (FILE *out);

  void inverse_matrix_GG (FILE *out);

  void assemble_ff (double *rhs);
  void evector (FILE *out);


  void qvector (double *q,double **ff);
  
  void feti_projection (double *v);


  //void get_jumps (double *jump);


  void scaling(double *invect,double *outvector,long n,FILE *out);
  void lscaling(double *invect,double *outvect,long ndom,FILE *out);
  
  void lumpedprec (long nd,double *dd,double *pp);
  void dirichletprec (long nd,double *dd,double *pp);

  void mpcg (FILE *out);
  void mpcg_old (gtopology *top,gmatrix *gm,double *w,double **ff,double *q,double *h,double *h1,FILE *out);
  void mprcg (gtopology *top,gmatrix *gm,
	      double *lambda,double **ff,double *e,double *g,double *g1,FILE *out);

  void nodalunknowns (double *lhs,FILE *out);
 
  void lagrmultnodunknowns (FILE *out);
  
  void matrices_assembl (gmatrix *gm,FILE *out);
  void vectors_assembl (double *rhs,FILE *out);

  void define_h (FILE *out);
  void define_b (FILE *out);

  void solve_system (gtopology *top,gmatrix *gm,double *lhs,double *rhs,FILE *out);




  ///  type of FETI implementation
  ///  fetiimpl=no_impl=0 - no implementation is defined
  ///  fetiimpl=boolean_matrices=1 - Boolean %matrix is assembled, it is used for tests only, it is not efficient implementation
  ///  fetiimpl=nonredundant=2 - nonredundant constraints are defined, %matrix B has linearly independent rows
  ///  fetiimpl=redundant=3 - redundant constraints are defined, %matrix B has linearly dependent rows
  fetiimplem fetiimpl;
  
  ///  indicator of the presence of discontinuities
  ///  discont=0 - the classical FETI method
  ///  discont=1 - the FETI method with interface discontinuities
  long discont;

  ///  type of storage of subdomain matrices
  ///  see galias.h
  storagetype smst;
  
  ///  indicator of assembled matrices
  ///  matassem=0 - matrices are not assembled
  ///  matassem=1 - matrices are assembled
  long matassem;
  ///  indicator of assembled right hand side vectors
  ///  vecassem=0 - vectors of right hand sides are not assembled
  ///  vecassem=1 - vectors of right hand sides are assembled
  long vecassem;
    
  
  ///  the number of subdomains
  long ns;

  ///  estimated number of rigid body modes (estimated dimension of the kernel)
  long ense;

  ///  threshold for kernel detection
  double thresh;
  
  ///  the maximum number of iterations in conjugate gradient method
  long nicg;
  ///  the number of performed iterations in conjugate gradient method
  long anicg;

  ///  required error
  double errcg;
  ///  attained error
  double aerrcg;

  ///  computer zero
  double zero;
  
  ///  type of preconditioner
  precondtype prec;


  ///  size of the %matrix G
  ///  it is determined in the function g_matrixsize
  long gsize;

  ///  the maximum number of degrees of freedom on one subdomain
  ///  it is determined in the function det_ndofmax
  long ndofmax;

  ///  the number of DOFs (unknowns) in coarse problem
  ///  it is defined in the function initiate
  ///  it is equal to the variable selnodfeti->tndofsn
  long ndofcp;

  ///  the numbers of DOFs on subdomains
  ///  it contains ns components
  ///  ndofmas[i]=j - the i-th subdomain contains j DOFs
  ///  array is assembled in the function assemble_subdom_unknowns
  long *ndofmas;
  
  ///  node-subdomain correspondence
  ///  it contains ns components
  ///  nsid[i]=j - the i-th node belongs to the j-th subdomain
  ///  array is assembled in the function assemble_subdom_unknowns
  long *nsid;
  
  ///  list of DOFs on subdomains
  ///  it is used in connection with preconditioning
  ///  it contains ns, ndofmas[i] components
  ///  cndom[i][j]=k - the j-th DOF on the i-th subdomain has global glued number k
  ///  array is assembled in the function assemble_subdom_unknowns
  long **cndom;
  
  ///  array of numbers of unknowns (DOFs) contributing to the coarse problem
  ///  ncdofd contains nproc components
  ///  ncdofd[i]=j - the i-th subdomains contributes to coarse problem by j contributions
  ///  ncdofd is a copy of selnodfeti->snndofmas
  ///  array is assembled in the function initiate
  long *ncdofd;
  
  ///  array containing code numbers contributing to the coarse problem
  ///  extracted values from subdomains to the coarse problem
  ///  it contains ns, ncdofd[i] components
  ///  edofs[i][j]=k - the j-th components contributing to the coarse problem from the i-th subdomains has number k
  ///  array is assembled in the function initiate
  long **edofs;



  ///  array containing numbers of RBM on subdomains
  ///  it contains ns components
  ///  nrbmdom[i]=j - the i-th subdomain contains j rigid body modes
  ///  array is assembled in the function kernel
  long *nrbmdom;

  ///  rigid body modes / kernel
  ///  array is assembled in the function kernel
  double **rbmdom;
  
  ///  array containing addresses of first RBM in coarse matrix
  ///  it contains ns+1 components
  ///  rbmadr[i]=j - rigid body modes of the i-th subdomains start from the index j
  ///  array is assembled in the function hmatrixsize
  long *rbmadr;

  ///  list of linearly dependent equations
  ///  se[i][j]=k - the j-th base %vector of the kernel of the i-th subdomain creates the k-th column of the %matrix
  ///  array is assembled in the function kernel
  long **se;

  ///  %matrix G
  ///  rigid body modes are stored in columns of the %matrix G
  ///  it contains ndofcp,gsize components
  ///  array is assembled in the function gmatrix
  double *g;


  ///  the number of contributions in the arrays booldatar, booldatac and booldata
  ///  it contains ns components
  ///  it is assembled in the function read_booldata
  long *ncbool;

  ///  array containing row indices for construction of the Boolean matrices
  ///  it contains ns,ncbool[i] components
  ///  booldatar[i][j] = k - the j-th interface unknown on the i-th subdomain contributes to the k-th row
  ///  it is assembled in the function read_booldata
  long **booldatar;

  ///  array containing column indices for construction of the Boolean matrices
  ///  it contains ns,ncbool[i] components
  ///  booldatac[i][j] = k - the j-th interface unknown on the i-th subdomain contributes to the k-th column
  ///  it is assembled in the function read_booldata
  long **booldatac;

  ///  array containing %matrix entries of the Boolean matrices
  ///  it contains ns,ncbool[i] components
  ///  booldata[i][j] = k - the j-th interface unknown on the i-th subdomain contributes to the Boolean matrix with the value k
  ///  it is assembled in the function read_booldata
  double **booldata;
  
  ///  array of coarse code numbers
  ///  it contains tnbn rows and ncdofd[i] columns
  ///  ccn[i][j]=k - the j-th contribution from the i-th subdomain goes to the k-th coarse unknown
  long **ccn;
  
  
  ///  array for the inverse %matrix to the %matrix G^T G
  ///  it contains gsize rows and columns
  double *invgg;
  
  ///  array for the right hand sides
  ///  it contains ns vectors, each vector contains ndofmas[i] entries
  ///  ff[i][j]
  double **ff;
  
  ///  array for the e %vector
  ///  it contains gsize entries
  double *e;
  
  ///  array of nodal variables on subdomains
  ///  it contains ns vectors, the vectors contain ndofmas[i] entries
  double **d;
  
  ///  array of Lagrange multipliers
  double *lambda;
    
  ///  subdomain matrices stored in the %skyline storage
  skyline *smsky;
  ///  subdomain matrices stored in the %dense format
  densemat *smdm;
  
  ///  array for the %vector b
  ///  it is a constant %vector containing prescribed discontinuities
  double *b;
  
  ///  array for the compliances in the %matrix H
  double *h;
  




  
  
  //  nize uvedene promenne nejsou zkontrolovany

  

  //long ndof;


  
  

  ///  numbers of nodes on subdomains
    //long *nnsd;
  ///  first node numbers on subdomains
    //long *fnnsd;
  
  


  
  
  
  
  
  ///  number of DOFs on subdomains used for preconditioning
  long *ndofprec;
  
  ///  code numbers for preconditioning
  ///  cnprec[i][j]=k - the j-th DOF on the i-th subdomain has global glued number k
  long **cnprec;
  
  ///  cpreccn[i][j]=k - the j-th DOF on the i-th subdomain has local number k
  long **cpreccn;

  double *wscalmat;
  long *nlwscalmat;
  double **lwscalmat;
  



  

  ///   subdomain matrices stored in the %skyline storage used for precondition
  skyline *smskyprec;
  
  
  ///  subdomain matrices for preconditioning
  ///  %matrix in the compressed row storage scheme
  symcomprow *smscr;
  ///  %matrix in the dense storage scheme (it contains Schur complements which are dense)
  densemat *psmdm;

  ///  array containing jumps
    //double **jum;
};

#endif
