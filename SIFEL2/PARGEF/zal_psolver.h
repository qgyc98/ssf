#ifndef PSOLVER_H
#define PSOLVER_H

#include <stdio.h>
#include "mpi.h"
#include "paral.h"
#include "../GEFEL/gtopology.h"
#include "../GEFEL/gmatrix.h"

typedef enum ldomdectype {lprimaldd=1,lfetidd=5};
typedef enum lredsystsolver {lpldl=1,lplu=2,lpgemp=3,lpcg=5};

class psolver
{
 public:
  psolver (int np,int mr,int nd);
  ~psolver ();
  void movedata (paral *plg);
  
  void redsys_parldl (paral *plg,double *condmat,double *condvect,FILE *out);
  void redsys_parlu (paral *plg,double *condmat,double *condvect,FILE *out);
  void redsys_pargemp (paral *plg,double *condmat,double *condvect,FILE *out);
  void redsys_parcg (paral *plg,double *condmat,double *condvect,FILE *out);

  void hmatrixsize (paral *plg,double *rbm,long &maxnrbm);
  void hmatrix (gtopology *top,paral *plg,double *h,double *rbm,long maxnrbm);
  void qvector (paral *plg,double *q,double *rbm,double *f,long maxnrbm);
  void feti_projection (double *v,double *h,double *h1);
  void mpcg (gtopology *top,paral *plg,gmatrix *gm,
	     double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out);
  void lagrmultdispl (gtopology *top,paral *plg,gmatrix *gm,
		      double *w,double *d,double *f,double *rbm,long *rbmi,
		      double *h,double *h1);
  void lumpedprec (gtopology *top,double *v,double *pv);
  
  void par_linear_solver (gtopology *top,paral *plg,gmatrix *gm,
			  double *lhs,double *rhs,FILE *out);
  
  
  //  DOMAIN DECOMPOSITION TYPE
  ldomdectype tdd;
  
  //  REDUCED SYSTEM SOLVER
  lredsystsolver rssol;

  //  SUBDOMAIN DESCRIPTION
  //  number of global unknowns
  long ngdof;
  //  total number of unknowns on subdomain
  long ndof;
  //  number of internal unknowns on subdomain
  long indof;
  //  dimension of kernel = number of rigid body modes
  long nrbm;
  //  estimate of kernel dimension
  long enrbm;
  //  threshold for linear independence
  double lithr;

  //  total number of processors
  int nproc;
  //  rank of processor
  int myrank;
  //  number of subdomain
  int ndom;

  //  CONJUGATE GRADIENT METHOD
  //  maximum number of iterations in conjugate gradient methods
  long nicg;
  //  number of performed iterations in conjugate gradient methods
  long anicg;
  //  required error in conjugate gradient methods
  double errcg;
  //  actual residual after iteration in conjugate gradient methods
  double aerrcg;
  
  //  size of corse problem in FETI method
  long hsize;
  
  //  threshold for rejection from compressed storages
  double limit;
  //  computer zero
  double zero;
};

#endif
