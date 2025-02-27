#ifndef DENSEMAT_H
#define DENSEMAT_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"

class precond;

/**
   class densemat
   
   this class serves for storage of dense matrices
   
   JK
*/

class densemat
{
 public:
  densemat (void);
  ~densemat (void);
  void alloc (long m);
  void dealloc(void);
  void copy (densemat *dm);
  void copy_dm (densemat &dm);
  double* status ();
  long decomp ();
  void changedecomp ();
  void setfact ();
  void setnotfact ();
  void nullmat ();
  void localize (matrix &b,long *cn);
  void localized (double *b,long *cn,long m);
  void glocalize (matrix &b,long *rcn,long *ccn);
  void mult_localize (long nm,long *ncn1,long *ncn2,long *mcn);
  void initiate (long ndof,long mespr);
  
  void mxv_dm (double *b,double *c);
  void addmat_dm (double c,densemat &dm);
  void scalmat_dm (double c);
  void maxmb(long nm, double *ma, double *mb, double *mc);
  void gemp (double *x,double *y,long m,double limit,long pivot);
  void gempkon (double *b,double *c,double *x,double *y,long m,double zero,long tc);
  void lu (double *x,double *y,double zero,long tc);
  void ll (double *x,double *y,double zero,long tc);
  void ill (double *x,double *y,double zero,double limit,long tc);
  void ker (double *r,long &dim,long *se,long ense,double limit);
  
  
  void cg (double *x,double *y,long ni,double res,long &ani,double &ares,double zero,long iv);
  void cg_prec (densemat &dm,double *x,double *y,long ni,double res,long &ani,double &ares,
		double zero,long iv,long tprec,double par);
 
  void cg_prec_new (precond &pr,double *x,double *y,long ni,double res,long &ani,double &ares,
		    double zero,long iv);

  void jacobi (double *x, double *y,long ni,double res,long &ani,double &ares,FILE *out);
  void gauss_seidel (double *x, double *y,long ni,double res,long &ani,double &ares,FILE *out);
  
  
  void printmat (FILE *out);
  void printdiag (FILE *out);
  double give_entry (long ri,long ci);
  void add_entry (double e,long ri,long ci);
  long give_negm ();

  void diag_scale (double *d);

  double estim_spect_radius ();
  
  void assemble_dense_from_scr(long *scradr,long *scrci,double *scra,long neq);
  
  ///  number of rows (columns)
  long n;
  ///  number of entries in the global %matrix
  long negm;
  ///  array containing global %matrix
  double *a;
  ///  decomposition (factorization) indicator
  long decompid;
  
  ///   array containing pseudoinverse (inverse) of global %matrix
  double *a_plus;
  ///   number of base vectors of null space of global %matrix
  long nbvns;
  ///   base vectors of null space of global %matrix
  double *bvns;
  ///   address of first entry of i-th base vector in bvns
  long *abvns;
  ///   array with number of row indices of dependent rows in global matix
  long *iirm;
  
};

#endif
