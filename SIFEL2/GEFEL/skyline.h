#ifndef SKYLINE_H
#define SKYLINE_H

#include "matrix.h"
#include "vector.h"
#include "gtopology.h"

/**
   class skyline
   
   serves for %matrix storage called skyline or profile
   skyline storage is available for symmetric real matrices
   columns are stored from diagonal entry to the farthest non-zero off-diagonal entry
   
   basic data
   double array a where %matrix entries are stored
   long array adr where addresses of diagonal entries are stored
   long number n stands for number of rows (columns) of the %matrix
   long number negm stands for number of entries of %matrix
   
   basic relations
   array a has negm components
   array adr has n+1 components
   adr[n] is equal to negm
   
   matrices stored in skyline format can be eliminated
   
   JK
*/
class skyline
{
 public:
  skyline (void);
  ~skyline (void);
  void allocadr (long m);
  double* status ();
  long decomp ();
  void changedecomp ();
  void setfact ();
  void setnotfact ();

  void column_lengths_elem (long *cn,long ndofe);
  void column_lengths_mult (long *ncn1,long *ncn2,long *mcn,long nm);
  void column_lengths (gtopology *top);
  void addresses ();
  void neglobmat ();
  void diagaddresses ();
  void allocglomat ();
  void localize (matrix &b,long *cn);
  void localized (double *b,long *cn,long m);
  void glocalize (matrix &b,long *rcn,long *ccn);
  void mult_localize (long nm,long *ncn1,long *ncn2,long *mcn);
  void nullsky ();
  void copy_sky (skyline &sky);
  void initiate (gtopology *top,long ndof,long mespr);
  
  //  functions for cr->sky transfer deals with equal matrices stored in different formats
  void column_lengths_cr (long *cradr,long *ci);
  void mat_entries (long *cradr,long *ci,double *cra);
  void assemble_from_cr (long crn,long *cradr,long *crci,double *cra);
  
  //  functions for scr->sky transfer select some entries from whole matrix stored in scr
  void column_lengths_scr (long *scradr,long *scrci,long *se);
  void mat_entries_scr (long *scradr,long *scrci,double *scra,long *se);
  void assemble_from_scr (long *scradr,long *scrci,double *scra,long neq,long *se);
  
  
  void ldl_sky (double *x,double *y,double zero,long tc);
  //void ldl_sky_new (double *x,double *y,double zero,long tc);
  long auxmax(long i,long j);
  long auxmin(long i,long j);
  void ldl_sky3 (double *x,double *y,long tc);
  void ldl_sky4 (double *x,double *y,long tc);
  void ldl_sky_10 (double *x,double *y,double zero,long tc);
  void eliminuj_4i_rev(double *a,long n,long s);
  void napln_a(long i1,long i2,long band,double *pole,double *a,long *adr);
  void uloz_a(long i1,long i2,long band,double *pole,double *a,long *adr);
  void napln_b(long i1,long i2,long band,double *pole,double *b,long *adr);
  void uloz_b(long i1,long i2,long band,double *pole,double *b,long *adr);
  void faze1_block3(double *a,double *b,long n,long s,long n2,long s2);
  void faze2_block3(double *a,double *b,long n,long s,long n2,long s2);
  void blokove2(double *pole,long *adr,long n, long band);
 
  
  void mxv_sky (double *b,double *c);
  void utv (double *b,double *c);
  void ltv (double *b,double *c);
  void ldlmxv_sky (double *b,double *c);
  void addmat_sky (double c,skyline &sky);
  void addmat_gm (double c,gmatrix &gm);
  void scalmat_sky (double c);
  void ker (double *r,long &nse,long *se,long ense,double limit,long tc);
  void ldlkon_sky (double *b,double *c,double *x,double *y,long m,long tc);

  double ldlkoncount_sky (double *b,double *c,double *x,double *y,long m,long tc);

  void ldl_feti_sky (double *x,double *y,long nse,long *se,double zero);
  void ldl_a12block (double *block,long nrdof);
  void printmat (FILE *out);
  void printdiag (FILE *out);
  
  double give_entry (long ri,long ci);
  void add_entry (double e,long ri,long ci);
  long give_negm ();
  
  void diag_scale (double *d);
  void select_submatrix(skyline *smsky,long nsdof,long *sdof,FILE *out);
  
  void diag_check (double thr);
  

  ///  number of rows of the %matrix
  long n;
  ///  number of entries in the skyline, it is equal to the size of the array a
  long negm;
  ///  addresses of diagonal entries
  long *adr;
  ///  stored %matrix
  double *a;
  /**  decomposition (factorization) indicator
       decompid=0 - %matrix is not decomposed (factorized)
       decompid=1 - %matrix is decomposed (factorized)
  */
  long decompid;

};

#endif
