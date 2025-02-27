#ifndef SHEFPLAST_H
#define SHEFPLAST_H

#include "iotools.h"
#include "alias.h"
#include "strretalg.h"
#include "matrix.h"
//struct matrix;
struct vector;

/**
   This class defines Sheffield plastic material model

   14.1.2004
   Frederic Dufour, Tomas Koudelka
*/

class shefplast
{
 public:
  shefplast (void);
  ~shefplast (void);
  void read (XFILE *in);

  void compute_khat (vector &q);
  void compute_rhoc ();
  void compute_r ();
  void compute_classvar (vector &sig, vector &q);
  double yieldfunction (vector &sig, vector &q);
  void deryieldfsigma (long ipp,vector &sig, vector &q, vector &dfds,long ido);
  void dderyieldfsigma (matrix &ddfds, strastrestate ssst);
  void derqdsigma (vector &sig, vector &drds);
  void derpotsigma (long ipp, vector &sig, vector &q, vector &dgds,long ido);
  void deryieldfq (vector &sig, vector &q, vector &dq);
  void derkhatds (long ipp, vector &sig, vector &q, vector &dkhatds,long ido);
  //  void der_q_gamma(long ipp, matrix &sig, vector &eps, vector &qtr, vector &dqdg);
  //  double plasmodscalar(long ipp, matrix &sig, vector &eps, vector &qtr);
  void updateq(long ipp, vector &eps, vector &sig, vector &q,long ido);
  //  void updateq(vector &epsp, strastrestate ssst, vector &q);
  void matstiff (matrix &d, long ipp, long ido);
  void tangentstiff (matrix &d, matrix &td, long ipp,long ido);
  void nlstresses (long ipp,long ido);
  void updateval (long ipp, long ido);
  double give_consparam (long ipp, long ido);

  void stress_return (long ipp,double &lambda,vector &k,vector &eps,vector &epsp, long ido);
  void compzeta();
  double hardening(long ipp,vector &sigtens,vector &q,long ido);
  void numdiff_dfdsds(long ipp,vector &sigtens,vector &q,matrix &dfdsds,long ido);
  void numdiff_dfdsdsc(long ipp,vector &sigtens,vector &q,matrix &dfdsds,long ido);
  void numdiff_dfdsdk(long ipp,vector &sigtens,vector &q,vector &dfdsdk,long ido);
  void numdiff_dfdsdkc(long ipp,vector &sigtens,vector &q,vector &dfdsdk,long ido);
  void numdiff_dhds(long ipp,vector &sigtens,vector &q,vector &dhds,long ido);
  void numdiff_dhdsc(long ipp,vector &sigtens,vector &q,vector &dhds,long ido);
  void numdiff_dhdk(long ipp,vector &sigtens,vector &q,double &dhdk,long ido);
  void numdiff_dhdkc(long ipp,vector &sigtens,vector &q,double &dhdk,long ido);

  double maxim (double a,double b);

  /// 
  double rc;
  ///
  double rt;
  ///
  double gamma;
  ///
  double p;
  ///
  double a;
  ///
  double k0;
  /// 
  double alpha;
  ///
  double ah, bh, ch;
  /// initial value of hardening parameter
  double kh0;

  // Internal temporary variables

  /// stress deviator
  vector dev;
  /// the first invariant of the stress tensor
  double i1s;
  /// the second invariant of the stress deviator
  double j2s;
  /// the third invariant of the stress deviator
  double j3s;
  double khat, rhoc, r, m, xi, ft;
  double b0, b1, c0, c1, c2, c3, c4, c5, theta, d0, d1, d2;
  
  
  double zeta;

  ///  stress return algorithm
  strretalg sra;
};

#endif
