#ifndef CONSOL_H
#define CONSOL_H

#include "alias.h"
#include "xfile.h"
struct matrix;
struct vector;


class consol
{
 public:
  consol (void);
  ~consol (void);
  void read (XFILE *in);
  long numberOfConsol ();
  void nlstressesincr (long ipp);
  void nlstresses (long ipp);
  void phase1 (long ipp);
  void phase2 (long ipp);
  void matstiff (matrix &d, long ipp);
  void seps_time (matrix &screep,vector &sig, long nc, long ncc, long ipp);
  void get_hc (long nc, long ipp);
  void updateval (long ipp, long im, long ido);
  
  void give_dstresses_eqother (vector &dsigma,long ipp,long ido);
  

  //  array containing numbers of components of stress and strain tensors
  long *cncomp;
  //  total number of components of stress and strain tensors
  long tncomp;
  //  number of approximated functions on the element
  long napfun;
  //  stress/strain state
  strastrestate ssst;
  long nc;
  //  number of stresses on the element
  long ncc;
  //  stress/strain position for computing
  double vlivTCSum;
  double ccTime;
  double ddTime;
  double t0;
  double timemat;
  long nRetTime;
  long imat;
  long iwink;
  double cn;
  double gama;
  //  depht of basic line foundation
  double vv;
  double ni;
  double e0;
  double m;
  //  active depht
  double hpodl;
  //  max. active depht
  double hpodlmax;
  //  radius of area for stress from type of consolidations
  double rr;
  //  part of loading for type of consolidations
  double nf;
  //  coeff. depend on active depht
  double nh;
  double nc1;
  double nc2;
  double ncn;
};

#endif
