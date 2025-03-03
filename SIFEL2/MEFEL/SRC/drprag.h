#ifndef DRPRAG_H
#define DRPRAG_H

#include "iotools.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   This class defines Drucker-Prager plasticity model

   structure of stored data
   plastic strains, consistency parameter, hardening parametr

   28.3.2002
*/

class drprag
{
 public:
  drprag (void);
  ~drprag (void);
  void read (XFILE *in);
  void print (FILE *out);

  double yieldfunction (vector &sig, vector &q);
  void dfdsigma (vector &sig, vector &dfds, vector &q);
  void dfdsigmadsigma (vector &sig, matrix &dfdsds);
  void dgdsigma (vector &sig, vector &dgds, vector &q);
  void dgdsigmadsigma (vector &sig,matrix &dgdsds, vector &q);

  void dderpotsigma (vector &sig,matrix &ddgdds, vector &q);
  void dderpotsigma (matrix &sig,matrix &ddgdds);
  void deryieldfq(vector &qtr, vector &dfq);
  void der_q_gamma(vector &dqdg);

  double cohesion(vector &qtr);
  double plasmodscalar (vector &qtr);
  void updateq(long ipp, vector &epsp, vector &q);
  void matstiff (matrix &d, long ipp, long ido);
  void tangentstiff (matrix &tde,long ipp,long ido);
  void nlstresses (long ipp,long im,long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void hardvect (vector &hv);
  void updateval (long ipp, long im,long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  double give_consparam (long ipp, long ido);
  void changeparam (atsel &atm,vector &val);
  

  ///  friction angle
  double phi;
  ///  cohesion
  double c;
  ///  dilatation
  double psi;
  /// angle of linear hardening/softening
  double theta;
  /// limit cohesion
  double clim;


  ///  material constant alpha
  double alpha;
  ///  material constant alpha1
  double alpha1;
  /// material constant beta
  double beta;

  ///  stress return algorithm
  strretalg sra;
};

#endif
