#ifndef CHEN_H
#define CHEN_H

#include <stdio.h>
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   class chen defines material model of plasticity based on Chen model
   
   structure of stored data; ncompother=ncompstr+2
   plastic strains, consistency parameter, hardening parametr
   
   MS,PS,JK
*/

class chen
{
 public:
  chen (void);
  ~chen (void);
  void read (FILE *in);
  void matstiff (matrix &d,long ipp,long ido);
  void tangentstiff (matrix &d,matrix &td,long ipp,long ido);

  double yieldfunction (matrix &sig,vector &q);
  void deryieldfdsigma (matrix &sig,vector &q,matrix &dfds);
  void deryieldfdsigma_old (matrix &sig,matrix &dfds);
  void deryieldfdsigmadsigma (matrix &sig,matrix &dfdsds);
  void deryieldfdsigma_old_old (matrix &sig,matrix &dfds);
  void deryieldfdq (matrix &sig,vector &q,vector &dfdq);
  void deryieldfdqdq (matrix &sig,vector &q,matrix &dfdqdq);
  void deryieldfdsigmadq (matrix &sig,vector &q,matrix &dfdsdq);
  
  void plasmod (long ipp,vector &epsp,matrix &sig,matrix &h);
  double plasmodscalar (vector &qtr);
  void updateq(long ipp,double dgamma,vector &epsp,matrix &sig,vector &q);
  void nlstresses (long ipp, long im, long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  void changeparam (atsel &atm,vector &val);
  double give_consparam (long ipp, long ido);
  long give_num_interparam ();
  void give_interparam (long ipp,long ido,vector &q);

  
  ///  compression yield stress
  double fyc;
  ///  tension yield stress
  double fyt;
  ///  bi-axial compression yield stress
  double fybc;
  ///  compression ultimate stress
  double fc;
  ///  tension ultimate stress
  double ft;
  ///  bi-axial compression ultimate stress
  double fbc;
  ///  limit strain
  double epsu;

  ///  coefficient A
  double ay;
  ///  coefficient A_u
  double au;
  ///  coefficient k
  double ky;
  ///  coefficient k_u
  double ku;
  ///  coefficient alpha
  double alpha;
  ///  coefficient beta
  double beta;

  ///  stress return algorithm
  strretalg sra;
  
  ///  auxiliary parameter
  long hard;

};

#endif
