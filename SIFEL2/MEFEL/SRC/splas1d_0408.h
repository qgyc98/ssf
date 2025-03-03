#ifndef SPLAS1D_H
#define SPLAS1D_H

#include <stdio.h>
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   class splas1d defines onedimensional model of plasticity
   it is used especially for program debugging
   
   JK
*/

class splas1d
{
 public:
  splas1d (void);
  ~splas1d (void);
  void read (FILE *in);
  double yieldfunction (matrix &sig,vector &q);
  void deryieldfsigma (matrix &sig,matrix &dfds);
  void deryieldfq (vector &dq);
  void plasmod (matrix &h);
  double plasmodscalar (vector &qtr);
  void   updateq(double dgamma, vector &q);
  void nlstresses (long ipp, long im,long ido);
  void nonloc_nlstresses (long ipp, long im,long ido);
  void matstiff (matrix &d,long ipp,long ido);
  void updateval (long ipp, long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  void changeparam (atsel &atm,vector &val);
  double give_consparam (long ipp,long ido);
  long give_num_interparam ();
  void give_interparam (long ipp,long ido,vector &q);
 
  ///  flow stress
  double fs;
  ///  plastic modulus
  double k;
  ///  stress return algorithm
  strretalg sra;
};

#endif
