#ifndef SPLAS1D_H
#define SPLAS1D_H

//#include <stdio.h>
#include "iotools.h"
#include "strretalg.h"
#include "hardsoft.h"
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

  void read (XFILE *in);
  void print(FILE *out);
  
  double yieldfunction (vector &sig,vector &q);
  
  void dfdsigma (vector &sig,vector &dfds);
  void dfdsigmadsigma(matrix &dfdsds);
  void dfdqpar (vector &dq);
  void dfdsigmadq (matrix &dfdsdq);
  void dfdqpardqpar (matrix &dfdqdq);
  //void dhdgamma (long ipp,vector &epsp, vector &sig,vector &dhdc);
  void hardvect (vector &hv);
  
  
  
  
  void plasmod (matrix &h);
  double plasmodscalar (vector &sig,vector &epsp,vector &qtr,double gamma);
  void updateq(long ipp,double dgamma, vector &epsp, vector &q);
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
  ///  hardening/softening
  hardsoft hs;

};

#endif
