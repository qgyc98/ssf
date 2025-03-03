#ifndef CHEN_H
#define CHEN_H

#include "xfile.h"
#include "strretalg.h"
#include "hardsoft.h"
struct matrix;
struct vector;

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
  void read (XFILE *in);
  void matstiff (matrix &d,long ipp,long ido);
  void tangentstiff (matrix &d,matrix &td,long ipp,long ido);

  double yieldfunction (vector &sig,vector &q);
  void deryieldfdsigma (vector &sig,vector &q,vector &dfds);
  void deryieldfdsigmadsigma (vector &sig,matrix &dfdsds);
  void deryieldfdq (vector &sig,vector &q,vector &dfdq);
  void deryieldfdsigmadq (vector &sig,vector &q,matrix &dfdsdq);
  void dhdsigma (vector &sigt,vector &q,vector &dhds);
  void dhdqpar (vector &sigt,vector &q,vector &dhdq);
  void dhdgamma (vector &dhdg);
  void hvalues (vector &sigt,vector &q,vector &h);
  

  void nlstresses (long ipp, long im, long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  void changeparam (atsel &atm,vector &val);
  double give_consparam (long ipp, long ido);
  long give_num_interparam ();
  void give_interparam (long ipp,long ido,vector &q);

  void hardsoftfunction (long ipp,vector &epsp,vector &q);
  
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

  ///  coefficient A
  double ay;
  ///  coefficient A_u
  double au;
  ///  coefficient k
  double ky;
  ///  coefficient k_u
  double ku;
  
  ///  stress return algorithm
  strretalg sra;
  ///  hardening/softening
  hardsoft hs;
  
  ///  auxiliary parameter
  long state;
  
  //  provizorium
  double epslim;
};

#endif
