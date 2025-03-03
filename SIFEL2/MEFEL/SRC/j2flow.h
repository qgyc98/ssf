#ifndef J2FLOW_H
#define J2FLOW_H

#include "iotools.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   class j2flow defines material model of plasticity based on J2 flow
   
   structure of stored data
   plastic strains, consistency parameter, hardening variable
   
   JK
*/

class j2flow
{
 public:
  j2flow (void);
  ~j2flow (void);
  void read (XFILE *in);
  void print (FILE *out);

  void matstiff (matrix &d,long ipp,long ido);
  void tangentstiff (matrix &d, matrix &td, long ipp,long ido);

  double yieldfunction (vector &sig,vector &q);
  void dfdsigma (vector &sig,vector &dfds);
  void dfdqpar (vector &dq);
  void dfdsigmadsigma (vector &sig,matrix &dfdsds);
  void dfdsigmadq (matrix &dfdsdq);
  void dfdqpardqpar (matrix &dfdqdq);
  void hardvect (vector &hv);

  double plasmodscalar ();
  void updateq(double dgamma, vector &q);
  void nlstresses (long ipp, long im, long ido);
  void nonloc_nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  void changeparam (atsel &atm,vector &val);
  double give_consparam (long ipp, long ido);
  long give_num_interparam ();
  void give_interparam (long ipp,long ido,vector &q);
  
  ///  flow stress
  double fs;
  ///  hardening parameter
  double k;
  ///  stress return algorithm
  strretalg sra;
};

#endif
