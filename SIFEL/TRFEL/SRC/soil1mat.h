#ifndef SOIL1MAT_H
#define SOIL1MAT_H

#include "genfile.h"

class soil1mat
{
 public:
  soil1mat();    //constructor
  ~soil1mat();   //destructor

  double _dg();
  double _alpha();
  double _ks();
  double _kt();
  double _emod();
  double _nu();
  double _betas();
  double _rhocp();
  double _lambdaa();

  double sat(double pc,double t);
  double dsat_dpc(double pc,double t);
  double dsat_dt(double pc,double t);

  double _krg(double pc,double t);
  double _krw(double pc, double t);
  double _phi();
  double _kintr();
  double _cps();
  double _rhos();

  void read(XFILE *in);
  void print(FILE *out);

 private:
  double mw;
  double ma;
  double gasr;
  
  double t0;
  double p0;
  double tcr;
};  

#endif
