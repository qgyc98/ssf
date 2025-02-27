#ifndef CREEPBBEAM_H
#define CREEPBBEAM_H

#include "alias.h"
struct matrix;
struct vector;


class creepbbeam
{
 public:
  creepbbeam (void);
  ~creepbbeam (void);
  void creepinit (long mie);
  void read (FILE *in);
  double approx (vector &areacoord,vector &nodval);
  void inv_sym (matrix &a);
  void nlstresses (long ipp);
  void phase1 (long ipp);
  void phase2 (long ipp);
  void get_h (long ipp);
  void get_temp (long ipp);
  void matstiff (matrix &d, long ipp);
  void seps_time (matrix &screep,vector &sig);
  void b2_law (double &jt,double &esht, double t0, double t);
  void b3_law (double &jt,double &esht, double t0, double t);


  // =1 constant h, 
  long type_h;
  // =1 constant temperature, 
  long type_temp;
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
  vector retTime;
  vector ert;
  double esht;
  double deps0;
  double ccTime;
  double ddTime;
  double t0;
  double timemat;
  double timeMax;
  long nRetTime;
  long imat;
  double e0;
  double mi;
  double alfa;
  double tb;
  double t_w;
  double fc;
  double wc;
  double sc;
  double gc;
  double c_s;
  double a1;
  double h_s;
  double h_slast;
  double temp_s;
  double templast;
  double k_s;
  double r_s;
  double ts;
  double k_d;
};

#endif
