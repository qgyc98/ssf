#ifndef CREEPB_H
#define CREEPB_H

#include "alias.h"
#include <stdio.h>
struct matrix;
struct vector;

/**
   class contains B3 model of creep
   
   from t0  to t 
   K_s  shape factor slab=1.0, cylinder=1.15, sguare prism.=1.25, sphere=1.3, cube=1.55 
   tb   from concrete starts
   t_w  age when drying begins 
   (fc') is 28 day average cilinder strenght fc' [ksi] ksi=1000psi=6.895 MPa(f.e.6.454=44.5MPa)***6.381
   (w/c) is water-cement ratio of the mix by weight   ***0.43
   (s/c) is send-cement ratio of the mix by weight    ***3.4
   (g/c) is gravel-cement ratio of the mix by weight g/c=a/c-s/c     ***1.98
   (a/c) is aggregate-cement ratio of the mix by weight a/c=g/c+s/c   
   (a1)  is coef. for cements of type I,II a1=1.00, III a1=0.93, IV a1=1.05   ***1.05
   (ro)  is mass of concrete in [lb/ft3] =16.03 kg/m3 ***156
   (k_d)  effective cross section thickness       D=2*vs_s
   cs  cement content in m3  .. kg/m3
   E0=(0.09+1/(1.7*(0.5*ro*ro*fc*1e-4)*(0.5*ro*ro*fc*1e-4)))
   Et=E0*sqrt(t/(4+0.85*t))          podle ACI Commite 209/II
   
   components of other array:
   previous total strains (nc components)
   internal variables describing history (7 x nc components for this model)
   shrinkage and thermal strain (1 component)
   previous moisture (1 component)
   previous temperature (1 component)
   
*/
class creepb
{
 public:
  creepb (void);
  ~creepb (void);
  void creepinit (long ipp,double val,nonmechquant nmq);
  void read (FILE *in);
  double approx (vector &areacoord,vector &nodval);
  void inv_sym (matrix &a);
  void updateval();
  void nlstresses (long ipp);
  void phase1 (long ipp);
  void phase2 (long ipp);
  void get_h (long ipp);
  void get_temp (long ipp);
  void matstiff (matrix &d, long ipp);
  void seps_time (matrix &screep,vector &sig);
  void get_desht (double &des_hn, double t0, double t);
  void b3_law (double &jt, double t0, double t);


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
  strastrestate ss;
  //  number of stresses on the element
  long nc;

  vector retTime;
  vector ert;
  double desht;
  double ccTime;
  double ddTime;
  double t0;
  double timemat;
  double timeMax;
  long nRetTime;
  long imat;
  //  Young's modulus of elasticity
  double e0;
  //  Poisson's ratio
  double mi;
  //  coefficient of thermal dilatancy
  double alfa;
  //  time of casting
  double tb;
  
  double t_w;
  //  compression strength in MPa
  double fc;
  //  water-cement ratio
  double wc;
  //  sound/cement ratio
  double sc;
  //  gravel-cement ratio
  double gc;
  double c_s;
  //  coefficient of shape of structure
  double a1;
  //  relative moisture
  double h_s;
  double h_slast;
  //  temperature
  double temp_s;
  double temp_slast;
  //  
  double k_s;
  //
  double r_s;
  //
  double ts;
  //
  double k_d;
  
  matrix apom;
};

#endif
