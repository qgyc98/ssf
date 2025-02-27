#ifndef CREEPEFFYM_H
#define CREEPEFFYM_H

#include "alias.h"
#include "iotools.h"
#include "gfunct.h"
#include <stdio.h>
struct matrix;
struct vector;

/**
  This class defines material model for simple creep modelling
  which uses effective Young modulus approach. The model is
  reasonable only with elasttime material model.
  
  Order of internal variables in the other array :
  -----------------------------------------------
  0          : Young's modulus from the previous time step
  1-ncompstr : computed strains due to ageing
*/
class creep_effym
{
 public:
  creep_effym (void);
  ~creep_effym (void);
  void read (XFILE *in);
  void print (FILE *out);
  void matstiff (matrix &d,long ipp,long im,long ido);
  void nlstressesincr (long ipp, long im, long ido);
  void nlstresses (long ipp, long im, long ido);
  void matstiff (matrix &d,strastrestate ssst);
  void elmatstiff (matrix &d,strastrestate ssst);
  void matstiff (matrix &d,strastrestate ssst,  double ea);
  void matstiff_bar (matrix &d, double ea);
  void matstiff_plbeam (matrix &d, double ea);
  void matstiff_spacebeam (matrix &d, double ea);
  void matstiff_plstress (matrix &d, double ea);
  void matstiff_plstrain (matrix &d, double ea);
  void matstiff_axi (matrix &d, double ea);
  void matstiff_plate (matrix &d, double ea);
  void matstiff_spacestr (matrix &d, double ea);
  void matcompl (matrix &c,strastrestate ssst);
  void matcompl_bar (matrix &c, double ea);
  void matcompl_plbeam (matrix &c, double ea);
  void matcompl_plstress (matrix &c, double ea);
  void matcompl_plstrain (matrix &c, double ea);
  void matcompl_axi (matrix &c, double ea);
  void matcompl_spacestr (matrix &c, double ea);
  double actual_modulus();
  void updateval (long ipp, long im, long ido);
  void giveirrstrains (long ipp, long im, long ido, vector &epscr);

  double give_actual_ym ();
  double give_initial_ym ();
  double give_actual_nu ();

  double give_actual_ft (long ipp, long im, long ido);

  ///  initial value of Young's modulus
  double e0;
  ///  Poisson's number
  double nu;
  /// evolution function of Young modulus
  ym_evolfunc evf;

  //material parameters for b3 law:
 
  // = 1 measured Young's modulus E_28 (28 day) psi 
  long type_e;

  //according to Bazant's notation:
  // t  = t  ... age of concrete in days
  // tl = t' ... age at loading [days]
  double tl;

  //parameters
  double q1,q2,q3,q4,q5;

  //  e28 is 28 day Young's modulus
  double e28;
  //  fc' is 28 day average cylinder strength fc' [psi] 1000psi=6.895 MPa(f.e.6.454=44.5MPa)     = 6381.0
  double fc;
  //  w/c is water-cement ratio of the mix by weight                                             = 0.43
  double wc;
  //  s/c is sand-cement ratio of the mix by weight                                              = 3.4
  double sc;
  //  g/c is gravel-cement ratio of the mix by weight g/c=a/c-s/c                                = 1.98
  double gc;
  //  cs  cement content in m3  .. kg/m3
  double cs;
  //  coefficient of shape of structure
  double a1;
  //  coefficient for curing
  double a2;
  //  k_s shape factor slab=1.0, cylinder=1.15, square prism.=1.25, sphere=1.3, cube=1.55 
  double ks;
  //  k_d effective cross section thickness D=2*vs_s in inches (inch = 25.4mm)
  double kd;
  //age of concrete on the computation begining in seconds
  double tb_time;


  // material parameters for builtin double power law - spatne je nutne opravit??!!
  double qs;
  double psi;
  double n;
  double m;
  double alpha;
  double t0;
  /// user defined evolution function
  gfunct gf;
};

#endif
