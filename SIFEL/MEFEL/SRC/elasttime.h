#ifndef ELASTTIME_H
#define ELASTTIME_H

#include "alias.h"
#include "iotools.h"
#include "gfunct.h"
#include "itablefunct.h"
struct matrix;
struct vector;
struct atsel;

/**
   artificial material model for computation of elastic analysis by mtsolver
   mtsolver solves slow mechanical problems (inertial forces are neglected)
   this material model works only with one elastic material model
   
   remark: in time dependent problems the array stress contains increments
   of stress, not stress; stress components are located in the array eqother!
   
   structure of eqother array of a elastic time material:
   total stress components, previous total strains,
   
   JK, 28.9.2004					
*/
class elasttime
{
 public:
  elasttime (void);
  ~elasttime (void);
  void read (XFILE *in);
  void print (FILE *out);

  void matstiff (matrix &d, long ipp);
  void matstiff (matrix &d, long ipp, double time);

  void elmatstiff (matrix &d, long ipp);
  void elmatstiff (matrix &d, long ipp, double time);

  void matstiff_bar (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_plbeam (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_spacebeam (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_plstress (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_plstrain (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_axi (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_plate (matrix &d, long ipp, double time, stiffmatrix smt);
  void matstiff_spacestr (matrix &d, long ipp, double time, stiffmatrix smt);
  
  void matcompl (matrix &c, long ipp);
  void matcompl (matrix &c, long ipp, double time);
  void matcompl_bar (matrix &c, long ipp, double time, stiffmatrix smt);
  void matcompl_plbeam (matrix &c, long ipp, double time, stiffmatrix smt);
  void matcompl_plstress (matrix &c, long ipp, double time, stiffmatrix smt);
  void matcompl_plstrain (matrix &c, long ipp, double time, stiffmatrix smt);
  void matcompl_axi (matrix &c, long ipp, double time, stiffmatrix smt);
  void matcompl_spacestr (matrix &c, long ipp, double time, stiffmatrix smt);
  
  double actual_modulus(long ipp);
  double actual_modulus(long ipp, double time);
  double initial_modulus();

  void nlstresses (long ipp, long ido);
  void changeparam (atsel &atm,vector &val);
  void updateval(long ipp, long im, long ido);
  
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

  /** table for switching between time dependent stiffness and strain dependent stiffness
      according to actual time:
      itab.getval(t) = 1 => time dependent stiffness function (udstiff[0]) is used 
      itab.getval(t) = 2 => strain dependent stiffness function (udstiff[1]) is used
  */
  itablefunct itab;

  /** user defined stiffness function according to itab
    udstiff[0] = pointer to the time dependent stiffness function
    udstiff[1] = pointer to the strain dependent stiffness function
  */
  gfunct **udstiff;
};

#endif
