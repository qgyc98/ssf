#ifndef ORTODAMROT_H
#define ORTODAMROT_H

#include "iotools.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
  This class defines ortoropic damage material model with crack rotation.
  The different type of norms for the computing parameters of
  the damage function can be used.
  The order of internal variables is following:
   0    - actual value of param At (actualized for the quasibrittle evolution law only)
   1    - actual value of param Ac (actualized for the quasibrittle evolution law only)
   2 ..  7 - attained values of damage tensor for tension (order of values corresponds to stress vector)
   8 .. 13 - attained values of damage tensor for compression (order of values corresponds to stress vector)
  14 .. 19 - attained values of inelastic strain tensor caused by damage (actualized for fatigue only)
  20 .. 22 - crack widths in the principal directions
*/
class ortodamrot
{
 public:
  ortodamrot (void);
  ~ortodamrot (void);
  
  void   read (XFILE *in);
  double damdrvforce_vol(double nu, vector &peps);
  void   damdrvforce_dev(double nu, vector &peps, vector &pyt, vector &pyc);
  void   loadfunc(long ipp, double nu, vector &peps, vector &damt, vector &damc, double aat, double aac, vector &lft, vector &lfc);
  double brittle_damage(long ipp, double y, double e, double f, double uf, double omegao);
  double qbezier_damage(long ipp, double y, double e, double f, double uf, double ul, double omegao,long pq);
  void   princ_dam(long ipp, vector &peps, double aat, double aac, vector &pdamt, vector &pdamc);
  void   give_actual_param_a(long ipp, long ido, double &aat, double &aac);
  double give_actual_ft(long ipp);
  double give_actual_fc(long ipp);
  void   initvalues(long ipp, long ido);
  void   matstiff (matrix &d,long ipp,long ido);
  void   pelmatstiff (long ipp, matrix &d,double e, double nu);
  void   elmatstiff (matrix &d,long ipp);
  void   nlstresses (long ipp, long im, long ido);
  void   updateval (long ipp,long im,long ido);


  /// type of damage evolution function
  dam_evolfunc damevf;
  /// tensile strength
  double ft;
  /// initial crack opening u_ft or eps_ft - controls initial slope of softening branch for tension
  double uft;
  /// limit crack opening u_lt or eps_lt - limit value for zero stresses in tension
  double ult;
  /// compressive strength
  double fc;
  /// limit crack opening u_lc or eps_lc - limit value for zero stresses in compression
  double ulc;
  /// initial crack opening u_fc or eps_fc - controls initial slope of softening branch for compression
  double ufc;
  /// parameters for Newton tangent method for damage parameter computation for brittle evolution fucntion and correction of dissipated energy
  strretalg sra;
  /// correction of disipated energy switch
  corr_disip_en cde;
  /// material parameter At for tension
  double at;  
  /// material parameter Bt for tension
  double bt;  
  /// initial treshold for dimensionless damage driving force Y0 for tension
  double y0t;  
  /// material parameter Ac for compression
  double ac;  
  /// material parameter Bc for compression
  double bc;  
  /// initial treshold for dimensionless damage driving force Y0 for compression
  double y0c;  
  /// fracture energy of damage for tension
  double gft;
  /// fracture energy of damage for compression
  double gfc;  
  /// indicator for singularity in quadratic Bezier evolution function for tension (pure quadratic function will be used instead Bezier one)
  long pqt;
  /// indicator for singularity in quadratic Bezier evolution function for compression (pure quadratic function will be used instead Bezier one)
  long pqc;
  /// flag for counting with fatigue
  fatigue_flag fat;
  /// coefficient of inelastic strains computed from damage parameter for tension
  double betat;
  /// coefficient of inelastic strains computed from damage parameter for compression
  double betac;
};

#endif
