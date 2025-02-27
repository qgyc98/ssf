#ifndef ORTODAM2_H
#define ORTODAM2_H

#include "iotools.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
  This class defines scalar isotropic damage material model.
  The different type of norms for the computing parameters of
  the damage function can be used.
  The order of internal variables is following :
  0..2 - reached values of principal damage parameters for tension
  3..5 - reached values of principal damage parameters for compression
  6..11 - damage tensor for tension
  12..15 - cracks openings computed from tensile damage tensor
*/
class ortodam2
{
 public:
  ortodam2 (void);
  ~ortodam2 (void);
  void   read (XFILE *in);
  double dam_eq_strain(long i,vector &peps, vector &psig, matrix &de);
  double brittle_damage(long ipp, double y, double e, double f, double uf, double omegao);
  double qbezier_damage(long ipp, double y, double e, double f, double uf, double ul, double omegao,long pq);
  void   princ_dam(long ipp, vector &peps, matrix &pde, vector &pdamt, vector &pdamc);
  double give_actual_ft(long ipp);
  double give_actual_fc(long ipp);
  void   matstiff (matrix &d, long ipp, long ido);
  void   elmatstiff (matrix &d, long ipp);
  void   nlstresses (long ipp, long im, long ido);
  void   updateval (long ipp, long im, long ido);


  /// type of damage evolution function
  dam_evolfunc damevf;
  /// type of equivalent strain norm
  paramf_type  dameqstr;
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
