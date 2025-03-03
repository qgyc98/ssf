#ifndef ANISODAMROT_H
#define ANISODAMROT_H

#include "xfile.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;

/**
  This class defines scalar isotropic damage material model.
  The different type of norms for the computing parameters of
  the damage function can be used.
  The order of internal variables is following :
  0 - actual value of param a
  1 - actual value of param At
  2 - actual value of param Ac
  3 - reached value of volumetric damage parameter
  4..9 - reached values of deviatoric damage parameters for tension
  10..15 - reached values of deviatoric damage parameters for compression
  16..18 - reached principal values of damage tensor for tension
  19..21 - reached principal values of damage tensor for compression
*/
class anisodamrot
{
 public:
  anisodamrot (void);
  ~anisodamrot (void);
  void   read (XFILE *in);
  double damdrvforce_vol(double nu, vector &peps);
  void   damdrvforce_dev(double nu, vector &peps, vector &pyt, vector &pyc);
  double loadfuncvol(long ipp, double nu, vector &peps, double d, double aa);
  void   loadfuncdev(long ipp, double nu, vector &peps, vector &damt, vector &damc, double aat, double aac, vector &lft, vector &lfc);
  double dam_vol(long ipp, double y, double aa, double dvo);
  void   pdam_dev(long ipp, vector &pyc, vector &pyt, double aat, double aac, vector &pdamt, vector &pdamc);
  double brittle_damage(long ipp, double y, double e, double f, double uf, double omegao);
  void   give_actual_param_a(long ipp, long ido, double &aa, double &aat, double &aac);
  double give_actual_ft(long ipp);
  double give_actual_fc(long ipp);
  void   initvalues(long ipp, long ido);
  void   matstiff (matrix &d, long ipp, long ido);
  void   elmatstiff (matrix &d, long ipp, long ido);
  void   nlstresses (long ipp, long im, long ido);
  void   updateval (long ipp, long im, long ido);


  /// type of damage evolution function
  dam_evolfunc damevf;
  /// tensile strength
  double ft;
  /// initial crack opening u_ft or eps_ft - controls initial slope of softening branch for tension
  double uft;
  /// compressive strength
  double fc;
  /// initial crack opening u_fc or eps_fc - controls initial slope of softening branch for compression
  double ufc;
  /// parameters for Newton tangent method for damage parameter computation for brittle evolution fucntion and correction of dissipated energy
  strretalg sra;
  /// correction of disipated energy switch
  corr_disip_en cde;
  /// material parameter A for volumetric damage
  double a;  
  /// material parameter B for volumetric damage
  double b;  
  /// initial treshold for damage driving force Y0 for volumetric damage
  double y0;  
  /// material parameter At for tension
  double at;  
  /// material parameter Bt for tension
  double bt;  
  /// initial treshold for damage driving force Y0 for tension
  double y0t;  
  /// material parameter Ac for compression
  double ac;  
  /// material parameter Bc for compression
  double bc;  
  /// initial treshold for damage driving force Y0 for compression
  double y0c;  
  /// fracture energy of volumetric damage 
  double gf;
  /// fracture energy of deviatoric damage for tension
  double gft;
  /// fracture energy of deviatoric damage for compression
  double gfc;
  
};

#endif
