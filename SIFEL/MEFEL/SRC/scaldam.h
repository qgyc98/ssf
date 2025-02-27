#ifndef SCALDAM_H
#define SCALDAM_H

#include "iotools.h"
#include "alias.h"
#include "gfunct.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
  This class defines scalar isotropic damage material model.
  The different type of norms for the computing parameters of
  the damage function can be used.
  The order of internal variables is following :
  0 - actual value of equivalent strain
  1 - actual value of damage parameter
  2 - actual value of crack opening
*/
class scaldam
{
 public:
  scaldam (void);
  ~scaldam (void);
  void read (XFILE *in);
  void print (FILE *out);
  void damfuncpar(long ipp, vector &eps, vector &kappa);
  void derdamfuncpar(long ipp, vector &eps, vector &dkappa);
  double damfunction(long ipp, vector &kappa, vector &omegao);
  void matstiff (matrix &d,long ipp,long ido);
  void elmatstiff (matrix &d,long ipp);
  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp,long im,long ido);
  double give_actual_ft (long ipp, long im, long ido);
  double epsefunction (long ipp);
  double givedamage (long ipp, long ido);
  double give_proczonelength (long ipp, long ido);
  double give_crackwidth (long ipp, long ido);
  void changeparam (atsel &atm,vector &val);


  /// tensile strength for room temperature
  double ft;
  /// determines the softening -> corresponds to crack opening (not strain) when tension stress vanishes
  double uf;
  /// parameter for von Mises norm - ratio of compression and tensile strength
  double k;
  /// type of function for evaluation damage funtion parameter
  paramf_type ftype;
  /// correction of disipated energy switch
  corr_disip_en cde;
  ///  stress return algorithm
  strretalg sra;
  /// switch for effects of temperature on tensile strength 
  long  cftt; 
  /// function describing dependency of tensile strength on temperature
  gfunct ft_temp;
  /// the minimum argument of exponential function
  double min_exp_arg;  
};

#endif
