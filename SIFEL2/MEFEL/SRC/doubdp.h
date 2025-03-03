#ifndef DOUBDP_H
#define DOUBDP_H

#include "alias.h"
#include "iotools.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   Class doubdp defines plasticity model for concrete
   consisting of two Drucker-Prager criterions.
   
   Josef Fiedler, 01/2015
*/

class doubdp
{
 public:
  doubdp (void);
  ~doubdp (void);
  void read (XFILE *in);

  void   params         (long ipp, vector &q);
  double yieldfunction  (vector &sig, vector &q);
  void   deryieldfsigma (vector &sig, vector &dfds, vector &q);
  void   matstiff       (matrix &d, long ipp,long ido);
  void   tangentstiff   (matrix &d,matrix &td,long ipp,long ido);
  void   nlstresses     (long ipp, long im, long ido);
  void   updateval      (long ipp, long im, long ido);
  double plasmodscalar  (vector &sig, vector &q);
  void   updateq        (long ipp, double dgamma, vector &q);

  /// strength in single compression
  double fc;
  /// strength in single tension
  double ft;
  /// strength in double compression
  double fb;
  /// alpha parameter for tensile DP
  double alphat;
  /// alpha parameter for compressive DP
  double alphac;
  /// tau parameter for tensile DP
  double taut;
  /// tau parameter for compressive DP
  double tauc;
  /// I1 value of intersection point
  double icr;
  /// sqrt(J2) value of the intersection point
  double jcr;
  /// I1 value of top point of cone
  double ivrch;
  /// sqrt(J2) value of the top point of cone
  double jvrch;
  /// index of stress return areas
  long v;
  /// concrete energy of cracking in tension
  double gtf;
  /// gammma for ultimate cracking
  double gamult;
  /// switcher to (not) consider softening, 0 -> off, 1 -> on
  long sft;

  ///  stress return algorithm
  strretalg sra;
  
  

};

#endif
