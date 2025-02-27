#ifndef CAMCLAY_H
#define CAMCLAY_H

#include <stdio.h>

#include "xfile.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;

/**
   This class defines modified cam-clay material model
   the material is defined by three material constants :
   m - is frictional constant and it depends on the frictional angle phi
   lambda - slope of the normal consolidationa line
   kappa  - slope of the swelling line
   
   In addition to these material constants, several initial values 
   have to be specified at each integration point in the following order :
   v_kappa1   - initial specific volume at the reference pressure p1 after 
                unloading from initial consolidation pressure p_c0 (v_kappa0)
   p1         - reference pressure (p_1), it must be negative
   p_c0       - initial consolidation pressure, it must be negative
   iepsp_1    -
      .        \
      .         \ components of initial plastic strains
      .         /
   iepsp_ncomp -  
   
   Order of the eqother array components
   id                | description
   ------------------+------------------------------------------------
   <0;  ncompstr-1>  | plastic strains 
    ncompstr         | gamma - consistency parameter 
    ncompstr+1       | hardening parameter p_c - consolidation pressure
    ncompstr+2       | v_lambda1 - initial value of specific volume on the NCL at reference pressure p1
    ncompstr+3       | v_ini - initial specific volume
        .            | i1s - mean stress  
        .            | j2s - the second invariant of stress deviator
        .            | eps_v - total volumetric strain 
    ncompstr+7       | eps_vp - plastic volumetric strain  
    ncompstr+8       | degree of saturation Sr
    ncompstr+9       | depsv/dt - volumtric strain rate
    ncompstr+10      | actual porosity e
   

   18.1.2005
   TKo
*/

class camclay
{
 public:
  camclay (void);
  ~camclay (void);
  void read (XFILE *in);
  void print (FILE *out);

  void initval(long ipp, long ido);
  //double cohesion(vector &qtr);
  double yieldfunction (vector &sig, vector &q);
  void deryieldfsigma (vector &sig, vector &q, vector &dfds);
  void dderyieldfsigma (matrix &ddfds);
  void derpotsigma (vector &sig, vector &q, vector &dgds);
  void deryieldfq(vector &sig, vector &dq);
  void deryieldfdqdq(matrix &dfdqdq);
  void deryieldfdsdq(matrix &dfdsdqt);
  void dqdsigma(vector &sigt, vector &qt, vector &dqds);
  void dhardfdq(long ipp, long ido, double dgamma, vector &qt, vector &dqdq);
  void der_q_gamma(long ipp, long ido, vector &sig, vector &qtr, vector &dqdg);
  double plasmodscalar(long ipp, long ido, vector &sig, vector &qtr);
  void updateq(long ipp, long ido, vector &eps, vector &epsp, vector &sig, vector &q);
//  void updateq(vector &epsp, strastrestate ssst, vector &q);

  double give_actual_ym(long ipp, long im, long ido);
  double give_actual_nu(long ipp, long im, long ido);
  void matstiff (matrix &d, long im, long ipp, long ido);
  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  double give_consparam (long ipp, long ido);
  double give_saturation_degree(long ipp, long ido);
  double give_preconspress(long ipp, long ido);
  double give_virgporosity(long ipp, long ido);
  double give_iniporosity(long ipp, long ido);
  double give_porosity(long ipp, long ido);
  double give_strain_vol_rate(long ipp, long ido);
  void changeparam (atsel &atm, vector &val);

  /// frictional constant
  double m;

  /// slope of normal consolidation line
  double lambda;

  /// slope of swelling line
  double kappa;

  ///  stress return algorithm
  strretalg sra;
};

#endif
