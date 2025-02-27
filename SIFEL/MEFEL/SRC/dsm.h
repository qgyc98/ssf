#ifndef DSMMAT_H
#define DSMMAT_H

#include "galias.h"
#include "gfunct.h"
#include "xfile.h"
#include "alias.h"
#include "strretalg.h"
struct matrix;
struct vector;
struct atsel;
class nonlinman;


/**
   This class defines double structure plasticity model for expansive clays based on cam-clay material model coupled with moisture transfer 
   and on a double structure concept according to:
   
   Sanchez, M., Gens, A., Do Nascimento Guimares, L., and Olivella, S., A double structure generalized plasticity model for expansive materials,
   International Journal for Numerical and Analytical Methods in Geomechanics, vol. 29, no. 8, pp. 751-787, 2005. doi:10.1002/nag.434.

   - this concept is not working for the curent MEFEL version, which assumes total values (not their rates) of strains and stresses

   - it must be reformulated and derived 

   - not completed!! - not working!!


   Plasticity model - modified cam-clay material model:

   The material is defined by 
   three usual camclay material constants :
   m - is frictional constant and it depends on the frictional angle phi
   lambda - slope of the normal consolidationa line
   kappa  - slope of the swelling line

   Additionally, there are three parameters which controls
   evolution of hardening with respect to degree of saturation and suction
   c1_tilde, c2 - controls the ratio of specific volumes v/v_sat (partially saturated/fully saturated)
   ks - controls apparent adhesion (intersection of yield surface with p axis)
   
   In addition to these material constants, several initial values 
   have to be specified at each integration point in the following order :
   v_kappa1   - initial specific volume at the reference pressure p1 after 
                unloading from initial consolidation pressure p_c0 (v_kappa0)
   p1         - reference pressure (p_1), it must be negative
   p_c0       - initial consolidation pressure, it must be negative
   sig0_1    -
      .        \
      .         \ components of eigenstresses
      .         /
   sig0_ncomp -  
   
   Order of the eqother array components
   id                | description
   ------------------+------------------------------------------------
   <0;  ncompstr-1>  | plastic strains 
    ncompstr         | gamma - consistency parameter 
    ncompstr+1       | hardening parameter p_c - consolidation pressure
    ncompstr+2       | hardening parameter p_s - apparent adhesion
    ncompstr+3       | v_lambda1 - initial value of specific volume on the NCL at reference pressure p1
        .            | v_ini - initial specific volume
        .            | i1s - mean stress  
        .            | j2s - the second invariant of stress deviator
        .            | epsv - total volume strain 
    ncompstr+8       | epsvp - plastic volume strain  
    ncompstr+9       | degree of saturation
    ncompstr+10      | depsv/dt - volumtric strain rate
    ncompstr+11      | actual porosity e
   

   10.2012
   TKo
*/
class dsmmat
{
 public:
  dsmmat (void);
  ~dsmmat (void);
  void read (XFILE *in);
  void initval(long ipp, long im, long ido);
  double cohesion(vector &qtr);
  double yieldfunction (vector &sig, vector &q);
  void deryieldfsigma (vector &sig, vector &q, vector &dfds);
  void dderyieldfsigma (matrix &ddfds);
  void derpotsigma (vector &sig, vector &q, vector &dgds);
  void deryieldfq(vector &sig, vector &q, vector &dq);
  void deryieldfdqdq(matrix &ddfdq);
  void deryieldfdsdq(matrix &dfdsdqt);
  void der_q_gamma(long ipp, long ido, vector &sig, vector &qtr, vector &epsp, vector &dqdg);
  void dqpardgammadsigma(long ipp, long ido, vector &sig, vector &qtr, vector &epsp, matrix &dqdgds);
  void dqpardgammadqpar(long ipp, long ido, vector &sig, vector &qtr, vector &epsp, matrix &dqdgdq);
  double plasmodscalar(long ipp, long ido, vector &sig, vector &qtr, vector &epsp);
  void updateq(long ipp, long ido, vector &eps, vector &epsp, vector &sig, vector &q);
//  void updateq(vector &epsp, strastrestate ssst, vector &q);
  void matstiff (matrix &d, long ipp,long ido);
  void depmatstiff (matrix &d, long ipp, long ido);
  void nlstressesincr (long ipp,long im,long ido);
  void givestressincr (long lcid, long ipp, long im, long ido, long fi, vector &sig);
  void nlstresses (long ipp,long im,long ido);
  long actualize_stiffmat(const nonlinman &nlman, long istep, long jstep);
  void updateval (long ipp, long im, long ido);
  double comp_actual_ym(long ipp, long im, long ido);
  double give_actual_ym(long ipp, long im, long ido);
  double give_actual_nu(long ipp, long im, long ido);
  double dstep_red(long ipp, long im, long ido);
  void giveirrstrains (long ipp, long ido, vector &epsp);
  double give_consparam (long ipp, long ido);
  double give_preconspress(long ipp, long ido);
  double give_saturation_degree(long ipp, long ido);
  double give_virgporosity(long ipp, long ido);
  double give_iniporosity(long ipp, long ido);
  double give_porosity(long ipp, long ido);
  double give_strain_vol_rate(long ipp, long ido);
  void give_reqnmq(long *anmq);
  void changeparam (atsel &atm, vector &val);

  /// frictional constant
  double m;

  /// slope of normal consolidation line
  double lambda;

  /// slope of swelling line
  double kappa;

  /// first parameter for control of the v/v_sat ratio evolution
  double c1_tilde;

  /// second parameter for control of the v/v_sat ratio evolution
  double c2;

  /// coefficent of apparent adhesion due to suction pressure
  double ks;

  /// soil swelling index reference value
  double kappa_s0;

  /// parameters for the strain rate caculation caused by the suction change
  double aks1,aks2;
  /// parameters for the thermal strain rate caculation
  double a0,a2;

  /// flag for suction definition (yes=suction evolution is prescribed by suc_fn, no=suction is taken from TRFEL)
    //answertype pr_suc_f;
  /// function of prescribed suction evolution
    //gfunct suc_fn;
  /// flag for temperature definition (yes=temperature is prescribed by tempr_fn, no=temperature is taken from TRFEL)
    //answertype pr_tempr_f;
  /// function of prescribed temperature evolution
    //gfunct tempr_fn;

  ///  stress return algorithm
  strretalg sra;
};

#endif
