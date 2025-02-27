#ifndef HYPOPLAST_UNSAT_EXP_THERM_H
#define HYPOPLAST_UNSAT_EXP_THERM_H

#include "galias.h"
#include "gfunct.h"
#include "generalmod.h"
#include "strretalg.h"

struct XFILE;
struct matrix;
struct vector;
struct atsel;


/**
  The class defines hypoplastic model for unsaturated expansive soils
  with double porosity structure and thermal effects according to D. Masin.

  Ordering of eqother/other array:
  --------------------------------

  History of qunatities handeled by MEFEL:
  - previous strain vector               /eqother[0] - eqother[ncompstr-1]/,
  - previous stress vector               /eqother[ncompstr] - eqother[2*ncompstr-1]/,
  - initial increments of 
    normal stress components 
    due to pore pressure change          /eqother[2*ncompstr - eqother[2*ncompstr+2],

  Internal state variables required and handeled by hypoplasticity model (nstatev=10):
  - total void ratio                    e      /eqother[2*ncompstr+3]/,
  - previous/actual suction             s      /eqother[2*ncompstr+3+1]/,
  - total degree of saturation          S_r    /eqother[2*ncompstr+3+2]/,
  - temperatutre                        T      /eqother[2*ncompstr+3+3]/,
  - void ratio on microlevel            em     /eqother[2*ncompstr+3+4]/,
  - void ratio on macrolevel            eM     /eqother[2*ncompstr+3+5]/,
  - degree of saturation on microlevel SrM     /eqother[2*ncompstr+3+6]/,
  -                                    ascan   /eqother[2*ncompstr+3+7]/,
  -                                    re      /eqother[2*ncompstr+3+8]/,
  - swelling indicator                         /eqother[2*ncompstr+3+9]/, 

  Auxiliary quantities related to the substepping and testing purposes
  - actual value of the substep length for RKF     dtsub      /eqother[2*ncompstr+3+nstatev]/,
  - actual value of time coefficient for solver    pnewdt     /eqother[2*ncompstr+3+nstatev+1]/,
  - volumetric strain                              eps_v      /eqother[2*ncompstr+3+nstatev+2]/,
  - the third component of strain deviator vector  e_ax       /eqother[2*ncompstr+3+nstatev+3]/,
  - the third component of stress deviator vector  q_ax       /eqother[2*ncompstr+3+nstatev+4]/,
  - the number of model evaluations                neval      /eqother[2*ncompstr+3+nstatev+5]/,
  - the number of performed RKF steps              nstep      /eqother[2*ncompstr+3+nstatev+6]/,
  - the minimum step RKF length                    mindt      /eqother[2*ncompstr+3+nstatev+7]/,
  - der. of degree of stauration w. r. suction     dSr/ds     /eqother[2*ncompstr+3+nstatev+8]/,
  - der. of degree of stauration w. r. temperature dSr/dT     /eqother[2*ncompstr+3+nstatev+9]/,
  - rate of the volumetric strain                  deps_v     /eqother[2*ncompstr+3+nstatev+10]/
  - der. of degree of stauration w. r. vol. strain dSr/depsv  /eqother[2*ncompstr+3+nstatev+11]/,

  'Consistent' generalized material stiffness matrix obtianed by numerical integration RKF
  - Generalized material stiffness matrix components stored by rows, 
    dimension of the matrix is (ngstress, ngstrain)           /eqother[2*ncompstr+3+nstatev+12 - 
                                                               eqother[2*ncompstr+3+nstatev+12 + ngstress*ngstrain -1]/,

  The total number of eqother components ncompeqother= 2*ncompstr+3+nstatev+12+ngstress*ngstrain
  Values of nstatev, ngstress and ngstrain are defined in the generalmod.h (actual values are ngstress=7, ngstrain=8 and nstatev=10)
  ncompstr represents the number of stress/strain vector components and it depends on the problem solved

  Generalized stress vector is given according to D. Masin in genralmod.cpp as follows:
  sig_g = {sig_x, sig_y, sig_z, tau_xy, tau_xz, tau_yz, Sr}^T (the ABACUS ordering of stress components is being used)
  Generalized strain vector is given according to D. Masin in genralmod.cpp as follows:
  eps_g = {eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz, s, T}^T (the ABACUS ordering of strain components is being used)  
*/
class hypoplast_unsat_exp_therm
{
 public:
  hypoplast_unsat_exp_therm();
  ~hypoplast_unsat_exp_therm();
  void read(XFILE *in);
  void print(FILE *out);
  void initval (long ipp,long ido, bool rinit);
  void matstiff (matrix &d,long ipp,long ido);
  void nlstressesincr (long ipp, long im, long ido);
  void givestressincr (long lcid, long ipp, long im, long ido, long fi, vector &sig);
  void nlstresses (long ipp, long im, long ido);
  long adfwdeuler(vector &eps, vector &sig, vector &statev, vector &deps, double &dtsub, long &neval, long &nstep, double &mindt);
  long rkf23(long ipp, vector &eps, vector &sig, vector &statev, vector &deps, double &dtsub, long &neval, long &nstep, double &mindt);
  long rkf23bs(vector &eps, vector &sig, vector &statev, vector &deps, double &dtsub, long &neval, long &nstep, double &mindt);
  long rkf34(vector &eps, vector &sig, vector &statev, vector &deps, double &dtsub, long &neval, long &nstep, double &mindt);
  long rkf45(vector &eps, vector &sig, vector &statev, vector &deps, double &dtsub, long &neval, long &nstep, double &mindt);
  long rkf_redstep(double rc, double &h, double hmin);
  double dstep_red(long ipp, long im, long ido);
  void updateval (long ipp,long im,long ido);
  void give_reqnmq(long *anmq);
  long givencompeqother(long ipp);
  double give_sr(long ipp, long ido);
  double give_dsr_ds(long ipp, long ido);
  double give_dsr_dtemp(long ipp, long ido);
  double give_dsr_depsv(long ipp, long ido);
  double give_porosity(long ipp, long ido);
  double give_void_ratio(long ipp, long ido);
  double give_strain_vol_rate(long ipp, long ido);
  double give_bulk_modulus(long ipp, long ido);
  void changeparam (atsel &atm,vector &val);
  
  double phideg;      /// critical state friction angle
  double lam_star;    /** slope of the normal compression lines (isotropic/oedometric/critical state line) of a saturated soil 
                          in the ln(p^M/p_r) vs ln(1+e) plane */
  double kap_star;    /// slope of the macrostructural isotropic unloading line of a saturated soil in the ln(p^M/p_r) vs ln(1+e) plane (saturated state)
  double n_star;      /// position of the isotropic normal compression line, i.e. the value of ln(1+e) for p^M = p_r = 1kPa (saturated state)
  double nu;          /// 
  double n;           /// defines position of the normal compression line with respect to suction  N(s) = n_star + n*<ln(s/s_e)>
  double l;           /// defines slope of the normal compression line with respect to suction  lambda(s) = lam_star + l*<ln(s/s_e)>
  double nt;          /// defines position of the normal compression line with respect to temperature ???
  double lt;          /// defines slope of the normal compression line with respect to temperature ???
  double m;           /** controls the dependency of wetting-induced collapse on the overconsolidation ratio and the macropore occlusion 
                          by microporosity on relative void ratio */
  double alpha_s;     /// thermal dilatancy coefficient ???
  double kap_m;       /// dependency of microstructural swelling/shrinkage on the microstructural effective stress (saturated effective stress)
  double sm_star;     /// suction corresponding to em_star
  double em_star;     /// void ratio of micropores for fully saturated micropores (S^m_r) and dry macropores (S^M_r)
  double csh;         /// 
  double s_airentry0; /// the air entry value of suction for the reference macrostructural void ratio e_0^M (fully saturated microstructure)
  double em0;         /// reference macrostructural void ratio e_0
  double tr;
  double at;
  double bt;
  double lambdap0;
  double aer;
  double p_t;  

  /// flag for suction definition (yes=suction evolution is prescribed by suc_fn, no=suction is taken from TRFEL)
  answertype pr_suc_fl;
  /// function of prescribed suction evolution
  gfunct suc_fn;
  /// flag for temperature definition (yes=temperature is prescribed by tempr_fn, no=temperature is taken from TRFEL)
  answertype pr_tempr_fl;
  /// function of prescribed temperature evolution
  gfunct tempr_fn;
  /// integration method parameters, i.e. RKF parametrs
  strretalg sra;
  /// maximum value of 'zero' level of state variables for RKF time integration
  double sv_rkf_zero;

  /// number of parameters
  long nparams;
  /// number of state variables
  long nstatev;
  /// number of generalised strain tensor components
  long ngstrain;
  /// number of generalised stress tensor components
  long ngstress;
  /// number of RKF method parameters
  long nrkfparams; 
  /// number of RKF state variables used for the step length handling and performance measurements
  long nrkfstatev; 
  /// array of the model parameters for the passing to the generalmod proceudres
  double *params;
  /// array of Runge-Kutta method parameters for passing to the generalmod proceudres
  double *rkfparams;


  Hypoplasti_unsat_expansive_thermal huet;
  //  Hypoplasti_unsat_expansive_thermal_original huetor;
};

#endif
