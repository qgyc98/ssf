#ifndef HYPOPLAST_H
#define HYPOPLAST_H

#include "galias.h"
#include "gfunct.h"

struct XFILE;
struct matrix;
struct vector;
struct atsel;


/**
  The class defines hypoplastic model for unsaturated soils
  according to D. Masin.

  MASIN HYPO - Masin hypoplastic model with intergranular strains
  D. Masin (2005) A hypoplastic constitutive model for clays. IJNAMG, 29:311-336
 
  Implementation based on:
  Fellin, W. and Ostermann, A. (2002):
  Consistent tangent operators for constitutive rate equations.
  International Journal for Numerical and Analytical Methods in Geomechanics

  Ordering of eqother/other array:
  --------------------------------

  History of qunatities handeled by MEFEL:
  previous strain vector               /eqother[0] - eqother[ncompstr-1]/,
  previous stress vector               /eqother[ncompstr] - eqother[2*ncompstr-1]/,

  Internal variables required and handeled by hypoplasticity model:
  actual void ratio       e            /eqother[2*ncompstr]/,
  previous/actual suction s            /eqother[2*ncompstr+1]/,
  degree of saturation    S_r          /eqother[2*ncompstr+2]/,
  number of function evaluation nfev   /eqother[2*ncompstr+3]/,
  phi_mob in degrees                   /eqother[2*ncompstr+4]/,
  suggested substep size dstub         /eqother[2*ncompstr+5]/,

  Internal variables handeled by hypoplasticity model:
  required time step size coef. pnewdt /eqother[2*ncompstr+6]/,
*/
class hypoplast
{
 public:
  hypoplast();
  ~hypoplast();
  void read(XFILE *in);
  void initval (long ipp,long ido);
  void matstiff (matrix &d,long ipp,long ido);
  void givestressincr (long ipp, long im, long ido, vector &sig);
  void nlstresses (long ipp, long im, long ido);
  double dstep_red(long ipp, long im, long ido);
  void updateval (long ipp,long im,long ido);
  void give_reqnmq(long *anmq);
  double give_virgporosity(long ipp, long ido);
  double give_iniporosity(long ipp, long ido);
  void changeparam (atsel &atm,vector &val);
  void umatunsat(long ipp, vector &stress, vector &statev, matrix &ddsdde,
                 vector &stran, vector &dstran, double dtime, 
                 double &pnewdt, vector &unsatvar);
  long check_rkf(vector &y);
  void calc_khalili_stress(vector &sig_net, vector &sig_ef, double suction);
  long rkf23_update(long ipp, vector &y, long nasv, double &dtsub, double err_tol, 
                    long maxnint, double dtmin, vector &deps_np1,
                    double dsuction, long &nfev, double dtime);
  void rhs(vector &y, long nasv, vector &deps, vector &krk, long &nfev, double dsuction);
  void get_f_sig_q(vector &sig, vector &q, vector &deps, vector &f_sig, vector &f_q, double dsuction);
  void get_tan(vector &deps, vector &sig, vector &q, matrix &hh,
               matrix &ll, vector &nn, vector &hunsat);
  void inv_sig(vector &sig, double &pp, double &qq, double &cos3t, double &i1, double &i2, double &i3);
  double norm_res(vector &y_til, vector &y_hat, long nasv);
  void perturbate(vector &y_n, long nasv, vector &deps_np1, matrix &dd);
  void calc_elasti(vector &y, vector &deps_np1, matrix &ddtan, double youngel, double nuel);
  void solout(vector &stress, vector &asv, matrix &ddsdde, vector &y, matrix &dd);
  void calc_statev(vector &stress, vector &statev, long nasv);
  
  double phi;
  double phideg;
  double p_t;
  double lam_star;
  double kap_star;
  double n_star;
  double rr;
  double n;
  double l;
  double m;
  double s_e0;
  double e_0;
  long nprops;
  double *props;
  long nstatv;
  answertype pr_suc_fl;
  gfunct suc_fn;
};

#endif
