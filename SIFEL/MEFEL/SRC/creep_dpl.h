#ifndef DPL_H
#define DPL_H

#include "xfile.h"
#include "alias.h"
#include <stdio.h>

struct vector;
struct matrix;

/** 
    Concrete creep double-power-law model according to Bazant.
*/

class dplmat
{
 public:
  dplmat (void);
  ~dplmat (void);
  void read (XFILE *in);
  void print (FILE *out);
  void compute_ages (long ipp,long ido);
  void give_ages (double &t_age_dt,double &t_age,double &tl_age,double &th_age,double &dt,double &maxtime,long &napptime,long ipp);
  long give_nceqother (long ipp);
  long give_ncother (long ipp);
  long give_nret_time (void);

  void give_rettimes (vector &rettimes,long n_ret_times,long ipp);

  void store_emu_eqother(vector &e_mu,long n_ret_times,long ipp,long ido);
  void give_emu_eqother(vector &e_mu,long n_ret_times,long ipp,long ido);

  void store_ym_eqother(double ym,long ipp,long ido);
  double give_ym_eqother(long ipp,long ido);

  void store_ym_old_eqother(double ym,long ipp,long ido);
  double give_ym_old_eqother(long ipp,long ido);

  void give_hidden_strains_eqother(matrix &gamma_mu,long ipp,long ido);
  void store_hidden_strains_eqother(matrix &gamma_mu,long ipp,long ido);

  void give_stresses_eqother(vector &sigma,long ipp,long ido);
  void store_stresses_eqother(vector &sigma,long ipp,long ido);

  void give_stresses_other(vector &sigma,long ipp,long ido);
  void store_stresses_other(vector &sigma,long ipp,long ido);
 
  void give_dstresses_eqother(vector &dsigma,long ipp,long ido);
  void store_dstresses_eqother(vector &dsigma,long ipp,long ido);
 
  void give_strains_eqother(vector &eps,long ipp,long ido);
  void store_strains_eqother(vector &eps,long ipp,long ido);
 
  void give_creepdstrains_eqother(vector &deps_cr,long ipp,long ido);
  void store_creepdstrains_eqother(vector &deps_cr,long ipp,long ido);
  
  void give_irrdstrains_eqother(vector &deps_sh,long ipp,long ido);
  void store_irrdstrains_eqother(vector &deps_sh,long ipp,long ido);

  long give_napproxtime(long ipp);

  double double_power_law (double tt, double t,long ipp,long ido);
  double give_q1();

  double creep_give_actual_ft (long ipp,long im,long ido);

  void addirrstrains_eqother (vector &deps,long ipp, long ido);
  void giveirrstrains_eqother (long ipp, long ido, vector &epscr);

  void store_agstrains_eqother(long ipp,long ido,vector &eps_ag);
  void give_agstrains_eqother(long ipp,long ido,vector &eps_ag);

  double get_othervalue(long compother,long ipp);
  void print_othervalue_name(FILE *out,long compother);

  //according to Bazant's notation:
  // t  = t  ... age of concrete in days
  // tt = t' ... age at loading [days]

  //parameters
  double q1;

  // = 1 measured Young's modulus E_28 (28 day) psi 
  long type_e;

  //  density .. kg/m3
  double ro;
  //  e28 is 28 day Young's modulus
  double e28;
  //  fc' is 28 day average cylinder strength fc' [psi] 1000psi=6.895 MPa(f.e.6.454=44.5MPa)     = 6381.0
  double fc;
  //  w/c is water-cement ratio of the mix by weight                                             = 0.43
  double wc;
  //  s/c is send-cement ratio of the mix by weight                                              = 3.4
  double sc;
  //  g/c is gravel-cement ratio of the mix by weight g/c=a/c-s/c                                = 1.98
  double gc;
  //  coefficient of shape of structure
  double a1;

  //time when structure is finished (concrete casting)
  double tb_time;
  // time when temperature and humidity start to change
  double th_time;
  //number of times of approximation
  long napproxtime;
  //number of retardation times
  long nRetTime;
  // type of reading of retardation times 0 = computing; 1 = reading
  long type_rt;
  //flag for tensile strength 0 = constatnt, 1 = time dependent
  int ft_flag;
  // ratio for tensile strength evolution
  double ft_ratio;
  //  ft' is 28 day average cylinder strength
  double ft;

  //  asymptotic Young's modulus
  double e0;

  //previous time
  double previoustime;
  //actual time
  double actualtime;
  //time step
  double dtb;
  //ages
  double tb_age_dt,tb_age,tbl_age,tbh_age,maxtimeb;
  long napptimeb;

  double *retTime;
  double *emu;
  double timeMax;
};

#endif
