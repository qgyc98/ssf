#ifndef B3_H
#define B3_H

#include "alias.h"
#include "xfile.h"
struct matrix;
struct vector;

/**
   class contains B3 model of creep
*/

class b3mat
{
 public:
  b3mat (void);
  ~b3mat (void);
  void read (XFILE *in);
  void print (FILE *out);
  void compute_ages (long ipp, long ido);
  void give_ages (double &tb_age_dt,double &tb_age,double &tl_age,double &th_age,double &dt,double &maxtime,long &napptime,long ipp);
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

  double compute_beta_t(long ipp);
  void store_beta_t_eqother(double beta,long ipp,long ido);
  double give_beta_t_eqother(long ipp,long ido);

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
  
  void give_stressirrdstrains_eqother(vector &deps_ss,long ipp,long ido);
  void store_stressirrdstrains_eqother(vector &deps_ss,long ipp,long ido);

  double give_shrinkage_eqother(long ipp,long ido);
  void store_shrinkage_eqother(double eps_sh,long ipp,long ido);

  void store_hum_eqother(long ipp,long ido);
  void store_temp_eqother(long ipp,long ido);

  double give_tb_time(long ipp);
  double give_th_time(long ipp);
  long give_napproxtime(long ipp);
  void initvalues (long ipp, long ido);
  void updatevalues (long ipp, long ido);

  double b3_law (double t0, double tl, double t,long ipp,long ido);
  double give_q1();

  void give_deps_free (vector &deps_sh, double t0, double t_dt, double t, double dt, long ipp,long im,long ido);
  void give_deps_stressinduced (vector &deps_ss, double t0, double t_dt, double t, vector &sigma, long ipp,long im,long ido);
  double creep_give_actual_ft (long ipp,long im,long ido);

  void addirrstrains_eqother (vector &deps,long ipp, long ido);
  void giveirrstrains_eqother (long ipp, long ido, vector &epscr);

  void store_agstrains_eqother(long ipp,long ido,vector &eps_ag);
  void give_agstrains_eqother(long ipp,long ido,vector &eps_ag);

  double get_othervalue(long compother,long ipp);
  void print_othervalue_name(FILE *out,long compother);
  void give_reqnmq(long *anmq);

  //previous time
  double previoustime;
 
 private:
  // reading parameters type, 
  long type_b3;
  // = 0 constant h, 
  long type_h;
  // = 0 constant temperature, 
  long type_temp;

  // = 1 measured Young's modulus E_28 (28 day) psi 
  long type_e;

  // environmental humidity and temperature (if they are not changing)
  double hum_env,temp_env;

  //according to Bazant's notation:
  // t  = t  ... age of concrete in days
  // tl = t' ... age at loading [days]
  // t0 = t0 ... age when drying begins [days] (only t0 <= tl is considered)

  //flags: drying shrinkage, drying shrinkage, thermal shrinkage
  int flag_drshr,flag_shr,flag_temp;

  //flag for tensile strength 0 = constatnt, 1 = time dependent
  int ft_flag;
  
  //temperature effect on creep
  double kappa,beta_t,beta_h,beta_s;

  //parameters
  // q1 = instantaneous straind due to unit stress = elastic compliance for asymptotic modulus
  // q2 = ageing viscoelastic compliance
  // q3 = non-ageing viscoelastic compliance
  // q4 = flow compliance
  // q5 = compliance for creep at drying
  double q1,q2,q3,q4,q5;
  double q1_ini,q2_ini,q3_ini,q4_ini,q5_ini;

  //  coefficient of thermal dilatancy
  double alpha;
  //  e28 is 28 day Young's modulus
  double e28;
  //  fc' is 28 day average cylinder strength fc' [psi] 1000psi=6.895 MPa(f.e.6.454=44.5MPa)     = 6381.0
  double fc;
  //  ft' is 28 day average cylinder strength ft' [psi] 1000psi=6.895 MPa(f.e.6.454=44.5MPa)     = 267.36
  double ft;
  // ft_ratio = actual tensile strength/actual young modulus
  double ft_ratio;
  //  w/c is water-cement ratio of the mix by weight                                             = 0.43
  double wc;
  //  s/c is sand-cement ratio of the mix by weight                                              = 3.4
  double sc;
  //  g/c is gravel-cement ratio of the mix by weight g/c=a/c-s/c                                = 1.98
  double gc;
  //  cs  cement content in m3  .. kg/m3
  double cs;
  //  coefficient of shape of structure
  double a1;
  //  coefficient for curing
  double a2;
  //  k_s shape factor slab=1.0, cylinder=1.15, square prism.=1.25, sphere=1.3, cube=1.55 
  double ks;
  //  k_d effective cross section thickness D=2*vs_s in inches (inch = 25.4mm)
  double kd;

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

  //  asymptotic Young's modulus
  double e0;

  //actual time
  double actualtime;
  //time step
  double dtb;
  //ages
  double tb_age_dt,tb_age,tbl_age,tbh_age,maxtimeb;
  long napptimeb;

  double *retTime;
  double timeMax;
  
  //  total autogenous shrinkage for concrete
  double eps_ainf;
};

#endif
