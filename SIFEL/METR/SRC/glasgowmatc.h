#ifndef GLASGOWMATC_H
#define GLASGOWMATC_H

#include <stdio.h>
#include "genfile.h"

class glasgowmatc
{
 public:
  glasgowmatc (void);    //constructor
  ~glasgowmatc (void);   //destructor
  
  double give_rhov (long nn);
  double give_pg (long nn);
  double give_temp (long nn);

  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond1d_u (matrix &d,long ri,long ci,long ipp);
  void matcond2d_u (matrix &d,long ri,long ci,long ipp);
  void matcond3d_u (matrix &d,long ri,long ci,long ipp);
  void matcap1d_u (matrix &d,long ri,long ci,long ipp);
  void matcap2d_u (matrix &d,long ri,long ci,long ipp);
  void matcap3d_u (matrix &d,long ri,long ci,long ipp);

  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void matcond1d_l (matrix &d,long ri,long ci,long ipp);
  void matcond2d_l (matrix &d,long ri,long ci,long ipp);
  void matcond3d_l (matrix &d,long ri,long ci,long ipp);
  void matcap1d_l (matrix &d,long ri,long ci,long ipp);
  void matcap2d_l (matrix &d,long ri,long ci,long ipp);
  void matcap3d_l (matrix &d,long ri,long ci,long ipp);

  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs1d_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d_u1 (matrix &d,long ri,long ci,long ipp);
  
  void read (XFILE *in);

  double ktt (double t,double pg,double rhov);
  double ktp (double t,double pg,double rhov);
  double ktv (double t,double pg,double rhov);

  double kat (double t,double pg,double rhov);
  double kap (double t,double pg,double rhov);
  double kav (double t,double pg,double rhov);

  double kmt (double t,double pg,double rhov);
  double kmp (double t,double pg,double rhov);
  double kmv (double t,double pg,double rhov);

  double ctt (double t,double pg,double rhov);
  double ctp ();
  double ctv (double t,double pg,double rhov);

  double cat (double t,double pg,double rhov);
  double cap (double t,double pg,double rhov);
  double cav (double t,double pg,double rhov);

  double cmt (double t,double pg,double rhov);
  double cmp ();
  double cmv (double t,double pg,double rhov);

  //added by Tomas Krejci 28.1.2005
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  
  double get_transmission_transcoeff_mv(double t,double pg,double rhov,long bc,long ipp);
  double get_transmission_nodval_mv(double bv,double t,double pg,double rhov,long bc,long ipp);
  double get_transmission_flux_mv(double bv,double t,double pg,double rhov,long bc,long ipp);

  double get_transmission_transcoeff_tt(double t,double pg,double rhov,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double t,double pg,double rhov,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double t,double pg,double rhov,long bc,long ipp);
  
  double get_transmission_transcoeff_mt(double t,double pg,double rhov,long bc,long ipp);
  double get_transmission_nodval_mt(double bv,double t,double pg,double rhov,long bc,long ipp);
  double get_transmission_flux_mt(double bv,double t,double pg,double rhov,long bc,long ipp);
  
  double get_othervalue(long compother,long ipp, double *r);
  void print_othervalue_name(FILE *out,long compother);

  double f_h (double rhov, double t);

  //new
  double kuu (double t,double pg,double rhov);
  double kut (double t,double pg,double rhov);
  double kut2 (double t,double pg,double rhov);
  double kup (double t,double pg,double rhov);
  double kuv (double t,double pg,double rhov);
  double ktu ();
  double kpu ();
  double kvu ();

  double ctu ();
  double cpu ();
  double cvu ();
  double cut ();
  double cup ();
  double cuv ();

  double fuv1(double t,double pg,double rhov);
  double fup1(double t,double pg,double rhov);
  double fut1(double t,double pg,double rhov);
  double emod ();

  //Prescribed Values (stated in glasgowmat.ccp or read from .in file)
  double rhoc;          //  density of cement
  double rhos;          //  density of skeleton
  double por0;          //  initial porosity
  double k0;            //  initial permeability
  double rhol0;         //  initial water content
  double ra;            //  ideal gas constant for dry air
  double rv;            //  ideal gas constant for water vapour
  double sssp;          //  solid saturation point
  double t0;            //  initial temperature
  double rhov0;         //  initial vapour content
  double pginf;         //  external gas pressure
  double emmi;          //  emmisivity
  double stef;          //  Stefan-Boltzman
  double alph;          //  thermal diffusivity of dry air
  double hq;            //  convective heat transfer coefficient
  double crhoair;       //  heat capacity of dry air
  double tfirestart;    //  time at which fire begins

  long model;           //  flag for model type (1=original tenchev, 2=modified tenchev)

  private:

  double f_rhovinf();
  double f_tinf (double time);
  double f_pv (double rhov, double t);
  double f_psat(double t);
  double f_dpsatdt (double t);
  double f_rhow0 (double t);
  double f_rhow (double t);
  double f_drhowdt (double t);
  double f_por (double t);
  double f_dpordt (double t);
  double f_fracl0 (double t);
  double f_m (double t);
  double f_mi (double t);
  double f_fracl (double rhov, double t);
  double f_dfracldt (double rhov, double t);
  double f_dfracldrhov (double rhov, double t);
  double f_dmdt (double t);
  double f_dmidt (double t);
  double f_mcbwrel (double t);
  double f_dmcbwreldt (double t);
  double f_fracd (double t);
  double f_dfracddt(double t);
  double f_fracg(double rhov, double t);
  double f_s(double rhov, double t);
  double f_pl(long ipp, double rhov, double pg, double t);
  double f_pc(double rhov, double t);
  double f_dpcdt (double rhov, double t);  
  double f_dpcdrhov (double rhov, double t);
  double f_sln (double rhov, double t);
  double f_ppore (double rhov, double pg, double t, double pginf);
  double f_sb (double rhov, double t);
  double f_pa (double rhov, double pg, double t);
  double f_rhoa (double rhov, double pg, double t);
  double f_rhog (double rhov, double pg, double t);
  double f_davex (double pg, double tinf);
  double f_dav (double rhov, double pg, double t);
  double f_db(double rhov, double t);
  double f_kk (double t);
  double f_kg (double rhov, double t);
  double f_kl (double rhov, double t);
  double f_mul (double t);
  double f_muv (double t);
  double f_mua (double t);
  double f_mug (double rhov, double pg, double t);
  double f_cl (double t);
  double f_cv (double t);
  double f_ca (double  t);
  double f_cs (double t);
  double f_keff (double t);
  double f_fracs (double rhov, double t);
  double f_crho (double rhov, double pg, double t);
  double f_le (double t);
  double f_ld ();
  double f_hrad(double t, double tinf);
  double f_hqr(double t, double tinf);
  double f_beta (double pg, double tinf);
  //new
  double f_alpha(double t);

};

#endif

