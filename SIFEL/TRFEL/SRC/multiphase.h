#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "genfile.h"

class multiph
{
 public:
  multiph();    //constructor
  ~multiph();   //destructor

  double give_pg (long nn);
  double give_pc (long nn);
  double give_temp (long nn);

  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void matcond2 (matrix &d,long ri,long ci,long ipp);
  
  void matcond1d_2 (matrix &d,long ri,long ci,long ipp);
  void matcond2d_2 (matrix &d,long ri,long ci,long ipp);
  void matcond3d_2 (matrix &d,long ri,long ci,long ipp);

  void rhs_volume (matrix &d,long ri,long ci,long ipp);
  void rhs1d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d1 (matrix &d,long ri,long ci,long ipp);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  void gaspress_check(double pc,double &pg,double t,long ipp);
  void cappress_check(double &pc,double pg,double t,long ipp);
  void cappress_stop(double &pc,double pg,double t,long ipp);

  void values_correction (vector &nv);

  double get_kcc(double pc,double pg,double t,long ipp);
  double get_capcc(double pc,double pg,double t,long ipp);
  double get_kcg(double pc,double pg,double t,long ipp);
  double get_capcg(double pc,double pg,double t,long ipp);
  double get_kct(double pc,double pg,double t,long ipp);
  double get_capct(double pc,double pg,double t,long ipp);

  double get_kgg(double pc,double pg,double t,long ipp);
  double get_capgg(double pc,double pg,double t,long ipp);
  double get_kgc(double pc,double pg,double t,long ipp);
  double get_capgc(double pc,double pg,double t,long ipp);
  double get_kgt(double pc,double pg,double t,long ipp);
  double get_capgt(double pc,double pg,double t,long ipp);

  double get_ktt1(double pc,double pg,double t,long ipp);
  double get_ktt2(double pc,double pg,double t,long ipp);
  double get_captt(double pc,double pg,double t,long ipp);
  double get_ktg(double pc,double pg,double t,long ipp);
  double get_captg(double pc,double pg,double t,long ipp);
  double get_ktc(double pc,double pg,double t,long ipp);
  double get_captc(double pc,double pg,double t,long ipp);

  double get_ktt2a(double pc,double pg,double t,long ipp);
  double get_ktt2b(double pc,double pg,double t,long ipp);
  double get_ktt2c(double pc,double pg,double t,long ipp);
  double get_ktt2d(double pc,double pg,double t,long ipp);

  double get_fc1(double pc,double pg,double t,long ipp);
  double get_fg(double pc,double pg,double t,long ipp);
  double get_ft1(double pc,double pg,double t,long ipp);

  double get_transmission_transcoeff_cc(double pc,double pg,double t,long bc,long ipp);
  double get_transmission_transcoeff_cc(double pc,double pg,double t,long bc,long ipp,int flag);
  double get_transmission_nodval_cc(double bv,double pc,double pg,double t,long bc,long ipp);
  double get_transmission_flux_cc(double bv,double pc,double pg,double t,long bc,long ipp);

  double get_transmission_transcoeff_tt(double pc,double pg,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double trr,double pc,double pg,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double trr,double pc,double pg,double t,long bc,long ipp);

  double heat_rate(double stime,double time);
  double iso_fire (double time, double t0, double tfirestart);

  double get_othervalue(long compother,long ipp,double *r);
  void print_othervalue_name(FILE *out,long compother);

 private: 
  
  double ma;//molar mass of dry air
  double mw;//molar mass of water
  double gasr;//universal gas constant

  double scale_pc,scale_pg,scale_t;
};  

#endif
