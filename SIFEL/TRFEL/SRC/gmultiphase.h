#ifndef GMULTIPHASE_H
#define GMULTIPHASE_H

#include "genfile.h"

class gmultiph
{
 public:
  gmultiph();    //constructor
  ~gmultiph();   //destructor
  
  double give_pw (long nn);
  double give_pg (long nn);
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
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  void gaspress_check(double pc,double &pg,double t,long ipp);
  void cappress_check(double &pc,double pg,double t,long ipp);
  void cappress_stop(double &pc,double pg,double t,long ipp);

  void values_correction (vector &nv);

  double get_kww(double pw,double pg,double t,long ipp);
  double get_capww(double pw,double pg,double t,long ipp);
  double get_kwg(double pw,double pg,double t,long ipp);
  double get_capwg(double pw,double pg,double t,long ipp);
  double get_kwt(double pw,double pg,double t,long ipp);
  double get_capwt(double pw,double pg,double t,long ipp);

  double get_kgg(double pw,double pg,double t,long ipp);
  double get_capgg(double pw,double pg,double t,long ipp);
  double get_kgw(double pw,double pg,double t,long ipp);
  double get_capgw(double pw,double pg,double t,long ipp);
  double get_kgt(double pw,double pg,double t,long ipp);
  double get_capgt(double pw,double pg,double t,long ipp);

  double get_ktt1(double pw,double pg,double t,long ipp);
  double get_ktt2(double pw,double pg,double t,long ipp);
  double get_captt(double pw,double pg,double t,long ipp);
  double get_ktg(double pw,double pg,double t,long ipp);
  double get_captg(double pw,double pg,double t,long ipp);
  double get_ktw(double pw,double pg,double t,long ipp);
  double get_captw(double pw,double pg,double t,long ipp);

  double get_ktt2a(double pw,double pg,double t,long ipp);
  double get_ktt2b(double pw,double pg,double t,long ipp);
  double get_ktt2c(double pw,double pg,double t,long ipp);
  double get_ktt2d(double pw,double pg,double t,long ipp);

  double get_fc1(double pw,double pg,double t,long ipp);
  double get_fg(double pw,double pg,double t,long ipp);
  double get_ft1(double pw,double pg,double t,long ipp);

  double get_transmission_transcoeff_ww(double pw,double pg,double t,long bc,long ipp);
  double get_transmission_nodval_ww(double bv,double pw,double pg,double t,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double pw,double pg,double t,long bc,long ipp);

  double get_transmission_transcoeff_tt(double pw,double pg,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp);

  double heat_rate(double stime,double time);

  double get_othervalue(long compother,long ipp,double *r);
  void print_othervalue_name(FILE *out,long compother);

 private: 
  
  double ma;//molar mass of dry air
  double mw;//molar mass of water
  double gasr;//universal gas constant

  double scale_pw,scale_pg,scale_t;
};  

#endif
