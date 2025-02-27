#ifndef CONSOLWF2MATC_H
#define CONSOLWF2MATC_H

#include "genfile.h"
#include "lewis_schrefler.h"
#include "van_genuchten.h"

class con_wf2matc
{
 public:
  con_wf2matc();    //constructor
  ~con_wf2matc();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond1d_u (matrix &d,long ri,long ci,long ipp);
  void matcond2d_u (matrix &d,long ri,long ci,long ipp);
  void matcond2d_ax_u (matrix &d,long ri,long ci,long ipp);
  void matcond3d_u (matrix &d,long ri,long ci,long ipp);
  void matcap1d_u (matrix &d,long ri,long ci,long ipp);
  void matcap2d_u (matrix &d,long ri,long ci,long ipp);
  void matcap2d_ax_u (matrix &d,long ri,long ci,long ipp);
  void matcap3d_u (matrix &d,long ri,long ci,long ipp);

  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void matcond1d_l (matrix &d,long ri,long ci,long ipp);
  void matcond2d_l (matrix &d,long ri,long ci,long ipp);
  void matcond2d_ax_l (matrix &d,long ri,long ci,long ipp);
  void matcond3d_l (matrix &d,long ri,long ci,long ipp);
  void matcap1d_l (matrix &d,long ri,long ci,long ipp);
  void matcap2d_l (matrix &d,long ri,long ci,long ipp);
  void matcap2d_ax_l (matrix &d,long ri,long ci,long ipp);
  void matcap3d_l (matrix &d,long ri,long ci,long ipp);

  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs1d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_ax_1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d1 (matrix &d,long ri,long ci,long ipp);

  void rhs_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs1d2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_ax_2 (matrix &d,long ri,long ci,long ipp);
  void rhs3d2 (matrix &d,long ri,long ci,long ipp);

  double get_sw(double pw, double pg,long ipp);
  double get_porosity(long ipp);

  double get_kuw(double pw, double pg,long ipp);
  double get_kwu(double pw, double pg,long ipp);

  double get_kug(double pw, double pg,long ipp);
  double get_kgu(double pw, double pg,long ipp);

  double get_capuw(double pw, double pg,long ipp);
  double get_capwu(double pw, double pg,long ipp);

  double get_capug(double pw, double pg,long ipp);
  double get_capgu(double pw, double pg,long ipp);

  double get_fu1(double pw, double pg,long ipp);
  double get_fu2(double pw, double pg,long ipp);

  /// Lewis and Schrefler's retention curve:
  lewis_reten lewis_ret;
  /// van Genuchten's retention curve:
  van_genuchten_reten van_genuchten_ret;

 private:
  
  airwaterflowmechtype model_type;
  int compress;
  double alpha,ks,phi0,kw,rhow,muw0,kintr,mug0;
  double p_atm,rhos0;
  
  double mefel_units; //basic units for pressures = Pa (Pascals)
};  

#endif
