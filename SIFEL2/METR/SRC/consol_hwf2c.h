#ifndef CONSOLHWF2MATC_H
#define CONSOLHWF2MATC_H

#include "genfile.h"
#include "lewis_schrefler.h"

class con_hwf2matc
{
 public:
  con_hwf2matc();    //constructor
  ~con_hwf2matc();   //destructor

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

  double get_sw(double pw, double t, long ipp);
  double get_porosity(long ipp);
  double get_alpha(double pw, double t,long ipp);

  double get_kuw(double pw, double t,long ipp);
  double get_kwu(double pw, double t,long ipp);

  double get_kut(double pw, double t,long ipp);
  double get_ktu(double pw, double t,long ipp);

  double get_capuw(double pw, double t,long ipp);
  double get_capwu(double pw, double t,long ipp);

  double get_caput(double pw, double t,long ipp);
  double get_captu(double pw, double t,long ipp);

  double get_fuw1(double pw, double t,long ipp);
  double get_fut2(double pw, double t, long ipp);

  double get_betas(long ipp);
  double get_e(long ipp);
  double get_nu(long ipp);

  /// Lewis and Schrefler's retention curve:
  lewis_reten lewis_ret;

 private:
  
  heatwaterflowmechtype model_type;
  double alpha0,ks0,phi0,kintr0,betas0,rhos0,cps0,lambda0,rhow0,emod,nu;

  
  double mefel_units; //basic units for pressures = Pa (Pascals)
};  

#endif
