#ifndef CON_AWF1MATC_H
#define CON_AWF1MATC_H

#include "genfile.h"
#include "aliast.h"
#include "aliasc.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "gardner.h"
#include "potts.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

/**
   This class defines model for water flow in soils (deformable porous medium)
*/

class con_awf1matc
{
 public:
  con_awf1matc();    //constructor
  ~con_awf1matc();   //destructor

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

  double get_sw(double pw,long ipp);
  double get_xi(double pw, long ipp);
  double get_phi(double pw,long ipp);
  double get_rhos(double pw,long ipp);

  double get_kuw(double pw,long ipp);
  double get_kwu(double pw,long ipp);

  double get_capuw(double pw,long ipp);
  double get_capwu(double pw,long ipp);

  double get_fu1(double pw,long ipp);

  double get_rhogw(double pw);
  double get_pgw(double pw);
  double get_pgws(double t);

  double nu;

  /// Baroghel retention curve:
  baroghel_reten baroghel_ret;
  /// Bazant retention curve:
  bazant_reten bazant_ret;
  /// Lewis and Schrefler's retention curve:
  lewis_reten lewis_ret;
  /// Gardner's retention curve:
  gardner_reten gardner_ret;
  /// Potts' retention curve:
  potts_reten potts_ret;
  /// van Genuchten's retention curve:
  van_genuchten_reten van_genuchten_ret;
  /// Masin's retention curve for bentonite
  masin_reten masin_ret;
  ///  general function for retention curve given by set of data
  gfunct data;
  /// FEBEX retention curve
  febex_granit_reten febex_granit_ret; 
  
 private:
  
  waterflowmechtype model_type;
  int vol_strain_effectc,pore_press_effectc,wrc_vol_strain_effectc,sr_type,xi_type;
  double alpha,phi0,rhos0,rhow0,t0;
  double gamma,lambda0,s_entry;
  double c8,c9,c10,c11,c12,c13;
  double mw,gasr;
};  

#endif
