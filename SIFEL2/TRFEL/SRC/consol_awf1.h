#ifndef CON_AWF1MAT_H
#define CON_AWF1MAT_H

#include "genfile.h"
#include "aliast.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

/**
   This class defines model for water flow in soils.
*/

class con_awf1mat
{
 public:
  con_awf1mat();    //constructor
  ~con_awf1mat();   //destructor

  void read(XFILE *in);
  void print(FILE *out);
 
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &cc,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void rhs_volume (matrix &d,long ri,long ci,long ipp);
  void rhs_volume2 (double &cc,long ri,long ci,long ipp);
  void rhs1d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d1 (matrix &d,long ri,long ci,long ipp);

  double get_sw(double pw, long ipp);
  double get_dsw_dpw(double pw, long ipp);
  double get_alpha();
  double get_kw();
  double get_muw();

  double get_kww(double pw,long ipp);
  double get_capww(double pw,long ipp);
  double get_fw1(double pw,long ipp);
  double get_fwu(double pw, long ipp);
    
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_transcoeff_ww(double pw,long bc,long ipp);
  double get_transmission_transcoeff_ww(double pw,long bc,long ipp,int flag);
  double get_transmission_nodval_ww(double bv,double pw,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double pw,long bc,long ipp);

  double get_othervalue(long compother,double pw, long ipp);
  void print_othervalue_name(FILE *out,long compother);
  void values_correction (vector &nv, long ipp);
  void waterpress_check(double &pw,long ipp);
  void updateval (long ipp);
  void initval(long ipp);

  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);

  double give_water_pressure(long ipp);
  double give_pore_pressure(long ipp);
  double give_suction(long ipp);
  double give_saturation_degree(long ipp);
  double get_w(double pw,long ipp);
  double get_xi(double pw, long ipp);
  double get_kintr(double pw, long ipp);
  double get_krw(double pw, long ipp);
  double get_porosity(long ipp);
  double get_dv(double pw,long ipp);
  double get_rhogw(double pw);
  double get_pgw(double pw);
  double get_drhogw_dpw(double pw);
  double get_dpgw_dpc(double pw);
  double get_pgws(double t);
  double give_effective_pore_pressure(long ipp);

  /// marks required non-transport quantities
  void give_reqntq(long *antq);

  /// Baroghel retention curve:
  baroghel_reten baroghel_ret;
  /// Bazant retention curve:
  bazant_reten bazant_ret;
  /// Lewis and Schrefler's retention curve:
  lewis_reten lewis_ret;
  /// van Genuchten's retention curve:
  van_genuchten_reten van_genuchten_ret;
  /// Masin's retention curve for bentonite
  masin_reten masin_ret;   
  ///  general function for retention curve given by set of data
  gfunct data;
  /// FEBEX retention curve
  febex_granit_reten febex_granit_ret; 

 private:
  
  waterflowtype model_type;
  int compress;
  int vol_strain_effect,krw_type,por_type,kintr_type,sr_type,xi_type,deff_type;
  double alpha0,phi0,ks;
  double mw,gasr,kw,rhow,muw0,t0;
  double kintr0,rhos0,krw0;
  double gamaw,k0,e0,lam0,lam,m,sigma_m_bar;
  double kappa,sigma_m_eff_bar,sigma_m_tot_bar,sigma_m_eff_init;
  double pw_bc;
  double mefel_units; //basic units for pressures = Pa (Pascals)
  double b1,bb1,phi01;
  double sirr,ssat,lambda_krw,lambda0,gamma,s_entry,beta_krw;
  double dv0,c8,c9,c10,c11,c12,c13,tau0;
};  

#endif
