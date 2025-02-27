#ifndef CONSOLHWF2MAT_H
#define CONSOLHWF2MAT_H

#include "aliast.h"
#include "genfile.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "gardner.h"
#include "potts.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

class con_hwf2mat
{
 public:
  con_hwf2mat();    //constructor
  ~con_hwf2mat();   //destructor

  void read(XFILE *in);
  void print(FILE *out);
  double get_betas(long ipp);
  double get_cps(double pw, double t, long ipp);
  double get_rhos(double t);
  double get_xi(double pw, double t, long ipp);
  double get_sw(double pw, double t, long ipp);
  double get_dsw_dpw(double pw, double t, long ipp);
  double get_dsw_dt(double pw, double t, long ipp);
  double get_kintr(double pw, double t, long ipp);
  double get_krw(double pw, double t, long ipp);
  double get_porosity(long ipp);

  double get_alpha(double pw, double t,long ipp);
  double get_ks(double pw,double t,long ipp);
  double get_betasw(double pw,double t,long ipp);
  double get_rhow(double t);
  double get_kw(double pw,double t,long ipp);

  double get_muw(double t);
  double get_cpw();
  double get_betaw(double t);
  double get_rhocp(double pw,double t,long ipp);
  double get_lambdaeff(double pw,double t,long ipp);

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
  void rhs_volume2 (double &cc,long ri,long ci,long ipp);
  void rhs1d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d1 (matrix &d,long ri,long ci,long ipp);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  void waterpress_check(double &pw,double t,long ipp);

  void values_correction (vector &nv,long ipp);

  double get_kww(double pw,double t,long ipp);
  double get_capww(double pw,double t,long ipp);
  double get_kwt(double pw,double t,long ipp);
  double get_capwt(double pw,double t,long ipp);

  double get_ktt1(double pw,double t,long ipp);
  double get_ktt2(double pw,double t,long ipp);
  double get_captt(double pw,double t,long ipp);
  double get_ktw(double pw,double t,long ipp);
  double get_captw(double pw,double t,long ipp);

  double get_ktt2a(double pw,double t,long ipp);
  double get_ktt2b(double pw,double t,long ipp);
  double get_ktt2c(double pw,double t,long ipp);

  double get_fw1(double pw,double t,long ipp);
  double get_ft1(double pw,double t,long ipp);
  double get_fwu(double pw, double t, long ipp);

  double get_transmission_transcoeff_ww(double pw,double t,long bc,long ipp);
  double get_transmission_transcoeff_ww(double pw,double t,long bc,long ipp,int flag);
  double get_transmission_nodval_ww(double bv,double pw,double t,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double pw,double t,long bc,long ipp);

  double get_transmission_transcoeff_tt(double pw,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double trr,double pw,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double trr,double pw,double t,long bc,long ipp);

  double get_othervalue(long compother,long ipp,double *r);
  void print_othervalue_name(FILE *out,long compother);
  void updateval (long ipp);
  void initval(long ipp);
  
  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);

  double give_temperature(long ipp);
  double give_water_pressure(long ipp);
  double give_effective_pore_pressure(long ipp);
  double give_pore_pressure(long ipp);
  double give_suction(long ipp);
  double give_saturation_degree(long ipp);
  double get_rh(double pw,double t);
  double get_w(double pw,double t,long ipp);

  /// marks required non-transport quantities
  void give_reqntq(long *antq);
  double get_drhogw_dpw(double pw,double t);
  double get_drhogw_dt(double pw,double t);
  double get_rhogw(double pw,double t);
  double get_dv(double pw,double t,long ipp);
  double get_pgws(double t);
  double get_dpgws_dt(double t);
  double get_pgw(double pw, double t);
  double get_dpgw_dpc(double pw, double t);
  double get_dpgw_dt(double pw,double t);
  double get_dhvap(double t);

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
  
  heatwaterflowtype model_type;
  int compress;
  int vol_strain_effect,por_type,kintr_type,krw_type,sr_type,xi_type,lambda_type,cps_type,betas_type,deff_type;

  double alpha0; //initial Boit's constant
  double ks0;    //inital bulk modulus of solid phase
  double kt0;    //inital bulk modulus of porous medium
  double phi0;   //inital porosity
  double kintr0; //intial intrinsic permeability
  double betas0; //inital cubic thermal dilatation coefficient
  double betas_dry;//cubic thermal dilatation coefficient of dry soil
  double betas_wet;//cubic thermal dilatation coefficient of saturated soil
  double rhos0;  //inital volume density of soil skeleton
  double cps0;   //inital specific heat of soil skeleton
  double cps_dry;//specific heat of dry soil
  double cps_wet;//specific heat of saturated soil
  double cps_lin; //linear part of specific heat of soil skeleton
  double lambda_eff0;//effective thermal conductivity of soil skeleton
  double lambda_dry;//effective thermal conductivity of dry soil
  double lambda_wet;//effective thermal conductivity of saturated soil
  double sr_dry; //saturation degree according to lambda_dry
  double sr_wet; //saturation degree according to lambda_wet
  double hvap0;
  
  double kw0,rhow0,muw0,cpw0;

  double pw_bc;
  double mefel_units; //basic units for pressures = Pa (Pascals)
  double sirr,ssat,lambda_krw,beta_krw;
  double bb1,phi01; //permeability parameters
  double tcr,mw,gasr;
  double c8,c9,c10,c11,c12,c13;
  double dv0,tau0,cpgw0;
  double gamma,lambda0,s_entry;

  //parameters for the arificial material:
  double kww0,kwt0,ktw0,ktt0;
  double capww0,capwt0,captw0,captt0;
};  

#endif
