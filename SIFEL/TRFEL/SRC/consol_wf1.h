#ifndef CON_WF1MAT_H
#define CON_WF1MAT_H

#include "genfile.h"
#include "aliast.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "gardner.h"
#include "potts.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

//debug version:
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "errno.h"

/**
   This class defines model for water flow in soils.
*/

class con_wf1mat
{
 public:
  con_wf1mat();    //constructor
  ~con_wf1mat();   //destructor

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

  double get_xi(double pw, long ipp);
  double get_sw(double pw, long ipp);
  double get_dsw_dpw(double pw, long ipp);
  double get_krw(double pw, long ipp);
  double get_porosity(long ipp);
  double get_kintr(long ipp);
  double get_k(double pw,long ipp);
  double get_alpha();
  double get_ks(long ipp);
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

  double give_capillary_pressure(long ipp);
  double give_water_pressure(long ipp);
  double give_effective_pore_pressure(long ipp);
  double give_pore_pressure(long ipp);
  double give_suction(long ipp);
  double give_saturation_degree(long ipp);
  double get_w(double pw,long ipp);
  double masin_hypopl_psi(double suction,double dsuction, double e, long ipp);

  /// marks required non-transport quantities
  void give_reqntq(long *antq);

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
  
  waterflowtype model_type;
  int compress;
  int vol_strain_effect,wrc_vol_strain_effect,por_type,kintr_type,krw_type,sr_type,xi_type;
  double mefel_units; //basic units for pressures = Pa (Pascals)
  double alpha,phi0,ks0;
  double kw,rhow0,muw0,t0;
  double kintr0,rhos0;
  double krw0,sirr,ssat,lambda_krw,beta_krw;
  double b1,phi01;
  double pw_bc;
  double gamma,lambda0,s_entry;
  //debug version:
  double sairentry0,aer,a_scan,eM0,kappam,smstar,emstar,lambdap0,em;

};  

#endif
