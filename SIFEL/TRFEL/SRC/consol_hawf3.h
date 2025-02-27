#ifndef CONSOLHAWF3MAT_H
#define CONSOLHAWF3MAT_H

#include "aliast.h"
#include "genfile.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

class con_hawf3mat
{
 public:
  con_hawf3mat();    //constructor
  ~con_hawf3mat();   //destructor

  void read(XFILE *in);
  void print(FILE *out);

  double get_betas(long ipp);
  double get_cps(double pw, double pg, double t, long ipp);
  double get_rhos(double t);
  double get_cdiff(double pw,double pg,double t);
  double get_dg(double pc,double pg,double t,long ipp);
  double get_xi(double pw, double pg, double t, long ipp);
  double get_sw(double pw, double pg, double t, long ipp);
  double get_ds_dpc(double pw, double pg, double t, long ipp);
  double get_ds_dt(double pw, double pg, double t, long ipp);
  double get_porosity(long ipp);
  double get_kintr(double pw, double pg, double t, long ipp);
  double get_kintrg(double pw, double pg, double t, long ipp);
  double get_kintrw(double pw, double pg, double t, long ipp);
  double get_krw(double pw, double pg, double t, long ipp);
  double get_krg(double pw, double pg, double t, long ipp);

  double get_alpha(double pw, double pg, double t,long ipp);
  double get_ks(double pw,double pg,double t,long ipp);
  double get_kt(double pw,double pg,double t,long ipp);

  double get_rhogw(double pw,double pg,double t);
  double get_cpgw();
  double get_pgw(double pw,double pg,double t);
  double get_dpgw_dpc(double pw,double pg,double t);
  double get_dpgw_dt(double pw,double pg,double t);
  double get_pgws(double t);
  double get_dpgws_dt(double t);
  double get_rhow(double t);
  double get_kw(double pw,double pg,double t,long ipp);
  double get_betasg(double pw,double pg,double t,long ipp);
  double get_betasw(double pw,double pg,double t,long ipp);
  double get_betasgw(double pw,double pg,double t,long ipp);
  double get_drhow_dt(double t);
  double get_dhvap(double t);
  double get_muw(double t);
  double get_cpw();
  double get_lambdaw(double t);
  double get_betaw(double t);
  double get_mg(double pw,double pg,double t);
  double get_rhog(double pw,double pg,double t);
  double get_mug(double pw,double pg,double t);
  double get_muga(double t);
  double get_mugw(double t);
  double get_rhoga(double pw,double pg,double t);
  double get_cpga();
  double get_rhocp(double pw,double pg,double t,long ipp);
  double get_cpg(double pw,double pg,double t);
  double get_lambdaeff(double pw,double pg,double t,long ipp);

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

  double get_fwu(double pw, double pg, double t, long ipp);
  double get_fgu(double pw, double pg, double t, long ipp);
  double get_ftu(double pw, double pg, double t, long ipp);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  void gaspress_check(double pw,double &pg,double t,long ipp);
  void waterpress_check(double &pw,double pg,double t,long ipp);

  void values_correction (vector &nv, long ipp);

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

  double get_fw1(double pw,double pg,double t,long ipp);
  double get_fg(double pw,double pg,double t,long ipp);
  double get_ft1(double pw,double pg,double t,long ipp);

  double get_transmission_transcoeff_ww(double pw,double pg,double t,long bc,long ipp);
  double get_transmission_transcoeff_ww(double pw,double pg,double t,long bc,long ipp,int flag);
  double get_transmission_nodval_ww(double bv,double pw,double pg,double t,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double pw,double pg,double t,long bc,long ipp);

  double get_transmission_transcoeff_tt(double pw,double pg,double t,long bc,long ipp);
  double get_transmission_nodval_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp);
  double get_transmission_flux_tt(double bv,double trr,double pw,double pg,double t,long bc,long ipp);

  double get_othervalue(long compother,long ipp,double *r);
  void print_othervalue_name(FILE *out,long compother);
  void updateval (long ipp);
  void initval(long ipp);
  
  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);

  double give_temperature(long ipp);
  double give_water_pressure(long ipp);
  double give_gas_pressure(long ipp);
  double give_effective_pore_pressure(long ipp);
  double give_pore_pressure(long ipp);
  double give_suction(long ipp);
  double give_saturation_degree(long ipp);
  double get_w(double pw,double pg,double t,long ipp);

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
  
  heatairwaterflowtype model_type;
  int compress;
  int vol_strain_effect,wrc_vol_strain_effect,por_type,kintr_type,krw_type,krg_type,sr_type,xi_type,lambda_type,cps_type,betas_type,deff_type;
  int thermal_capacity_type;
  double ma;//molar mass of dry air
  double mw;//molar mass of water
  double gasr;//universal gas constant

  double t0;//[K] reference temperature 
  double p0;//Pa reference athmospheric pressure

  // PHYSICAL PROPERTIES OF DRY AIR
  double muga0,alphaa,betaa;

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  // from D.Gawin, F.Pesavento (PRVAP.f90)
  double dv0,bv0;
  double c8,c9,c10,c11,c12,c13;
  double mugw0,alphaw;

  // PHYSICAL PROPERTIES OF WATER
  // from Dariusz Gawin (WATPROP.f90)
  double rhow0;//density at refernce temperature and pressure
  double tcr;//[K] critical temperature of water
  double cwat;
  double betawat;
  double hvap0;//latent heat of vaporization at reference temperature
  double a0,a1,a2;
  double a3,a4,a5;
  double b0,b1,b2;
  double b3,b4,b5;
  double pr1,prif;
  double muw0;//reference dynamic viscosity of water
  double muw0const;//constant dynamic viscosity of water in 20 deg.C
  double conb;
  double conc;
  double cpw0;//water specific heat at const. pressure = 4181 J/kg.K
  double lambdaw;//water heat conductivity = 0.6 W/(m.K)
  double kw0;//compresibility coefficient of water

  //PHYSICAL PROPERTIES OF SOIL
  double alpha0; //initial Boit's constant
  double ks0;    //inital bulk modulus of solid phase
  double kt0;    //inital bulk modulus of porous medium
  double phi0;   //inital porosity
  double kintr0; //general intial intrinsic permeability (for both phases)
  double kintrw0; //intial intrinsic permeability for liquid water, if it is separated from gas
  double kintrg0; //intial intrinsic permeability for gas, if it is separated from water
  double betas0; //inital cubic thermal dilatation coefficient
  double betas_dry;//cubic thermal dilatation coefficient of dry soil
  double betas_wet;//cubic thermal dilatation coefficient of saturated soil
  double deff0;  //initial effective diffusivity
  double cdiff0; //initial water vapour diffusivity in air
  double rhos0;  //inital volume density of soil skeleton
  double cps0;   //inital specific heat of soil skeleton
  double cps_dry;//specific heat of dry soil
  double cps_wet;//specific heat of saturated soil
  double cps_0;  //constant part of specific heat of soil skeleton
  double cps_lin; //linear part of specific heat of soil skeleton
  double lambda_eff0;//effective thermal conductivity of soil skeleton
  double lambda_dry;//effective thermal conductivity of dry soil
  double lambda_wet;//effective thermal conductivity of saturated soil
  double sr_dry; //saturation degree according to lambda_dry
  double sr_wet; //saturation degree according to lambda_wet
  double rhocp_dry; //thermal capacity according to lambda_dry
  double rhocp_wet; //thermal capacity according to lambda_wet  
  double tau0; //tortuosity factor 0.4 to 0.6

  double mefel_units; //basic units for pressures = Pa (Pascals)
  double pw_bc;
  double scale_pw,scale_pg,scale_t;
  int rel_gas_press;  //relative gas pressure according to ambient air; yes=1=relative pg, no=0=absolute pg
  double krw0,sirr,ssat,lambda_krw,beta_krw;
  double bb1,phi01,ng; //intrinsic permeability parameters
  double krg0,s_crit,ag; //gas relative permeability parameters

  //parameters for the arificial material:
  double kww0,kwg0,kwt0,kgw0,kgg0,kgt0,ktw0,ktg0,ktt0;
  double capww0,capwg0,capwt0,capgw0,capgg0,capgt0,captw0,captg0,captt0;
  double gamma,lambda0,s_entry;
  double kggmin;
  double kg0,kgn;
};  

#endif
