#ifndef CONSOLAWF2MAT_H
#define CONSOLAWF2MAT_H

#include "genfile.h"
#include "aliast.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

class con_awf2mat
{
 public:
  con_awf2mat();    //constructor
  ~con_awf2mat();   //destructor

  void read(XFILE *in);
  void print(FILE *out);
  
  void matcond (matrix &d,long ri,long ci,long ipp);
  void matcap (double &c,long ri,long ci,long ipp);
  
  void matcond1d (matrix &d,long ri,long ci,long ipp);
  void matcond2d (matrix &d,long ri,long ci,long ipp);
  void matcond3d (matrix &d,long ri,long ci,long ipp);

  void rhs_volume (matrix &d,long ri,long ci,long ipp);
  void rhs_volume2 (double &cc,long ri,long ci,long ipp);
  void rhs1d1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d1 (matrix &d,long ri,long ci,long ipp);

  double get_sw(double pw, double pg, long ipp);
  double get_ds_dpc(double pw, double pg, long ipp);
  double get_krw(double pw, double pg, long ipp);
  double get_krg(double pw, double pg, long ipp);
  double get_kintr(double pw, double pg, long ipp);

  double get_kww(double pw, double pg, long ipp);
  double get_kwg(double pw, double pg, long ipp);
  double get_kgw(double pw, double pg, long ipp);
  double get_kgg(double pw, double pg, long ipp);

  double get_capww(double pw, double pg, long ipp);
  double get_capwg(double pw, double pg, long ipp);
  double get_capgw(double pw, double pg, long ipp);
  double get_capgg(double pw, double pg, long ipp);

  double get_fw1(double pw, double pg, long ipp);
  double get_fg1(double pw, double pg, long ipp);
  double get_fwu(double pw, double pg, long ipp);
  double get_fgu(double pw, double pg, long ipp);

  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp);
  double transmission_transcoeff(double trc,long ri,long ci,long nn,long bc,long ipp,int flag);
  double transmission_nodval(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);
  double transmission_flux(double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);

  double get_transmission_transcoeff_ww(double pw,double pg,long bc,long ipp);
  double get_transmission_transcoeff_ww(double pw,double pg,long bc,long ipp,int flag);
  double get_transmission_nodval_ww(double bv,double pw,double pg,long bc,long ipp);
  double get_transmission_flux_ww(double bv,double pw,double pg,long bc,long ipp);

  double get_transmission_transcoeff_gg(double pw,double pg,long bc,long ipp);
  double get_transmission_nodval_gg(double bv,double trr,double pw,double pg,long bc,long ipp);
  double get_transmission_flux_gg(double bv,double trr,double pw,double pg,long bc,long ipp);


  double get_alpha(double pw, double pg,long ipp);
  double get_ks(double pw,double pg,long ipp);
  double get_kt(double pw,double pg,long ipp);
  double get_rhos(double t);
  double get_cdiff(double pw,double pg,double t);
  double get_dg(double pc,double pg,double t,long ipp);
  double get_rhogw(double pw,double pg,double t);
  double get_pgw(double pw,double pg,double t);
  double get_dpgw_dpc(double pw,double pg,double t);
  double get_xi(double pw, double pg, long ipp);
  double get_pgws(double t);
  double get_rhow(double t);
  double get_kw(double pw,double pg,long ipp);
  double get_muw(double t);
  double get_mg(double pw,double pg,double t);
  double get_rhog(double pw,double pg,double t);
  double get_mug(double pw,double pg,double t);
  double get_muga(double t);
  double get_mugw(double t);
  double get_rhoga(double pw,double pg,double t);

  double get_othervalue(long compother,double pw, double pg, long ipp);
  void print_othervalue_name(FILE *out,long compother);
  void values_correction (vector &nv, long ipp);
  void gaspress_check(double pw,double &pg,long ipp);
  void waterpress_check(double &pw,double pg,long ipp);
  void updateval (long ipp);
  void initval(long ipp);
  
  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);
  double give_effective_pore_pressure(long ipp);
  double give_water_pressure(long ipp);
  double give_gas_pressure(long ipp);
  double give_pore_pressure(long ipp);
  double give_suction(long ipp);
  double give_saturation_degree(long ipp);
  double get_porosity(long ipp);
  double get_w(double pw,double pg,long ipp);

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
  
  airwaterflowtype model_type;
  int compress;  //compressible grains yes=1, no=0
  int vol_strain_effect,por_type,kintr_type,krw_type,krg_type,sr_type,xi_type,deff_type;

  double p_atm;//athmospheric pressure [Pa]
  double p_atm_kpa;//athmospheric pressure [kPa]

  double ma;//molar mass of dry air
  double mw;//molar mass of water
  double gasr;//universal gas constant

  double t0;//[K] reference temperature 
  double p0;//Pa reference athmospheric pressure

  // PHYSICAL PROPERTIES OF DRY AIR
  double muga0;//reference dynamic viscosity of dry air

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  // from D.Gawin, F.Pesavento (PRVAP.f90)
  double dv0,bv0;
  double c8,c9,c10,c11,c12,c13;
  double mugw0;

  // PHYSICAL PROPERTIES OF WATER
  double rhow0;//density at refernce temperature and pressure
  double tcr;//[K] critical temperature of water
  double muw0;//reference dynamic viscosity of water
  double kw0;//compresibility coefficient of water
  double mug0;//reference dynamic viscosity of moist air


  //PHYSICAL PROPERTIES OF SOIL
  double alpha0; //initial Boit's constant
  double ks0;    //inital bulk modulus of solid phase
  double kt0;    //inital bulk modulus of porous medium
  double phi0;   //inital porosity
  double kintr0; //intial intrinsic permeability
  double betas0; //inital cubic thermal dilatation coefficient
  double deff0;  //initial effective diffusivity
  double rhos0;  //inital volume density of soil skeleton
  double cdiff0; //initial water vapour diffusivity in air
  double tau0; //tortuosity factor 0.4 to 0.6

  double mefel_units; //basic units for pressures = Pa (Pascals)
  double pw_bc;//water vapour pressure for switching of free surface b.c.
  int rel_gas_press;  //relative gas pressure according to ambient air; yes=1=relative pg, no=0=absolute pg
  double krw0,sirr,ssat,lambda_krw,beta_krw;
  double bb1,phi01; //permeability parameters
  double krg0,s_crit,ag; //gas relative permeability parameters
  double kg0,kgn;
  double gamma,lambda0,s_entry;

};  

#endif
