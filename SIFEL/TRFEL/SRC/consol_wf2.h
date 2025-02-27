#ifndef CONSOLWF2MAT_H
#define CONSOLWF2MAT_H

#include "genfile.h"
#include "aliast.h"
#include "lewis_schrefler.h"
#include "van_genuchten.h"

class con_wf2mat
{
 public:
  con_wf2mat();    //constructor
  ~con_wf2mat();   //destructor

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
  double get_cs(double pw, double pg, long ipp);
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

  double get_othervalue(long compother,double pw, double pg, long ipp);
  void print_othervalue_name(FILE *out,long compother);
  void values_correction (vector &nv, long ipp);
  void gaspress_check(double pw,double &pg,long ipp);
  void waterpress_check(double &pw,double pg,long ipp);
  void updateval (long ipp);
  void initval(long ipp);
  
  /// returns ordered dof names
  void give_dof_names(namevart *dofname, long ntm);

  double give_water_pressure(long ipp);
  double give_gas_pressure(long ipp);
  double give_pore_pressure(long ipp);
  double give_suction(long ipp);
  double get_suction(double pw, double pg, long ipp);
  double give_saturation_degree(long ipp);
  double get_porosity(long ipp);
  double get_w(double pw,double pg,long ipp);

  /// marks required non-transport quantities
  void give_reqntq(long *antq);

  /// Lewis and Schrefler's retention curve:
  lewis_reten lewis_ret;
  /// van Genuchten's retention curve:
  van_genuchten_reten van_genuchten_ret;
 
  ///  general function for retention curve given by set of data
  gfunct data;

 private:
  
  airwaterflowtype model_type;//it is the same as in consol_awf2.h
  int compress;  //compressible grains yes=1, no=0
  int vol_strain_effect,por_type,kintr_type,krw_type,krg_type,sr_type;

  double p_atm;//athmospheric pressure [Pa]
  double p_atm_kpa;//athmospheric pressure [kPa]

  // PHYSICAL PROPERTIES OF WATER
  double rhow0;//density at refernce temperature and pressure
  double tcr;//[K] critical temperature of water
  double muw0;//reference dynamic viscosity of water
  double kw;//compresibility coefficient of water
  double mug0;//reference dynamic viscosity of moist air


  //PHYSICAL PROPERTIES OF SOIL
  double alpha;  //initial Boit's constant
  double ks;    //inital bulk modulus of solid phase
  double kt0;    //inital bulk modulus of porous medium
  double phi0;   //inital porosity
  double kintr0; //intial intrinsic permeability
  double rhos0;  //inital volume density of soil skeleton

  double mefel_units; //basic units for pressures = Pa (Pascals)
  double pw_bc;//water vapour pressure for switching of free surface b.c.
  int rel_gas_press;  //relative gas pressure according to ambient air; yes=1=relative pg, no=0=absolute pg
  double sirr,ssat,lambda_krw;
  double bb1,phi01; //permeability parameters
  double t0;
};  

#endif
