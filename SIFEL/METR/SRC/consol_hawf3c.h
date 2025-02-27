#ifndef CONSOLHAWF3MATC_H
#define CONSOLHAWF3MATC_H

#include "genfile.h"
#include "baroghel_reten.h"
#include "bazant_reten.h"
#include "lewis_schrefler.h"
#include "gardner.h"
#include "potts.h"
#include "van_genuchten.h"
#include "masin_reten.h"
#include "febex_granit_reten.h"

class con_hawf3matc
{
 public:
  con_hawf3matc();    //constructor
  ~con_hawf3matc();   //destructor

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

  void gaspress_check(double pw,double &pg,double t,long ipp);
  void waterpress_check(double &pw,double pg,double t,long ipp);

  void values_correction (vector &nv);

  double get_kuw(double pw, double pg, double t,long ipp);
  double get_kwu(double pw, double pg, double t,long ipp);

  double get_kug(double pw, double pg, double t,long ipp);
  double get_kgu(double pw, double pg, double t,long ipp);

  double get_capuw(double pw, double pg, double t,long ipp);
  double get_capwu(double pw, double pg, double t,long ipp);

  double get_capug(double pw, double pg, double t,long ipp);
  double get_capgu(double pw, double pg, double t,long ipp);

  double get_kut(double pw, double pg, double t,long ipp);
  double get_ktu(double pw, double pg, double t,long ipp);

  double get_caput(double pw, double pg, double t,long ipp);
  double get_captu(double pw, double pg, double t,long ipp);


  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs1d_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_ax_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs3d_u1 (matrix &d,long ri,long ci,long ipp);

  double get_fuw1(double pw,double pg,double t,long ipp);
  double get_fug1(double pw,double pg,double t,long ipp);
  double get_fut1(double pw,double pg,double t,long ipp);

  void rhs_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs1d_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs2d_ax_u2 (matrix &d,long ri,long ci,long ipp);
  void rhs3d_u2 (matrix &d,long ri,long ci,long ipp);

  double get_fut2(double pw,double pg,double t,long ipp);

  double get_sw(double pw, double pg, double t,long ipp);
  double get_xi(double pw, double pg, double t, long ipp);
  double get_alpha(double pw, double pg, double t,long ipp);

  double get_rhogw(double pw,double pg,double t);
  double get_pgw(double pw,double pg,double t);
  double get_pgws(double t);
  double get_rhow(double t);
  double get_rhos(double t);
  double get_porosity(double pw, double pg,double t,long ipp);
  double get_rhoga(double pw,double pg,double t);
  double get_rhog(double pw,double pg,double t);
  double get_betas(long ipp);
  double get_e(long ipp);
  double get_nu(long ipp);
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
  

  heatairwaterflowmechtype model_type;
  int compress, vol_strain_effectc, wrc_vol_strain_effectc, pore_press_effectc, temper_effectc, sr_type, xi_type, betas_type;
  
  double ma;//molar mass of dry air
  double mw;//molar mass of water
  double gasr;//universal gas constant

  double t0;//[K] reference temperature 
  double p0;//Pa reference athmospheric pressure

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  // from D.Gawin, F.Pesavento (PRVAP.f90)
  double c8,c9,c10,c11,c12,c13;

  // PHYSICAL PROPERTIES OF WATER
  // from Dariusz Gawin (WATPROP.f90)
  double rhow0;//density at refernce temperature and pressure
  double tcr;//[K] critical temperature of water
  double hvap0;//latent heat of vaporization at reference temperature
  double a0,a1,a2;
  double a3,a4,a5;
  double b0,b1,b2;
  double b3,b4,b5;
  double pr1,prif;

  //PHYSICAL PROPERTIES OF SOIL
  double alpha0; //initial Boit's constant
  double phi0;   //inital porosity
  double betas0; //inital cubic thermal dilatation coefficient
  double betas_dry;//cubic thermal dilatation coefficient of dry soil
  double betas_wet;//cubic thermal dilatation coefficient of saturated soil
  double rhos0;  //inital volume density of soil skeleton
  double gamma,lambda0,s_entry;
  double emod,nu;
};  

#endif
