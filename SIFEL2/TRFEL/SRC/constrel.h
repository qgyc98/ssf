#ifndef CONSTREL_H
#define CONSTREL_H

#include "genfile.h"

class state_eq
{
 public:
  state_eq();//constructor
  ~state_eq();   //destructor
  
  //CAPILLARY PRESSURE
  double get_pc(double pg,double pw);
  //DEGREE OF SATURATION AND DERIVATIVES
  double get_s(double pc,double t,long ipp);
  double get_ds_dpc(double pc,double t,long ipp);
  double get_ds_dt(double pc,double t,long ipp);
  double get_ssp(long ipp);
  double get_ddbw(double pc,double pg,double t,long ipp);
  //GAS - DRY AIR+WATER VAPOR; DENSITY AND DERIVATIVES
  double get_rh (double pc,double t);//RELATIVE HUMIDITY
  double get_pcrh(double rh,double t);//inverse (RELATIVE HUMIDITY)
  double get_w(double pc,double pg,double t,long ipp);
  double get_drh_dpc (double pc,double t);
  double get_drh_dt (double pc,double t);
  double get_pg(double pga,double pgw,double t);
  double get_rhog(double pc,double pg,double t);
  double get_drhog_dpc(double pc,double t);
  double get_drhog_dpg(double t);
  double get_drhog_dt(double pc,double pg,double t);
  double get_mug(double pc,double pg,double t);
  double get_mg(double pc,double pg,double t);
  double get_rhocpg(double pc,double pg,double t);
  double get_cpg(double pc,double pg,double t);
  //DRY AIR; DENSITY AND DERIVATIVES
  double get_pga(double pc,double pg,double t);
  double get_rhoga(double pc,double pg,double t);
  double get_drhoga_dpg(double pc,double pg,double t);
  double get_drhoga_dpc(double pc,double pg,double t);
  double get_drhoga_dt(double pc,double pg,double t);
  double get_muga(double t);
  double get_cpga();
  //WATER VAPOR; DENSITY AND DERIVATIVES
  double get_cdiff(double pc,double pg,double t);
  double get_pgw(double pc,double t);
  double get_pcpgw(double pgw,double t);
  double get_dpgw_dpc(double pc,double t);
  double get_dpgw_dt(double pc,double t);
  double get_pgws(double t);
  double get_dpgws_dt(double t);
  double get_rhogw(double pc,double t);
  double get_pcrhogw(double rhogw,double t);
  double get_drhogw_dpc(double pc,double t);
  double get_drhogw_dt(double pc,double t);
  double get_mugw(double t);
  double get_cpgw();
  //WATER; DENSITY AND DERIVATIVES
  double get_pw(double pc,double pg,double t);
  double get_rhow(double t);
  double get_drhow_dt(double pc,double t);
  double get_dhvap(double t);
  double get_muw(double t);
  double get_cpw();
  double get_lambdaw(double t);
  double get_betaw(double t);
  double get_kw();
  //AVERAGED DENSITY
  double get_rho(double pc,double pg,double t,long ipp);
  // BIOT'S CONSTANT
  double get_alpha(double pc,double pg,double t,long ipp);
  //READ INPUT DATA
  //void read(FILE *in);
  //MATERIAL PROPERTIES
  double get_rhos(double t,long ipp);
  double get_kt(double pc,double pg,double t,long ipp);
  double get_ks(double pc,double pg,double t,long ipp);
  double get_krg(double pc,double t,long ipp);
  double get_krw(double pc,double t,long ipp);
  double get_phi(double t,long ipp);
  double get_dphi_dt(double pc, double pg,double t,long ipp);
  double get_dg(double pc,double pg,double t,long ipp);
  double get_deff(double pc,double pg,double t,long ipp);
  double get_betas(long ipp);
  double get_kintr(double pc,double pg,double t,long ipp);
  double get_cps(double t,long ipp);
  double get_rhocp(double pc,double pg,double t,long ipp);
  double get_cp(double pc,double pg,double t,long ipp);
  double get_lambdaeff(double pc,double pg,double t,long ipp);
  double get_betaswg(double pc,double pg,double t,long ipp);
  double get_betaswg_c(double pc,double pg,double t,long ipp);
  double get_betasw(double pc,double pg,double t,long ipp);
  double get_betasw_c(double pc,double pg,double t,long ipp);
  double get_betasg_c(double pc,double pg,double t,long ipp);

  /* tato cast je jenom pro beton za vysoke teploty */
  double get_dehydw_dt(double pc,double pg,double t,long ipp);
  double get_hydren(double pc,double pg,double t,long ipp);
  double get_fste(double pc,double pg,double t,long ipp);

  double pcmin; //minimal value of capillary pressure

 private: 

  double gasr;//universal gas constant
  double ma;//molar mass of dry air
  double mw;//molar mass of water
  
  double t0;//[K] reference temperature 
  double p0;//Pa reference athmospheric pressure

  // PHYSICAL PROPERTIES OF DRY AIR
  double muga0,alphaa,betaa;

  // PHYSICAL PROPERTIES OF WATER VAPOUR
  // from D.Gawin, F.Pesavento (PRVAP.f90)
  double dv0,bv;
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
  double conb;
  double conc;
  double cpw;//water specific heat at const. pressure = 4181 J/kg.K
  double lambdaw;//water heat conductivity = 0.6 W/(m.K)
  double kw0;//compresibility coefficient of water
};  

#endif
