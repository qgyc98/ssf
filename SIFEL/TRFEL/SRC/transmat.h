#ifndef TRANSMAT_H
#define TRANSMAT_H

#include <stdio.h>
#include "iotools.h"
#include "aliast.h"
#include "selection.h"
#include "intpointst.h"
#include "isotrmat.h"
#include "tdisotrmat.h"
#include "sejtkrmat.h"
#include "nlisotrmat.h"
#include "cerny_concrete.h"
#include "damisotrmat.h"
#include "bazped.h"
#include "baroghelB.h"
#include "concreteB.h"
#include "pedersen.h"
#include "carbmat1.h"
#include "carbmat1.h"
#include "dampermeability.h"
#include "kunzel.h"
#include "kunzel2.h"
#include "glasgowmat.h"
#include "genfile.h"
#include "C60baroghel.h"
#include "C30baroghel.h"
#include "o30bazant.h"
#include "C60bazant.h"
#include "C30bazant.h"
#include "saltmat1.h"
#include "saltmat2.h"
#include "saltmat3.h"
#include "saltmat4.h"
#include "grunewaldmat.h"
#include "aepointst.h"
#include "soil1mat.h"
#include "discisotrmat.h"
#include "cemhydmat.h"
#include "devriesmat.h"
#include "lincoupmat.h"
#include "discmat.h"
#include "consol_awf1.h"
#include "consol_wf1.h"
#include "consol_wf2.h"
#include "consol_awf2.h"
#include "consol_hawf3.h"
#include "consol_hwf2.h"
#include "millymat.h"
#include "homogmat.h"
#include "radiationmat.h"
#include "richards.h"
#include "moistheat.h"
#include "interfacematt.h"

class ipmap;

class transmat
{
 public:
  transmat (void);
  ~transmat (void);


  // *********************************************
  // Functions for alocation and reading materials
  // *********************************************

  /// determines number of integration points and sets integration point pointers on elements
  long intpnum (void);

  /// obsolate function for reading material properties on individual integration points
  void readip (FILE *in);

  /// allocation of array of integration points
  void intpointalloc ();

  /// allocation of array of non-transport quantities
  void alloc_nontransq(long n);

  /// allocation of array of auxiliary integration points used for the quantity transfer among different meshes
  void alloc_aux_intp(long n, ipmap *ipm);

  /// zero stage of int. point intialization - copying material ids from elements to int. points
  void intpointinit ();

  /// zero stage of auxiliary int. point intialization - copying material ids from elements to int. points
  void aipinit (long n, ipmap *ipm);

  /// assembles array with dof names and establishes the dof order  
  void assemble_dof_nameord(namevart *dofname, long ntm, long *var_dofid, long tnkv);

  /// returns dof names used in the material of given integration point
  void give_dof_names(long ipp, namevart *dofname, long ntm);

  /// reading of materials and their parameters
  void read (XFILE *in);

  /// printing of materials
  void print (FILE *out);

  /// reading of materials and their parameters
  void readmatchar (XFILE *in);

  /// reading of material parameters for the given material type
  void readmattype(XFILE *in, mattypet mtype, long numt);

  /// printing of individual materials and their parameters
  void printmatchar (FILE *out, mattypet mt, long numinst);

  /// initailizes material models for all int. points before begining of main computation procedure
  void initmaterialmodels (void);

  /// initializes material models on all auxiliary int. points before begining of main computation procedure (used in coupled problems)
  void aip_initmaterialmodels (void);

  /// initailizes material models for the given int. point
  void initvalues (long ipp,long im,long ido);

  /// computes initial condition values on auxiliary integration points
  void inic_aipval(long n, ipmap *ipm);

  /// updates internal variables for all regular int. points for the attained equilibrium state 
  void updateipval (void);

  /// updates internal variables for all auxiliary int. points for the attained equilibrium state 
  void update_aipval (void);

  void updateipvalmat (long ipp,long im,long ido);

  // ***********************************************************
  // Functions for conductivity and capacity matrices and fluxes
  // ***********************************************************

  void matcond (matrix &d,long ipp,long ri,long ci);


  void matcond2 (matrix &d,long ipp,long ri,long ci);

  ///  function returns capacity coefficient
  double capcoeff (long ipp,long ri,long ci);
  ///  function returns reaction coefficient for reaction-diffusion problems
  double reactcoeff (long ipp,long ri,long ci);


  void volume_rhs (matrix &d,long ipp,long ri,long ci,long ncomp);

  double volume_rhs2 (long ipp,long ri,long ci);

  void computenlfluxes (long lcid,long ipp);
  void computenlcapacities (long lcid,long ipp);

  void flux_contributions (long ipp);


  void transmission_transcoeff(double &new_trc,double trc,long ri,long ci,long nn,long bc,long ipp);


  void transmission_transcoeff(double &new_trc,double trc,long ri,long ci,long nn,long bc,long ipp,int flag);


  void transmission_nodval(double &new_nodval,double nodval,double trc2,long ri,long ci,long nid,long bc,long ipp);


  void transmission_flux(double &flux,double nodval,double trc2,long ri,long ci,long nn,long bc,long ipp);


  void values_correction_puc ();

  void values_correction_ipp (long ipp);

  void values_correction (vector &nv,long ipp);
  
  
  void values_transformation (mattypet mti,long idmi,double *inv,
			      double *iinv,mattypet mto,long idmo,double *outv,double *ioutv,double *jum);
  
  void fluxcomparing (long ipp);


  void aux_values (long lcid);
  
  void mat_aux_values (long ipp);

  
  double give_extinction_coeff (mattypet mt,long matid);


  // ************************************************
  // Functions for gradients and fluxes, other values
  // ************************************************

  void storegrad (long lcid,long ipp,vector &gr);


  void givegrad (long lcid,long ipp,vector &gr);


  void storeflux (long lcid,long ipp,vector &fl);


  void givefluxes (long lcid,long ipp,vector &fl);

  
  /// storage of given other variable components on given int. point
  void storeother (long ipp, long fi, long ncomp, double *comp);

  /// storage of id-th component to the other array on the given int. point
  void storeother (long ipp, long id, double comp);

  /// storage of given eqother variable components on given int. point
  void storeeqother (long ipp, long fi, long ncompeq, double *compeq);

  /// returns given other variable components on given int. point
  void giveother (long ipp, long fi, long ncomp, double *comp);

  /// returns given eqother variable components on given int. point
  void giveeqother (long ipp, long fi, long ncompeq, double *compeq);


  long givencompother ();


  /// returns the number of components of ipp's eqother array at the given ip
  long givencompeqother (long ipp,long im);

  /// returns the number of components of ipp's eqother array of the given matrial type
  long givencompeqothermat (mattypet tm);

  double givecompother(long compother,long ipp,double *r);


  void give_othervalue_name(FILE *out,long ipp,long compother);


  void give_eqothervalue_name(FILE *out,long ipp,long compeqother);
  
  ///  function returns the temperature in integration point
  double give_temperature (long ipp);

  ///  function returns the initial temperature in integration point
  double give_inittemperature (long ipp);

  ///  function returns the relative humidity in integration point
  double give_rel_hum (long ipp);

  /// returns required transport quantity for passed to MEFEL
  double givetransq(nonmechquant qt, long ipp);

  /// returns water pressure
  double give_water_pressure(long ipp);

  /// returns air pressure
  double give_gas_pressure(long ipp);

  /// returns effective pore pressure
  double give_effective_pore_pressure(long ipp);

  /// returns pore pressure
  double give_pore_pressure(long ipp);

  /// returns suction
  double give_suction(long ipp);

  /// returns saturation degree
  double give_saturation_degree(long ipp);

  /// returns volume change
  double give_volume_change (long ipp);

  /// returns volumetric moisture content
  double give_vol_moist_cont(long ipp);
  
  ///  returns volumetric moisture content stored in nodes
  double give_nodal_vol_moist_cont (long nid,long mattype);

  ///  returns saturated volumetric moisture content stored in nodes
  double give_nodal_sat_vol_moist_cont (long nid,long mattype);

  ///  function returns relative humidity for individual material models
  double give_nodal_rel_hum (long nid,long mattype);

  // **********************************************************************
  // Functions for retrieving of non-transport quantities from int. points
  // **********************************************************************

  /// searches for the number and types of required non-transport quantities at all material models
  long search_reqntq(nontransquant* &rntq);

  /// marks required types of non-transport quantities by the given material model in the given int. point
  void give_reqntq(long ipp, long *antq);

  /// returns given non-transport quantity from nontransq array on given integration point
  double givenontransq (nontransquant qt, long ipp);

  /// stores given non-transport quantity from nontransq array on given integration point
  void storenontransq (nontransquant qt, long ipp, double val);

  /// returns status of non-transport quantity (stored in nontransq 1, otherwise 0)
  long givestatusntq (nontransquant qt);

  /// returns status of non-transport quantity (stored in nontransq 1, otherwise 0)
  long givenontransqid (nontransquant qt);


  void advection_velocities (XFILE *in);

  // *****************************************
  // Funtions for backup of int. point content
  // *****************************************

  /// save content of int. points to one huge text file
  void save_intpointst_txt (FILE *aux, sel &selelems, sel *selother);

  /// save content of int. points to several text files according to stored quantity
  void save_intpointst_txt (long ni, sel &selelems, sel *selother);

  /// save content of int. points to one huge binary file
  void save_intpointst_bin (FILE *aux, sel &selelems, sel *selother);

  /// save content of int. points to several binary files according to stored quantity
  void save_intpointst_bin (long ni, sel &selelems, sel *selother);

  /// restore content of int. points to one huge text file
  void restore_intpointst_txt (FILE *aux, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several text files according to stored quantity
  void restore_intpointst_txt (sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to one huge binary file
  void restore_intpointst_bin (FILE *aux, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several binary files according to stored quantity
  void restore_intpointst_bin (sel &selelems, sel *selother, long **selid);

  // ***************************************************
  // Funtions for backup of auxiliary int. point content
  // ***************************************************

  /// save content of int. points to one huge text file
  void save_auxintpointst_txt (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother);

  /// save content of int. points to several text files according to stored quantity
  void save_auxintpointst_txt (long ni, long n, ipmap *ipm, sel &selelems, sel *selother);

  /// save content of int. points to one huge binary file
  void save_auxintpointst_bin (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother);

  /// save content of int. points to several binary files according to stored quantity
  void save_auxintpointst_bin (long ni, long n, ipmap *ipm, sel &selelems, sel *selother);

  /// restore content of int. points to one huge text file
  void restore_auxintpointst_txt (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several text files according to stored quantity
  void restore_auxintpointst_txt (long n, ipmap *ipm, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to one huge binary file
  void restore_auxintpointst_bin (FILE *aux, long n, ipmap *ipm, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several binary files according to stored quantity
  void restore_auxintpointst_bin (long n, ipmap *ipm, sel &selelems, sel *selother, long **selid);
  



  long cycle_detection (double *r,double *pr,double *ppr,long ipp);

  void freezing_thawing ();


  
  ///  number of material types
  long nmt;
  ///  material types
  mattypet *mattype;
  ///  numbers of particular material types
  long *numtype;
  long tnip;
  ///  pointer to integration points
  intpointst *ip;
  ///  integration point - element correspondation
  ///  elip[number of int point] = number of element
  long *elip;
  
  ///  total number of auxiliary integration points (coupled problems - quantity transfer among meshes)
  long tnaip;
  ///  auxiliary integration points
  intpointst *aip;
  /**  auxiliary integration point - element correspondation
       elaip[number of aux int point] = number of element */
  long *elaip;

  // ******************
  //  material models
  // ******************

  // one-phase medium
  //  isotropic material
  isotrmat *itrm;
  //  one-phase medium - heat transfer
  //  nonlinear isotropic material
  nlisotrmat *nlitrm;
  //  artificial material obtained from homogenization
  homogmat *hommat;
  // one-phase medium - moisture transfer
  //  isotropic material with optional influence of damage
  damisotrmat *damitrm;
  //  time dependent isotropic diffusion 
  tdisotrmat *tditrm;
  //  material for contact elements
  interfacematt *ifacemat;
  
  //  saturated-unsaturated one-phase flow in deforming medium
  sejtkrmat *sejtkrm;
  con_awf1mat *consol_awf1;
  con_wf1mat *consol_wf1;

  ///  one-phase medium - heat transfer
  ///  isotropic material with jumps
  discisotrmat *ditrm;

  carbmat1 *carb1;

  // concrete from Cerny
  cernymat *cernym;
  
  // heat-moisture transport in concrete
  //two-phase medium
  kunmat *kun;
  kunmat2 *kun2;
  grunewaldmat *grunw;
  bazpedmat *bazped;
  pedmat *ped;
  devriesmat *dvries;
  discmat *sdmat;
  con_wf2mat *consol_wf2;
  con_awf2mat *consol_awf2;
  con_hwf2mat *consol_hwf2;
  millymat *mill;
  moistheatmat *moisth;

  //three-phase medium
  baroghelmat *baroghel;
  concreteBmat *concrete;
  C60barmat  *C60baroghel;
  C30barmat  *C30baroghel;
  o30bazmat  *o30bazant;
  C60bazmat  *C60bazant;
  C30bazmat  *C30bazant;
  glasgowmat *tench;
  
  dampermeability *damper;

  con_hawf3mat *consol_hawf3;

  // moisture-salt transport in porous mat.
  //two-phase medium
  saltmat1 *salt1;
  
  //  three-phase
  saltmat2 *salt2;
  saltmat3 *salt3;
  //  four-phase materials
  saltmat4 *salt4;
  
  radiationmat *radmat;
  richards *richar;
  
  soil1mat *soil1;

  //cemhyd model
  cemhydmat *cemhydr;
  
  //  linear coupled model (for arbitrary number of matters)
  lincoupmat *lcmat;
  
  
  ///  strain/stress points
  aepointst grad,flux;
  
  
  ///  array containing initial values of selected unknown at regular integration points
  ///  initval[number of integ. point] - initial value at the required integration point
  double *initval;
  ///  array containing initial values of selected unknown at auxiliary integration points
  ///  aip_initval[number of integ. point] - initial value at the required integration point
  double *aip_initval;
  ///  array containing nontransport quantities used for regular integration points used in coupled problems
  ///  nontransq[quantity id][integ. point number] - nontransport quantity at the required integration point
  double **nontransq;
  ///  array containing nontransport quantities used for auxiliary integration points used in coupled problems
  ///  nontransq[quantity id][integ. aux point number] - nontransport quantity at the required integration point
  double **aip_nontransq;
  ///  number of nontransport quantities
  long nntq;
  /// total number of known(=implemented) non-transport quantities
  static const long tnkntq = sizeof(nontransquantstr)/sizeof(*nontransquantstr);
  /// array of indices in nontransq for particular non-transport quantities
  long ntqid[tnkntq];
  /// array with ordering of used non-transport quantities
  nontransquant *ntqo;
};

#endif
