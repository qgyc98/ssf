#ifndef MECHMAT_H
#define MECHMAT_H


#include <stdio.h>
#include "alias.h"
#include "aepoints.h"
#include "iotools.h"
#include "selection.h"
#include "intpoints.h"
#include "homogmatm.h"

class elastisomat;
class elastortomat;
class elastgmat3d;
class elastgmat2d;
class homogmatm;
class splas1d;
class j2flow;
class microM4;
class microSIM;
class microfiber;
class mohrcoulomb;
class mohrcoulombpar;
class boermat;
class drprag;
class doubdp;
class drprag2;
class bbmmat;
class dsmmat;
class camclay;
class camclaycoup;
class shefplast;
class hissplas;
class glasgmech;
class glasgowdam;
class creepdam;
class timeswmat;
class effstress;
class svisplas;
class simviscous;
class isoviscous;
class lemaitre;
class scaldam;
class fixortodam;
class scaldamcc;
class ortodam;
class ortodam2;
class ortodamrot;
class anisodam;
class anisodamrot;
class aci78;
class cebfip78;
class graphmat;
class varelastisomat;
class elastisopdmat;
class b3mat;
class rspecmat;
class dplmat;
class creepb;
class creep_effym;
class consol;
class winpast;
class nonlocmicroM4;
class therisomat;
class thervolisomat;
class nonlocplast;
class nonlocdamg;
class damplast;
class geoelastmat;
class visplast;
class viselast;
class chen;
class lenjonesmat;
class contactmat;
class cebfipcontactmat;
class damplastifacemat;
class plastifacemat;
class layplate;
class relaxeuroc;
class elasttime;
class shrinkmat;
class hypoplast;
class hypoplast_unsat_exp_therm;
class cusatismaterial;
struct matrix;
struct vector;
class sel;
class ipmap;


/**
  Class mechmat:
   
  It is one of the 5 most important classes of the program.
  (probdesc, mechtop, mechmat, mechbclc, mechcrsec)
   
  Class mechmat contains data about materials in problem.
   
  basic data:
   - nip (long) - stands for number of all integration points for stiffness matrix computation in problem
   - nmt (long) - stands for number of material types
   
  Class mechmat creates nip objects of class intpoints.
   
  Created by JK, TKo
*/
class mechmat
{
 public:
  mechmat ();
  ~mechmat ();

  /// determines number of integration points and sets integration point pointers on elements
  long intpnum (void);

  /// allocation of array of integration points
  void intpointalloc ();

  /// allocation of array of intial conditions
  void alloc_ic ();

  /// allocation of array of non-mechanical quantities
  void alloc_nonmechq(long n);

  /// allocation of array of auxiliary integration points used for the quantity transfer among different meshes
  void alloc_aux_intp(long n, ipmap *ipm);

  /// obsolate function for reading material properties on individual integration points
  void readip (XFILE *in);

  /// zero stage of int. point intialization - copying material ids from elements to int. points
  void intpointinit ();

  /// zero stage of auxiliary int. point intialization - copying material ids from elements to int. points
  void aipinit (long n, ipmap *ipm);

  /// function computes initial values from the given nodal values for auxiliary integration points
  void inic_aipval(long n, ipmap *ipm);

  /// update bit flag array hmt according to im-th material model
  void update_hmt (long ipp, long im);

  /// initialize eigenstrain/eigenstress arrays on auxiliary int. points
  void init_aipeigstr(long n, ipmap *ipm, double time);
  
  /// the first stage of int. point intialization - setting of ncompstr, ssst according to used elements
    //void init_ip_1 (void);

  /// the second stage of int. point intialization - allocation of strain/stress/eqother/other on int. points
  void init_ip_2 (void);

  /// reading of materials and their parameters
  void read (XFILE *in);

  /// reading of materials and their parameters
  void readmatchar (XFILE *in);

  /// reading of numinst material parameters for the given material type
  void readmattype (XFILE *in, mattype mtype, long numinst);

  /// printing of materials
  void print (FILE *out);

  /// printing of individual materials and their parameters
  void printmatchar (FILE *out, mattype mt, long numinst);

  /// allocation of eigstrid array on the given array of int. points
  void alloceigstrid (long** &estrid, intpoints *iparray, long nip);
  
  /// allocates array for eigenstrain storage at regular int. points
  void alloceigstrains ();

  /// allocates array for eigenstrain storage
  void alloceigstrains (double** &estra, intpoints *iparray, long nip);

  /// stes eigenstrains to zero
  void nulleigstrains ();

  /// allocates array for thermal strain storage at regular int. points
  void alloctempstrains ();

  /// allocates array for thermal strain storage
  void alloctempstrains (double** &tstra, intpoints *iparray, long nip);


  /// sets thermal strains to zero
  void nulltempstrains ();

  /// sets thermal strains to zero
  void nullaiptempstrains ();

  void nullstrains ();
  
  void nullstresses ();

  /// allocates array for eigenstress storage set it to zero
  void alloceigstresses ();

  /// allocates array for eigenstress storage set it to zero
  void alloceigstresses (double** &estre, intpoints *iparray, long nip);

  /// allocates array for macro-stresses and set it to zero (used in strain based approach to homogenization)
  void allocmacrostresses ();

  /// initailizes material models on all regular int. points before begining of main computation procedure
  void initmaterialmodels (long lcid, bool rinit);

  /// initializes material models on all auxiliary int. points before begining of main computation procedure (used in coupled problems)
  void aip_initmaterialmodels (long lcid, bool rinit);

  /// initailizes material models for the given int. point
  void initvalues (long lcid, long ipp,long im,long ido, bool rinit);

  /// returns the number of the normal components in the stress/strain state of the given int. point
  long give_num_norm_comp(long ipp);


  // *****************************************************************************
  // Function used for access to stress/strain/eqother/other arrays on int. points
  // *****************************************************************************

  /// storage of strains on given int. point for given load case
  void storestrain (long lcid,long ipp,vector &eps);

  /// storage of strain components from given index on int. point for given load case
  void storestrain (long lcid,long ipp,long fi,vector &eps);

  /// storage of given number of strain components from given index on int. point for given load case
  void storestrain (long lcid,long ipp,long fi,long ncomp,vector &eps);

  /// storage of given number of strain components from given index on int. point for given load case
  void storestrain (long lcid, long ipp, long fi, long ncomp, double *eps);

  /// returns strains on given int. point for given load case
  void givestrain (long lcid,long ipp,vector &eps);

  /// returns strain components from given index on int. point for given load case
  void givestrain (long lcid,long ipp,long fi,vector &eps);

  /// retruns given number of strain components from given index on int. point for given load case
  void givestrain (long lcid,long ipp,long fi,long ncomp,vector &eps);

  /// sets strains to zero
  void cleanstrain ();

  /// storage of stresses on given int. point for given load case
  void storestress (long lcid,long ipp,vector &sig);

  /// storage of stress components from given index on int. point for given load case
  void storestress (long lcid,long ipp,long fi,vector &sig);

  /// storage of stress components from given index on int. point for given load case
  void storestress (long lcid,long ipp,long fi,long ncomp,vector &sig);

  /// storage of stress components from given index on int. point for given load case
  void storestress (long lcid,long ipp,long fi,long ncomp,double *sig);

  /// returns stresses on given int. point for given load case
  void givestress (long lcid,long ipp,vector &sig);

  /// returns stress components from given index on int. point for given load case
  void givestress (long lcid,long ipp,long fi,vector &sig);

  /// retruns given number of stress components from given index on int. point for given load case
  void givestress (long lcid,long ipp,long fi,long ncomp,vector &sig);

  /// storage of eigenstrains for given int. point
  void storeeigstrain (long ipp,vector &eps);

  /// store given eigenstrain components for given int. point
  void storeeigstrain (long ipp,long fi,long ncomp,vector &eps);

  /// retruns eigenstrains for given int. point
  void giveeigstrain (long ipp,vector &eps);

  /// retruns given eigenstrain components for given int. point
  void giveeigstrain (long ipp,long fi,long ncomp,vector &eps);

  /// storage of eigenstresses for given int. point
  void storeeigstress (long ipp,vector &sig);

  /// storage of given eigenstress components on given int. point
  void storeeigstress (long ipp,long fi,long ncomp,vector &sig);

  /// returns eigenstresses for given int. point
  void giveeigstress (long ipp,vector &sig);

  /// returns eigenstress components from given id on given int. point
  void giveeigstress (long ipp,long fi,vector &sig);

  /// returns given eigenstress components on given int. point
  void giveeigstress (long ipp,long fi,long ncomp,vector &sig);
  
  /// storage of given internal variable components on given int. point
  void storeeqother (long ipp,long fi,long ncomp,double *comp);

  /// returns given internal variable components from the last equilibrium step on given int. point
  void giveeqother (long ipp,long fi,long ncomp,double *comp);

  /// returns given internal variable components from the last equilibrium step on given int. point
  void giveeqother (long ipp,long fi,long ncomp, vector &comp);

  /// returns given actual internal variable components on given int. point
  void giveother (long ipp,long fi,long ncomp,double *comp);

  /// returns given actual internal variable components on given int. point
  void giveother (long ipp,long fi,long ncomp, vector &comp);

  /// returns given actual internal variable components on given int. point
  void giveother (long ipp, sel &selcomp, vector &comp);

  /// storage of given internal variable components on given int. point
  void storeother (long ipp,long fi,long ncomp,double *comp);

  /// returns the number of other array components on given int. point (num. of internal variables used by material models)
  long givencompother (long ipp,long im);

  /// returns the number of eqother array components on given int. point (num. of internal variables used by material models)
  long givencompeqother (long ipp,long im);

  /// storage of given nonlocal values on given integration point
  void storenonloc (long ipp,long fi,long ncomp,double *comp);

  /// returns given nonlocal values from integration point
  void givenonloc (long ipp,long fi,long ncomp,double *comp);



  // ******************************************************************
  // Functions for retrieving of mechanical quantities from int. points
  // ******************************************************************

  /// returns tensor quantity for integral \int B^T \mbf{q} dV on element (nodal force computation)
  void givequantity (integratedquant iq,long lcid,long ipp,long fi,vector &ipv);

  /// returns actual Young modulus for given integration point
  double give_actual_ym(long ipp, long im=0, long ido=0);

  /// returns initial Young modulus for given integration point  
  double give_initial_ym(long ipp, long im=0, long ido=0);

  /// returns actual Poissons ratio for given integration point  
  double give_actual_nu(long ipp, long im=0, long ido=0);

  /// returns initial Poissons ratio for given integration point  
  double give_initial_nu(long ipp, long im=0, long ido=0);

  /// returns actual tensile strength (for damage models)
  double give_actual_ft (long ipp, long im=0, long ido=0);

  /// returns actual compression strength (for damage models)
  double give_actual_fc (long ipp, long im=0, long ido=0);

  /// returns quantity representation (scalar/vector/tensor) and for the tensor quantities, it returns quantity storage notation
  /// (Voigt reduced notation/Voigt full notation/Full matrix 3x3) and tensor type identfier (stress/strain/other)
  quantrep give_quant_rep(mechquant mq, mechquant refmq, long &ncomp, tensqnot &tn, strastre &tti);
  
  /// returns value of the required scalar quantity in the given integration point
  void give_quant(long ipp, mechquant mqn, mechquant reftensq, long lcid, sel &selrefcomp, strastre ref_strastre, double &qv);

  /// returns value of the required %vector quantity in the given integration point
  void give_quant(long ipp, mechquant mqn, long lcid, sel &selcomp, vector &qv);

  /// returns value of the required %tensor quantity in the given integration point
  void give_quant(long ipp, mechquant mqn, mechquant reftensq, long lcid, sel &selcomp, strastre array_strastre, matrix &qv);

  /// returns value of the required quantity defined in the given integration point
  void give_quant(long ipp, mechquant mq, mechquant reftensq, long lcid, quantrep qr, sel &selcomp, strastre array_strastre,
                  double &sv, vector &vv, matrix &tv);

  /// returns required mechanical quantity for passed to TRFEL
  double givemechq(nontransquant qt, long ipp);

  /// returns preconsolidation pressure for Cam-Clay type models
  double give_preconspress(long ipp, long im=0, long ido=0);

  /// returns saturation degree for hypo-plasticity models
  double give_saturation_degree(long ipp, long im=0, long ido=0);

  /// returns derivative of saturation degree with respect to suction
  double give_der_saturation_degree(long ipp, long im=0, long ido=0);

  /// returns derivative of saturation degree with respect to temperature
  double give_der_saturation_degree_dtemp(long ipp, long im=0, long ido=0);

  /// returns derivative of saturation degree with respect to volumetric strain
  double give_der_saturation_degree_depsv(long ipp, long im=0, long ido=0);
  
  /// returns derivative of saturation degree with respect to volumetric strain
  double give_bulk_modulus(long ipp, long im=0, long ido=0);
  
  /// returns virgin porostity for Cam-Clay type models
  double give_virgporosity(long ipp, long im=0, long ido=0);

  /// returns initial porosity for Cam-Clay type models
  double give_iniporosity(long ipp, long im=0, long ido=0);

  /// returns actual porosity for soil models
  double give_porosity(long ipp, long im=0, long ido=0);

  /// returns actual void_ratio for soil models
  double give_void_ratio(long ipp, long im=0, long ido=0);

  /// returns actual rate of the volumetric strain
  double give_strain_vol_rate(long ipp, long im=0, long ido=0);

  /// returns consistency parameter for plasticity models
  double give_consparam (long ipp, long im=0, long ido=0);

  /// returns damage parameter omega
  double give_dampar (long ipp,long im=0,long ido=0);

  /// returns length of process zone
  double give_proczonelength (long ipp,long im=0,long ido=0);

  /// returns crack width
  double give_crackwidth (long ipp,long im=0,long ido=0);

  /// returns equivalent strain for damage models
  void give_kappa(long ipp,long im,long ido, vector &kappa);

  /// return irreversible strains for nonlinear material models
  void giveirrstrains(long ipp, long im, long ido, vector &epsirr);




  // **********************************************************************
  // Functions for retrieving of non-mechanical quantities from int. points
  // **********************************************************************
  
  /// searches for the number and types of required non-mechanical quantities at all material models
  long search_reqnmq(nonmechquant* &rnmq);

  /// marks required types of non-mechanical quantities by the given material model in the given int. point
  void give_reqnmq(long ipp, long im, long *anmq);
  
  /// returns given non-mechanical quantity from nonmechq array on given integration point
  double givenonmechq (nonmechquant qt, long ipp);

  /// returns components of non-mechanical quantity array nonmechq on given integration point
  void give_nonmechqcomp (long ipp, sel selcomp, vector &nmqcomp);

  /// stores given non-mechanical quantity from nonmechq array on given integration point
  void   storenonmechq (nonmechquant qt, long ipp, double val);

  /// returns status of non-mechanical quantity (stored in nonmechq 1, otherwise 0)
  long   givestatusnmq (nonmechquant qt);

  /// returns status of non-mechanical quantity (stored in nonmechq 1, otherwise 0)
  long   givenonmechqid (nonmechquant qt);




  // *******************************************
  // Functions connected with stress calculation
  // *******************************************
  
  /// computes the stiffness matrix for the given integration point and material
  void matstiff (matrix &d, long ipp, long im=0, long ido=0);

  /// computes elastic stiffness matrix
  void elmatstiff (matrix &d, long ipp, long ido=0);

  /// computes elastic stiffness matrix
  void elmatstiff (matrix &d, long ipp, long im, long ido, strastrestate ssst);

  /// computes elastic stiffness matrix in the local coordinate system
  void loc_elmatstiff (matrix &d, long ipp);

  /// returns transformation %matrix to the material local coordinate system (x_g = T x_l)
  void loc_transf_mat (matrix &tmat, long ipp);

  /// computes elastic complience matrix
  void elmatcompl (matrix &c, long ipp, long ido=0);
  
  ///  computes damping matrix of material
  void matdamp (matrix &d,long ipp,long im=0,long ido=0);

  /// computes eigenstress and corresponding eigenstrains defined by relaxation material models
  void eigenstresses (vector &sig, long ipp, long id);

  /// computes stresses with respect to material model and store them in the given int. point
    void computenlstresses (long ipp, intpoints &aip, long im=0, long ido=0);

  /// computes stress increments with respect to material model and store them in the given int. point
  void computenlstressesincr (long ipp, long im=0, long ido=0);

  /// computes stress with respect to non-local material model and store them in the given int. point
  void compnonloc_nlstresses (long ipp, long im=0, long ido=0);

  /// compute averaged values of strains //???!!! candidate for removal
  void compute_averstrains ();

  /// compute strains due to temperature changes (thermal strains)
  void temprstrains ();

  /// compute strains due to temperature changes and accumulates them in Mm->tempstrains (thermal strains) on all int. points
  void cumultemprstrains ();

  /// compute strains due to temperature changes and accumulates them in Mm->tempstrains (thermal strains) at one int. point
  void cumultemprstrainsmat (long ipp);

  /// subtracts eigenstrains and thermal strains from total strains on integration points
  void eigstrmod ();

  /// subtracts eigenstrains and thermal strains from total strains at auxiliary integration points
  void aip_eigstrmod ();

  /// function adds macro-strain components to the strain %vector which contains fluctuation strains
  void add_macro_strains (long lcid, vector &macrostrains);

  /// reassambles total strains - adds eigenstrains and thermal strains to strains on integration points
  void totstrains ();
  
  /// computes thermal dilatancy matrix
  void matdilat (matrix &d, long ipp);

  /// stores appropriate volume of element belonging to given int. point
  void storeipvol (long ipp, double vol);
  
  /// computes averaged local values for non-loacl material models
  void nonlocaverage (long ipp, long im=0, long ido=0);

  /// returns the number of averaged local values.
  long give_num_averq(long ipp, long im);

  /// returns vector of averaged quantities
  void give_aver_quantv(long ipp, long im, long ido, vector &qv);

  /// returns radius of averaged neighbourhood for non-local material models
  double nonlocradius (long ipp, long im);

  /// returns id of the first non-local material model defined in the given int. point
  long givenonlocid(long ipp);

  /// updates internal variables for all regular int. points for the attained equilibrium state 
  void updateipval (void);

  /// updates internal variables for all auxiliary int. points for the attained equilibrium state 
  void update_aipval (void);

  /// updates internal variables for given int. point after attaining of equilibrium state 
  void updateipvalmat (long ipp,long im,long ido);

  /// copies internal variables for all regular int. points for the attained equilibrium state (eqother array) to the other array
  void updateother (void);

  /// copies internal variables at given int. point for the attained equilibrium state (eqother array) to the other array
  void updateothermat (long ipp,long im,long ido);

  /// updates actual material index and optionaly calls initialization of new activated material
    long update_actual_mat_id(long lcid, long init, bool rinit);

  /// sets all arrays on int. points to zero
  void clean_ip ();

  /// check time step size requirements which originate in material models for all integration points
  double dstep_red_mat();

  /// check time step size requirements which originate in material models for the given ip
  double dstep_red_ip(long ipp, long im, long ido);

  
  // **********************************************************
  // Functions used for handling of stress return in plasticity
  // **********************************************************

  /// retruns actual value of yield function for plasticity materials
  double yieldfunction (long ipp, long idpm, vector &sig, vector &q);

  /// retruns derivative of yield function with respect to stress tensor
  void dfdsigma (long ipp, long idpm, vector &sig, vector &q, vector &dfds);

  /// retruns derivative of plastic potential with respect to stress tensor
  void dgdsigma (long ipp, long idpm, vector &sig, vector &q, vector &dgds);

  /// retruns derivative of yield function with respect to %vector of hardening parameters
  void dfdqpar (long ipp, long idpm,vector &sig, vector &q, vector &dq);

  /// retruns the second derivative of yield function with respect to stress tensor
  void dfdsigmadsigma (long ipp, long idpm, vector &sig, matrix &dfdsds);

  /// retruns the second derivative of plastic potential with respect to stress tensor
  void dgdsigmadsigma (long ipp, long idpm, vector &sig, vector &q, matrix &dgdsds);

  /// retruns the second derivative of yield function with respect to stress tensor and %vector of hardening parameters
  void dfdsigmadqpar (long ipp, long idpm, vector &sig, vector &q,matrix &dfdsdq);

  /// returns the second derivative of plastic potential function with respect to stress tensor and hardening parameters
  void dgdsigmadqpar (long ipp, long idpm, vector &sig, vector &q, matrix &dgdsdq);
    
  /// the derivative of hardening parameters with respect to consistency parameter gamma
  void dqpardgamma (long ipp, long ido, long idpm, vector &sig, vector &q, vector &epsp, vector &dqdg);

  /// the second derivative of the hardening parameters with respect to gamma and stresses
  void dqpardgammadsigma (long ipp, long ido, long idpm, vector &sig, vector &q, vector &epsp, matrix &dqdgds);

  /// the second deritavives of the hardening paramaters with respect to gamma and hardening parameters
  void dqpardgammadqpar (long ipp, long ido, long idpm, vector &sig, vector &q, vector &epsp, matrix &dqdgdq);
  

  /// retruns derivative of plastic potential with respect to %vector of hardening parameters
    //void dgdqpar (long ipp, long idpm, vector &sig, vector &q, vector &dq);


  /// derivative of hardening function with respect to stress tensor
  void dhdsigma (long ipp,long idpm, vector &sig,vector &q,matrix &dhds);

  /// derivative of hardening function with respect to consistency parameter 
    //void dhdgamma (long ipp,long idpm,vector &epsp,vector &sig,vector &q,vector &dhdc);
  void hardvect (long ipp,long idpm,vector &hv);

  /// derivative of hardening function with respect to %vector of hardening parameters
  void dhdqpar (long ipp, long idpm, long ido, vector &sigt, double dgamma, vector &q, matrix &dhdq);

  /// retruns the actual value of  hardening/softening function
  void hardsoftfunction (long ipp,long idpm, vector &sig, vector &q, vector &h);

  /// returns %matrix of plastic moduli
  void plasticmoduli (long ipp,long idpm,vector &epsp,vector &sig,matrix &h);

  /// returns contribution to denominator from hardening in cutting plane stress return algorithm
  void plasmodscalar(double &denom,long ipp, long idpm, long ido, vector &sig, vector &eps, vector &epsp, vector &qtr,double gamma);

  /// returns the number of internal(hardening) parameters for plasticity models  
  long give_num_interparam (long ipp,long idpm);

  /// returns vector of hardening parameters
  void give_interparam (long ipp,long im,long ido,vector &q);

  /// updates of hardening parameters 
  void updateq (long ipp, long idpm, long ido, double dgamma, vector &eps, vector &epsp, vector &sig, vector &q);
  
  /// cutting plane stress return algorithm
  void cutting_plane (long ipp, long im, long ido, double &gamma, vector &epsn, vector &epsp, vector &q, long ni, double err);
  
  /// cutting plane stress return algorithm with algorithmic consistent matrix evaluation
  long cutting_plane (long ipp, long im, long ido, double &gamma, vector &epsn, vector &epsp, vector &q, long ni, double err,
                      long dcomp, matrix &d);

  /// cutting plane stress return algorithm with general modification for plane stress
  void cutting_plane3 (long ipp,long im,long ido,double &gamma,vector &epsn,vector &epsp,vector &q,long ni,double err);

  /// multi-surface cutting plane stress return algorithm
  void mult_surf_cutting_plane (long ipp,vector &epsn,vector &epsp,vector &gamma,vector &q);

  /// closest point projection stress return algorithm
  void newton_stress_return (long ipp,long im,long ido,double &gamma,vector &epsn,vector &epsp,vector &q,long ni,double err);

  /// closest point projection stress return algorithm
  void newton_stress_return_2 (long ipp,long im,long ido,double &gamma,vector &epsn,vector &epsp,vector &q,long ni,double err);

  /// detection of plasticity for adaptivity computations
  void refresh_plast (long im=0, long ido=0);
  //void dfunctdqdq (long ipp, long idpm, matrix &sig, vector &q, matrix &dfdqdq);
  //void store_interparam (long ipp,long im,long ido,vector &q);



  // *************************************************
  // Functions for handling of scalar damage materials
  // *************************************************

  /// computation of equivalent strain
  void damfuncpar(long ipp, long im, vector &eps, vector &kappa);

  /// computation of new damage parameter omega
  double damfunction(long ipp, long im, long ido,vector &kappa, vector &eps, vector &omegao);

  /// computation of intial elsatic strain for damage models
  double epsefunction(long ipp, long idem);

  /// solver for scalar damage models
  double scal_dam_sol (long ipp,long im,long ido,vector &eps,vector &kappa,vector &omegao,vector &sigma,matrix &d);
  double scal_dam_sol2 (long ipp,long im,long ido,vector &eps,vector &kappa,vector &omegao,vector &sigma,matrix &d);
  


  // *************************************************
  // Functions for handling of visco-plastic materials
  // *************************************************

  /// returns stress increments due to visco-plastic strains on int. point for given load case
  void givestressincr (long lcid,long ipp,long im,long ido,long fi,vector &sig);

  //void storestressincr (long lcid,long ipp,long im,long ido,vector &sig);

  /// computation of stresses Perzyna type visco-plastic materials
  void stiff_deps_vispl (long ipp,long im,long ido,vector &sig,vector &q,double dt);

  /// computation of contribution to right hand side for creep models - obsolate, candidate for removal
  void stiff_eps (long ipp);



  // *************************************************
  // Functions connected with Len-Jons potential model
  // *************************************************
 
  /// returns the first derivative of lenjones potential
  double give_first_derivative (long ipp,double r);

  /// returns the second derivative of lenjones potential
  double give_second_derivative (long ipp,double r);

  // *****************************************
  // Funtions for backup of int. point content
  // *****************************************

  /// save content of int. points to one huge text file
  void save_intpoints_txt (FILE *aux, sel &selelems, sel *selother);

  /// save content of int. points to several text files according to stored quantity
  void save_intpoints_txt (long ni, sel &selelems, sel *selother);

  /// save content of int. points to one huge binary file
  void save_intpoints_bin (FILE *aux, sel &selelems, sel *selother);

  /// save content of int. points to several binary files according to stored quantity
  void save_intpoints_bin (long ni, sel &selelems, sel *selother);

  /// restore content of int. points to one huge text file
  void restore_intpoints_txt (FILE *aux, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several text files according to stored quantity
  void restore_intpoints_txt (sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to one huge binary file
  void restore_intpoints_bin (FILE *aux, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several binary files according to stored quantity
  void restore_intpoints_bin (sel &selelems, sel *selother, long **selid);

  // ***************************************************
  // Funtions for backup of auxiliary int. point content
  // ***************************************************

  /// save content of int. points to one huge text file
  void save_auxintpoints_txt (FILE *aux, sel &selelems, sel *selother);

  /// save content of int. points to several text files according to stored quantity
  void save_auxintpoints_txt (long ni, sel &selelems, sel *selother);

  /// save content of int. points to one huge binary file
  void save_auxintpoints_bin (FILE *aux, sel &selelems, sel *selother);

  /// save content of int. points to several binary files according to stored quantity
  void save_auxintpoints_bin (long ni, sel &selelems, sel *selother);

  /// restore content of int. points to one huge text file
  void restore_auxintpoints_txt (FILE *aux, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several text files according to stored quantity
  void restore_auxintpoints_txt (sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to one huge binary file
  void restore_auxintpoints_bin (FILE *aux, sel &selelems, sel *selother, long **selid);

  /// restore content of int. points to several binary files according to stored quantity
  void restore_auxintpoints_bin (sel &selelems, sel *selother, long **selid);

  
  

  
  ///  indicator of plastic regime in domain
  long plast;
  
  ///  type of eigstrain
  eigstraintype est;
  
  ///  number of material types
  long nmt;
  ///  actual total number of all integration points (regular or auxiliary - used in the int. point loops)
  long tnip;
  ///  total number of regular integration points (used for integration on elements)
  long tnrip;
  ///  integration points
  intpoints *ip;
  /**  integration point - element correspondation
       elip[number of int point] = number of element */
  long *elip;

  ///  total number of auxiliary integration points (coupled problems - quantity transfer among meshes)
  long tnaip;
  ///  auxiliary integration points
  intpoints *aip;
  /**  auxiliary integration point - element correspondation
       elaip[number of aux int point] = number of element */
  long *elaip;

  ///  pointers to local coordinate systems in integration points
  long *lcsip;
  ///  array containing local coordinate system basis for integration points
  double **lcsbase;

  /// array with initial conditions at regular integration points (allocated at mechbclc::readinic)
  double **ic;
  
  /// array with initial conditions at auxiliary integration points (allocated at mechbclc::aip_inic)
  double **aip_ic;

  ///  array containing volume appropriate to integration points
  double *ipv;
  /// maximum number of component of strain/stress arrays at nodes
  long max_ncompstrn;
  /// maximum number of component of strain/stress array on elements
  long max_ncompstre;
  /// maximum number of component of other array at nodes
  long max_ncompothern;
  /// maximum number of component of eqother array on elements
  long max_ncompothere;

  ///  array containing eigenstrains id
  ///  the array is allocated in function mechbclc::read
  long **eigstrid;
  ///  array containing eigenstrains
  ///  the array is allocated only for several problems
  ///  eigstrains[number of int. point][number of component]
  ///  the array is allocated in function mechbclc::read
  double **eigstrains;
  ///  array containing eigenstresses
  ///  the array is allocated only for several problems
  ///  eigstresses[number of int. point][number of component]
  ///  the array is allocated in function mechbclc::read
  double **eigstresses;
  ///  array containing strains caused by temperature
  ///  tempstrains[number of int. point][number of component]
  double **tempstrains;

  ///  array containing eigenstrains id at auxiliary int. points
  ///  the array is allocated in function aipinit
  long **aip_eigstrid;
  ///  array containing eigenstrains at auxiliary int. points
  ///  the array is allocated only for several problems
  ///  aip_eigstrains[number of aux. int. point][number of component]
  ///  the array is allocated in function mechbclc::read
  double **aip_eigstrains;
  ///  array containing eigenstresses at auxiliary int. points
  ///  the array is allocated only for several problems
  ///  aip_eigstresses[number of aux int. point][number of component]
  ///  the array is allocated in function mechbclc::read
  double **aip_eigstresses;
  ///  array containing strains caused by temperature at auxiliary int. points
  ///  aip_tempstrains[number of aux. int. point][number of component]
  double **aip_tempstrains;

  /// array of macro-stresses (used in macro-strain based approach in homogenization - Mp->homog==4)
  double *mstress;
    

  ///  array containing non-mechanical quantities for regular integration points used in coupled problems
  ///  nonmechq[quantity no.][number of integ. point] - nonmechanical quantity at the required integration point
  double **nonmechq;

  ///  array containing non-mechanical quantities for auxiliary integration points used in coupled problems
  ///  nonmechq[quantity no.][number of aux. integ. point] - nonmechanical quantity at the required integration point
  double **aip_nonmechq;

  /// number of actually used non-mechanical quantities (it is set by METR or according to temperature in the input file)
  long nnmq;
  /// array with ordering of used non-mechanical quantities (its size is nnmq)
  nonmechquant *nmqo;
  /// total number of known(=implemented) non-mechanical quantities
  static const long tnknmq = sizeof(nonmechquantstr)/sizeof(*nonmechquantstr);
  /// array of indices in nonmechq for particular non-mechanical quantities
  long nmqid[tnknmq];

  ///  elastic isotropic linear material
  elastisomat *eliso;
  ///  elastic orthotropic linear material
  elastortomat *elorto;
  ///  elastic general 3D material
  elastgmat3d *elgm3d;

  ///  elastic general 2D material
  elastgmat2d *elgm2d;
  
  //  artificial material obtained from homogenization
  homogmatm *hommatm;

  ///  simple plastic 1D material
  splas1d *spl1d;
  ///  J2 flow material
  j2flow *j2f;
  ///  Mohr-Coulomb material
  mohrcoulomb *mcoul;
  ///  Parabolic Mohr-Coulomb material
  mohrcoulombpar *mcpar;
  ///  Boer model of plasticity
  boermat *boerm;
  ///  Drucker-Prager model of plasticity
  drprag *drprm;
  /// Double Drucker-Prager model
  doubdp *ddpm;
  ///  Drucker-Prager model of plasticity prepare by Standa Sysala from IGN
  drprag2 *drprm2;
  ///  double structure model
  dsmmat *dsm;
  ///  modified barcelona basic model
  bbmmat *bbm;
  ///  modified cam-clay model
  camclay *cclay;
  ///  modified cam-clay model coupled with moisture transfer
  camclaycoup *cclayc;
  ///  Sheffield plasticity model
  shefplast *shpl;
  ///  HISS plasticity model (Desai model)
  hissplas *hisspl;
  /// Glasgow scalar damage with a thermal influence  
  glasgmech *glasgmat;
  /// Glasgow damage model
  glasgowdam *glasgdam;
  ///  artificial model for combination of damage and creep model
  creepdam *crdam;
  ///  artificial model for combination of several time dependent material models switched on in the given time 
  timeswmat *tswmat;
  /// number of artificial model for combination of several time dependent material models
  long      ntswmat;
  /// artificial material for effective stress concept   
  effstress *effstr;
  ///  microplane M4
  microM4 *mpM4;
  /// simplifed version of microplane
  microSIM *mpSIM;
  ///
  microfiber *mpfib;
  
  ///  simple model of viscosity
  simviscous *svis;
  ///  viscous model for isotropic material
  isoviscous *isovis;
  ///  simple visco-plastic material
  svisplas *svipl;
  ///  Lemaitre visco-plastic model
  lemaitre *lmtr;

  ///  isotropic damage material model
  scaldam *scdam;

  ///  orthotropic damage material model
  fixortodam *fixodam;

  ///  simple damage 3D material with crack closure
  scaldamcc *scdamcc;

  ///  ortotropic damage material model
  ortodam *ortdam;
  ///  ortotropic damage material model with orthotropic elastic material
  ortodam2 *ortdam2;
  ///  ortotropic damage material model with crack rotation
  ortodamrot *ortdamrot;
  ///  anisotropic damage material model
  anisodam *anidam;
  ///  anisotropic damage material model with crack rotation
  anisodamrot *anidamrot;
  ///  aci78 model for creep and shrinkage
  aci78 *aci78mod;
  ///  cebfip78 model for creep and shrinkage
  cebfip78 *cebfip78mod;

  /// material with explicit prescribed diagram force/displacement
  graphmat *grmat;

  /// material with variable elastic parameters E and nu
  varelastisomat *veliso;
  ///  elastic isotropic linear material with pressure dependent initial stiffness
  elastisopdmat *elisopd;

  /// combination of plasticity and special kind of elasticity (different loading and unloading modulus)
  geoelastmat *geoel;

  ///  B3 model of creep based on Bazant theory
  b3mat *crb3;
  ///  retardation spectrum
  rspecmat *crrs;
  ///  Double power law model of creep based on Bazant theory
  dplmat *crdpl;
  ///  model of creep based on Bazant theory
  creepb *crbaz;
  ///  model of creep based on effective Young modulus
  creep_effym *creffym;
  ///  consolidation model
  consol *csol;
  ///  Winkler-Pasternak material model of subsoil
  winpast *wpast;

  nonlocmicroM4 *nmpM4;

  ///  thermal isotropic dilatancy model
  therisomat *tidilat;
  ///  thermal and volume isotropic dilatancy model
  thervolisomat *tvidilat;
  
  ///  model of stress relaxation based on Eurocode 2
  relaxeuroc *relaxec;
  
  /// nonlocal plasticity models
  nonlocplast *nlplast;

  /// nonlocal damage models
  nonlocdamg *nldamg;

  ///  artificial model for combination of damage and plasticity model
  damplast *dampl;
  /// artificail model for combination of plasticity and viscous material models
  visplast *visplas;
  /**  artificial model for elastic and viscous material model
       the model is intended for definition of damping */
  viselast *viselas;
  
  ///  model with time dependent elastic modulus for time dependent problems
  elasttime *eltimemat;

  ///  Chen model of plasticity of concrete
  chen *chplast;

  ///  strain/stress points
  aepoints stra,stre;
  
  ///  Lennard-Jones interatomic potential
  lenjonesmat *lenjon;

  ///  material model for contact problems
  contactmat *conmat;

  ///  material model for contact problems according to MC90
  cebfipcontactmat *cebfipconmat;

  ///  damage plasticity material model for interface elements
  damplastifacemat *damplifm;

  ///  plasticity material model for interface elements
  plastifacemat *plastifm;

  /// layered model for plates
  layplate *lplate;

  /// hypoplastic material for unsaturated soils
  hypoplast *hypopl;

  /// hypoplastic material for unsaturated soils with thermal effects
  hypoplast_unsat_exp_therm *hypoplustherm;

  ///  material types
  mattype *mtype;
  ///  numbers of particular material types
  long *numtype;
  
  ///  material for shrinkage strain computation
  shrinkmat *shmat;
  
  ///  material for the Cusatis lattice discrete particle models
  cusatismaterial *cusmat;

};

#endif
