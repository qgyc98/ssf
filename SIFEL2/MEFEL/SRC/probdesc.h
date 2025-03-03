#ifndef PROBDESC_H
#define PROBDESC_H

#include <stdio.h>
#include "alias.h"
#include "galias.h"
#include "gfunct.h"
#include "timecontr.h"
#include "hdbcontr.h"
#include "slesolv.h"
#include "nonlinman.h"
#include "eigvalsol.h"
#include "xfile.h"


class auxdatout;
class gmatrix;

/**
   class probdesc
   
   it is one of the five most important classes of the program
   (probdesc, mechtop, mechmat, mechbclc, mechcrsec)
   
   
   JK, TK
*/

class probdesc
{
 public:
  probdesc ();
  ~probdesc ();
  void read (XFILE *in);
  void print (FILE *out);
  
  //  for open dx
  long  detnodstrain;

  /// keyword processing option
  pkwd_sw kwdsw;

  ///  problem name
  char name[1001];

  ///  type of problem
  problemtype tprob;

  ///  STORAGE TYPES
  ///  type of storage of stiffness matrix
  storagetype tstorsm;
  ///  type of storage of mass matrix
  storagetype tstormm;
  ///  type of storage of damping matrix
  storagetype tstordm;
  ///  type of storage of initial stress (geometric stiffness) matrix
  storagetype tstorgm;

  
  ///  DAMPING
  ///  type of damping
  damping damp;
  ///  mass matrix coefficient
  double dmass;
  ///  stiffness matrix coefficient
  double dstiff;

  ///  DIAGONALIZATION
  /// diagonalization of mass matrix
  long diagmass;
  

  ///  time
  double time;
  /// actual time increment
  double dtime;
  ///  time controller
  timecontr timecon;
  /// load coefficient
  double lambda;
  /// actual load coeficient increment
  double dlambda;
  /// index of actual load step (arclength,. ..)
  long istep;
  /// index of actual iteration step (Newton-Raphson)
  long jstep;
  
  
  ///  number of subdomains in hemivariational inequalities
    ///long nsub;

  
  ///  SOLVER TYPES
  ///  type of eigenvibration solver
  eigvalsol eigsol;
  ///  type of forced vibration solver
  forcedsolver tforvib;
  /// type of earth pressure solver
  epsolvertype tepsol;
  
  ///  data about solver of system of linear equations
  slesolv *ssle;
  ///  data about solver of system of nonlinear equations
  nonlinman *nlman;
  ///  data about solver of eigenvalues and eigenvectors
  //eigvalsol eig;
  

  ///  threshold for rejection from compressed storages
  double limit;
  ///  computer zero
  double zero;

  

  
  ///  FORCED VIBRATION SOLVER
  ///  coefficients of the Newmark method
  double alphafvn,deltafvn;
  
  

  //  VISCO-PROBLEM SOLVER
  
  /**  Flag for non-local matrial models
       phase=1 - evaluation of local variables
       phase=2 - evaluation of non-local variables */
  long nonlocphase;

  /**  variable for computation of internal forces at nodes
       nodeintfor=0 - internal forces are not computed at nodes
       nodeintfor=1 - internal forces are  computed at nodes */
  long nodeintfor;

  /**  Stress computation flag
       strcomp=0 - stresses are not computed during evaluation of internal forces
                   this is used for initiation of static nonlinear problems
       strcomp=1 - stresses are computed during evaluation of internal forces
                   except of the previous case, this is always used */
  long strcomp;


  /**  Eigenstress computation flag
       strcomp=0 - Eigenstresses are not computed during evaluation of eigenstress forces
                   this is used for computation of nodal forces due to pore pressures in partialy coupled problems
       strcomp=1 - Eigenstresses are computed during evaluation of eigenstress forces (due to prescribed eigenstrains or 
                   due to thermal starins. Except of the previous case, this is always used */
  long eigstrcomp;
  
  
  //  EARTH-PRESSURE SOLVER
  ///  maximum number of increments
  long niep;
  ///  required norm of vector of unbalanced forces
  double errep;
  /// load step size for the load iteration algorithm
  double stepep;
  /// required error of lambda value in the arclength method
  double rlerrep;
  /// number of selected nodes for earth pressure problems, where will be variable supports
  long nselnodep;
  /// selected nodes for earth pressure problems, where will be variable supports
  long *selnodep;
  /// coeficient for static load, which simulates excaving material; number of coef. equals to niep
  double *loadcoef;



  ///  type of material model (local x non-local)
  materialmodel matmodel;
  
  ///  compute strains
  ///  straincomp=0 - strains are not computed
  ///  straincomp=1 - strains are computed
  long straincomp;
  ///  strain averaging
  ///  strainaver=0 - values are not averaged
  ///  strainaver=1 - values are averaged, arithmetic average is used
  ///  strainaver=2 - values are averaged, volume average is used
  long strainaver;
  ///  state of strain computation
  ///  strainstate = 0 - strains at integration points have not been computed
  ///  strainstate = 1 - strains at integration points have been computed
  long strainstate;
  ///  position where strains are computed
  ///  strainpos = 1 - integration points
  ///  strainpos = 2 - nodes, values are copied from the closest integration point
  ///  strainpos = 3 - nodes, values are computed at nodes
  long strainpos;

  ///  compute stresses
  ///  stresscomp=0 - stresses are not computed
  ///  stresscomp=1 - stresses are computed
  long stresscomp;
  ///  stress averaging
  ///  stressaver=0 - values are not averaged
  ///  stressaver=1 - values are averaged, arithmetic average is used
  ///  stressaver=2 - values are averaged, volume average is used
  long stressaver;
  ///  state of stress computation
  ///  stressstate = 0 - stresses at integration points have not been computed
  ///  stressstate = 1 - stresses at integration points have been computed
  long stressstate;
  ///  position where stresses are computed
  ///  stresspos = 1 - integration points
  ///  stresspos = 2 - nodes, values are copied from the closest integration point
  ///  stresspos = 3 - nodes, values are computed at nodes
  long stresspos;

  ///  compute other values
  ///  othercomp=0 - components of other array are not computed
  ///  othercomp=1 - components of other array are computed
  long othercomp;
  ///  other values averaging
  ///  otheraver=0 - components of other array are not averaged
  ///  otheraver=1 - components of other array are averaged, arithmetic average is used
  ///  otheraver=2 - components of other array are averaged, volume average is used
  long otheraver;
  ///  state of other computation
  ///  otherstate = 0 - other components have not been computed
  ///  otherstate = 1 - other components have been computed
  long otherstate;
  ///  position where components of array other are computed
  ///  otherpos = 1 - integration points
  ///  otherpos = 2 - nodes, values are copied from the closest integration point
  long otherpos;

  ///  compute reactions
  ///  reactcomp=0 - reactions are not computed
  ///  reactcomp=1 - reactions are computed
  long reactcomp;
  
  ///  adaptivity calculation indicator
  ///  adaptivity == 0 - no adaptivity calculation
  ///  adaptivity == 1 - adaptivity calculation - SPR type
  ///  adaptivity == 2 - adaptivity calculation - Z2 type
  long adaptivityflag;
  
  ///  stochastic calculation indicator
  ///  stochasticcalc=0 - deterministic computation
  ///  stochasticcalc=1 - stochastic/fuzzy computation, data are read all at once
  ///  stochasticcalc=2 - stochastic/fuzzy computation, data are read sequentially
  ///  stochasticcalc=3 - stochastic/fuzzy computation, data are generated in the code
  long stochasticcalc;
  
  ///  detection of hanging nodes
  ///  hangnodesdetect = 0 - no detection
  ///  hangnodesdetect = 1 - detection is performed
  answertype hangnodesdetect;
  sel barlist;
  
  
  ///  indicator of eigenstrains
  long eigstrains;
  
  ///  indicator of temperature
  ///  temperature=0 - problem without influence of temperature
  ///  temperature=1 - problem with influence of defined temperature
  ///  temperature=2 - problem with influence of computed temperature
  long temperature;
  
  /**  indicator of pore pressures
       pore_press=0 - problem without influence of pore pressures
       pore_press=1 - problem with influence of defined pore pressure, partially coupled approach, set by MEFEL
       pore_press=2 - problem with influence of computed pore pressure, partially coupled approach, set by METR
       pore_press=3 - problem with influence of computed pore pressure, fully coupled approach, Set by METR 
  */
  long pore_press;

  ///  type of layered problem
  long tlm;

  ///  number of simulation
  long ns;
  
  ///  definition of homogenization
  ///  homog = 0 - homogenization is switched off
  ///  homog = 1 - homogenization is switched on
  ///  homog = 2 - homogenization is based on Wang's tiles
  ///  homog = 3 - homogenization is based on stress approach
  ///  homog = 4 - homogenization is based on strain approach
  ///  homog = 5 - homogenization is based on stress approach of trabicular structure = prescribed forces on RVE boundary
  ///  homog = 6 - homogenization is based on strain approach of trabicular structure = prescribed displacements on RVE boundary
  long homog;

  ///  stress-strain state for homogenization problem (option homog=3 or homog=4)
  strastrestate hstrastre;
  ///  direction of unit load for homogenization for homog=5 and homog=6
  long homogdir;


  ///  type of macro-micro problem correspondence
  ///  mami=1 - elements are connected with microproblems
  ///  mami=2 - aggregates of elements are connected with microproblems, each element from the aggregate sends its own data
  ///  mami=3 - aggregates of elements are connected with microproblems, values from elements are averaged and one packet of data is sent to microproblem
  macromicrotype mami;
 
  ///  the number of Wang's tiles (the number of
  ///  different tiles, e.g. 8)
  long ntiletypes;

  
  ///  number of particles in one cell describing molecular structure
  ///  it has similar meaning as number nodes on one finite element
  long mcnne;
  
  ///  dimension of solved problem
  long dim;
  
  ///  estimated number of rigid body motions
  long ense;
  
  ///  number of iterations in method of variable stiffness
  long nivsm;
  
  
  //  number of increments in floating subdomain problems
  long nincr;

  // Decomposed filename of the output(input) filename
  char *path;      /// path of the output(input) filename
  char *filename;  /// filename of the output(input) filename
  char *suffix;    /// suffix of the output(input) filename
    
  //  NEWTON-RAPHSON METHOD
  ///  maximum number of iterations in inner loop
  ///  maximum number of iterations in one increment
  long niilnr;
  ///  required norm of vector of unbalanced forces
  double errnr;

  // CONSTRUCTION PHASE SOLVER
  // These parameters controls smooth transition between two contruction phases
  // i.e. element group addition or removal. In the case of cpsmooth=yes=1, the traction 
  // forces (due to element addition/removal) on interface nodes between new and old parts 
  // will be applied gradually in load steps defined by cpincrnr and cpminincrnr

  /** indicator of smoothed transition between contruction phases
      cpsmooth=yes=1 - gradual application of traction on interfaces between new and old structure part (suitable for complex soil model)
      cpsmooth=no=0 - direct application of traction on interfaces between new and old structure part */
  answertype cpsmooth;
  /** default load step length for smoothed application of traction on interfaces between new and old structre part 
      It is equvalent to the initial value of the load coefficient increment in the Newton-Raphson procedure. */
  double cpincrnr;
  /** minimum load step length for smoothed application of traction on interfaces between new and old structre part 
      It is equvalent to the required minimum value of the load coefficient increment in the Newton-Raphson procedure. */
  double cpminincrnr;

  /**  type of fully coupled solver
       nrsolv = 1 - full newton method (matrices are recomputed in all internal loops)
       nrsolv = 2 - modified Newton method (matrices are recomputed only in new increment) */
  nonlintimesolvertype nrsolv;

  /** flag fro computing of initial displacements of new structure parts used in growing problems
      comp_inidispl = 0 (no - default value) - initial displacements are taken only from interface nodes
      comp_inidispl = 1 (yes) - initial displacements due to attained displacements on interface nodes are calculated */
  answertype comp_inidispl;

  /** prescribed rotation of selected nodes at the initial displacement computation
      in case of growing problem
      rot_inidispl = 0 (no - default value) - initial displacements are calculated according to comp_inidispl setup or 
                       prescribed initial displacements are taken from the initial condition section 
      rot_inidispl = 1 (yes) - prescribed initial displacements according to the rotation of selected nodes initial displacements are taken into account */
  answertype rot_inidispl;

  ///  minimum time increment
  double dtmin;
 

  ///  type of time printing
  timetypeprin tpr;
  
  /**  type of algorithm for load balancing method
       lbtype=1 - program stops after negm is known
       lbtype=0 - program goes on after negm is known */
  long lbtype;

  /// backup controler
  hdbcontr hdbcont;

  /** array used in connection with important times of the time controller
      hdbtime[i] = 0 (no)  - the backup will not be performed in the i-th important time
      hdbtime[i] = 1 (yes) - the backup will performed in the i-th important time */
  answertype *hdbtime;

  /** Array used in the connection with important times of the time controller
      outdrv_time[i] = 0 (no)  - the output of the outdriver will not be performed in the i-th important time
      outdrv_time[i] = 1 (yes) - the output of the oputdriver will performed in the i-th important time
      The array is taken into account only in the case that the time step selection is set to sel_impvalues in the outdriver
      otherwise is ignored. */
  answertype *outdrv_time;
};

#endif
