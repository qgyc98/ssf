#ifndef PROBDESCT_H
#define PROBDESCT_H

#include <stdio.h>
#include "aliast.h"
#include "galias.h"
#include "gfunct.h"
#include "timecontr.h"
#include "hdbcontr.h"
#include "slesolv.h"
#include "nonlinmant.h"
#include "iotools.h"
#include "xfile.h"

/**
   class probdesc
   
   it is one of the 5 most important classes of the program
   (probdesc, transtop, transmat, transbclc, transcrsec)
   
   class contains general information about solved problem
*/

class probdesct
{
 public:
  probdesct (void);
  ~probdesct (void);
  void read (XFILE *in);
  void print (FILE *out);
  long give_var_dofid(namevart var);
  
  

  char *path;
  char *filename;
  char *suffix;
  
  /// keyword processing option
  pkwd_sw kwdsw;

  ///  problem name
  char name[1001];

  ///  type of problem
  problemtypet tprob;
  
  ///  type of transported matter
  transmattert tmatt;
  ///  number of transported matters=number of row(column) blocks in matrices
  long ntm;

  /// total number of known(=implemented) primary variables(= nodal unknowns)
  static const long tnkv = sizeof(namevartstr)/sizeof(*namevartstr);

  /** array of dof indices for individual known(=implemented) primary variables(= nodal unknowns)
      var_dofid[int(namevart)] = <0;ntm-1> <=> dof id at node for the given name(=namevart) of nodal unknown 
                               = -1 if the given name(=namevart) is not defined in the problem*/
  long var_dofid[tnkv];

  /**  array with ordered names of used primary variables(= nodal unknowns)
       dofname[i] = name of i-th nodal unknown*/
  namevart *dofname;

  ///  names of transported media
  mednamest mednam;

  /// DOF vs. media name map
  char **dof_mednam;
  
  /// string with media names used separated by \0
  char *dof_mednam_str;

  /// array of scales of transported media
  double *scale;
  
  ///  type of solver of nonstationary problems
  nonstatsolver tnpsolver;

  ///  data about solver of system of nonlinear equations
  nonlinmant nlman;

  ///  gravity acceleration is taken into account
  gravityaccelerationt tgravity;
  double gr[3];

  ///  geometrical dimension
  ///  it is determined in elementt.cpp
  long gdim;
  ///  savemode - deals with the number of integration points
  ///  savemode=0 - greater number of integration points
  ///  savemode=1 - smaller number of integration points
  long savemode;
  
  ///  type of transport solver
  ///  trsolv = 1 - full newton method (matrices are recomputed in all internal loops)
  ///  trsolv = 2 - modified Newton method (matrices are recomputed only in new increment)
  transpsolver trsolv;

  ///  type of residuum computation
  ///  trestype = 1 - residuum is computed from internal fluxes
  ///  trestype = 2 - residuum is computed as difference between right and left hand side
  transpresiduumtype trestype;
  
  ///  convergence control for full newton-raphson
  ///  yes = 1
  ///  no  = 0 
  answertype convergcontrolt;

  //  STORAGE TYPE
  ///  type of storage of conductivity matrix
  storagetype tstorkm;
  ///  type of storage of capacity matrix
  storagetype tstorcm;

  /// diagonalization of capacity matrix
  long diagcap;
  /// diagonalization of reaction matrix
  long diagreact;
  long react_capac;

  ///  data about solver of system of linear equations
  slesolv *ssle;

  ///  time controller
  timecontr timecont;
  ///  time
  double time;
  /// index of actual time step
  long istep;
  /// index of actual iteration step (Newton-Raphson)
  long jstep;
  
  ///  threshold for rejection from compressed storages
  double limit;
  ///  computer zero
  double zero;
  ///  threshold for the size of the right hand side %vector
  double *threshrhs;
  
  ///  definition of homogenization
  ///  homogt = 0 - homogenization is switched off
  ///  homogt = 1 - homogenization is switched on, homogenization is computed on a single processor
  ///  homogt = 2 - homogenization is performed on a parallel computer
  long homogt;
  
  ///  type of macro-micro problem correspondence
  ///  mami=1 - elements are connected with microproblems
  ///  mami=2 - aggregates of elements are connected with microproblems, each element from the aggregate sends its own data
  ///  mami=3 - aggregates of elements are connected with microproblems, values from elements are averaged and one packet of data is sent to microproblem
  macromicrotype mami;
  
  ///  advection
  ///  advect=0 - no advection
  ///  advect=1 - there is advection described by velocity vector
  long advect;
    
  ///  reaction
  ///  react=0 - no reaction term
  ///  react=1 - there is reaction term
  long react;
    
  /***************************/
  /*  non-stationary solver  */
  /***************************/

  ///  coefficient alpha
  ///  alpha=0.0 - explicit algorithm
  ///  alpha=1/2 - Crank-Nicholson scheme
  ///  alpha=1.0 - implicit algorithm
  double alpha;
    
  ///  required error in Newton-Raphson method
  double err;
  double *errarr;
  ///  maximum number of inner iterations in Newton-Rapson method
  long nii;

  ///  stochastic calculation indicator
  long stochasticcalc;
  
  /// status flag of approximation of nodal values to integration points
  /// ipvcomp=0 - approximated nodal values in Tm->nv have not been actualized in the given step and must be computed
  /// ipvcomp=1 - approximated nodal values in Tm->nv have already been actualized in the given step and no computation is needed
  long ipvcomp;
  
  ///  compute gradients
  ///  gradcomp=0 - gradients are not computed
  ///  gradcomp=1 - gradients are computed
  long gradcomp;
  ///  compute fluxes
  ///  fluxcomp=0 - fluxes are not computed
  ///  fluxcomp=1 - fluxes are computed
  long fluxcomp;
  ///  compute other values
  ///  othercomp=0 - values from array other are not computed
  ///  othercomp=1 - values from array other are computed
  long othercomp;
  ///  compute eq_other values 
  ///  eqothercomp=0 - values from array eqother are not computed
  ///  eqothercomp=1 - values from array eqother are computed
  long eqothercomp;
  
  ///  gradient averaging
  ///  gradaver=0 - components of other array are not averaged
  ///  gradaver=1 - components of other array are averaged, arithmetic average is used
  ///  gradaver=2 - components of other array are averaged, volume average is used
  long gradaver;
  ///  flux averaging
  ///  fluxaver=0 - components of other array are not averaged
  ///  fluxaver=1 - components of other array are averaged, arithmetic average is used
  ///  fluxaver=2 - components of other array are averaged, volume average is used
  long fluxaver;
  ///  other values averaging
  ///  otheraver=0 - components of other array are not averaged
  ///  otheraver=1 - components of other array are averaged, arithmetic average is used
  ///  otheraver=2 - components of other array are averaged, volume average is used
  long otheraver;
  ///  eqother values averaging
  ///  eqotheraver=0 - components of other array are not averaged
  ///  eqotheraver=1 - components of other array are averaged, arithmetic average is used
  ///  eqotheraver=2 - components of other array are averaged, volume average is used
  long eqotheraver;
  
  ///  position where gradients are computed
  ///  gradpos = 1 - integration points
  ///  gradpos = 2 - nodes, values are copied from the closest integration point
  ///  gradpos = 3 - nodes, values are computed at nodes
  long gradpos;
  
  ///  position where fluxes are computed
  ///  fluxpos = 1 - integration points
  ///  fluxpos = 2 - nodes, values are copied from the closest integration point
  ///  fluxpos = 3 - nodes, values are computed at nodes
  long fluxpos;
  
  ///  position where other are computed
  ///  otherpos = 1 - integration points
  ///  otherpos = 2 - nodes, values are copied from the closest integration point
  ///  otherpos = 3 - nodes, values are computed at nodes
  long otherpos;
  
  ///  position where eqother are computed
  ///  eqotherpos = 1 - integration points
  ///  eqotherpos = 2 - nodes, values are copied from the closest integration point
  ///  eqotherpos = 3 - nodes, values are computed at nodes
  long eqotherpos;
  
  ///  adaptivity calculation indicator
  ///  adaptivity == 0 - no adaptivity calculation
  ///  adaptivity == 1 - adaptivity calculation - SPR type
  ///  adaptivity == 2 - adaptivity calculation - Z2 type
  long adaptivityflag;
  
  ///  storage of actual nodal values on nodes
  ///  nvs=0 - actual nodal values are not stored
  ///  nvs=1 - actual nodal values are stored
  long nvs;
  ///  storage of initial nodal values on nodes
  ///  invs=0 - initial nodal values are not stored
  ///  invs=1 - initial nodal values are stored
  long invs;
  ///  storage of nodal values from previous time step on nodes
  ///  pnvs=0 - previous nodal values are not stored
  ///  pnvs=1 - previous nodal values are stored
  long pnvs;
  ///  storage of time derivative of  nodal values on nodes
  ///  tdnvs=0 - time derivative of nodal values are not stored
  ///  tdnvs=1 - time derivative of nodal values are stored
  long tdnvs;
  
  ///  type of time printing
  timetypeprint tprt;
  

  ///  number of simulation
  long ns;

  /// backup controler
  hdbcontr hdbcont;
};

#endif
