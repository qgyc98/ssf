#ifndef NONLINMAN_H
#define NONLINMAN_H

#include <stdio.h>
#include "alias.h"
#include "galias.h"
#include "iotools.h"

class gmatrix;

/**
  Class nonlinman defines type of solver of systems of nonlinear algebraic equations.
   
  Created by JK, 16.6.2004
*/

class nonlinman
{
 public:
  
  nonlinman (void);
  ~nonlinman (void);
  void read (XFILE *in,long mespr);
  void print (FILE *out);

  void read_arclength(XFILE *in, nonlinsolvertype nlst, long mespr);
  void read_newtonraphson(XFILE *in, nonlinsolvertype nlst, long mespr);
  void read_dissipincr(XFILE *in, long mespr);
  
  void read_displnorm(XFILE *in, displacementnorm dn);
  void read_seldof (XFILE *in);
  void read_seldofcoord (XFILE *in);
  void read_selmstr (XFILE *in);
  void read_selnod (XFILE *in);
  void read_selnoddistincr (XFILE *in);

  void print_arclength(FILE *out, nonlinsolvertype nlst);
  void print_newtonraphson(FILE *out, nonlinsolvertype nlst);
  void print_dissipincr(FILE *out);
  
  void print_displnorm(FILE *out, displacementnorm dn);
  void print_seldof (FILE *out);
  void print_seldofcoord (FILE *out);
  void print_selmstr(FILE *out);
  void print_selnod (FILE *out);
  void print_selnoddistincr (FILE *out);
  
  void allocavec(long n);
  void deallocavec();
  void initiate (nonlinman &nm);
  
  
  ///  type of solver of nonlinear algebraic equation system
  nonlinsolvertype tnlinsol;

  /**  TYPE OF STIFFNESS MATRIX
       stmat=1 - initial elastic stiffness
       stmat=2 - tangent stiffness
       stmat=3 - secant stiffness %matrix
       stmat=5 - incr_tangent_stiff - the stiffness %matrix is recalculated after one increment of the outer loop
       stmat=6 - ijth_tangent_stiff - the stiffness %matrix is recalculated after i-th step of the outer loop and after j-th step in the inner loop */
  stiffmatrix stmat;
  ///  the stiffness %matrix is recalculated after ithstiffchange steps in outer loop
  long ithstiffchange;
  ///  the stiffness %matrix is recalculated after jthstiffchange steps in inner loop
  long jthstiffchange;
  
  //  ARC-LENGTH METHOD
  ///  norm measure of displacement increments
  displacementnorm displnorm;
  ///  determination of increment of load parameter
  detlambda dlam;
  ///  hard-disc backup
  long hdbackupal;
  ///  maximum number of increments
  long nial;
  /**  maximum number of iterations in inner loop
       maximum number of iterations in one increment */
  long niilal;
  ///  required norm of vector of unbalanced forces
  double erral;
  ///  length of arc
  double dlal;
  ///  minimum length of arc
  double dlminal;
  ///  maximum length of arc
  double dlmaxal;
  ///  displacement-loading driving switch
  double psial;
  ///  number of selected nodes (for norm computation)
  long nsnal;
  ///  selected nodes (for norm computations)
  long *selnodal;
  ///  number of selected DOFs (for norm computation)
  long nsdofal;
  ///  selected degrees of freedom
  long *seldofal;
  ///  coordinates of selected nodes
  double *selnodcoord;
  ///  problem dimension (for increment of two node distance)
  long probdimal;
  ///  unit vector between selected nodes
  double nxal,nyal,nzal;
  /// the number of selected macro-strain/stress components, it is dimension of mstrastre and mstrid arrays
  long nmstrcomp;
  /// array of macro-strain/stress component types
  strastre *mstrastre;
  /// array of id of selected macro-strain/stress components
  long *mstrid;

  /*
    Setup of nonliniear solver controlled by dissipation increment
  */
  /// initial value of dissipation increment
  double tau_ini;
  /// limit dissipation increment where the switching to the common arclength will be performed
  double tau_lim;
  /// maximum allowable dissipation increment
  double tau_max;
  /// actual value of the dissipation increment
  double tau;
  ///  type of alternative solver of nonlinear algebraic equation system for solver controlled by the dissipation increment
  nonlinsolvertype altnlst;
  
  
  /*  the stiffness %matrix can be modified in outer as well as in the inner loops
      for this purpose, two variables (ithstiffchange, jthstiffchange) are defined
      
  */

  /// flag for checking of maximum total(cumulative) performed length of the arc
  flagsw check_tot_al;
  /// maximum total(cumulative) performed length of the arc
  double max_tot_al;

  /// flag for checking of lambdar
  flagsw check_lambdar;
  /// required load coefficient for the proportional %vector
  double lambdar;

  /// required error for lambdar tresholds
  double errl;

  /// flag for detection of divergency of inner iteration loop
  flagsw check_div;
  /** the minimum number of performed inner iteration steps 
      before the loop can be interrupted due to divergency 
      (div_min_steps has to be greater or equal to 3)*/
  long div_min_steps;
  /// array of stored step id (length of array = div_min_steps)
  double *divc_step;
  /// array of stored errors (length of array = div_min_steps)
  double *divc_err; 

  //  NEWTON-RAPHSON METHOD
  ///  maximum number of increments
  long ninr;
  /**  maximum number of iterations in inner loop
       maximum number of iterations in one increment */
  long niilnr;
  ///  required norm of vector of unbalanced forces
  double errnr;
  ///  magnitude of increment of loading
  double incrnr;
  ///  minimum magnitude of increment of loading
  double minincrnr;
  ///  maximum magnitude of increment of loading
  double maxincrnr;
  /// type of computation of residual vector norm (rel_load_norm=1, rel_react_norm=2, rel_loadreact_norm=3, absol_norm=4)
  resnormt rnormtnr;

  // the following data are used for Hinkley analysis only
  /// number of steps in the initial loop for equilibrium
  long nienr;
  /// harddisk backup restart indicator
  long hdbr;
  /// index of harddisk backup in the harddisk backup file
  long hdbid;
  /// backup filename
  char backupfname[1025];
};


#endif
