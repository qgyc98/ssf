#ifndef NONLINMANT_H
#define NONLINMANT_H

#include <stdio.h>
#include "aliast.h"
#include "galias.h"
#include "iotools.h"

/**
   class nonlinmant defines type of solver of systems of nonlinear algebraic equations
   
   JK, 18.10.2007
*/

class nonlinmant
{
 public:
  
  nonlinmant (void);
  ~nonlinmant (void);
  void read (XFILE *in,long mespr);
  void print (FILE *out);
  
  /*
  void read_selnod (XFILE *in);
  void read_seldof (XFILE *in);
  void read_seldofcoord (XFILE *in);
  void print_selnod (FILE *out);
  void print_seldof (FILE *out);
  void print_seldofcoord (FILE *out);
  
  void initiate (nonlinman &nm);
  */
  
  ///  type of solver of nonlinear algebraic equation system
  nlsolvertype tnlinsol;
  
  /*
  //  ARC-LENGTH METHOD
  ///  norm measure of displacement increments
  displacementnorm displnorm;
  ///  hard-disc backup
  long hdbackupal;
  ///  maximum number of increments
  long nial;
  ///  maximum number of iterations in inner loop
  ///  maximum number of iterations in one increment
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
  */

  //  NEWTON-RAPHSON METHOD
  ///  maximum number of increments
  long ninr;
  ///  maximum number of iterations in inner loop
  ///  maximum number of iterations in one increment
  long niilnr;
  ///  required norm of vector of unbalanced forces
  double errnr;
  ///  magnitude of increment of loading
  double incrnr;
  ///  minimum magnitude of increment of loading
  double minincrnr;
  ///  maximum magnitude of increment of loading
  double maxincrnr;
  /// required value of lambda
  double rvlam;
  /// required error of reached required lambda
  double errl;
  
  /*
  /// number of steps in the initial loop for equilibrium
  long nienr;
  /// harddisk backup restart indicator
  long hdbr;
  /// index of harddisk backup in the harddisk backup file
  long hdbid;
  /// backup filename
  char backupfname[1025];
  */
};


#endif
