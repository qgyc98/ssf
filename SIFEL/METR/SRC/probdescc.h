#ifndef PROBDESCC_H
#define PROBDESCC_H

#include <stdio.h>
#include "aliasc.h"
#include "galias.h"
#include "timecontr.h"
#include "slesolv.h"
#include "iotools.h"
#include "xfile.h"

class gmatrix;

/**
   class probdescc
   
   it is one of the 5 most important classes of the program
   (probdesc, mechtop, mechmat, mechbclc, mechcrsec)
   
   of thermo-hydro-mechanical problem
*/
class probdescc
{
 public:
  probdescc (void);
  ~probdescc (void);
  void read (XFILE *in);
  void print (FILE *out);

  
  /// keyword processing option
  pkwd_sw kwdsw;

  char *path;
  char *filename;
  char *suffix;
  
  ///  problem name
  char name[1001];
  ///  MEFEL input file name
  char minfile[1001];
  ///  TRFEL input file name
  char tinfile[1001];
  ///  output file name
  char outfile[1001];
  ///  auxiliary output file name
  char auxfile[1001];
  
  ///  type of problem
  problemtypec tprob;
  ///  type of transported matter
  transmatterc tmatt;
  ///  number of transported matters=number of row(column) blocks in matrices
  //long ntm;
  ///  names of transported media
  mednamesc mednam;
  //  STORAGE TYPE
  ///  type of storage of zero-order matrix
  storagetype tstord0;
  ///  type of storage of first-order matrix
  storagetype tstord1;
  
  ///  data about solver of system of linear equations
  slesolv *ssle;
  
  ///  gravity acceleration is taken into account
  gravityaccelerationc cgravity;
  double gr[3];

  ///  type of solver of nonlinear algebraic equation system
  nonlinsolvertypec tnlinsol;
  
  ///  Matrices cleaning for more memory
  ///  cleanmatrix = 0 - no
  ///  cleanmatrix = 1 - yes
  coupcleanmatrices cleanmatrix;

  ///  Babuska-Brezzi condition
  ///  bb=1 - linear app. in mechanics, linear app. in transport
  ///  bb=2 - quadratic app. in mechanics, linear app. in transport
  ///  bb=3 - quadratic app. in mechanics, quadratic app. in transport
  babuskabrezzi bb;
  
  ///  savemode - deals with the number of integration points
  ///  savemode=0 - greater number of integration points
  ///  savemode=1 - smaller number of integration points
  long savemode;
  
  ///  type of fully coupled solver
  ///  fcsolv = 0 - linear solver (no update of system matrices is required)
  ///  fcsolv = 1 - full newton method (matrices are recomputed in all internal loops)
  ///  fcsolv = 2 - modified Newton method (matrices are recomputed only in new increment)
  coupsolver fcsolv;
  
  ///  type of residuum computation
  ///  restype = 1 - residuum is computed from internal fluxes and forces
  ///  restype = 2 - residuum is computed as difference between right and left hand side
  residuumtype restype;

  ///  threshold for rejection from compressed storages
  double limit;
  ///  computer zero
  double zero;

  
  ///  number of printed unknowns in non-linear solver
  long npun;
  ///  array containing nodes and code numbers of printed values
  long *requn;

  ///  stiffness of material
  ///  stmat=0 - initial elastic stiffness
  ///  stmat=1 - tangent stiffness
  long stmat;
  
  
  
  //  NEWTON-RAPHSON METHOD
  ///  maximum number of iterations in inner loop
  ///  maximum number of iterations in one increment
  long niilnr;
  ///  required norms of residual vector of flux resultants
  double errnr;
  
  ///  coefficient in numerical solver of system of ODE
  ///  alpha=0 - forward Euler method
  ///  alpha=1 - backward Euler method
  double alpha;
  ///  actual time
  double time;
  
  //  step number in Newton-Raphson method
  long incr_step;
  long inner_step;
  
  // type of the passing data between modules (TRFEL, MEFEL) in coupled problems
  // pass_by_closest_ip=1 -> data are copied to nodes from the closest int. point, requires the same number of elements in both modules - minimum additional memory requirements, fastest way, may be inaccurate for certain combination of elements
  // pass_by_nodes_comp=2 -> data are calculated at nodes with the help of state varibales from the closest int. point, requires the first n elements in both modules to be the same in shape and same ordering of the first m nodes,
  //                         where n is the minimum number of elements in modules, m is the minimum number of nodes on elements - minimum additional memory requirements, slower than type 1, better accuracy than type 1
  // pass_by_aux_ip = 3 -> data are calculated in auxiliary int. points, the number of nodes and elements in meshes of particular modules can be independent - may leed to large memory requirements, slower than type 1, exact.
  datapasstype dpt;
  // the maximum number of itartions in the computation of auxiliary int. point coordinates
  long aip_ni;
  // the maximum distance of two points at which they may be assumed to be identical in the computation of auxiliary int. point coordinates
  double aip_err;

  ///  time controller
  timecontr timecon;

  ///  indicator of eigenstrains
  long eigstrains;

};

#endif
