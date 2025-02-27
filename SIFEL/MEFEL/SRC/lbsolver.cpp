#include "lbsolver.h"
#include "global.h"
#include "probdesc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"
#include "node.h"
#include "elemswitch.h"
#include <string.h>

/**
   function solves problems of load balancing
   
   function contains the following steps:
   1. stiffness matrix assembling
   2. right hand side assembling (load vectors assembling)
   3. solution of equation systems by choosen method
   4. strain computation (eligible)
   5. stress computation (eligible)
   6. reaction computation (eligible)
   7. results printing
   
   JK, 22.11.2007
*/
void solve_load_balancing ()
{
  long i,nrdof;
  double *lhs,*rhs,*condmat,*condvect;
  
  double no;
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  

  //Smat->printmat (Out);
  //abort ();



  //  array is allocated in Lsrs
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }
  
  nrdof=Gtm->nbdof;
  condmat = new double [nrdof*nrdof];
  condvect = new double [nrdof];
  if (Mespr==1)
    fprintf (stdout,"\n number of interface unknowns  %ld",nrdof);
  
  //  solution of equation system
  for (i=0;i<Lsrs->nlc;i++){
    Smat->condense (Gtm,condmat,condvect,lhs,rhs,nrdof,1,Out);
    //no = Smat->sky->ldlkoncount_sky (condmat,condvect,lhs,rhs,nrdof,1);
  }
  
  Smat->printmat (Out);
  abort ();

  fprintf (stdout,"\n");
  fprintf (stdout,"\n number of operations %le",no);
  fprintf (stdout,"\n");
  
  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();
  
}

