#include "hvisolver.h"
#include "global.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "gmatrix.h"
#include "mechprint.h"
#include <stdlib.h>


/**
   function solves problems with hemivariational inequalities
   
   JK, 22.10.2007
*/
void solve_hemivariational_inequalities ()
{
  long i,j,lcid,nbdof,nidof;
  double *lhs,*rhs,*condmat,*condvec;
  
  //  load case must be equal to zero in this type of problems
  lcid=0;
  
  //  stiffness matrix assembling
  stiffness_matrix (lcid);
  
  //  arrays are allocated in Lsrs
  lhs = Lsrs->give_lhs (lcid);
  rhs = Lsrs->give_rhs (lcid);
  
  //  right hand side of system
  mefel_right_hand_side (lcid,rhs);
  
  //  solution of equation system
  //Mp->ssle->solve_system (Gtm,Smat,lhs,rhs,Out);
  
  //  number of boundary/interface unknowns
  nbdof=Gtm->nbdof;
  //  number of internal unknowns
  nidof=Gtm->nidof;
  
  condmat = new double [nbdof*nbdof];
  condvec = new double [nbdof];
  
  //Smat->condense (Gtm,condmat,lhs,rhs,nbdof,1,Out);
  
  j=0;
  for (i=nidof;i<Ndofm;i++){
    condvec[j]=rhs[i];  j++;
  }
  




  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();
  
  delete [] condvec;
  delete [] condmat;
}


