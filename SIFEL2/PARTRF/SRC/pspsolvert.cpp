#include "mpi.h"
#include "pspsolvert.h"
#include "pglobalt.h"
#include "seqfilest.h"
#include <string.h>
#include <math.h>

void par_solve_stationary_problem ()
{
  long i;
  double *lhs,*rhs;
  
  //  actual load case
  i=0;
  
  //  conductivity matrix assembling
  conductivity_matrix (0);
  
  rhs = Lsrst->give_rhs (0);
  lhs = Lsrst->give_lhs (0);
  
  //Tb->lc[i].assemble (i,rhs+i*Ndoft,Ndoft);
  //Tb->lc[i].assemble (i,rhs+i*Ndoft);
  trfel_right_hand_side (0,rhs,Ndoft);
  

  //  solution of equation system
  Psolt->par_linear_solver (Gtt,Kmat,lhs,rhs,Outt,Mesprt);
  
  print_initt(-1, "wt");    
  for (i=0;i<Lsrst->nlc;i++){
    //  computes and prints required quantities
    print_stept(i, 0, 0.0, NULL);    
  }
  print_closet();
  
}

