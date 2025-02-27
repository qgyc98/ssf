#include "llssolver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"

/**
   function solves layered linear static problems
   
   function contains the following steps:
   1. stiffness matrix assembling
   2. right hand side assembling (load vectors assembling)
   3. solution of equation systems by choosen method
   4. strain computation (eligible)
   5. stress computation (eligible)
   6. reaction computation (eligible)
   7. results printing
   
   JK
*/
void solve_layered_linear_statics ()
{
  long i;
  double *rhs;
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //Smat->printmat (Out);

  rhs = Lsrs->give_rhs (0);
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    Mb->lc[i].assemble (i,rhs+i*Ndofm,NULL,1.0);
  }
  

  fprintf (Out,"\n\n right hand side before solution");
  for (i=0;i<Ndofm;i++){
    fprintf (Out,"\n %4ld  %20.10f",i,Lsrs->rhs[i]);
  }


  //  solution of system of linear equations
  for (i=0;i<Lsrs->nlc;i++){
    //Smat->solve_system (Gtm,Lsrs->give_lhs(i),Lsrs->give_rhs(i));
    Mp->ssle->solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
  }
  
  Smat->printmat(Out);

  fprintf (Out,"\n\n right hand side after solution");
  for (i=0;i<Ndofm;i++){
    fprintf (Out,"\n %4ld  %20.10f",i,Lsrs->rhs[i]);
  }

  //Smat->printdiag (Out);
  
  fprintf (Out,"\n\n solution of the system");
  for (i=0;i<Ndofm;i++){
    fprintf (Out,"\n %4ld  %20.10f",i,Lsrs->lhs[i]);
  }


  /*
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    
    //  strains
    if (Mp->straincomp==1){
      Mm->computestrains (i);
      Mm->stra.transformvalues (1);
    }
    
    //  stresses
    if (Mp->stresscomp==1){
      Mm->computestresses (i);
      Mm->stre.transformvalues(0);
    }
    
    //  reactions
    if (Mp->reactcomp==1){
      Mb->lc[i].compute_reactions (i);
    }
    
    //  output data
    Lsrs->output (Out,i);
    print_multipliers (Out);
  }
  */
  
  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();
  
}

