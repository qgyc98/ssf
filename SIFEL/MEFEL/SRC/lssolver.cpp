#include "lssolver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "lhsrhs.h"
#include "outdriverm.h"
#include "adaptivity.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"
#include "node.h"
#include "elemswitch.h"
#include <string.h>

/**
   function solves linear statics problems
   
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
void solve_linear_statics ()
{
  long i,n,aux;
  double *rhs,*fl,*fi;
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  number of mechanical degrees of freedom
  n=Ndofm;
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  fl   = new double [n];
  nullv (fl,n);
  //  vector of resulting internal forces
  fi   = new double [n];
  nullv (fi,n);

  rhs = Lsrs->give_rhs (0);

  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }

  /*
  FILE *ma;
  ma=fopen ("prava_strana.txt","w");
  fprintf (ma,"%ld\n",n);
  //Smat->printmat (ma);
  for (i=0;i<n;i++){
    fprintf (ma,"% 16.12le\n",rhs[i]);
  }
  fclose (ma);
  fprintf (stdout,"\n\n\n RHS %le\n\n",scprd(rhs,rhs,n));
  */

  /* fprintf (Out,"\n");
     fprintf (Out,"\n");
     for (i=0;i<n;i++){
     fprintf (Out,"rhs[%ld] %16.12le\n",i,rhs[i]);
     }
  */

  //Smat->printmat (Out);

  //  solution of equation system
  for (i=0;i<Lsrs->nlc;i++){
    Mp->ssle->solve_system (Gtm,Smat,Lsrs->give_lhs(i),Lsrs->give_rhs(i),Out);
  }
  

  /* fprintf (Out,"\n");
     fprintf (Out,"\n");
     for (i=0;i<n;i++){
     fprintf (Out,"lhs[%ld] %16.12le\n",i,Lsrs->lhs[i]);
     }
  */
  /*
  ma=fopen ("reseni.txt","w");
  fprintf (ma,"pocet slozek %ld\n",n);
  for (i=0;i<n;i++){
    fprintf (ma,"% 16.12le\n",Lsrs->lhs[i]);
  }
  fclose (ma);
  */
  
  print_init(-1, "wt");    
  // clean temperature strains because they are stored in one array used for all load cases
  Mm->nulltempstrains();
  rhs = Lsrs->give_rhs(0);
  for (i=0;i<Lsrs->nlc;i++){
    //  indicator of strain computation
    //  it means, strains at integration points have not been computed
    Mp->strainstate=0;
    //  indicator of stress computation
    //  it means, stresses at integration points have not been computed
    Mp->stressstate=0;
    //  indicator of computation of other array
    //  it means, stresses at integration points have not been computed
    Mp->otherstate=0;
    // rhs was rewritten by the solver call -> rebuild load vector
    aux = Mp->homog;
    Mp->homog = 0;
    mefel_right_hand_side(i,rhs,fl);
    Mp->homog = aux;
    //  computes required quantities
    compute_req_val (i); 
    // in artificial time 0.0 the nodal forces contains only the load vector fl
    print_step(i, 0, 0.0, fl);

    if (Mp->homog == 9){
      long j, ncomp = Mt->max_ncompstr;
      vector tifor(ASTCKVEC(ncomp)), ifor(ASTCKVEC(ncomp));

      for(j=0; j<Mm->tnip; j++)
        Mm->computenlstresses(j,Mm->ip[j]);
      for (j=0; j<Mt->ne; j++){
        nullv(ifor);
        elem_volintegration_quant(j, locstress, i, ifor);
        addv(tifor, ifor, tifor);
      }
      fprintf(stdout, "\nMacro-stress components:\n");
      printv(stdout, tifor);
    }

    // if the nodal forces are required for graphical output
    // reaction + load vector are printed in the artificial time 1.0
    if (Outdm->nog.selnforce.st != sel_no)
    {
      //  indicator of strain computation
      //  it means, strains at integration points have not been computed
      Mp->strainstate=0;
      //  indicator of stress computation
      //  it means, stresses at integration points have not been computed
      Mp->stressstate=0;
      //  indicator of computation of other array
      //  it means, stresses at integration points have not been computed
      Mp->otherstate=0;
      // calculate nodal forces and reactions
      internal_forces(i, fi);
      //  computes required quantities
      compute_req_val (i);
      //  prints required quantities    
      // in artificial time 1.0 the nodal forces contains the load vector + reactions
      Outdm->print_graphics(Outdm->outgr, i, 1.0, 1, fi);
    }
    // clean temperature strains because they are stored in one array used for all load cases
    Mm->nulltempstrains();
  }
  print_close();

  
  if (Mp->adaptivityflag)
    Ada->run (2,2,0);
  
  delete [] fl;  
  delete [] fi;

}

