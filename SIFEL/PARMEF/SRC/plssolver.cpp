#include "mpi.h"
#include "plssolver.h"
#include "pglobal.h"
#include "seqfilesm.h"
#include "genfile.h"
#include "vector.h"
#include <string.h>

void par_solve_linear_statics ()
{
  long i;
  double *lhs,*rhs;
  
  //  actual load case
  i=0;
  
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);
  

  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  Psolm->computation_stat (Gtm,Smat);
  
  rhs = Lsrs->give_rhs (0);
  
  
  
  Mb->lc[i].assemble (i,rhs+i*Ndofm,NULL,1.0);

  //fprintf (Out,"\n\n kontrola RHS ");
  //for (i=0;i<Ndofm;i++){
  //fprintf (Out,"\n rhs %ld    %le",i,rhs[i]);
  //}
  
  //for(i = 0; i < Ndofm; i++){
  //rhs[i]/=double(Psolm->ptop->dofmultip[i]);
  //}
  
  
  
  /*
  double norm=0.0;
  for (i=0;i<Ndofm;i++){
    norm+=rhs[i]*rhs[i];
  }
  printf ("\n norm of the RHS %lf",norm);
  */
  
  if (Psolm->tdd==layered_plate){
    double *th;
    
    //  array containing thicknesses
    th = new double [Mt->nn];
    for (i=0;i<Mt->nn;i++){
      th[i]=Mc->give_onethickness (Mt->nodes[i].crst,Mt->nodes[i].idcs);
    }
    
    //  constraint matrix assembling
    Psolm->constr_mat (th,Out);
    
    delete [] th;
  }
  
  
  fflush(stdout);
  fflush(stderr);

  /*
  if (Myrank==5){
    long i,j;
    FILE *p;
    
    p = fopen ("matice.txt","w");

    fprintf (p,"%ld  %ld\n",Smat->cr->n,Smat->cr->negm);
    for (i=0;i<Smat->cr->n;i++){
      for (j=Smat->cr->adr[i];j<Smat->cr->adr[i+1];j++){
	fprintf (p,"%8ld %8ld  % 20.15le\n",i,Smat->cr->ci[j],Smat->cr->a[j]);
      }
    }
    fclose (p);
  }
  */

  //fprintf (Out,"\n\n kontrola matice podoblasti \n");
  //Smat->printmat (Out);
  //fprintf (Out,"\n\n konec kontroly matice podoblasti \n");


  //  solution of equation system
  Psolm->par_linear_solver (Gtm,Smat,lhs,rhs,Out,Mespr);


  // kontrola
  
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
  }
  
  */
  
  //  output data
  //Lsrs->output (Out,0);


  print_init(-1, "wt",Pmp->fni,Pmp->fei);
  print_flush ();
  for (i=0;i<Lsrs->nlc;i++){
    print_step(i, 0, 0.0, NULL);    
    print_flush ();
  }
  print_close();


  if (Mp->adaptivityflag)
    Ada->run (2,0,0); // 2 - for testing

  Psolm->computation_stat_print (Out);

}

