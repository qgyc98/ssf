#include "mdsolver.h"
#include "global.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"
#include "node.h"
#include "lhsrhs.h"
#include "mechmat.h"
#include <string.h>

/**
   function solves problems of molecular dynamics
   
   JK, 19.6.2005
*/
void solve_molecular_dynamics ()
{






















  long i;
  double *lhs,*rhs;
  
  //  stiffness matrix assembling
  stiffness_matrix (0);
  
  //  array is allocated in Lsrs
  lhs = Lsrs->give_lhs (0);
  rhs = Lsrs->give_rhs (0);
  
  //  loop over load cases
  for (i=0;i<Lsrs->nlc;i++){
    //  right hand side of system
    mefel_right_hand_side (i,rhs);
  }
  
  /*
  Smat->printmat(Out);
  FILE *p;
  p = fopen ("vektor.txt","w");
  fprintf (p,"%ld\n",Ndofm);
  for (i=0;i<Ndofm;i++){
    fprintf (p,"%le\n",rhs[i]);
  }
  fclose (p);
  */

  //  solution of equation system

  for (i=0;i<Lsrs->nlc;i++){
    Smat->solve_system (Lsrs->give_lhs(i),Lsrs->give_rhs(i));
  }
  
  //  zacatek Risoviny
  /*
    long j;
  long n=6;
  double *condmat,*condvect,*condsolv;
  condsolv = new double [n];
  condvect = new double [n];
  condmat = new double [n*n];
  nullvr (condmat,n*n);
  nullvr (condvect,n);
  nullvr (condsolv,n);

  //Smat->condense (Gtm,condmat,lhs,rhs,3,1);
  Smat->condense (Gtm,condmat,lhs,rhs,n,1);

  fprintf (Out,"\n\n kondenzovany vektor");
  j=0;
  for (i=Ndofm-n;i<Ndofm;i++){
    condvect[j]=rhs[i];
    fprintf (Out,"\n %e",condvect[j]);
    j++;
  }
  
  fprintf (Out,"\n\n kondenzovana matice");
  for (i=0;i<n;i++){
    fprintf (Out,"\n");
    for (j=0;j<n;j++){
      fprintf (Out," %e",condmat[i*n+j]);
    }
  }

  gemp (condmat,condsolv,condvect,n,1,1.0e-15,1);
  
  fprintf (Out,"\n\n vyresena cast vysledku");
  j=0;
  for (i=Ndofm-n;i<Ndofm;i++){
    lhs[i]=condsolv[j];
    fprintf (Out,"\n %e",condsolv[j]);
    j++;
  }
  
  //  back substitution on subdomains
  //Smat->condense (Gtm,condmat,lhs,rhs,3,2);
  Smat->condense (Gtm,condmat,lhs,rhs,n,2);
  
  delete [] condmat;
  delete [] condvect;
  delete [] condsolv;
  */
  //  konec Risoviny
  

  /*
  //  zacatek JK testu
  long nse=6,dim;
  long *se;
  double *rbm;
  
  se = new long [nse];
  rbm = new double [nse*Ndofm];
  
  Smat->kernel (rbm,dim,se,nse,1.0e-3,3);
  
  fprintf (Out,"\n\n\n RIGID BODY MOTIONS \n");
  fprintf (Out,"\n pocet RBM  %ld\n",dim);
  for (i=0;i<dim;i++){
    fprintf (Out," %ld",se[i]);
  }
  for (i=0;i<dim;i++){
    fprintf (Out,"\n\n RBM number %ld",i);
    for (j=0;j<Ndofm;j++){
      fprintf (Out,"\n %4ld %4ld  %e",i,j,rbm[i*Ndofm+j]);
    }
  }
  //  konec JK testu
  */
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
    

  }
*/
  
  print_init(-1, "wt");    
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    Mm->allipstrains (i);
    Mm->allipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();

  
  /*
  for (i=0;i<Ndofm;i++){
    printf ("\n posun %20.10lf",Lsrs->lhs[i]);
  }
  */

}

