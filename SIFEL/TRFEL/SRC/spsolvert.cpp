#include "spsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "transprint.h"
#include <string.h>
#include "npsolvert.h"

/**
   function solves linear stationary transport problem

   JK, 20.12.2002
*/
void solve_stationary_problem ()
{
  long i;
  double *rhs;
  
  //  conductivity matrix assembling
  conductivity_matrix (0);
  
  rhs=Lsrst->give_rhs (0);
  
  //  right hand side of the system
  //for (i=0;i<Lsrst->nlc;i++){
  i = 0;
  trfel_right_hand_side (i,rhs,Ndoft);
  //}
  
  /*
  for (i=0;i<Ndoft;i++){
    fprintf (stdout,"\n rhs %20.15le",rhs[i]);
  }
  */
  /*
  long j;
  FILE *p;
  p = fopen ("matice.txt","w");
  fprintf (p,"%ld %ld\n",Ndoft,Kmat->cr->negm);
  for (i=0;i<Ndoft;i++){
    for (j=Kmat->cr->adr[i];j<Kmat->cr->adr[i+1];j++){
      fprintf (p,"%ld %ld  %20.15le\n",i,Kmat->cr->ci[j],Kmat->cr->a[j]);
    }
  }
  fclose (p);
  p = fopen ("vektor.txt","w");
  fprintf (p,"%ld\n",Ndoft);
  for (i=0;i<Ndoft;i++){
    fprintf (p,"%20.15le\n",rhs[i]);
  }
  fclose (p);

  for (i=0;i<Ndoft;i++){
    fprintf (stdout,"\n rhs %20.15le",rhs[i]);
  }
  
  
  abort ();
  */
  
  /*
  for (i=0;i<Ndoft;i++){
    Lsrst->rhs[i]=12.0;
  }
  */
  
  //  solution of equation system
  //for (i=0;i<Lsrst->nlc;i++){
  //Kmat->solve_system (Gtt,Lsrst->give_lhs(i),Lsrst->give_rhs(i));
  i = 0;
  Tp->ssle->solve_system (Gtt,Kmat,Lsrst->give_lhs(i),Lsrst->give_rhs(i),Outt);
  
  /*
  long nrdof=2;
  double *condvect,*condvect2;
  densemat condmat;
  condmat.alloc (nrdof);
  condvect = new double [nrdof];
  condvect2 = new double [nrdof];
  Kmat->condense (Gtt,condmat.a,condvect,Lsrst->give_lhs(i),Lsrst->give_rhs(i),nrdof,1,Outt);
  condmat.gemp (condvect2,condvect,1,1.0e-15,1);
  Kmat->condense (Gtt,condmat.a,condvect2,Lsrst->give_lhs(i),Lsrst->give_rhs(i),nrdof,2,Outt);
  */
  
  //}

  //for (i=0;i<Ndoft;i++){
  //fprintf (stdout,"\n  %20.15le",Lsrst->lhs[i]);
  //}
  //fprintf (stdout,"\n");
  
  if(Tp->homogt == 0)//in case of no homogenization
    compute_req_valt (0);
  
/*
  FILE *p;
  //double zdroj;
  p = fopen ("elpole.txt","w");
  for (i=0;i<Tt->nn;i++){
    //zdroj = (Tt->nodes[i].gradient[0][0]*Tt->nodes[i].gradient[0][0] + Tt->nodes[i].gradient[0][1]*Tt->nodes[i].gradient[0][1])*0.0026;
    //fprintf (p,"%10.8le %10.8le   %12.10le\n",Gtt->gnodes[i].x,Gtt->gnodes[i].y,zdroj);
    fprintf (p,"%10.8le %10.8le   %12.9le  %12.9le\n",Gtt->gnodes[i].x,Gtt->gnodes[i].y,Tt->nodes[i].gradient[0][0],Tt->nodes[i].gradient[0][1]);
  }
  fclose (p);
*/
  
  
  surface_fluxes (Outt);
  
  
  print_initt(-1, "wt");    
  //for (i=0;i<Lsrst->nlc;i++){
  //  computes and prints required quantities
  print_stept(i, 0, 0.0, rhs);  //print_stept(i, 0, 0.0, NULL);
  //}
  print_closet();
  
  /// adaptivity
  if (Tp->adaptivityflag)
    Adat->run (2, false);
}


/**
   function solves linear stationary transport problem

   JK, 26.7.2011
*/
void solve_radiation_stationary_problem ()
{
  long i,ni,lcid;
  double *rhs;
  
  //  conductivity matrix assembling
  conductivity_matrix (0);
  
  rhs=Lsrst->give_rhs (0);
  
  lcid = 0;
  
  ni=100;
  
  print_initt(-1, "wt");    
  for (i=0;i<ni;i++){
    
    nullv (rhs,Ndoft);
    trfel_right_hand_side (lcid,rhs,Ndoft);
    
    Tt->edge_temperature ();
    Tt->heat_fluxes (rhs,Outt);

    fprintf (Outt,"\n\n\n i %ld",i);
    //for (j=0;j<Ndoft;j++){
    //fprintf (Outt,"\n rhs %20.15le",rhs[j]);
    //}

    fprintf (Outt,"\n\n\n matice vodivosti pred eliminaci, krok  %ld",i);
    Kmat->printmat (Outt);
    Tp->ssle->solve_system (Gtt,Kmat,Lsrst->give_lhs(0),Lsrst->give_rhs(0),Outt);
    fprintf (Outt,"\n\n\n matice vodivosti po eliminaci, krok  %ld",i);
    Kmat->printmat (Outt);
    
    print_stept(lcid, 0, 0.0, rhs);  //print_stept(i, 0, 0.0, NULL);
    
  }
  
  
  //for (i=0;i<Lsrst->nlc;i++){
  //  computes and prints required quantities
  //}
  print_closet();
}

