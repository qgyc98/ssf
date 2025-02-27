#include "nspsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "elemswitcht.h"
#include "transprint.h"
#include <string.h>
#include "npsolvert.h"




/**
   function solves non-linear stationary transport problem

   TKr, 16/03/2010
*/
void solve_nonlinear_stationary_problem_pokus ()
{
  long j,k,n,lcid,ni,ini;
  double zero,err,lambda,dlambda,dlambdamin,dlambdamax,norf;
  double *r,*rb,*dr,*f,*fi,*fb,*fp,*rhs,*jk;
  
  //  load case id must be equal to zero in this case
  lcid=0;
  //  number of rows of the matrix
  n = Ndoft;
  //  number of increments in the Newton-Raphson method
  ni = Tp->nlman.ninr;
  //  maximum number of iterations in one increment
  ini = Tp->nlman.niilnr;
  //  computer zero
  zero = Tp->zero;
  //  required error in inner loop
  err = Tp->nlman.errnr;
  //  increment size
  dlambda=Tp->nlman.incrnr;
  //  minimum increment size
  dlambdamin=Tp->nlman.minincrnr;
  //  maximum increment size
  dlambdamax=Tp->nlman.maxincrnr;

  rhs=Lsrst->give_rhs (0);
  //lhsi=Lsrst->give_lhsi (0);
  

  rb  = new double [n];
  f  = new double [n];
  fb = new double [n];
  fi = new double [n];
  dr = new double [n];
  jk = new double [n];
  memset (rb,0,n*sizeof(double));
  memset (f,0,n*sizeof(double));
  memset (fb,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));
  memset (dr,0,n*sizeof(double));
 
  //  initiation of transport material models
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation (); 

  //  right hand side of the system
  trfel_right_hand_side (lcid,rhs,n);
  
  //  nodal values
  r  = Lsrst->give_lhs (lcid);
  //  vector of proportional fluxes
  fp = Lsrst->give_rhs (lcid*2);
  //fc = Lsrst->give_rhs (lcid*2);
  //  vector of constant fluxes
  //fc = Lsrst->give_rhs (lcid*2+1);
  //fp = Lsrst->give_rhs (lcid*2+1);
  //  proportional variable
  lambda=0.0;
  
  //  conductivity matrix assembling
  conductivity_matrix (lcid);
  
  copyv (fp,fb,n);
  
  //  solution of equation system
  Tp->ssle->solve_system (Gtt,Kmat,r,fb,Outt);
  
  nullv (fb,n);
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  print_initt(-1, "wt");
  
  //  internal iteration loop
  for (j=0;j<ini;j++){
    
    if(Tp->homogt != 1)//29/03/2011 this must be checked!
      solution_correction ();

    compute_req_valt (0);
    approximation ();
    
    // full newton-raphson
    conductivity_matrix (0);
    
    if (Tp->trestype==lrhst){
      //  right hand side of transport part
      trfel_right_hand_side (0,rhs,n);
      // Solver computes residuum from system of equations
      Kmat->gmxv (r,fb);
      // Res. = F - K*r
      for (k=0;k<n;k++){
	fi[k] = rhs[k] - fb[k];
      }
    }
    
    if (Tp->trestype==fluxest){
      internal_fluxes (fi,n);
    }
    
    //  norm of vector of unbalanced fluxes
    norf=normv(fi,n);
    //norf /= norfa;
    
    fprintf (stdout,"\n inner loop j %ld     norf=%20.16le  dlambda %e",j,norf,dlambda);
    
    if (norf<err)
      break;
    
    //  solution of K(r).v=F
    Tp->ssle->solve_system (Gtt,Kmat,dr,fi,Outt);
    
    for (k=0;k<n;k++){
      r[k]+=dr[k];
    }
    
  }
  
  print_stept(lcid, 0, lambda, f);
  print_flusht();
  print_closet();

}




/**
   function solves linear stationary transport problem

   JK, 20.12.2002
*/
void solve_nonlinear_stationary_problem ()
{
  long i,j,k,n,lcid,ni,ini;
  double zero,err,lambda,dlambda,dlambdamin,dlambdamax,norf;
  double *r,*rb,*dr,*f,*fi,*fb,*fp,*rhs,*lhsi,*jk;
  
  //  load case id must be equal to zero in this case
  lcid=0;
  //  number of rows of the matrix
  n = Ndoft;
  //  number of increments in the Newton-Raphson method
  ni = Tp->nlman.ninr;
  //  maximum number of iterations in one increment
  ini = Tp->nlman.niilnr;
  //  computer zero
  zero = Tp->zero;
  //  required error in inner loop
  err = Tp->nlman.errnr;
  //  increment size
  dlambda=Tp->nlman.incrnr;
  //  minimum increment size
  dlambdamin=Tp->nlman.minincrnr;
  //  maximum increment size
  dlambdamax=Tp->nlman.maxincrnr;
  i = 0;

  rhs=Lsrst->give_rhs (0);
  lhsi=Lsrst->give_lhsi (0);
  

  rb  = new double [n];
  f  = new double [n];
  fb = new double [n];
  fi = new double [n];
  dr = new double [n];
  jk = new double [n];
  memset (rb,0,n*sizeof(double));
  memset (f,0,n*sizeof(double));
  memset (fb,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));
  memset (dr,0,n*sizeof(double));

  //  right hand side of the system
  trfel_right_hand_side (lcid,rhs,n);
  
  //  nodal values
  r  = Lsrst->give_lhs (lcid);
  //  vector of proportional fluxes
  fp = Lsrst->give_rhs (lcid*2);
  //fc = Lsrst->give_rhs (lcid*2);
  //  vector of constant fluxes
  //fc = Lsrst->give_rhs (lcid*2+1);
  //fp = Lsrst->give_rhs (lcid*2+1);
  //  proportional variable
  lambda=0.0;
  
  
  
  //  conductivity matrix assembling
  conductivity_matrix (lcid);
  
  copyv (fp,fb,n);
  copyv (rhs,fi,n);

  //  solution of equation system
  Tp->ssle->solve_system (Gtt,Kmat,r,fb,Outt);

  //  conductivity matrix assembling
  conductivity_matrix (lcid);
  //Kmat->gmxv (r,fc);
  //subv(fi,fc,fi);
  
  //nullv (fc,n);
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  print_initt(-1, "wt");
  print_stept(lcid, 0, 0.0, f);
  print_flusht();

  solution_correction ();
  compute_req_valt (0);
  approximation ();
  
  
  //internal_fluxes (fi,n);
  
  
  //  norm of vector of unbalanced fluxes
  norf=normv(fi,n);
  //norf /= norfa;
  
  if (norf<err){
    
    print_stept(lcid, i+1, lambda, f);
    print_flusht();
    
  }
  else{
    
    //  internal iteration loop
    for (j=0;j<ini;j++){
      
      conductivity_matrix (lcid);
      //  solution of K(r).v=F
      Tp->ssle->solve_system (Gtt,Kmat,dr,fi,Outt);
      
      for (k=0;k<n;k++){
	r[k]+=dr[k];
      }
      
      conductivity_matrix (lcid);
      //Kmat->gmxv (r,fc);
      //subv(fi,fc,fi);

      solution_correction ();
      compute_req_valt (0);
      approximation ();
      
      //internal_fluxes (fi,n);
      
      //  norm of vector of unbalanced fluxes
      norf=normv(fi,n);
      //norf /= norfa;
      
      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf=%20.16le  dlambda %e",i,j,norf,dlambda);
      //      fprintf (Out,"\n j=%ld     norf=%e",j,norf);
      
      if (norf<err)
	break;
    }
  }
  
  
  print_stept(lcid, i+1, lambda, f);
  print_flusht();
  print_closet();

}



/**
   function solves linear stationary transport problem

   JK, 20.12.2002
*/
void solve_nonlinear_stationary_problem_old ()
{
  long i,j,k,n,lcid,ni,ini,modif;
  double zero,err,lambda,blambda,dlambda,dlambdamin,dlambdamax,norf;
  double *d,*dd,*db,*fc,*fp,*fb,*f,*rhs;
  
  //  load case id must be equal to zero in this case
  lcid=0;
  //  the number of rows of the matrix
  n = Ndoft;
  //  the number of increments in the Newton-Raphson method
  ni = Tp->nlman.ninr;
  //  maximum number of iterations in one increment
  ini = Tp->nlman.niilnr;
  //  computer zero
  zero = Tp->zero;
  //  required error in inner loop
  err = Tp->nlman.errnr;
  //  increment size
  dlambda=Tp->nlman.incrnr;
  //  minimum increment size
  dlambdamin=Tp->nlman.minincrnr;
  //  maximum increment size
  dlambdamax=Tp->nlman.maxincrnr;
  
  
  dd = new double [n];
  memset (dd,0,n*sizeof(double));
  db = new double [n];
  memset (db,0,n*sizeof(double));
  fc = new double [n];
  memset (fc,0,n*sizeof(double));
  fp = new double [n];
  memset (fp,0,n*sizeof(double));
  fb = new double [n];
  memset (fb,0,n*sizeof(double));
  f  = new double [n];
  memset (f,0,n*sizeof(double));

  //dr = new double [n];
  //memset (dr,0,n*sizeof(double));
  
  rhs=Lsrst->give_rhs (0);

  //  right hand side of the system
  trfel_right_hand_side (lcid,rhs,n);
  
  //  nodal values
  d  = Lsrst->give_lhs (lcid);
  //  vector of proportional fluxes
  fp = Lsrst->give_rhs (lcid*2);
  //  vector of constant fluxes
  fc = Lsrst->give_rhs (lcid*2+1);
  
  //  proportional variable
  lambda=0.0;
  
  
  //  conductivity matrix assembling
  conductivity_matrix (lcid);
  
  copyv (fc,fb,n);

  //  solution of equation system
  Tp->ssle->solve_system (Gtt,Kmat,d,fb,Outt);

  // ***************************
  //  main iteration loop  ****
  // ***************************
  print_initt(-1, "wt");
  print_stept(lcid, 0, 0.0, f);
  print_flusht();

  for (i=0;i<ni;i++){
    
    //  backup of left hand side vector
    for (j=0;j<n;j++){
      db[j]=d[j];
    }
    //  backup of reached lambda parameter
    blambda=lambda;
   
    fprintf (stdout,"\n increment %ld   lambda %e,    dlambda %e",i,lambda,dlambda);
 
    //  vector of maximum load and vector of load increment
    for (j=0;j<n;j++){
      f[j]=fc[j]+(lambda+dlambda)*fp[j];
      fb[j]=dlambda*fp[j];
    }
    
    //  conductivity matrix assembling
    conductivity_matrix (lcid);
    
    //  solution of equation system
    Tp->ssle->solve_system (Gtt,Kmat,dd,fb,Outt);
    
    for (j=0;j<n;j++){
      d[j]+=dd[j];
    }
    //  conductivity matrix assembling
    conductivity_matrix (lcid);

    solution_correction ();
    compute_req_valt (0);
    approximation ();
    
    //  computation of internal fluxes form actual values
    //internal_fluxes (fi,n);
    Kmat->gmxv (d,fb);
    
    //  vector of unbalanced fluxes
    for (j=0;j<n;j++){
      fb[j]=f[j]-fb[j];
    }
    //  norm of vector of unbalanced fluxes
    norf=normv(fb,n);
    
    if (norf<err){
      
      lambda+=dlambda;
      print_stept(lcid, i+1, lambda, f);
      print_flusht();

      if (lambda>1.00)
	break;
      
      modif++;
      if (modif>1){
	dlambda*=2.0;
	if (dlambda>dlambdamax)
	  dlambda=dlambdamax;

	if (Mesprt==1){
	  fprintf (stdout,"\n increment must be modified (dlambda increase) dlambda=%e",dlambda);
	}
	modif=0;
      }
      
      continue;
    }
    
    //  internal iteration loop
    for (j=0;j<ini;j++){
      
      //  solution of K(r).v=F
      Tp->ssle->solve_system (Gtt,Kmat,dd,fb,Outt);
    
      for (k=0;k<n;k++){
	d[k]+=dd[k];
      }

      //  conductivity matrix assembling
      conductivity_matrix (lcid);
      Kmat->gmxv (d,fb);
      
      solution_correction ();
      compute_req_valt (0);
      approximation ();
      
      //  computation of internal fluxes
      //internal_fluxes (fi,n);
      
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
	fb[k]=f[k]-fb[k];
      }
      
      //  norm of vector of unbalanced fluxes
      norf=normv(fb,n);
      //norf /= norfa;
      
      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf=%e  dlambda %e",i,j,norf,dlambda);
      //     fprintf (Out,"\n j=%ld     norf=%e",j,norf);
      
      if (norf<err)  break;
    }//  end of the loop for (j=0;j<ini;j++){
    
    modif=0;

    if (j==ini || norf>err){
      dlambda/=2.0;
      if (dlambda<dlambdamin){
	dlambda=dlambdamin;
	break;
      }

      if (Mesprt==1){
	fprintf (stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
//	fprintf (Out,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
      }
      
      for (j=0;j<n;j++){
      	d[j]=db[j];
      }
      lambda=blambda;
    }
    else{
      lambda+=dlambda;
      print_stept(lcid, i+1, lambda, f);
      print_flusht();
      
      if (lambda>1.0)
	break;
    }
    
  }//  end of the loop for (i=0;i<ni;i++){
  
  print_closet();
}
