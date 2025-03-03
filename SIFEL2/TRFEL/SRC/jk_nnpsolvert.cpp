#include "nnpsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "transprint.h"
#include <string.h>

void solve_nonlinear_nonstationary_problem ()
{
  //  solution of the problem
  //nonlinear_nonstat_solv_pokus ();
  nonlinear_nonstat_solv_new ();
  //nonlinear_nonstat_solv_linesearch ();
}


/**
   function solves nonlinear nonstationary transport problem
   + selective norm
   
   JK, 23.12.2002, 12.5.2005
*/
void nonlinear_nonstat_solv ()
{
  long i,j,k,n,ini,stop;
  double dt,end_time,alpha,norfb,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;

  n=Ndoft;
    
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  arhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);
  nullv (arhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (fb,n);
  nullv (fi,n);
  
  //  nodes - integration points interpolation
  approximation ();
  
  alpha=Tp->alpha;
  err=Tp->err;
  ini=Tp->nii;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);

    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    Kmat->solve_system (tdlhs,fb);
    
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    solution_correction ();
    approximation ();
    internal_fluxes (fi,n);
    
    //  vector of unbalanced fluxes
    for (j=0;j<n;j++){
      fb[j]=fi[j];
      //fprintf (Outt,"\nfb %ld  %e",j,fb[j]);
    }
    
    //norfb = ss (fb,fb,n);
    
    stop = selectivenorm (fb,err,i,0);

    //if (Mesprt==1)  fprintf (stdout,"\n%e %e",norfb,err);
    //if (Mesprt==1)  fprintf (Outt,"\nnorfb %e   err  %e",norfb,err);
    
    //if (norfb<err){
    if (stop==Tp->ntm){
      //  time increment
      Tp->time = Tp->timecont.newtime ();
      dt = Tp->timecont.actualbacktimeincr ();
      
      i++;
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      continue;

    }
    
    for (j=0;j<ini;j++){
      
      Kmat->solve_system (p,fb);
      
      for (k=0;k<n;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
      }
      
      solution_correction ();
      approximation ();
      internal_fluxes (fi,n);
    
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
	fb[k]= fi[k];
      }
      
      stop = selectivenorm (fb,err,i,j);
      //norfb = ss (fb,fb,n);
      
      //if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);
      
      //if (norfb<err){
      if (stop==Tp->ntm){
	break;
      }
    }
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d;

  print_closet();
}


/**
   function solves nonlinear nonstationary transport problem

   JK, 23.12.2002
*/
void nonlinear_nonstat_solv_old ()
{
  long i,j,k,n,ini;
  double dt,end_time,alpha,norfb,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;

  n=Ndoft;
    
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  arhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);
  nullv (arhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (fb,n);
  nullv (fi,n);
  
  //  nodes - integration points interpolation
  approximation ();
  
  alpha=Tp->alpha;
  err=Tp->err;
  ini=Tp->nii;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
          
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    Kmat->solve_system (tdlhs,fb);
    
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    solution_correction ();
    approximation ();
    internal_fluxes (fi,n);
    
    //  vector of unbalanced fluxes
    for (j=0;j<n;j++){
      fb[j]=fi[j];
      //fprintf (Outt,"\nfb %ld  %e",j,fb[j]);
    }
    
    norfb = ss (fb,fb,n);
    
    if (Mesprt==1)  fprintf (stdout,"\n%e %e",norfb,err);
    //if (Mesprt==1)  fprintf (Outt,"\nnorfb %e   err  %e",norfb,err);
    
    if (norfb<err){
      //  time increment
      Tp->time = Tp->timecont.newtime ();
      dt = Tp->timecont.actualbacktimeincr ();
      
      i++;
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      continue;

    }
    
    for (j=0;j<ini;j++){
      
      Kmat->solve_system (p,fb);
      
      for (k=0;k<n;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
      }
      
      solution_correction ();
      approximation ();
      internal_fluxes (fi,n);
    
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
	fb[k]= fi[k];
      }
      
      norfb = ss (fb,fb,n);
      
      if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);
      
      if (norfb<err){
	break;
      }
    }
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d;

  print_closet();
}




/**
   function solves nonlinear nonstationary transport problem

   TKr, 15.5.2005
*/
void nonlinear_nonstat_solv_pokus ()
{
  long i,j,k,n,ini;
  double dt,end_time,alpha,norfb,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;
  double norfb_last,*lhs_last;

  n=Ndoft;
    
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  lhs_last = new double [n];
  arhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (lhs_last,n);
  nullv (tdlhs,n);
  nullv (rhs,n);
  nullv (arhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (fb,n);
  nullv (fi,n);
  
  //  nodes - integration points interpolation
  approximation ();
  
  alpha=Tp->alpha;
  err=Tp->err;
  ini=Tp->nii;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
          
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    Kmat->solve_system (tdlhs,fb);
    
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    //solution_correction ();
    approximation ();
    internal_fluxes (fi,n);

    //  vector of unbalanced fluxes
    for (j=0;j<n;j++){
      fb[j]=fi[j];
      //fprintf (Outt,"\nfb %ld  %e",j,fb[j]);
    }
    
    norfb = ss (fb,fb,n);
    
    //tady nove??!!!
    for (k=0;k<n;k++){
      lhs_last[k] = lhs[k]; //storing results from previous inner step
    }
    norfb_last = norfb; //storing norfb from previous inner step      
    

    if (Mesprt==1)  fprintf (stdout,"\n%e %e",norfb,err);
    //if (Mesprt==1)  fprintf (Outt,"\nnorfb %e   err  %e",norfb,err);
    
    if (norfb<err){
      //  time increment
      Tp->time = Tp->timecont.newtime ();
      dt = Tp->timecont.actualbacktimeincr ();
      
      i++;
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      continue;

    }
    
    for (j=0;j<ini;j++){
      
      /*********************************************/
      //tady nove??!!!
      //znova se pocitaji matice
      
      //  capacity matrix
      capacity_matrix (0);
      
      //  conductivity matrix
      conductivity_matrix (0);
      //  auxiliary vector  K (d+(1-alpha)*dt*v)
      Kmat->gmxv (d,p);

      //  matrix of the system of equations
      //  C + alpha.dt.K
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      /*********************************************/

      Kmat->solve_system (p,fb);
      
      for (k=0;k<n;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
      }
      
      //solution_correction ();
      approximation ();
      internal_fluxes (fi,n);
    
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
	fb[k]= fi[k];
      }
      
      norfb = ss (fb,fb,n);
      
      if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);
      
      approximation ();//tady??!!

      if (norfb<err){
	break;
      }

      /*******************************************************/
      //tady nove??!!!
      //condition of decreesing of error

      if (norfb>norfb_last){
	for (k=0;k<n;k++){
	  lhs[k] = lhs_last[k]; //storing results from previous inner step
	}
	approximation ();
	break;
      }
      for (k=0;k<n;k++){
	lhs_last[k] = lhs[k]; //storing results from previous inner step
      }
      
      norfb_last = norfb; //storing norfb from previous inner step
      
      /*******************************************************/ 
    }
    

    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d; delete [] lhs_last;

  print_closet();
}


/**
   function solves nonlinear nonstationary transport problem

   TKr, 15.5.2005
*/
void nonlinear_nonstat_solv_new ()
{
  long i,j,k,n,ini;
  double dt,end_time,alpha,norfb,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;
  double norfb_last,*lhs_last;

  n=Ndoft;
    
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  lhs_last = new double [n];
  arhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (lhs_last,n);
  nullv (tdlhs,n);
  nullv (rhs,n);
  nullv (arhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (fb,n);
  nullv (fi,n);
  
  //  nodes - integration points interpolation
  approximation ();
  
  alpha=Tp->alpha;
  err=Tp->err;
  ini=Tp->nii;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
          
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    Kmat->solve_system (tdlhs,fb);
    
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }

    //solution_correction ();
    //approximation ();//jeste jednou??!!
    
    for (j=0;j<ini;j++){
      
      approximation ();//jeste jednou??!!
      
      /*********************************************/
      //znova se pocitaji matice + prava strana
      //  capacity matrix
      capacity_matrix (0);
      
      //  conductivity matrix
      conductivity_matrix (0);

      //  right hand side
      trfel_right_hand_side (0,rhs,n);
    
      //  auxiliary vector  K (d+(1-alpha)*dt*v)
      Kmat->gmxv (d,p);

      //  matrix of the system of equations
      //  C + alpha.dt.K
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      Kmat->gmxv (tdlhs,fi);//tady??!!

      //  Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
      for (k=0;k<n;k++){
	fb[k] = rhs[k] - p[k] - fi[k];
      }
      /*********************************************/
      
      norfb = ss (fb,fb,n);
      
      if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);

      //approximation ();//jeste jednou??!!

      if (norfb<err){
	break;
      }

      Kmat->solve_system (p,fb);
      
      for (k=0;k<n;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
      }
      
      //solution_correction ();
      //approximation ();//jeste jednou??!!
    }
    
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d; delete [] lhs_last;

  print_closet();
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

/**
   function solves nonlinear nonstationary transport problem

   JK, 21.2.2005
*/
void nonlinear_nonstat_solv_linesearch ()
{
  long i,j,k,l,n,ini,nils;
  double dt,end_time,alpha,norfb,err,eta,deta,a,errls,nom;
  double *d,*fb,*fi,*p,*pp,*lhs,*lhsi,*tdlhs,*rhs,*arhs;
  
  nils=10;
  errls=1.0e-6;

  n=Ndoft;
    
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  p = new double [n];
  pp = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (fb,n);
  nullv (fi,n);
  

  //  nodes - integration points interpolation
  approximation ();
  
  alpha=Tp->alpha;
  err=Tp->err;
  ini=Tp->nii;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    //  predictor
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    //  solution of the system of equations
    Kmat->solve_system (p,fb);
    
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      tdlhs[j]=p[j];
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    //solution_correction ();
    approximation ();
    
    //  vector of unbalanced fluxes
    internal_fluxes (fb,n);
    
    
    norfb = ss (fb,fb,n);
    
    if (Mesprt==1)  fprintf (stdout,"\n%e %e",norfb,err);
    //if (Mesprt==1)  fprintf (Outt,"\nnorfb %e   err  %e",norfb,err);
    
    if (norfb<err){
      //  time increment
      Tp->time = Tp->timecont.newtime ();
      dt = Tp->timecont.actualbacktimeincr ();
      
      i++;
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      continue;
      
    }
    
    
    for (j=0;j<ini;j++){
      
      //approximation ();
      //  capacity matrix
      capacity_matrix (0);
      //  conductivity matrix
      conductivity_matrix (0);
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      
      
      //  line search
      Kmat->gmxv (p,pp);
      a = ss (pp,p,n);
      eta =1.0;
      
      fprintf (stdout,"\n jsme ve vnitrni smycce");

      for (k=0;k<nils;k++){
	nom = ss (fb,p,n);
	deta = nom/a;
	eta-=deta;


	for (l=0;l<n;l++){
	  tdlhs[l]=eta*p[l];
	  lhs[l]=d[l]+alpha*dt*tdlhs[l];
	}
	
	//solution_correction ();
	approximation ();
	internal_fluxes (fb,n);
	
	norfb = ss (fb,fb,n);
	fprintf (stdout,"\n jsme v line search,  nom=%le      deta %le   eta %le  norfb %le",nom,deta,eta,norfb);

	if (fabs(nom)<errls)
	  break;
	
      }
      
      
      /*
      for (k=0;k<n;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*tdlhs[k];
      }
      
      approximation ();
      //  vector of unbalanced fluxes
      */

      
      norfb = ss (fb,fb,n);
      
      if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);
      
      if (norfb<err){
	break;
      }
      Kmat->solve_system (p,fb);
    }
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d;

  print_closet();
}


long selectivenorm (double *f,double err,long ni,long nii)
{
  long i,j,dof,stop;
  double norf;
  
  if (Mesprt==1)  fprintf (stdout,"\n\n number of iter. %ld   number of inner iter. %ld",ni,nii);

  stop=0;
  for (i=0;i<Tp->ntm;i++){
    norf=0.0;
    for (j=0;j<Tt->nn;j++){
      dof = Tt->give_dof (j,i);
      if (dof>0){
	dof--;
	norf+=f[dof]*f[dof];
      }
    }
    if (Mesprt==1)  fprintf (stdout,"\nvalue %ld  residuum  %e   req. error %e",i+1,norf,err);
    if (norf<err)  stop++;
  }
  
  return stop;
}
