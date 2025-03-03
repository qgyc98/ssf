#include "nnpsolvert.h"
#include "npsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "elemswitcht.h"
#include "transprint.h"
#include "backupsolt.h"
#include <string.h>

void solve_nonlinear_nonstationary_problem ()
{
  //  solution of the problem

  //nonlinear_nonstat_solv ();
  nonlinear_nonstat_solv_fnr_dform ();

  //nonlinear_nonstat_solv_dform ();
  //nonlinear_nonstat_solv_fnr_dform ();
  //nonlinear_nonstat_solv_nr_dform ();
  //nonlinear_nonstat_solv_oldd ();
  //nonlinear_nonstat_solv_vform ();
  //nonlinear_nonstat_solv_dform ();
  //nonlinear_nonstat_solv_old ();
  //nonlinear_nonstat_solv_new ();
  //nonlinear_nonstat_solv_pokus ();
  //nonlinear_nonstat_solv_linesearch ();
  //nonlinear_nonstat_solv_dform_dneska ();
}

/**
   Function solves system of non-linear TRFEL algebraic by Newton-Raphson method
   for time-dependent problems

   TKr 5.4.2007
*/
void nonlinear_nonstat_solv ()
{
  long i,j,k,nt,ini,nsts;
//  long stop;
  double zero,dt,end_time,alpha;
  double *p,*d,*lhs,*tdlhs,*rhs;
  double *lhsb, *tdlhsb;
  double *fbt,*fit,*lhs_last;
  //double s;
  double norf_last;
  double norf,*err,*thresh;  

  //  maximum number of iterations in inner loop
  ini = Tp->nii;
  //  required norm of vector of unbalanced forces
  err = Tp->errarr;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;
  //  number of transport degrees of freedom
  nt=Ndoft;

  //  left hand side
  lhs=Lsrst->give_lhs (0);
  tdlhs=Lsrst->give_tdlhs (0);
  //  right hand side
  rhs=Lsrst->give_rhs (0);

  //  array containing predictors
  d = new double [nt];
  //  auxiliary vector
  p = new double [nt];
  fbt = new double [nt];
  fit = new double [nt];
  lhs_last = new double [nt];
  lhsb = new double[nt];
  tdlhsb = new double[nt];


  //  initial values
  nullv (lhs,nt);
  nullv (tdlhs,nt);
  nullv(lhsb, nt);
  nullv(tdlhsb, nt);
  nullv (d,nt);
  nullv (p,nt);
  nullv (fbt,nt);
  nullv (fit,nt);
  
  //  coefficient of the trapezoidal rule  
  alpha=Tp->alpha;
  //  computer zero
  zero=Tp->zero;
  
  //  starting time
  Tp->time=Tp->timecont.starttime ();
  //  time increment
  dt=Tp->timecont.initialtimeincr ();
  //  end time
  end_time = Tp->timecont.endtime ();
  
  //  initiation of transport material models
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation (); 
  
  // vlozil JM 15.1.2008
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();

  // **************************
  //  main iteration loop  ***
  // **************************  
  i=0;
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (lhs,tdlhs,rhs,i,Tp->time,dt,Tp->timecont,nt);
    print_initt(-1, "at");
    //print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  }

  nsts = 0;
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;

    //  conductivity matrix
    conductivity_matrix (0);
    
    
    //  capacity matrix
    capacity_matrix (0);
    
    for (j=0;j<nt;j++){
      p[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side of transport part
    trfel_right_hand_side (0,rhs,nt);
    
    /*
    fprintf (Outt,"\n\n\n\n");
    for (j=0;j<nt;j++){
      fprintf (Outt,"\n rhs %4ld    %lf",j,rhs[j]);
    }
    fprintf (Outt,"\n\n\n\n");
    */

    //  auxiliary vector
    Kmat->gmxv (p,d);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    Kmat->copygm (*Cmat);
    
    for (j=0;j<nt;j++){
      rhs[j] = rhs[j] - d[j];
    }
    
    //  solution of the system of algebraic equations
    Tp->ssle->solve_system (Gtt,Cmat,tdlhs,rhs,Outt);

    for (j=0;j<nt;j++){
      //verze1:
      //s=(1.0-alpha)*d[j]+alpha*tdlhs[j];
      //lhs[j]+=dt*s;
      //verze2:
      lhs[j] = p[j] + alpha*dt*tdlhs[j];
    }
    
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into integration points
    approximation ();

    // vlozil JM 15.1.2008
    // nulleqother ();
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    norf_last = 1.0e20;
    //inner iteration loop//////////////////
    for (j=0;j<ini;j++){
      
      //  physically corrected solution
      //  JK, 25.11.2010
      //solution_correction ();
      //  approximation of nodal values into ontegration points
      // JK, 25.11.2010
      //approximation ();

      // vlozil JM 15.1.2008
      // nulleqother ();
      // JK, 25.11.2010
      //if (Tp->nvs==1 && Tp->pnvs==1)
      //actual_previous_nodval ();
      
      // full newton-raphson
      if (Tp->trsolv == fullnewtont){
	      // matrices are computing in each inner iteration
	      //  capacity matrix
	      capacity_matrix (0);
	
	      //  conductivity matrix
	      conductivity_matrix (0);
	
	      //  auxiliary vector  K (d+(1-alpha)*dt*v)
	      Kmat->gmxv (p,d);
	
	      //  matrix of the system of equations
	      //  C + alpha.dt.K
	      Kmat->scalgm (dt*alpha);
	      Kmat->addgm (1.0,*Cmat);
      }
      
      if (Tp->trestype==lrhst){
	      //  correct right hand side of transport part
	      trfel_right_hand_side (0,rhs,nt);
	      // Solver computes residuum from system of equations
	      Kmat->gmxv (tdlhs,fit);
	      // Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
	      for (k=0;k<nt;k++){
	        fbt[k] = rhs[k] - d[k] - fit[k];
	      }
      }
      
      if (Tp->trestype==fluxest){
	      // Solver computes unbalanced fluxes
	      internal_fluxes (fit,nt);
	      //  vector of unbalanced fluxes
	      for (k=0;k<nt;k++){
	        fbt[k]=fit[k];
	      }
      }
      
      norf = ss (fbt,fbt,nt);
      if (Mesprt==1)  fprintf (stdout,"\n inner iteration %ld   error %e",j,norf);
      if (norf<err[0]) break;

      //this is not working now:
      //stop = norm_computation_vec (rhs,fbt,err,thresh,2,1);      
      //if (stop==1)
      //break;

      Tp->ssle->solve_system (Gtt,Kmat,d,fbt,Outt);
      
      for (k=0;k<nt;k++){
	      tdlhs[k]+=d[k];
	      lhs[k]+=alpha*dt*d[k];
      }
      
      //  physically corrected solution
      solution_correction ();
      //  approximation of nodal values into ontegration points
      approximation ();
      
      // vlozil JM 15.1.2008
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	      actual_previous_nodval ();
      
      //  convergence control for newton-raphson
      //  condition of decreesing of error
      if (Tp->convergcontrolt==yes){
	      if (norf > norf_last){
	        for (k=0;k<nt;k++){
	          lhs[k] = lhs_last[k]; //storing results from previous inner step
	        }
	  
	        //  physically corrected solution
	        solution_correction ();
	        //  approximation of nodal values into ontegration points
	        approximation ();	  

	        // vlozil JM 15.1.2008
	        // nulleqother ();
	        if (Tp->nvs==1 && Tp->pnvs==1)
	          actual_previous_nodval ();
	  
	        if (Mesprt==1)  fprintf (stdout,"\n\n convergence control: inner iteration skiped %ld error %e\n\n",j,norf);	  
	  
	        break;
	      }
	      for (k=0;k<nt;k++){
	        lhs_last[k] = lhs[k]; //storing results from previous inner step
      	}
	
      	norf_last = norf; //storing norf from previous inner step
      }
    }
    ////////////////////////////////////////    

    if (norf < err[0]) {
      nsts++;
      if (nsts == 2)
      {
        dt *= 2.0;
        if (dt > Tp->timecont.dtmax)
          dt = Tp->timecont.dtmax;
        nsts = 0;
      }

      //  new time and new time increment
      Tp->time = Tp->timecont.newtime(dt);
      dt = Tp->timecont.actualbacktimeincr();

      // backup nodal values, lshsb = lhs
      copyv(lhs, lhsb, nt);
      // backup time derivatives of nodal values, tdlhsb = tdlhs
      copyv(tdlhs, tdlhsb, nt);

      //  printing of output and graphical informations
      compute_req_valt(0);
      print_stept(0, i, Tp->time, rhs);//print_stept(0,i,Tp->time,NULL);
      print_flusht();

      if ((Tp->timecont.isitimptime() == 1) && Tp->hdbcont.save_stat())
        solvert_save(lhs, tdlhs, rhs, i, Tp->time, dt, Tp->timecont, nt);
    }
    else {
      nsts = 0;
      dt *= 0.5;
      if (Tp->timecont.tct > 0) // for adaptive time controler
      {
        //restoring results from previous time step
        copyv(lhsb, lhs, nt);
        copyv(tdlhsb, tdlhs, nt);
      }
    }
  }while(Tp->time<=end_time);
  
  print_closet ();
  
  delete [] p;
  delete [] d;

  delete [] fit;
  delete [] fbt;
  delete [] lhs_last;
  delete [] lhsb;
  delete [] tdlhsb;
}





/**
   function solves nonlinear nonstationary transport problem
   + selective norm
   
   JK, 23.12.2002, 12.5.2005
*/
void nonlinear_nonstat_solv_oldd ()
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

  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);

    Tm->updateipval ();

    fprintf (Outt,"\n\n krok %ld",i);
    fprintf (Outt,"\n %le   %le   %le   %le",Tm->ip[0].av[0],Tm->ip[1].av[0],Tm->ip[4].av[0],Tm->ip[5].av[0]);

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
    
    //Kmat->solve_system (Gtt,tdlhs,fb);
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,fb,Outt);
    
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
    
    stop = norm_computation (fb,rhs,err,i,0);

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
      
      //Kmat->solve_system (Gtt,p,fb);
      Tp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);
	  
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
      
      //stop = selectivenorm (fb,rhs,err,i,j);
      norfb = ss (fb,fb,n);
      
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
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 23.12.2002,
   revised 24.8.2006
*/
void nonlinear_nonstat_solv_vform ()
{
  long i,j,k,n,ini,stop;
  double dt,end_time,alpha,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  allocation of auxiliary arrays
  d = new double [n];
  arhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initialization
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

  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    //  predictor
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v) = K dd
    Kmat->gmxv (d,p);

    //discont_contributions (p);

    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  right hand side f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    
    //  f - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    //  time derivatives v_{n+1} are solved
    //Kmat->solve_system (Gtt,tdlhs,fb);
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,fb,Outt);

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    //  computation of values at integration points from nodal values
    approximation ();
    
    //  computation of internal fluxes (conductivity multiplies by gradient)
    internal_fluxes (fi,n);
    
    stop = norm_computation (fi,rhs,err,1,1);
    
    if (stop==0){
      //  iteration for equilibrium is required
      for (j=0;j<ini;j++){
	
	//Kmat->solve_system (Gtt,p,fi);
	Tp->ssle->solve_system (Gtt,Kmat,p,fi,Outt);

	for (k=0;k<n;k++){
	  tdlhs[k]+=p[k];
	  lhs[k]+=alpha*dt*p[k];
	}
	
	approximation ();
	//  vector of unbalanced fluxes
	internal_fluxes (fi,n);
	
	if (Mesprt==1)  fprintf (stdout,"\n iteration %ld",j);
	stop = norm_computation (fi,rhs,err,1,1);
	
	if (stop==1){
	  break;
	}
      }
    }
    
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d;
  
  print_closet();
}

/**
   function solves nonlinear nonstationary transport problem
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 23.12.2002,
   revised 24.8.2006
*/
void nonlinear_nonstat_solv_dform ()
{
  long i,j,k,n,ini,stop;
  double dt,end_time,alpha,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*prhs,*plhs,*ptdlhs;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  allocation of auxiliary arrays
  d = new double [n];
  prhs = new double [n];
  plhs = new double [n];
  ptdlhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initialization
  nullv (lhs,n);
  nullv (plhs,n);
  nullv (tdlhs,n);
  nullv (ptdlhs,n);
  nullv (rhs,n);
  nullv (prhs,n);
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
  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    
    //nullv (fb,n);
    //discont_contributions (fb,alpha,dt);

    
    //Kmat->gmxv (lhs,d);
    Cmat->gmxv (lhs,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  right hand side f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    
    fprintf (Outt,"\n\n\n rhs %lf",rhs[0]);

    //  C.d + alpha dt f_n+1
    for (j=0;j<n;j++){
      fb[j] = rhs[j]*alpha*dt + p[j];
    }
    
    copyv (lhs,plhs,n);

    fprintf (Outt,"\n fb %lf\n",fb[0]);

    //  nodal values d_{n+1} are solved
    //Kmat->solve_system (Gtt,lhs,fb);
    Tp->ssle->solve_system (Gtt,Kmat,lhs,fb,Outt);

    if (alpha<Tp->zero){
      for (j=0;j<n;j++){
	tdlhs[j]=ptdlhs[j];
      }
    }
    else{
      for (j=0;j<n;j++){
	tdlhs[j]=(lhs[j]-plhs[j])/dt/alpha;
      }
    }
    
    //  computation of values at integration points from nodal values
    approximation ();
    
    //  computation of internal fluxes (conductivity multiplies by gradient)
    internal_fluxes (fi,n);
    
    stop = norm_computation (fi,rhs,err,1,1);
    
    for (j=0;j<n;j++){
      //  ??????????????????????? f_z  ???????????
      //fb[j]=rhs[j]-fi[j];
      fb[j]=0.0-fi[j];
    }
    
    if (stop==0){
      //  iteration for equilibrium is required
      for (j=0;j<ini;j++){
	
	//conductivity_matrix (0);
	
	//Kmat->solve_system (Gtt,p,fi);
	Tp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

	for (k=0;k<n;k++){
	  lhs[k]-=p[k];
	}
	
	if (alpha<Tp->zero){
	  for (k=0;k<n;k++){
	    tdlhs[j]=ptdlhs[j];
	  }
	}
	else{
	  for (k=0;k<n;k++){
	    tdlhs[k]=(lhs[k]-plhs[k])/dt/alpha;
	  }
	}
	
	approximation ();
	//  vector of unbalanced fluxes
	internal_fluxes (fi,n);
	
	for (k=0;k<n;k++){
	  //  ??????????????????????? f_z  ???????????
	  fb[k]=rhs[k]-fi[k];
	}

	if (Mesprt==1)  fprintf (stdout,"\n iteration %ld",j);
	stop = norm_computation (fi,rhs,err,1,1);
	
	
	
	if (stop==1){
	  break;
	}
      }
    }
    
    
    copyv (rhs,prhs,n);
    copyv (tdlhs,ptdlhs,n);
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] prhs;  delete [] tdlhs;  delete [] d;
  
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

  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
          
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
    
    //Kmat->solve_system (Gtt,tdlhs,fb);
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,fb,Outt);

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
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

      //Kmat->solve_system (Gtt,p,fb);
      Tp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

      for (k=0;k<n;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
      }
      
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
  long i,j,k,n,ini,stop;
  double dt,end_time,alpha,norfb,norrhs,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;
  double *lhs_last;

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

  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;
         
    //  capacity matrix
    capacity_matrix (0);
    
    
    //fprintf (Outt,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    //fprintf (Outt,"C75 %le    C76 %le   C77 %le",Cmat->give_entry (6,4),Cmat->give_entry (6,5),Cmat->give_entry (6,6));

    //  conductivity matrix
    conductivity_matrix (0);
    
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    
    norrhs = ss (rhs,rhs,n);

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
    
    
    
    //Kmat->printmat (Outt);
    
    //fprintf (Outt,"\nC75 %le    C76 %le   C77 %le",Kmat->give_entry (6,4),Kmat->give_entry (6,5),Kmat->give_entry (6,6));
    
    /*
    fprintf (Outt,"\n\n\n nova prava strana \n");
    for (j=0;j<n;j++){
      fprintf (Outt,"   %le",fb[j]);
    }
    */
    fprintf (Outt,"\n");
    
    
    
    //  solution of system of algebraic equations
    //Kmat->solve_system (Gtt,tdlhs,fb);
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,fb,Outt);

    //fprintf (Outt,"\n\n TDLHS krok %ld\n",i);

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      //fprintf (Outt,"%e\n",tdlhs[j]);

      lhs[j] = d[j] + alpha*dt*tdlhs[j];
      //fprintf (Outt,"\n  lhs  %ld     %lf  %lf",j,lhs[j],lhsi[j]);
    }
    
    
    
    for (j=0;j<ini;j++){
      
      //  zakomentovano v pondeli 3.10. kvuli JM - soli
      solution_correction ();
      approximation ();
      
      fprintf (Outt,"\n\n krok %ld",i);
      fprintf (Outt,"\n %le   %le   %le   %le",Tm->ip[0].av[0],Tm->ip[1].av[0],Tm->ip[4].av[0],Tm->ip[5].av[0]);
      

      // *********************************************
      //znova se pocitaji matice + prava strana
      
      if (Tp->trsolv == fullnewtont){
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
	
      }
      
      
      //  umisteno v pondeli 3.10. kvuli solim
      solution_correction ();
      approximation ();
      
      
      if (Tp->trestype==lrhst){
	Kmat->gmxv (tdlhs,fi);//tady??!!
	//  Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
	for (k=0;k<n;k++){
	  fb[k] = rhs[k] - p[k] - fi[k];
	}
      }
      if (Tp->trestype==fluxest){
	internal_fluxes (fi,n);
	//  vector of unbalanced fluxes
	for (k=0;k<n;k++){
	  fb[k]=fi[k];
	  //fb[k]=rhs[k]-fi[k];
	}
      }
      // *********************************************
      
      norfb = ss (fb,fb,n);

      stop = norm_computation (fb,rhs,err,i,j);
      
      //if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);

      //fprintf (Outt,"\n iteration %ld   error %e,  normed error %e",j,norfb,norfb/norrhs);
      
      //approximation ();//jeste jednou??!!
      
      /*
      if (norrhs>zero){
	if (norfb/norrhs<err)
	  break;
      }
      else{
	if (norfb<err){
	  break;
	}
      }
      */
      
      if (stop==Tp->ntm)
	break;
      
      //Kmat->solve_system (Gtt,p,fb);
      Tp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

      //fprintf (Outt,"\n\n P-> TDLHS krok %ld  iterace %ld\n",i,j);
      for (k=0;k<n;k++){
	//fprintf (Outt,"%e\n",p[k]);

	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
	//fprintf (Outt,"\n  lhs  %ld     %lf  %lf",k,lhs[k],lhsi[k]);
      }
      
      //approximation ();//jeste jednou??!!
    }
    
    //  pondeli 10.7.2006
    //  MADERA
    /*
    long *cn;
    double *r;
    cn = new long [4];
    r = new double [4];
    for (j=0;j<Tt->ne;j++){
      Tt->give_code_numbers (j,cn);
      nodalvalues (0,j,r,cn,4);
      fprintf (Outt,"\n prvek %ld   %le %le %le %le",j,r[0],r[1],r[2],r[3]);
    }
    delete [] cn;
    delete [] r;
    */
    //  konec upravy
    
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
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
  double *d,*fb,*fi,*p,*pp,*lhs,*lhsi,*tdlhs,*rhs;
  
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
  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();

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
    //Kmat->solve_system (Gtt,p,fb);
    Tp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      tdlhs[j]=p[j];
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
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
      //Kmat->solve_system (Gtt,p,fb);
      Tp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

    }
    
    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] d;

  print_closet();
}

/**
   function computes norm of the %vector of unbalanced fluxes
   it enables computation of Euclidean norm of the %vector
   as well as selective norm (appropriate components are collected
   togehter, e.g. heat fluxes are handled separately from moisture
   fluxes, both norms are then compared with required error)
   
   norm of unbalanced fluxes can be scaled with respect to the
   prescribed fluxes (right hand side of the system of equations)
   
   function returns two possible values
   0 means that equilibrium is not reached
   1 means that equilibrium is reached and unbalanced fluxes are smaller than required
   
   @param f - %vector of unbalanced fluxes
   @param rhs - %vector of right hand side (prescribed fluxes)
   @param err - required (prescribed) error
   
   @param nt - type of norm
               1 - Euclidean norm
	       2 - selective norm
   @param sc - scaling by norm of right hand side %vector
               0 - no scaling
	       1 - scaling by the norm of the right hand side %vector

   JK, 24.8.2006
*/
long norm_computation (double *f,double *rhs,double err,long nt,long sc)
{
  long i,j,dof,stop,pstop;
  double norf,norrhs;
  
  if (nt==1){
    // *****************************
    //  Euclidean norm is computed
    // *****************************
    norf=0.0;
    for (i=0;i<Ndoft;i++){
      norf+=f[i]*f[i];
    }
    norf=sqrt(norf);
    
    if (sc==0){
      //  norm is not scaled
      if (norf<err)
	stop=1;
      else
	stop=0;
      
      if (Mesprt==1)  fprintf (stdout,"\n norm of residuum  %e,    required error  %e",norf,err);
    }

    if (sc==1){
      //  norm is scaled
      norrhs=0.0;
      for (i=0;i<Ndoft;i++){
	norrhs+=rhs[i]*rhs[i];
      }
      norrhs=sqrt(norrhs);
      
      if (norrhs>0.0){
	if (norf/norrhs<err)
	  stop=1;
	else
	  stop=0;
	
	if (Mesprt==1)  fprintf (stdout,"\n norm of residuum  %e,    norm of rhs  %e,   error  %e,   required error  %e",norf,norrhs,norf/norrhs,err);
      }
      else{
	if (norf<err)
	  stop=1;
	else
	  stop=0;
	
	if (Mesprt==1)  fprintf (stdout,"\n norm of residuum  %e,    norm of rhs  %e,  required error  %e",norf,norrhs,err);
      }
    }
  }
  
  
  if (nt==2){
    // ********************************************
    //  selective norm is computed
    //  appropriate components are added togehter
    // ********************************************
    
    pstop=0;
    for (i=0;i<Tp->ntm;i++){
      norf=0.0;
      norrhs=0.0;
      for (j = 0; j < Tt->nn; j++) {
        if (Gtt->give_ndofn(j) > 0) { // i-th node is NOT hanging node
          dof = Tt->give_dof(j, i);
          if (dof > 0) {
            dof--;
            //norf+=f[dof]*f[dof];
            norf += (f[dof] - rhs[dof])*(f[dof] - rhs[dof]);
            norrhs += rhs[dof] * rhs[dof];
          }
        }
      }
      norf=sqrt(norf);
      norrhs=sqrt(norrhs);
      
      if (sc==0){
        //  norm is not scaled
        if (norf<err)
          pstop++;
	
        if (Mesprt==1)  fprintf (stdout,"\n phase %ld,  norm of residuum  %.15e,    required error  %.15e",i+1,norf,err);
      }
      
      if (sc==1){
        //  norm is scaled
        if (norrhs>0.0){
          if (norf/norrhs<err)
            pstop++;
	  
          if (Mesprt==1)  fprintf (stdout,"\n phase %ld,   norm of residuum  %.15e,    norm of rhs  %.15e,   error  %e,   required error  %e",i+1,norf,norrhs,norf/norrhs,err);
        }
        else{
          if (norf<err)
            pstop++;
	  
          if (Mesprt==1)  fprintf (stdout,"\n phase %ld,   norm of residuum  %.15e,    norm of rhs  %.15e,  required error  %e",i+1,norf,norrhs,err);
        }
      }
    }
    
    if (pstop==Tp->ntm)
      stop=1;
    else
      stop=0;
  }

  return stop;
}


/**
   function computes norm of the %vector of unbalanced fluxes
   it enables computation of Euclidean norm of the %vector
   as well as selective norm (appropriate components are collected
   togehter, e.g. heat fluxes are handled separately from moisture
   fluxes, both norms are then compared with required error)
   
   norm of unbalanced fluxes can be scaled with respect to the
   prescribed fluxes (right hand side of the system of equations)
   
   function returns two possible values
   0 means that equilibrium is not reached
   1 means that equilibrium is reached and unbalanced fluxes are smaller than required
   
   @param f - %vector of unbalanced fluxes
   @param rhs - %vector of right hand side (prescribed fluxes)
   @param err - array of required (prescribed) errors
   @param thresh - array of thresholds for particular variables
   @param nt - type of norm
               1 - Euclidean norm
	       2 - selective norm
   @param sc - scaling by norm of right hand side %vector
               0 - no scaling
	       1 - scaling by the norm of the right hand side %vector
   @param norfb - computed norm of the residual

   @retval 0 - if the norm(s) is(are) NOT in the given tolerance(s)
   @retval 1 - if the norm(s) is(are) in the given tolerance(s)

   JK, 24.8.2006
*/
long norm_computation_vec (double *f,double *rhs,double *err,double *thresh,long nt,long sc, double &norfb)
{
  long i,j,dof,stop,pstop;
  double norf,norrhs;
  stop = 0L;

  if (nt==1){
    // *****************************
    //  Euclidean norm is computed
    // *****************************
    norf=0.0;
    for (i=0;i<Ndoft;i++){
      norf+=f[i]*f[i];
    }
    norfb = norf = sqrt(norf);
    
    if (sc==0){
      //  norm is not scaled
      if (norf<err[0])
        stop=1;  // norm is in the given tolerance
      else
        stop=0; // norm is not in the given tolerance
      
      if (Mesprt==1)  fprintf (stdout,"\n norm of residuum  %e,    required error  %e",norf,err[0]);
    }
    
    if (sc==1){
      //  norm is scaled
      norrhs=0.0;
      for (i=0;i<Ndoft;i++){
        norrhs+=rhs[i]*rhs[i];
      }
      norrhs=sqrt(norrhs);
      
      if (norrhs>thresh[0]){
        if (norf/norrhs<err[0])
          stop=1;  // norm is in the given tolerance
        else
          stop=0; // norm is not in the given tolerance
        norfb /= norrhs;	
        if (Mesprt==1)  fprintf (stdout,"\n norm of residuum  %.15e,    norm of rhs  %.15e,   error  %e,   required error  %e",norf,norrhs,norf/norrhs,err[0]);
      }
      else{
        if (norf<err[0])
          stop=1;  // norm is in the given tolerance
        else
          stop=0; // norm is not in the given tolerance
        if (Mesprt==1)  fprintf (stdout,"\n norm of residuum  %.15e,    norm of rhs  %.15e,  required error  %e",norf,norrhs,err[0]);
      }
    }
  }
  
  
  if (nt==2){
    // ********************************************
    //  selective norm is computed
    //  appropriate components are added togehter
    // ********************************************
    
    pstop=0;
    norfb = 0.0;
    for (i=0;i<Tp->ntm;i++){
      norf=0.0;
      norrhs=0.0;
      for (j = 0; j < Tt->nn; j++) {
        if (Gtt->give_ndofn(j) > 0) { // i-th node is NOT hanging node
          dof = Tt->give_dof(j, i);
          if (dof > 0) {
            dof--;
            norf += f[dof] * f[dof];
            norrhs += rhs[dof] * rhs[dof];
          }
        }
      }
      norf=sqrt(norf);
      norrhs=sqrt(norrhs);
      
      if (sc==0){
        //  norm is not scaled
        if (norf<err[i])
          pstop++;
        
        norfb += sqr(norf);
	
        if (Mesprt==1)  fprintf (stdout,"\n phase %ld,  norm of residuum  %.15e,    required error  %e",i+1,norf,err[i]);
      }
      
      if (sc==1){
        //  norm is scaled
        if (norrhs>thresh[i]){
          if (norf/norrhs<err[i])
            pstop++;
          norfb += sqr(norf)/sqr(norrhs);
          if (Mesprt==1)  fprintf (stdout,"\n phase %ld,   norm of residuum  %.15e,    norm of rhs  %.15e,   error  %e,   required error  %e",i+1,norf,norrhs,norf/norrhs,err[i]);
        }
        else{
          if (norf<err[i])
            pstop++;
          norfb += sqr(norf);	  
          if (Mesprt==1)  fprintf (stdout,"\n phase %ld,   norm of residuum  %.15e,    norm of rhs  %.15e,  required error  %e",i+1,norf,norrhs,err[i]);
        }
      }
      //fprintf (stdout,"\n phase %ld, pstop %ld",i+1,pstop);
    }
    
    if (pstop==Tp->ntm)
      stop=1;   // all norms ar in the tolerances
    else
      stop=0;   // norms ar not in the given tolerances
    norfb = sqrt(norfb);
  }
  
  return stop;
}












































/**
   function solves nonlinear nonstationary transport problem
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 31.7.2007
*/
void nonlinear_nonstat_solv_dform_dneska ()
{
  long i,j,k,n,ini;
  double dt,end_time,alpha,err,norm;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*plhs,*pplhs,*rhs;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values /  vector of unknowns
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  nodal values from the previous time step
  plhs = new double [n];
  //  nodal values from the time step i-2
  pplhs = new double [n];
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  allocation of auxiliary arrays
  //  predictor
  d = new double [n];
  //  auxiliary vector
  p = new double [n];
  //  auxiliary vector
  fb = new double [n];
  //  auxiliary vector
  fi = new double [n];

  //  initialization
  nullv (lhs,n);
  nullv (plhs,n);
  nullv (rhs,n);
  nullv (d,n);
  nullv (fb,n);
  
  for (i=0;i<n;i++){
    plhs[i]=lhsi[i];
    pplhs[i]=lhsi[i];
  }

  //  nodes - integration points interpolation
  //approximation ();

  //actual_previous_change ();
  
  
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


  approximation ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_change ();

  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
/* fprintf (Outt,"\n\n\n CAPACITY MATRIX \n\n\n");
    Cmat->printmat (Outt);
    Cmat->printdiag (Outt);
    fprintf (Outt,"\n\n\n CONDUCTIVITY MATRIX \n\n\n");
    Kmat->printmat (Outt);
    Kmat->printdiag (Outt);*/

    nullv (fb,n);
    if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem){
      //discont_contrib (fb,alpha,dt);
    }
    
    //  predictor
    for (j=0;j<n;j++){
      d[j]=lhs[j]+(1.0-alpha)*dt*tdlhs[j];
    }
    
    //  1/alpha/dt C d = p
    Cmat->gmxv (d,p);
    cmulv(1.0/alpha/dt,p,n);
    
    
    //  matrix of the system of equations
    //  1/alpha/dt C + K
    Cmat->scalgm (1.0/dt/alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  right hand side f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    
    //  f_n+1 + 1/alpha/dt/ C d + fb
    for (j=0;j<n;j++){
      fb[j] += rhs[j] + p[j];
    }
    
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,fb,Outt);

    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }

    //  physically corrected solution
    solution_correction ();    
    //  computation of values at integration points from nodal values
    approximation ();
    
    
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);

    Cmat->gmxv (tdlhs,p);
    Kmat->gmxv (lhs,fi);
    
    nullv (fb,n);
    if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem){
      //discont_contrib (fb,alpha,dt);
    }
    
    norm=0.0;
    for (j=0;j<n;j++){
      fi[j] = rhs[j] + fb[j] - fi[j] - p[j];
      norm+=fi[j]*fi[j];
    }
    
    //fprintf (Outt,"\n norma rezidua v kroku  %ld (cas %lf) je %le",i,Tp->time,norm);
    
    //stop = norm_computation (fi,rhs,err,1,1);

    if (norm>err){
      //  iteration for equilibrium is required
      for (j=0;j<ini;j++){
	
	
	//  matrix of the system of equations
	//  1/alpha/dt C + K
	Cmat->scalgm (1.0/dt/alpha);
	Kmat->addgm (1.0,*Cmat);
	
	Tp->ssle->solve_system (Gtt,Kmat,p,fi,Outt);

	for (k=0;k<n;k++){
	  lhs[k]+=p[k];
	}
	
	for (k=0;k<n;k++){
	  tdlhs[k]=(lhs[k]-d[k])/dt/alpha;
	}
	
	//  physically corrected solution
    solution_correction ();    
	approximation ();


	//  capacity matrix
	capacity_matrix (0);
	
	//  conductivity matrix
	conductivity_matrix (0);

	Cmat->gmxv (tdlhs,p);
	Kmat->gmxv (lhs,fi);
	
	nullv (fb,n);
	if (Tp->tprob == discont_nonstat_problem || Tp->tprob == discont_nonlin_nonstat_problem){
	  //discont_contrib (fb,alpha,dt);
	}
	
	norm=0.0;
	for (k=0;k<n;k++){
	  fi[k] = rhs[k] + fb[k] - fi[k] - p[k];
	  norm+=fi[k]*fi[k];
	}
	
//	fprintf (Outt,"\n norma rezidua v kroku  %ld (cas %lf) vnitrni iterace %ld       je %le",i,Tp->time,j,norm);
	
	//if (Mesprt==1)  fprintf (stdout,"\n iteration %ld",j);
	//stop = norm_computation (fi,rhs,err,1,1);
	
	if (norm<err){
	  break;
	}
      }
    }

    
    compute_cycles (plhs,pplhs);

    if (i>1){
      for (j=0;j<n;j++){
	pplhs[j]=plhs[j];
      }
    }
    for (j=0;j<n;j++){
      plhs[j]=lhs[j];
    }
    
/*	fprintf(Outt, "\n\n\n pocty cyklu \n");
	for (i = 0; i <Tt->nn;i++){
		fprintf(Outt, "\n time %lf node %6ld  %ld",Tp->time, i,Tt->ncycl[i]);
	}*/

    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
  }while(Tp->time<end_time);
   
  fprintf(Outt, "\n\n\n pocty cyklu \n");
  for (i = 0; i <Tt->nn;i++){
	  fprintf(Outt, "\n node %6ld  %ld",i,Tt->ncycl[i]);
  }
  
  delete [] fb;  delete [] p;
  delete [] d;
  
  print_closet();
}








































/**
   function solves nonlinear nonstationary transport problem

   time discretization is based on the d-form version of
   the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   nonlinearity is solved by the full Newton-Raphson method
   
   JK, 23.4.2008
*/
void nonlinear_nonstat_solv_fnr_dform_old ()
{
  long i,j,k,n,ini,ani,nsts,lcid;
  double dt,dtdef,dtmin,end_time,time,alpha,err,norres,norrhs;
  double *d,*f,*p,*v,*lhs,*lhsi,*tdlhs,*rhs,*lhsb,*tdlhsb;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  //  load case id, it must be equal to zero for this type of computation
  lcid=0;

  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  allocation of auxiliary arrays
  //  predictor (d with tilde)
  d = new double [n];
  //  auxiliary vector
  p = new double [n];
  //  prescribed fluxes
  f = new double [n];
  //  auxiliary vector
  v = new double [n];
  //  backup of nodal values
  lhsb = new double [n];
  //  backup of time derivatives of nodal values
  tdlhsb = new double [n];
  

  //  initialization
  nullv (lhs,n);
  nullv (rhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (f,n);
  nullv (v,n);

  //  interpolation of nodal values to integration points
  approximation ();
  
  //  parameter of the trapezoidal method
  alpha=Tp->alpha;
  //  prescribed norm of residual
  err=Tp->err;
  //  maximum number of iterations in one time step
  ini=Tp->nii;
  //  minimum time increment
  dtmin=Tp->timecont.dtmin;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  time=Tp->time;
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  //  number of time step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n time %f,   time step = %f,    increment number %ld",Tp->time,dt,i);
    
    Tm->updateipval ();
    i++;
    
    //  backup of attained nodal values and their time derivatives
    for (j=0;j<n;j++){
      lhsb[j]=lhs[j];
      tdlhsb[j]=tdlhs[j];
    }
    
    //  capacity matrix
    capacity_matrix (lcid);
    
    //  conductivity matrix
    conductivity_matrix (lcid);
    
    //  predictor d
    for (j=0;j<n;j++){
      d[j] = lhs[j]+dt*(1.0-alpha)*tdlhs[j];
    }
    
    //  C.d
    Cmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed fluxes f_{n+1}
    trfel_right_hand_side (lcid,f,n);
    
    //  C.d + alpha dt f_n+1
    norrhs=0.0;
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
      norrhs+=f[j]*f[j];
    }
    norrhs=sqrt(norrhs);
    
    
    //  solution of the system of equations
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    
    //  computation of time derivatives of the nodal values
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    //  computation of values at integration points from nodal values
    approximation ();
    
    
    //  iteration for equilibrium
    for (j=0;j<ini;j++){
      
      //  capacity matrix
      capacity_matrix (lcid);
      
      //  conductivity matrix
      conductivity_matrix (lcid);
      
      
      //  matrix of the system of equations
      //  C + alpha.dt.K
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      
      //  (C + alpha.dt.K).d
      Kmat->gmxv (lhs,v);
      
      //  C.d
      Cmat->gmxv (d,p);

      //  computation of residuals
      norres=0.0;
      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k]-v[k];
	norres+=rhs[k]*rhs[k];
      }
      norres=sqrt(norres);
      
      if (Mesprt != 0)  
	fprintf (stdout,"\n nres   %15.10le, nrhs   %15.10le",norres,norrhs);
      
      if (norrhs<Tp->zero){
	if (norres<err){
	  break;
	}
      }
      else{
	if (Mesprt != 0)  
	  fprintf (stdout,"   nres/nrhs   %15.10le",norres/norrhs);
	if (norres/norrhs<err){
	  break;
	}
      }
      
      if (Mesprt != 0)  fprintf (stdout,"\n iteration number %ld",j);
      
      //  solution of the system of equations
      Tp->ssle->solve_system (Gtt,Kmat,v,rhs,Outt);
      
      //  correction of nodal values
      for (k=0;k<n;k++){
	lhs[k]+=v[k];
      }
      
      
      //  computation of time derivatives of the nodal values
      for (k=0;k<n;k++){
	tdlhs[k]=(lhs[k]-d[k])/dt/alpha;
      }
      
      //  computation of values at integration points from nodal values
      approximation ();
      
    }
    
    //  actual number of performed iterations
    ani=j;
    
    if (ani==0)
      nsts++;
    else
      nsts = 0;
    
    if (ani==ini){
      //  backup is used for actual values
      for (j=0;j<n;j++){
	lhs[j]=lhsb[j];
	tdlhs[j]=tdlhsb[j];
      }
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      i--;
      fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
	fprintf (stderr,"\n\n time increment is less than minimum time increment");
	fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
      }
    }else{
      //  time increment
      Tp->time = Tp->timecont.newtime ();
      dtdef = Tp->timecont.actualforwtimeincr ();
      
      if (nsts==2){
	dt*=2.0;
	nsts=0;
	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
      }
      if (dt>dtdef)
	dt=dtdef;
      
      time+=dt;
      Tp->time=time;
      Tp->timecont.time=time;
      
      //  time control
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
    }
    
  }while(Tp->time<end_time);
  
  delete [] v;
  delete [] f;
  delete [] p;
  delete [] d;
  delete [] lhsb;
  delete [] tdlhsb;
  
  print_closet();
}



/**
   function solves nonlinear nonstationary transport problem

   time discretization is based on the d-form version of
   the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   nonlinearity is solved by the full Newton-Raphson method
   
   JK, 23.4.2008
   JK, 8.4.2019
*/
void nonlinear_nonstat_solv_fnr_dform ()
{
  long i,j,k,n,ini,ani,nsts,lcid;
  long stop;
  double dt,dtmin,end_time,time,alpha,*err,norrhs,zero,*thresh,norfb;
  double *d,*f,*p,*v,*lhs,*lhsi,*tdlhs,*rhs,*lhsb,*tdlhsb;
  double *z,*diag;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  //  load case id, it must be equal to zero for this type of computation
  lcid=0;
  //  threshold for detection of zero diagonal matrix entries
  zero=1.0e-20;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  allocation of auxiliary arrays
  //  predictor (d with tilde)
  d = new double [n];
  //  auxiliary vector
  p = new double [n];
  //  prescribed fluxes
  f = new double [n];
  //  auxiliary vector
  v = new double [n];
  //  backup of nodal values
  lhsb = new double [n];
  //  backup of time derivatives of nodal values
  tdlhsb = new double [n];

  z = new double [n];
  diag = new double [n];
  

  //  initialization
  nullv (lhs,n);
  nullv (rhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (f,n);
  nullv (v,n);
  nullv (z,n);

  //  interpolation of nodal values to integration points
  approximation ();
  
  //  parameter of the trapezoidal method
  alpha=Tp->alpha;
  //  prescribed norm of residual
  err=Tp->errarr;
  //  maximum number of iterations in one time step
  ini=Tp->nii;
  //  minimum time increment
  dtmin=Tp->timecont.dtmin;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  time=Tp->time;
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();

  // vlozil JM 25.4.2008----*****
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation (); 
  
  
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  // -----*******
  
  

  // ***************************
  //  main iteration loop  ****
  // ***************************
  //  number of time step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  stop=0;


  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  
  
  do{
    
    Tp->time=Tp->timecont.newtime (dt);
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n time %le,   time step = %le,    increment number %ld",Tp->time,dt,i);
    
    Tm->updateipval ();
    i++;
    
    //  backup of attained nodal values and their time derivatives
    for (j=0;j<n;j++){
      lhsb[j]=lhs[j];
      tdlhsb[j]=tdlhs[j];
    }
    
    //  capacity matrix
    capacity_matrix (lcid);

    //fprintf (Outt,"\n matice kapacity");
    //Cmat->printmat (Outt);
    
    //  conductivity matrix
    conductivity_matrix (lcid);

    //fprintf (Outt,"\n matice vodivosti");
    //Kmat->printmat (Outt);
    
    //fprintf (Outt,"\n nejvetsi vl. cislo %le",Kmat->power_method (diag,1000,1.0e-6));
    //fprintf (Outt,"\n nejmensi vl. cislo %le",Kmat->inverse_iteration (diag,1000,1.0e-6));
    //conductivity_matrix (lcid);

    
    //  predictor d
    for (j=0;j<n;j++){
      d[j] = lhs[j]+dt*(1.0-alpha)*tdlhs[j];
    }
    
    //  C.d
    Cmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed fluxes f_{n+1}
    trfel_right_hand_side (lcid,f,n);
    
    //  C.d + alpha dt f_n+1
    norrhs=0.0;
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
	//  fprintf (stdout,"\n rhs   %15.10le   k je %d",rhs[j],j);
      //norrhs+=f[j]*f[j];
      //norrhs+=rhs[j]*rhs[j];
      //z[j]=rhs[j];
    }
    //norrhs=sqrt(norrhs);
    
    //Kmat->printmat (Outt);
    
    //  solution of the system of equations
    //  nodal values d_{n+1} are solved
    

    
    //Kmat->diag_check (zero,rhs);
    
    /*
    Kmat->diag_scale (diag);
    for (k=0;k<n;k++){
      rhs[k]*=diag[k];
    }
    */

    //Kmat->printmat (Outt);
    
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    
    /*
    for (k=0;k<n;k++){
      lhs[k]*=diag[k];
    }
    */

    //  computation of time derivatives of the nodal values
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    //  computation of values at integration points from nodal values
    // approximation ();
    // *************************
    // vlozil JM 25.4.2008
    
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();
    
    
    // nulleqother ();
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    // ***********************************
    
    //  iteration for equilibrium
    for (j=0;j<ini;j++){
      
      //  capacity matrix
      capacity_matrix (lcid);
      
      //  conductivity matrix
      conductivity_matrix (lcid);
      
      
      //  matrix of the system of equations
      //  C + alpha.dt.K
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      
      //  (C + alpha.dt.K).d
      Kmat->gmxv (lhs,v);
      
      //  C.d
      Cmat->gmxv (d,p);
      
      //  computation of residuals
      //norrhs=0.0;
      //norres=0.0;
      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k]-v[k];
	//	fprintf (stdout,"\n rhs   %15.10le   k je %d",rhs[k],k);
	//norres+=rhs[k]*rhs[k];
	z[k]=f[k]*alpha*dt+p[k];
	//norrhs+=z[k]*z[k];
	//fprintf (stdout,"\n rhs   %15.10le   z  %15.10le  v %15.10le",rhs[k],z[k],v[k]);
      }
      //norres=sqrt(norres);
      //norrhs=sqrt(norrhs);

      //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
      
      /*
      //  nodal energy
      nullv (z,n);
      nodal_energy (z,n,dt);
      double jkjk=0.0;
      for (k=0;k<n;k++){
	jkjk += z[i]*z[i];
      }
      fprintf (stdout,"\n jkjk   %15.10le",jkjk);      
      if (jkjk<1.0e-9){
	stop=1;
      }
      */
      
      stop = norm_computation_vec (rhs,z,err,thresh,2,1,norfb);

      //fprintf (stdout,"\n pauza \n");
      
      if (stop==1){
	break;
      }
      
      /*
      fprintf (stdout,"\n norres   %15.10le, nrhs   %15.10le,   norres/norrhs  %15.10le",norres,norrhs,norres/norrhs);
      
      if (norres/norrhs<err){
	break;
      }
      */

      /*
      if (Mesprt != 0)  
	fprintf (stdout,"\n nres   %15.10le, nrhs   %15.10le",norres,norrhs);
      
      if (norrhs<Tp->zero){
	if (norres<err){
	  break;
	}
      }
      else{
	if (Mesprt != 0)  
	  fprintf (stdout,"   nres/nrhs   %15.10le",norres/norrhs);
	if (norres/norrhs<err){
	  break;
	}
      }
      */
      

      if (Mesprt != 0)  fprintf (stdout,"\n iteration number %ld",j);
      

      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k];
	//z[k]=lhs[k];
      }

      
      //  solution of the system of equations
      //Tp->ssle->solve_system (Gtt,Kmat,v,rhs,Outt);
      
      
      Kmat->diag_check (zero,rhs);
      
      /*
      Kmat->diag_scale (diag);
      for (k=0;k<n;k++){
	rhs[k]*=diag[k];
      }
      */
      
      Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
      
      /*
      for (k=0;k<n;k++){
	lhs[k]*=diag[k];
      }
      */

      /*
      fprintf (Outt,"\n\n\n");
      for (k=0;k<n;k++){
	fprintf (Outt,"\n lhs %5ld   %20.15le",k,lhs[k]);
      }
      */
      
      //  correction of nodal values
      //for (k=0;k<n;k++){
      //lhs[k]+=v[k];
      //}
      
      
      //  computation of time derivatives of the nodal values
      for (k=0;k<n;k++){
	tdlhs[k]=(lhs[k]-d[k])/dt/alpha;
      }
      
      //  computation of values at integration points from nodal values
      // approximation ();
      
      // *************************
      // vlozil JM 25.4.2008
      
      //  physically corrected solution
      solution_correction ();    
      //  approximation of nodal values into ontegration points
      approximation ();
      
      
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	actual_previous_nodval ();
      // ***********************************
      
      //stop = norm_computation (lhs,z,err,2,1);
      
      //for (k=0;k<n;k++){
      //z[k]=lhs[k];
      //}

      //if (stop==1){
      //break;
      //}
    }
    
    //  actual number of perfromed iterations
    ani=j;
    
    if (ani==0)
      nsts++;
    else
      nsts=0;
    
    if (ani==ini){
      //  backup is used for actual values
      for (j=0;j<n;j++){
	lhs[j]=lhsb[j];
	tdlhs[j]=tdlhsb[j];
      }
      
      // *************************
      // vlozil JM 25.4.2008
      
      //  physically corrected solution
      solution_correction ();    
      //  approximation of nodal values into ontegration points
      approximation ();
      
      
      
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	actual_previous_nodval ();
      // ***********************************
      
      
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      i--;
      Tp->timecont.oldtime ();

      //fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium dt = %le",dt);
      //  fprintf (Outt2,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
	fprintf (stderr,"\n\n time increment is less than minimum time increment");
	fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
	abort ();
      }
    }else{
      //  time increment
      //Tp->time = Tp->timecont.newtime ();
      //dtdef = Tp->timecont.actualforwtimeincr ();
      
      if (nsts==2){
	dt*=2.0;
	//dt*=1.2;
	nsts=0;
	//	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	//fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le, dtdef = %le",dt,dtdef);
	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le",dt);
	//fprintf (Outt2,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt= %lf",dt);
      }
      //if (dt>dtdef)
      //dt=dtdef;
      
      //Tp->time = Tp->timecont.newtime (dt);
      
      
      //time+=dt;
      //Tp->time=time;
      //Tp->timecont.time=time;
      //Tp->timecont.fdt=dt;
      
      //  time control
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
    }
    
  }while(Tp->time<end_time);
  
  delete [] v;
  delete [] f;
  delete [] p;
  delete [] d;
  delete [] lhsb;
  delete [] tdlhsb;
  
  print_closet();
}



/**
   function solves nonlinear nonstationary transport problem

   time discretization is based on the d-form version of
   the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   nonlinearity is solved by the modified Newton-Raphson method
   
   JK, 21.10.2008
*/
void nonlinear_nonstat_solv_nr_dform ()
{
  long i,j,k,n,ini,ani,nsts,lcid;
  long stop;
  double dt,dtmin,end_time,time,alpha,*err,norrhs,zero,*thresh,norfb;
  double *d,*f,*p,*v,*lhs,*lhsi,*tdlhs,*rhs,*lhsb,*tdlhsb;
  double *z,*diag;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  //  load case id, it must be equal to zero for this type of computation
  lcid=0;
  //  threshold for detection of zero diagonal matrix entries
  zero=1.0e-20;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  allocation of auxiliary arrays
  //  predictor (d with tilde)
  d = new double [n];
  //  auxiliary vector
  p = new double [n];
  //  prescribed fluxes
  f = new double [n];
  //  auxiliary vector
  v = new double [n];
  //  backup of nodal values
  lhsb = new double [n];
  //  backup of time derivatives of nodal values
  tdlhsb = new double [n];

  z = new double [n];
  diag = new double [n];
  

  //  initialization
  nullv (lhs,n);
  nullv (rhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (f,n);
  nullv (v,n);
  nullv (z,n);

  //  interpolation of nodal values to integration points
  approximation ();
  
  //  parameter of the trapezoidal method
  alpha=Tp->alpha;
  //  prescribed norm of residual
  err=Tp->errarr;
  //  maximum number of iterations in one time step
  ini=Tp->nii;
  //  minimum time increment
  dtmin=Tp->timecont.dtmin;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  time=Tp->time;
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();

  // vlozil JM 25.4.2008----*****
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation (); 
  
  
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  // -----*******
  
  

  // ***************************
  //  main iteration loop  ****
  // ***************************
  //  number of time step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  stop=0;


  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  
  
  do{
    
    Tp->time=Tp->timecont.newtime (dt);
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n time %le,   time step = %le,    increment number %ld",Tp->time,dt,i);
    
    Tm->updateipval ();
    i++;
    
    //  backup of attained nodal values and their time derivatives
    for (j=0;j<n;j++){
      lhsb[j]=lhs[j];
      tdlhsb[j]=tdlhs[j];
    }
    
    //  capacity matrix
    capacity_matrix (lcid);

    //fprintf (Outt,"\n matice kapacity");
    //Cmat->printmat (Outt);
    
    //  conductivity matrix
    conductivity_matrix (lcid);

    //fprintf (Outt,"\n matice vodivosti");
    //Kmat->printmat (Outt);
    
    //fprintf (Outt,"\n nejvetsi vl. cislo %le",Kmat->power_method (diag,1000,1.0e-6));
    //fprintf (Outt,"\n nejmensi vl. cislo %le",Kmat->inverse_iteration (diag,1000,1.0e-6));
    //conductivity_matrix (lcid);

    
    //  predictor d
    for (j=0;j<n;j++){
      d[j] = lhs[j]+dt*(1.0-alpha)*tdlhs[j];
    }
    
    //  C.d
    Cmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed fluxes f_{n+1}
    trfel_right_hand_side (lcid,f,n);
    
    //  C.d + alpha dt f_n+1
    norrhs=0.0;
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
	//  fprintf (stdout,"\n rhs   %15.10le   k je %d",rhs[j],j);
      //norrhs+=f[j]*f[j];
      //norrhs+=rhs[j]*rhs[j];
      //z[j]=rhs[j];
    }
    //norrhs=sqrt(norrhs);
    
    //Kmat->printmat (Outt);
    
    //  solution of the system of equations
    //  nodal values d_{n+1} are solved
    

    
    Kmat->diag_check (zero,rhs);
    
    /*
    Kmat->diag_scale (diag);
    for (k=0;k<n;k++){
      rhs[k]*=diag[k];
    }
    */

    //Kmat->printmat (Outt);
    
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    
    /*
    for (k=0;k<n;k++){
      lhs[k]*=diag[k];
    }
    */

    //  computation of time derivatives of the nodal values
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    //  computation of values at integration points from nodal values
    // approximation ();
    // *************************
    // vlozil JM 25.4.2008
    
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();
    
    
    // nulleqother ();
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    // ***********************************
    
    //  iteration for equilibrium
    for (j=0;j<ini;j++){
      
      //  capacity matrix
      //capacity_matrix (lcid);
      
      //  conductivity matrix
      //conductivity_matrix (lcid);
      
      
      //  matrix of the system of equations
      //  C + alpha.dt.K
      //Kmat->scalgm (dt*alpha);
      //Kmat->addgm (1.0,*Cmat);
      
      //  (C + alpha.dt.K).d
      Kmat->gmxv (lhs,v);
      
      //  C.d
      //Cmat->gmxv (d,p);
      
      //  computation of residuals
      //norrhs=0.0;
      //norres=0.0;
      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k]-v[k];
	//	fprintf (stdout,"\n rhs   %15.10le   k je %d",rhs[k],k);
	//norres+=rhs[k]*rhs[k];
	z[k]=f[k]*alpha*dt+p[k];
	//norrhs+=z[k]*z[k];
	//fprintf (stdout,"\n rhs   %15.10le   z  %15.10le  v %15.10le",rhs[k],z[k],v[k]);
      }
      //norres=sqrt(norres);
      //norrhs=sqrt(norrhs);

      //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
      
      stop = norm_computation_vec (rhs,z,err,thresh,2,1,norfb);

      //fprintf (stdout,"\n pauza \n");
      
      if (stop==1){
	break;
      }
      
      /*
      fprintf (stdout,"\n norres   %15.10le, nrhs   %15.10le,   norres/norrhs  %15.10le",norres,norrhs,norres/norrhs);
      
      if (norres/norrhs<err){
	break;
      }
      */

      /*
      if (Mesprt != 0)  
	fprintf (stdout,"\n nres   %15.10le, nrhs   %15.10le",norres,norrhs);
      
      if (norrhs<Tp->zero){
	if (norres<err){
	  break;
	}
      }
      else{
	if (Mesprt != 0)  
	  fprintf (stdout,"   nres/nrhs   %15.10le",norres/norrhs);
	if (norres/norrhs<err){
	  break;
	}
      }
      */
      

      if (Mesprt != 0)  fprintf (stdout,"\n iteration number %ld",j);
      

      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k];
	//z[k]=lhs[k];
      }

      
      //  solution of the system of equations
      //Tp->ssle->solve_system (Gtt,Kmat,v,rhs,Outt);
      
      
      Kmat->diag_check (zero,rhs);
      
      /*
      Kmat->diag_scale (diag);
      for (k=0;k<n;k++){
	rhs[k]*=diag[k];
      }
      */
      
      Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
      
      /*
      for (k=0;k<n;k++){
	lhs[k]*=diag[k];
      }
      */

      /*
      fprintf (Outt,"\n\n\n");
      for (k=0;k<n;k++){
	fprintf (Outt,"\n lhs %5ld   %20.15le",k,lhs[k]);
      }
      */
      
      //  correction of nodal values
      //for (k=0;k<n;k++){
      //lhs[k]+=v[k];
      //}
      
      
      //  computation of time derivatives of the nodal values
      for (k=0;k<n;k++){
	tdlhs[k]=(lhs[k]-d[k])/dt/alpha;
      }
      
      //  computation of values at integration points from nodal values
      // approximation ();
      
      // *************************
      // vlozil JM 25.4.2008
      
      //  physically corrected solution
      solution_correction ();    
      //  approximation of nodal values into ontegration points
      approximation ();
      
      
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	actual_previous_nodval ();
      // ***********************************
      
      //stop = norm_computation (lhs,z,err,2,1);
      
      //for (k=0;k<n;k++){
      //z[k]=lhs[k];
      //}

      //if (stop==1){
      //break;
      //}
    }
    
    //  actual number of perfromed iterations
    ani=j;
    
    if (ani==0)
      nsts++;
    else
      nsts=0;
    
    if (ani==ini){
      //  backup is used for actual values
      for (j=0;j<n;j++){
	lhs[j]=lhsb[j];
	tdlhs[j]=tdlhsb[j];
      }
      
      // *************************
      // vlozil JM 25.4.2008
      
      //  physically corrected solution
      solution_correction ();    
      //  approximation of nodal values into ontegration points
      approximation ();
      
      
      
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	actual_previous_nodval ();
      // ***********************************
      
      
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      i--;
      Tp->timecont.oldtime ();

      //fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium dt = %le",dt);
      //  fprintf (Outt2,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
	fprintf (stderr,"\n\n time increment is less than minimum time increment");
	fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
	abort ();
      }
    }else{
      //  time increment
      //Tp->time = Tp->timecont.newtime ();
      //dtdef = Tp->timecont.actualforwtimeincr ();
      
      if (nsts==2){
	dt*=2.0;
	//dt*=1.2;
	//	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	//fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le, dtdef = %le",dt,dtdef);
	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le",dt);
	//fprintf (Outt2,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt= %lf",dt);
      }
      //if (dt>dtdef)
      //dt=dtdef;
      
      //Tp->time = Tp->timecont.newtime (dt);
      
      
      //time+=dt;
      //Tp->time=time;
      //Tp->timecont.time=time;
      //Tp->timecont.fdt=dt;
      
      //  time control
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
    }
    
  }while(Tp->time<end_time);
  
  delete [] v;
  delete [] f;
  delete [] p;
  delete [] d;
  delete [] lhsb;
  delete [] tdlhsb;
  
  print_closet();
}
