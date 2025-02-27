#include "mpi.h"
#include "npglvec.h"
#include "backupsolt.h"
#include "transprint.h"
#include "pnpsolvert.h"
#include "pglobalt.h"
#include "seqfilest.h"
#include <string.h>
#include <math.h>

/**
   function solves parallel nonstationary transport problem
   
   function selects a suitable numerical scheme for solution
   
   TKr, 19/02/2013
*/
void par_solve_nonstationary_problem ()
{
  if (Tp->tprob == discont_nonstat_problem){
    //  formulation based on nodal unknowns is used
    //par_linear_nonstat_solv_dform ();//not implemented
  }
  else{
    //  formulation based on unknown first nodal derivatives
    //  with respect to time is used
    par_nonstat_solv_vform_comp ();
    //par_linear_nonstat_solv_vform ();
    //par_linear_nonstat_solv_dform ();//not implemented
    //par_linear_nonstat_radiation_solv_dform ();//not implemented
  }
}


/**
   function solves parallel linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   according to JK, 20.12.2002 and TKo mtsolver.cpp
   revised 05/11/2012 by TKr
*/
void par_nonstat_solv_vform_comp ()
{
  long i,li,nsts,ret,n;
  double newtime,dt,dtmin,dtmax,dtdef,end_time, prev_time;
  double *lhsb,*tdlhsb;
  np_glob_vec np_gv;
  double *gv = NULL;
  double *grhs = NULL;
  long lcid = 0;
  long rest_calc;
  
  //
  //  initialization phase
  //
  par_nonstat_solver_init(lcid, rest_calc, np_gv);
  if(Tp->tprob == nonlinear_nonstationary_problem){
    // global vectors of coarse problem for used computation of selective norm of residual vector
    gv = new double[Psolt->schcom->ndofcp];
    grhs = new double[Psolt->schcom->ndofcp];
    nullv(gv, Psolt->schcom->ndofcp);
    nullv(grhs, Psolt->schcom->ndofcp);
  }
  
  // store initial value for attained time in the previous time step
  // (prev_time is exploited in the growing transport problems)
  prev_time = Tp->time;

  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  //  minimum time increment
  dtmin=Tp->timecont.dtmin;
  //  maximum time increment
  dtmax=Tp->timecont.dtmax;
  //  end time
  end_time = Tp->timecont.endtime ();
  
  //  number of step
  li = i = np_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;
  // return status from the one_step function
  ret = 0;
  
  // only for adaptive time increment
  if(Tp->timecont.tct > 0){
    n=Ndoft;
    lhsb   = new double [n];
    tdlhsb   = new double [n];
    nullv (lhsb,n);
    nullv (tdlhsb,n);
  }

  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    i++;
    
    newtime = Tp->time = Tp->timecont.newtime();
    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Tp->timecont.timefun.tfunc != constant) && Tp->timecont.tct == 0)
      dt = Tp->timecont.actualbacktimeincr ();
    
    if (((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)) &&
        (prev_time == Tp->timecont.starttime()))
      prev_time = Tp->time;

    if(Tp->tprob == nonstationary_problem)
      // linear algorithm
      // perform one time step with linear solver
      ret = par_one_step_linear(lcid, newtime, dt, i, li, np_gv);
    
    
    if(Tp->tprob == nonlinear_nonstationary_problem){
      // non-linear algorithm
      // perform one time step with non-linear solver
      ret = par_one_step_nonlinear(lcid, newtime, dt, prev_time, rest_calc, i, li, np_gv, gv, grhs);
    }
    
    if (ret >= 0){ // equilibrium was attained
	
      //  Backup of lhs is used in case that the Newton-Raphson procedure does not converge.
      copyv(np_gv.lhs, np_gv.lhsb, Ndoft);
      //  Backup of tdlhs is used in case that the Newton-Raphson procedure does not converge.
      copyv(np_gv.tdlhs, np_gv.tdlhsb, Ndoft);
	
      if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)) {
        // Following two commands must be performed here because in case that
        // elements or nodes are changed in the next step and the step would not converge,
        // the lhs and tdlhs could not be recovered correctly
        
        // backup of nodal values because code numbers may be changed
        Tt->lhs_save (np_gv.lhs,Lsrst->lhsi,np_gv.tdlhs);
        //  backup of nodal forces because code numbers may be changed//from mefel??!!
        //Mt->save_nodforce(mt_gv.f);//from mefel???!!!
        
        // new elements can be definitely changed to old ones
        Gtt->update_auxinf();
        prev_time = Tp->time;
      }
        
      //////////////////////////////////////////////////
      //  printing of output and graphical informations
      //  update of values stored at integration points
      //  approximation(); // it has been already called in the solver step procedures
      Tm->updateipval();
      
      compute_req_valt (lcid);    
      print_stept(lcid,Tp->istep,Tp->time,np_gv.rhs);
      print_flusht();
      
      if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat()) {
        if (Mesprt==1)
          fprintf (stdout,"\n Creating TRFEL backup file\n");
        solvert_save (np_gv.lhs,np_gv.tdlhs,np_gv.f,Tp->istep,Tp->time,dt,Tp->timecont,Ndoft);
      }
    }
    // handling with time controler for nonlinear solver
    if(Tp->tprob == nonlinear_nonstationary_problem && Tp->timecont.tct > 0){//only for adaptive time controler
      // only for adaptive time controler
      if (ret >= 0) {
        if (ret == 0)
          nsts++;      
        dtdef = Tp->timecont.actualforwtimeincr();
        if (nsts==2){
          dt*=2.0;
          nsts=0;
          if (Myrank==0 && Mesprt==1)  
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
        }
        if (dt<=dtdef)
          dt = dtdef;
	
        if (dt>dtmax)//maximum time increment
          dt = dtmax;
      }//  end of the if (ret >= 0)
      else{
        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
        dt/=2.0;
        Tp->timecont.oldtime ();
        i--;
        
        //restoring results from previous time step
        copyv(np_gv.lhsb, np_gv.lhs, Ndoft);
        copyv(np_gv.tdlhsb, np_gv.tdlhs, Ndoft);
	
        //  approximation of nodal values into ontegration points
        approximation ();
        Tm->updateipval();
        if (Myrank==0 && Mesprt==1)  
          fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
	
        if (dt<dtmin){
          if (Myrank==0 && Mesprt==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
          if (Myrank==0 && Mesprt==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
          if (Myrank==0 && Mesprt==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
          compute_req_valt (lcid);
          print_stept(lcid, i, Tp->time, np_gv.rhs);
          print_flusht();
          break;
        }
        // change status of elements and nodes according to the previous time step
        if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
          update_elemnod_prev_timet(prev_time, np_gv.ncd, np_gv.nce);
      }//  end of else
    }
    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
      prev_time = newtime;
    
  }while(Tp->time < end_time);
  
  print_closet();
  
  //Psolt->computation_stat_print (Outt);
  
  if(Tp->timecont.tct > 0){
    delete [] lhsb;
    delete [] tdlhsb;
  }
  delete [] gv;
  delete [] grhs;
}



/**
  Function allocates and initializes global vectors used in npsolver
  It is used in one_step concept solver.

  Created by tomas Krejci according to Tomas Koudelka, 02/2013
*/

void par_nonstat_solver_init (long lcid, long &rest_calc, np_glob_vec &np_gv)
{
  long n;
  double dt;

  //  number of unknowns
  n=Ndoft;
  
  np_gv.alloc(n);
  //  nodal values
  np_gv.lhs = Lsrst->give_lhs (lcid);
  //  time derivatives of nodal values
  np_gv.tdlhs = Lsrst->give_tdlhs (lcid);
  //  prescribed fluxes / right hand side
  np_gv.rhs = Lsrst->give_rhs (lcid);
  
  if(Tp->tprob == nonlinear_nonstationary_problem)
    np_gv.alloc_aux(n);//auxiliary vectors

  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();

  //  nodes - integration points interpolation - ???!!! is it necessary
  approximation();  
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (np_gv.lhs, np_gv.tdlhs, np_gv.f, np_gv.istep, Tp->time, dt, Tp->timecont, n);
    //  Backup of lhs is used in case that the Newton-Raphson procedure does not converge.
    copyv(np_gv.lhs, np_gv.lhsb, Ndoft);

    //  Backup of tdlhs is used in case that the Newton-Raphson procedure does not converge.
    copyv(np_gv.tdlhs, np_gv.tdlhsb, Ndoft);
    
    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)){
      update_elnod_stat_after_hdbrestt(np_gv);
    }

      //  initiation of transport material models
    Tm->initmaterialmodels();
    //  initiation of transport material models at aux. int. points
    Tm->aip_initmaterialmodels();
    compute_req_valt (lcid);
    // set indicator of calculation started from restored backup data
    rest_calc = 1;
    print_initt(-1, "at",Ptp->fni,Ptp->fei);
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    //  initiation of transport material models at aux. int. points
    Tm->aip_initmaterialmodels();
    compute_req_valt (lcid);
    print_initt(-1, "wt",Ptp->fni,Ptp->fei);
    print_stept(lcid, np_gv.istep, Tp->time, np_gv.rhs);
  }
  print_flusht();
}



/**
   function solves parallel linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   according to JK, 20.12.2002,
   revised 19/02/2013 by TKr
*/
long par_one_step_linear (long lcid,double time, double dt, long istep, long /*li*/, np_glob_vec &np_gv)
{
  long j,n;
  double alpha,*f,*d,*p,*lhs,*tdlhs,*rhs;
  
  //new time increment
  Tp->time = time;
  //  new step number
  Tp->istep = istep;
  
  //  nodal values
  lhs = np_gv.lhs;
  //  time derivatives of nodal values
  tdlhs = np_gv.tdlhs;
  //  prescribed fluxes / right hand side
  rhs = np_gv.rhs;
  //  vector of prescribed fluxes (right hand side)
  f = np_gv.f;
  //  predictor
  d = np_gv.d;
  //  auxiliary vector
  p = np_gv.p;
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;

  //  number of transport degrees of freedom
  n=Ndoft;

  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------");
  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n PTRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n -------------------------------------------------------------------\n");
  
  // update material properties and auxiliary values
  Tm->updateipval ();
  
  //  capacity matrix    
  capacity_matrix (lcid);
  
  //  conductivity matrix
  conductivity_matrix (lcid);
  
  if (istep == 0){
    Psolt->computation_stat (Gtt,Kmat);
  }

  //  predictor
  //  dd = d_n + (1-alpha) dt v_n
  for (j=0;j<n;j++){
    d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
  }
  
  //  auxiliary vector
  //  K dd
  Kmat->gmxv (d,p);
  
  //  matrix of the system of equations
  //  C + K alpha dt
  Kmat->scalgm (dt*alpha);
  Kmat->addgm (1.0,*Cmat);
  
  //  prescribed nodal fluxes f_{n+1}
  trfel_right_hand_side (0,rhs,n);
  
  //  computation of the right hand side vector
  //  f_{n+1} - K dd
  for (j=0;j<n;j++){
    f[j] = rhs[j] - p[j];
  }
  
  //  solution of the system of algebraic equations
  //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
  //  time derivatives v_{n+1} are solved
  Psolt->par_linear_solver (Gtt,Kmat,tdlhs,f,Outt,Mesprt);
  
  //  nodal values computed from nodal derivatives
  //  d_{n+1} = dd + alpha dt v_{n+1}
  for (j=0;j<n;j++){
    lhs[j] = d[j] + alpha*dt*tdlhs[j];
  }
    
  solution_correction ();
  //  nodes - integration points interpolation
  approximation ();
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  if (j == 0)
    return 0;
  
  return 1;
}

/**
   function solves parallel non-linear nonstationary transport problem  by Newton-Raphson method
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   according to TKo
   revised 05/11/2012 by TKr
*/
long par_one_step_nonlinear (long lcid,double time, double dt, double prev_time, long /*rest_calc*/, long istep, long /*li*/, np_glob_vec &np_gv, double *gv, double *grhs)
{
  long j, k, n, ini, tncd, tnce, stop;
  double alpha, *f, *d, *p, *lhs, *tdlhs, *rhs;
  double *fb, *fi;
  double norf_last;
  double zero, norfb, *err, *thresh;  
  
  //new time increment
  Tp->time = time;
  //  new step number
  Tp->istep = istep;
  Tp->jstep = -1;
  
  //  number of transport degrees of freedom
  n=Ndoft;

  //  nodal values
  lhs = np_gv.lhs;
  //  time derivatives of nodal values
  tdlhs = np_gv.tdlhs;
  //  prescribed fluxes / right hand side
  rhs = np_gv.rhs;
  //  vector of prescribed fluxes (right hand side)
  f = np_gv.f;
  //  predictor
  d = np_gv.d;
  //  auxiliary vector
  p = np_gv.p;

  // auxiliary vectors:
  fb = np_gv.fb;
  fi = np_gv.fi;
  
  //  initial values
  nullv (fb,n);
  nullv (fi,n);
  nullv (f,n);

  
  //  maximum number of iterations in inner loop
  ini = Tp->nii;
  //  required norm of vector of unbalanced forces
  err = Tp->errarr;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;

  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  //  computer zero
  zero=Tp->zero;

  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)){
    update_elnod_statt(lcid, istep, prev_time, lhs, tdlhs, tnce, tncd);
    np_gv.ncd = tncd;
    np_gv.nce = tnce;
    n = Ndoft;
  }

  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------");
  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n PTRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n -------------------------------------------------------------------\n");
  
  //  number of transport degrees of freedom
  n=Ndoft;
  
  // update material properties and auxiliary values
  Tm->updateipval ();
  
  //  capacity matrix    
  capacity_matrix (lcid);
  
  //  conductivity matrix
  conductivity_matrix (lcid);

  if (istep == 0){
    Psolt->computation_stat (Gtt,Kmat);
  }
  
  //  predictor
  //  dd = d_n + (1-alpha) dt v_n
  for (j=0;j<n;j++){
    d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
  }
  
  //  auxiliary vector
  //  K dd
  Kmat->gmxv (d,p);
  
  //  matrix of the system of equations
  //  C + K alpha dt
  // old:
  //Kmat->scalgm (dt*alpha);
  //Kmat->addgm (1.0,*Cmat);
  //Kmat->copygm (*Cmat);
  // new:
  Cmat->addgm (dt*alpha,*Kmat);

  //  prescribed nodal fluxes f_{n+1}
  trfel_right_hand_side (0,rhs,n);
  
  //mechanical influence of transport problem
  trfel_right_hand_side2 (0,rhs,n);

  //  computation of the right hand side vector
  //  f_{n+1} - K dd
  for (j=0;j<n;j++){
    f[j] = rhs[j] - p[j];
  }
  
  //  solution of the system of algebraic equations
  //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
  //  time derivatives v_{n+1} are solved
  Psolt->par_linear_solver (Gtt,Cmat,tdlhs,f,Outt,Mesprt);
  
  //  nodal values computed from nodal derivatives
  //  d_{n+1} = dd + alpha dt v_{n+1}
  for (j=0;j<n;j++){
    lhs[j] = d[j] + alpha*dt*tdlhs[j];
  }

  //  physically corrected solution
  solution_correction ();
  //  nodes - integration points interpolation
  approximation ();
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  //////////////////////////////////////////////////
  stop = 0;
  norf_last = 1.0e20;
  //inner iteration loop//////////////////
  for (j=0;j<ini;j++){
    Tp->jstep = j;    
    // full newton-raphson
    if (Tp->trsolv == fullnewtont){
      // matrices are computing in each inner iteration
      // capacity matrix
      capacity_matrix (0);
      
      //  conductivity matrix
      conductivity_matrix (0);
      
      //  auxiliary vector  K (d+(1-alpha)*dt*v)
      Kmat->gmxv (d,p);
      
      //  matrix of the system of equations
      //  C + alpha.dt.K
      //old:
      //Kmat->scalgm (dt*alpha);
      //Kmat->addgm (1.0,*Cmat);
      //Kmat->copygm (*Cmat);
      //new:
      Cmat->addgm (dt*alpha,*Kmat);
    }
    
    if (Tp->trestype==fluxest){
      // Solver computes unbalanced fluxes
      internal_fluxes (fi,n);//new fi vector
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
	fb[k]=-fi[k];
      }
    }

    if (Tp->trestype==lrhst){
      //  correct right hand side of transport part
      trfel_right_hand_side (lcid,rhs,n);
      //mechanical influence of transport problem
      trfel_right_hand_side2 (0,rhs,n);
      
      // Solver computes residuum from system of equations
      if (Tp->trsolv == modnewtont)
	capacity_matrix (0);
      // Solver computes residuum from system of equations
      // old:
      // Kmat->gmxv (tdlhs,fi);//new fi vector
      // new:
      Cmat->gmxv (tdlhs,fi);//new fi vector
      // Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
      for (k=0;k<n;k++){
	fb[k] = rhs[k] - p[k] - fi[k];
      }
    }
    
    // compute selective norm of residual scaled by rhs
    stop = par_norm_computation_vec (fb, rhs, err, thresh, 2, 1, gv, grhs, norfb);
    if (stop) break;
    
    if (Myrank==0 && Mesprt==1)   fprintf (stdout,"\n inner iteration %ld   error %e", j, norfb);
    
    Psolt->par_linear_solver (Gtt,Cmat,fi,fb,Outt,Mesprt);//fi is output now
    
    for (k=0;k<n;k++){
      tdlhs[k]+=fi[k];
      lhs[k]+=alpha*dt*fi[k];
    }
    
    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
      Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

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
      if (norfb > norf_last){
		
	if (Myrank==0 && Mesprt==1)  
	  fprintf (stdout,"\n\n convergence control: inner iteration skiped %ld error %e\n\n", j, norfb);	  
	
	return -2;//break;
      }
            
      norf_last = norfb; //storing norf from previous inner step
    }
  }  

  // old:
  //if(norf >= err[0])
  //  return -1;
  // new:
  if(stop == 0)
    return -1;

  if (Tp->jstep > 1)
    return 1;
  else
    return 0;

}


/**
   function solves linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   JK
*/
void par_linear_nonstat_solv_vform ()
{
  long i,j,n;
  double s,zero,alpha,*d,*p,*lhs,*lhsi,*tdlhs,*rhs;

  double dt,end_time;
  

  n=Ndoft;
  
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);

  d = new double [n];
  p = new double [n];
  //tdlhs = new double [n];

  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  


  
  alpha=Tp->alpha;

  zero=Tp->zero;
  



  // **************************
  //  main iteration loop  ***
  // **************************

  //  loop counter
  i=0;
  
  //  starting time
  Tp->time = Tp->timecont.starttime ();
  //  time increment
  dt = Tp->timecont.initialtimeincr ();
  //  end time
  end_time = Tp->timecont.endtime ();
  

  //  nodes - integration points interpolation
  approximation ();
  actual_previous_change ();

  Tm->initmaterialmodels();
  

  print_initt(-1, "wt",Ptp->fni,Ptp->fei);
  print_stept(0,i,Tp->time,NULL);
  print_flusht();


  do{
    

    if (Myrank==0) fprintf (stdout,"\n iteration number %ld",i);
    
    //fprintf (Outt,"\n\n\n\n\n iteration number %ld",i);
    
    //par_aux_nonstat_print (gr,lhsi,d,time);
    
    
    //fprintf (Out,"\n\n\n\n\n\n kontrola integracnich bodu");
    //for (long ijk=0;ijk<Tm->nip;ijk++){
    //fprintf (Out,"\n %ld  %lf %lf",ijk,Tm->ip[ijk].av[0],Tm->ip[ijk].av[1]);
    //}

    
    
    conductivity_matrix (0);

    capacity_matrix (0);
    
    if (i==0){
      Psolt->computation_stat (Gtt,Kmat);
    }
    
    //Kmat->printdiag (Outt);
    //fprintf (Outt,"\n\n MATICE K \n");
    //Kmat->printmat (Outt);
    //fprintf (Outt,"\n\n MATICE C \n");
    //Cmat->printmat (Outt);

    for (j=0;j<n;j++){
      p[j]=lhs[j]+(1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    Kmat->gmxv (p,d);
    
    

    //fprintf (Outt,"\n\n\n\n\n kontrola K.(d+dt*(1-a)*v)");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n %ld %le",ijk,d[ijk]);
    //}
    

    
    Kmat->scalgm (alpha*dt);
    Kmat->addgm (1.0,*Cmat);
    
    
    
    //fprintf (Out,"\n\n\n\n\n kontrola d+dt*(1-a)*v");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Out,"\n %ld %lf",ijk,p[ijk]);
    //}

    
    
    
    //fprintf (Out,"\n\n\n\n\n kontrola K.(d+dt*(1-a)*v)");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Out,"\n %ld %lf",ijk,r[ijk]);
    //}

    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //fprintf (Outt,"\n\n\n\n\n kontrola prave strany");
    for (long ijk=0;ijk<n;ijk++){
      fprintf (Outt,"\n rhs %ld %le",ijk,rhs[ijk]);
    }
    
    
    
    for (j=0;j<n;j++){
      rhs[j] = rhs[j] - d[j];
      d[j]=tdlhs[j];
    }
    
    
    
    //fprintf (Outt,"\n\n\n\n\n kontrola prave strany");
    for (long ijk=0;ijk<n;ijk++){
      fprintf (Outt,"\n rhs2 %ld %le",ijk,rhs[ijk]);
    }
    
 
    //fprintf (Outt,"\n\n MATICE SOUSTAVY \n");
    Kmat->printmat (Outt);




    //par_linear_solver ();
    Psolt->par_linear_solver (Gtt,Kmat,tdlhs,rhs,Outt,Mesprt);
    



    //fprintf (Outt,"\n\n MATICE SOUSTAVY PO RESENI \n");
    //Kmat->printmat (Outt);
    
    //fprintf (Outt,"\n\n\n\n\n kontrola reseni soustavy");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n tdlhs %ld      %le",ijk,tdlhs[ijk]);
    //}
    
    
    
    for (j=0;j<n;j++){
      s=(1.0-alpha)*d[j]+alpha*tdlhs[j];
      lhs[j]+=dt*s;
    }
    
    
    //fprintf (Outt,"\n\n\n\n\n kontrola obnovenych promennych");
    for (long ijk=0;ijk<n;ijk++){
      fprintf (Outt,"\n lhs %ld   %le",ijk,lhs[ijk]);
    }

    //  nodes - integration points interpolation
    approximation ();
    


    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

    printf ("\n Tp->time   %e",Tp->time);
    printf ("\n i          %ld",i);
    
  }while(Tp->time<end_time);



  delete [] p;  delete [] d;
  //fclose (gr); 

  print_closet();
  
  Psolt->computation_stat_print (Outt);
  
  
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
   @param gf[in] - auxiliary global vector for storage values of f of coarse problem, dimension must be ndofcp (see schurcompl)
   @param grhs[in] - auxiliary global vector for storage values of rhs of coarse problem, dimension must be ndofcp (see schurcompl)
   @param norfb - computed norm of the %vector lv, i.e. residual

   Created by T. Koudelka, 08.2021, according to JK, 24.8.2006, 
*/
long par_norm_computation_vec (double *f,double *rhs,double *err,double *thresh, long nt,long sc, double *gf, double *grhs, double &norfb)
{
  long i, stop, pstop;
  double norf,norrhs;
  stop = 0L;

  if (nt==1){
    // *****************************
    //  Euclidean norm is computed
    // *****************************
    norf = Psolt->pss(f, f, Outt);
    norfb = norf = sqrt(norf);
    
    if (sc==0){
      //  norm is not scaled
      if (norf<err[0])
        stop=1;
      else
        stop=0;
      
      if ((Mesprt==1) && (Myrank == 0))  fprintf (stdout,"\n norm of residuum  %e,    required error  %e",norf,err[0]);
    }
    
    if (sc==1){
      //  norm is scaled
      norrhs = Psolt->pss (rhs, rhs, Outt);
      norrhs=sqrt(norrhs);
      
      if (norrhs>thresh[0]){
        if (norf/norrhs<err[0])
          stop=1;
        else
          stop=0;
        norfb /= norrhs;	
        if ((Mesprt==1) && (Myrank == 0))  fprintf (stdout,"\n norm of residuum  %e,    norm of rhs  %e,   error  %e,   required error  %e",norf,norrhs,norf/norrhs,err[0]);
      }
      else{
        if (norf<err[0])
          stop=1;
        else
          stop=0;
        if ((Mesprt==1) && (Myrank == 0)) fprintf (stdout,"\n norm of residuum  %e,    norm of rhs  %e,  required error  %e",norf,norrhs,err[0]);
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
      pstop += Psolt->selected_norm_calculation(i, err[i], thresh[i], 1, Gtt, f, gf, rhs, grhs, Ndoft, Outt, norf);
      norfb += sqr(norf);      
    }
    
    if (pstop==Tp->ntm)
      stop=1;
    else
      stop=0;
    norfb = sqrt(norfb);
  }
  
  return stop;
}

