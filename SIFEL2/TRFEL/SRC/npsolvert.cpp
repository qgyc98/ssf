#include "npglvec.h"
#include "npsolvert.h"
#include "nnpsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "elemswitcht.h"
#include "transprint.h"
#include "saddpoint.h"
#include "backupsolt.h"
#include <string.h>
#include <errno.h>

/**
   function solves nonstationary transport problem
   
   function selects a suitable numerical scheme for solution
   
   JK
*/
void solve_nonstationary_problem ()
{
  if (Tp->tprob == discont_nonstat_problem){
    //  formulation based on nodal unknowns is used
    linear_nonstat_solv_dform ();
  }
  else{
    
    switch (Tp->tnpsolver){
    case trapezoidd:{
      linear_nonstat_solv_dform ();
      break;
    }
    case trapezoidv:{
      nonstat_solv_vform_comp ();  
      break;
    }
    case optim_trapezoidd:{
      optim_linear_nonstat_solv_dform ();
      break;
    }
    case forwardfindif:{
      explicit_difference_method (0);
      break;
    }
    default:{
      print_err("unknown type of nonstationary solver is required", __FILE__, __LINE__, __func__);      
    }
    }
    
  }
}

/**
   function solves linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   according to JK, 20.12.2002 and TKo mtsolver.cpp
   revised 05/11/2012 by TKr
*/
void nonstat_solv_vform_comp ()
{
  long i,li,nsts,ret,rest_calc;
  double newtime,dt,dtmin,dtmax,dtdef,end_time,prev_time;
  np_glob_vec np_gv;
  long lcid = 0;

  //
  //  initialization phase
  //
  //  nodes - integration points interpolation
  approximation();  
  nonstat_solver_init(lcid, rest_calc, np_gv);

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

  bool breakloop = false;
  
  //saving values for the first iteration
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    Tt->lhs_save (np_gv.lhs,Lsrst->lhsi,np_gv.tdlhs);

  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    i++;
  
    Tp->ipvcomp = 0;
    newtime = Tp->time = Tp->timecont.newtime(dt);
    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Tp->timecont.timefun.tfunc != constant) && Tp->timecont.tct == 0)
      dt = Tp->timecont.actualbacktimeincr ();
   
    if (((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)) &&
        (prev_time == Tp->timecont.starttime()))
      prev_time = Tp->time;

    if((Tp->tprob == nonstationary_problem) || (Tp->tprob == growing_np_problem))
      // linear algorithm
      // perform one time step with linear solver
      ret = one_step_linear(lcid, newtime, dt, prev_time, rest_calc, i, li, np_gv);
    
    if((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin)){
      // non-linear algorithm
      // perform one time step with non-linear solver
      ret = one_step_nonlinear(lcid, newtime, dt, prev_time, rest_calc, i, li, np_gv);
    }
    
    if (ret >= 0) { // equilibrium was attained
      /*
      //  calculate correct rates of lhs
      subv(np_gv.lhs, np_gv.lhsb, np_gv.tdlhs, Ndoft);
      cmulv(1.0/dt, np_gv.tdlhs, Ndoft);*/

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

    // handling with time controler for nonlinear solver, only for adaptive time increment
    if(((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin)) && Tp->timecont.tct > 0) {
      // only for adaptive time controler
      if (ret >= 0) {

        if (ret == 0)
          nsts++;      
        else
          nsts=0;

        dtdef = Tp->timecont.actualforwtimeincr();
        if (nsts==2){
          dt*=2.0;
          nsts=0;
		
          if (Mesprt==1)  
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
        }
        if (dt<=dtdef)//increase of time increment according to prescribed one
          dt = dtdef;
	    
        if (dt>dtmax)//maximum time increment
          dt = dtmax;

	//surface_fluxes (Outt); //printing in TRFEL part moved into outdriver
      }
      else {
        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
            nsts = 0;
        dt/=2.0;
        Tp->timecont.oldtime ();
        i--;
	    
        //restoring results from previous time step
        copyv(np_gv.lhsb, np_gv.lhs, Ndoft);
        copyv(np_gv.tdlhsb, np_gv.tdlhs, Ndoft);

        //  approximation of nodal values into ontegration points
        approximation ();
        Tm->updateipval();
	    
        if (Mesprt==1)  
          fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
	    
        if (dt<dtmin) {
          if (Mesprt==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
          if (Mesprt==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
          if (Mesprt==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
          compute_req_valt (lcid);
          print_stept_forced(lcid, i, Tp->time, np_gv.rhs);
          print_flusht();
          break;
        }

        // change status of elements and nodes according to the previous time step
        if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
          update_elemnod_prev_timet(prev_time, np_gv.ncd, np_gv.nce);
      }
    }
    
    // adaptivity
    if (Tp->adaptivityflag  &&  Tp->time < end_time)
      if (i%2 == 0)
	breakloop = Adat->run (2, true);
    
    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
      prev_time = newtime;
    
  }while(Tp->time < end_time  &&  breakloop == false);
    
  print_closet();
}


/**
  Function allocates and initializes global vectors used in npsolver
  It is used in one_step concept solver.

  @param lcid      - load case id
  @param rest_calc - indicator of calculation restorage from backup 
                     (0=no restorage was required, 1=restorage was performed)
  @param np_gv     - structure with global vectors used in the solver

  Created by Tomas Krejci according to Tomas Koudelka, 11/2012
  Modified by Tomas Krejci 02/2017
*/

void nonstat_solver_init (long lcid, long &rest_calc, np_glob_vec &np_gv)
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
  
  if((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin))
    np_gv.alloc_aux(n);//auxiliary vectors

  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();

  ///  termitovo
  if ((Adat != NULL) && (Adat->give_step()>0)) {
    Adat->statedata_restore();
    // timectrl ??
    // i = Tp->istep;
    print_initt(-1, "at");  // ??
  }
  else {
    //  inital printing of output and graphical informations
    if (Tp->hdbcont.restore_stat()){
      if (Mesprt==1)
	fprintf (stdout,"\n Reading of TRFEL backup file\n");
      solvert_restore (np_gv.lhs, np_gv.tdlhs, np_gv.f, np_gv.istep, Tp->time, dt, Tp->timecont, n);
 
      //  Backup of lhs is used in case that the Newton-Raphson procedure does not converge.
      copyv(np_gv.lhs, np_gv.lhsb, Ndoft);
      //  Backup of tdlhs is used in case that the Newton-Raphson procedure does not converge.
      copyv(np_gv.tdlhs, np_gv.tdlhsb, Ndoft);
      
      if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
	{
	  update_elnod_stat_after_hdbrestt(np_gv);
	}
      
      //  initiation of transport material models at regular int. points
      Tm->initmaterialmodels();
      //  initiation of transport material models at aux. int. points
      Tm->aip_initmaterialmodels();
      compute_req_valt (lcid);
      // set indicator of calculation started from restored backup data
      rest_calc = 1;
      print_initt(-1, "at");
    }
    else{
      //  initiation of transport material models at regular int. points
      Tm->initmaterialmodels();
      //  initiation of transport material models at aux. int. points
      Tm->aip_initmaterialmodels();
      compute_req_valt (lcid);
      print_initt(-1, "wt");
      print_stept(lcid, np_gv.istep, Tp->time, np_gv.rhs);
    }
    print_flusht();
  }
  
  /// termitovo kontrola transferu state variables  DOCASNE
  if (Adat != NULL){
    if (Adat->give_step()) {
      Adat->run (2, true);
      Adat->answer = 0;
    }
  }
}



/**
   function solves linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   @param lcid - load case id
   @param time - actual time
   @param dt   - actual time increment
   @param dtr   - the minimum ratio of required time step size to the actual one required by material models (output)
   @param prev_time - attained time in the previous time step (needed only for growing structures)
   @param rest_calc - inidicator whether step follows restorage from the backup immediately (=1) or no (=0)
   @param istep - time step id
   @param li    - time step id of the first performed time step
   @param np_gv - structure with pointers to %vectors of righthand and lefthand side
   
   @retval  0 - if the equilibrium was attained with NO inner loop
   
   according to JK, 20.12.2002,
   revised 05/11/2012 by TKr
   modified 02/2017 by TKr
*/
long one_step_linear (long lcid,double time, double dt, double prev_time, long /*rest_calc*/, long istep, long /*li*/, np_glob_vec &np_gv)
{
  long j,n, tncd, tnce;
  double alpha,*f,*d,*p,*lhs,*tdlhs,*rhs;
  
  //new time increment
  Tp->time = time;
  //  new step number
  Tp->istep = istep;
  Tp->jstep = -1;
  
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

  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    {
      update_elnod_statt(lcid, istep, prev_time, lhs, tdlhs, tnce, tncd);
      np_gv.ncd = tncd;
      np_gv.nce = tnce;
      n = Ndoft;
    }
  
  if (Mesprt==1){
    fprintf (stdout,"\n\n ------------------------------------------------------------------------");
    fprintf (stdout,"\n TRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
    fprintf (stdout,"\n ------------------------------------------------------------------------\n");
  }

  // update material properties and auxiliary values
  Tm->updateipval ();

  //  capacity matrix    
  capacity_matrix (lcid);

  //  conductivity matrix
  conductivity_matrix (lcid);
 
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
  //for debug printing conductivity matrix:

  Tp->ssle->solve_system (Gtt,Kmat,tdlhs,f,Outt);
 
  //  nodal values computed from nodal derivatives
  //  d_{n+1} = dd + alpha dt v_{n+1}
  for (j=0;j<n;j++){
    lhs[j] = d[j] + alpha*dt*tdlhs[j];
  }
 
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

  solution_correction ();
  //  nodes - integration points interpolation
  approximation ();
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();

  return 0;
}



/**
   function solves non-linear nonstationary transport problem  by Newton-Raphson method
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   according to TKo
   revised 05/11/2012 by TKr
*/
long one_step_nonlinear (long lcid,double time, double dt,  double prev_time, long /*rest_calc*/, long istep, long /*li*/, np_glob_vec &np_gv)
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

  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    {
      update_elnod_statt(lcid, istep, prev_time, lhs, tdlhs, tnce, tncd);
      np_gv.ncd = tncd;
      np_gv.nce = tnce;
      n = Ndoft;
    }

  if (Mesprt==1)  fprintf (stdout,"\n\n --------------------------------------------------------------");
  if (Mesprt==1)  fprintf (stdout,"\n TRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
  if (Mesprt==1)  fprintf (stdout,"\n --------------------------------------------------------------\n");
  
  //fprintf (Outt,"\n\n --------------------------------------------------------------");
  //fprintf (Outt,"\n TRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
  //fprintf (Outt,"\n --------------------------------------------------------------\n");
  //fflush(Outt);

  //  number of transport degrees of freedom
  n=Ndoft;
  
  // update material properties and auxiliary values
  Tm->updateipval ();
  
  //  capacity matrix    
  capacity_matrix (lcid);
  //fprintf(Outt, "\nCapacity matrix at istep=%ld, jstep=%ld:\n", Tp->istep, Tp->jstep);
  //Cmat->printmat(Outt);
  
  //  conductivity matrix
  conductivity_matrix (lcid);
  //fprintf(Outt, "\nConductivity matrix at istep=%ld, jstep=%ld:\n", Tp->istep, Tp->jstep);
  //Kmat->printmat(Outt);
  
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
  //old:
  //Kmat->scalgm (dt*alpha);
  //Kmat->addgm (1.0,*Cmat);
  //Kmat->copygm (*Cmat);
  //new:
  Cmat->addgm (dt*alpha,*Kmat);
  //fprintf(Outt, "\nOverall matrix at istep=%ld, jstep=%ld:\n", Tp->istep, Tp->jstep);
  //Cmat->printmat(Outt);

  //  prescribed nodal fluxes f_{n+1}
  trfel_right_hand_side (0,rhs,n);
  
  //fprintf(Outt, "\nRight-hand side vector:\n");
  //  fprintf(Outt, "%e",Tp->time);
  //for (j=0;j<n;j++){
  //fprintf(Outt, " %e",rhs[j]);
  //}
  //fprintf(Outt, "\n");
  //fflush(Outt);  

  //mechanical influence of transport problem
  trfel_right_hand_side2 (0,rhs,n);
  
  //fprintf(Outt, "\nRight-hand2 side vector:\n");
  //fprintf(Outt, "%e",Tp->time);
  //for (j=0;j<n;j++){
  //fprintf(Outt, " %e",rhs[j]);
  //}
  //fprintf(Outt, "\n");
  //fflush(Outt);  

  //  computation of the right hand side vector
  //  f_{n+1} - K dd
  for (j=0;j<n;j++){
    f[j] = rhs[j] - p[j];
  }

  //  solution of the system of algebraic equations
  //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
  //  time derivatives v_{n+1} are solved
  Tp->ssle->solve_system (Gtt,Cmat,tdlhs,f,Outt);

  check_math_err();//debug??!!
  
  //  nodal values computed from nodal derivatives
  //  d_{n+1} = dd + alpha dt v_{n+1}
  for (j=0;j<n;j++){
    lhs[j] = d[j] + alpha*dt*tdlhs[j];
  }

  // arrays of actual values are moved to arrays of previous values
  actual_previous_change ();

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
      Cmat->gmxv (tdlhs,fi);//new fi vector
      // Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
      for (k=0;k<n;k++){
	fb[k] = rhs[k] - p[k] - fi[k];
      }
    }

    // compute selective norm of residual scaled by rhs
    stop = norm_computation_vec (fb,rhs,err,thresh,2,1,norfb);
    if (stop) break;
    
    if (Mesprt==1)  fprintf (stdout,"\n inner iteration %ld   error %.15e", j, norfb);

    Tp->ssle->solve_system (Gtt,Cmat,fi,fb,Outt);//fi is output now
    
    for (k=0;k<n;k++){
      tdlhs[k]+=fi[k];
      lhs[k]+=alpha*dt*fi[k];
    }
    

    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
      Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

    //  physically corrected solution
    solution_correction ();
    //  approximation of nodal values into integration points
    approximation ();   
    
    // vlozil JM 15.1.2008
    // nulleqother ();
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    //  convergence control for newton-raphson
    //  condition of decreesing of error
    if (Tp->convergcontrolt==yes){
      if (norfb > norf_last){
	if (Mesprt==1)  fprintf (stdout,"\n\n divergence control: inner iteration skiped %ld error %e\n\n",j,norfb);	  
	return -2;//break;
      }
      norf_last = norfb; //storing norf from previous inner step
    }
    // end of convergence control for newton-raphson
  }

  if(stop == 0)
    return -1;

  //debug??!!
  if (Tp->jstep > 1)
    return 1;
  else
    return 0;
}

/**
  The function updates status of all elements and nodes according 
  to their time functions after restorage from the backup file.
  The function is being used in the growing mechanical problems.

  @param lcid - load case id

  Created by Tomas Krejci according to Tomas Koudelka 05/02/2017
*/
void update_elnod_stat_after_hdbrestt(np_glob_vec &np_gv)
{
  //  marking of elements switched on
  Gtt->update_elem(Tp->time);
  //  at Mp->time=prev_timem must be auxinf=leso
  Gtt->update_auxinf();
  //  marking of nodes switched on
  Gtt->update_nodes();
  //  marking of DOFs switched on
  Gtt->update_dofs(Tp->time);
  // new DOF numbering
  Ndoft = Gtt->codenum_generation(Outt);
  // save nodal values 
  Tt->lhs_save (np_gv.lhs,Lsrst->lhsi,np_gv.tdlhs);
}



/**
  Regenerate DOF numbers according to old time step if there were changes
  in nodes or elements.

  @param prev_time - attained time from previous time step
  @param ncd - the number of changed DOFs in the actual time step
  @param nce - the number of changed elements in the actual time step

  Created by Tomas Krejci according to Tomas Koudelka 05/02/2017
*/
void update_elemnod_prev_timet(double prev_time, long ncd, long nce)
{
  //    
  // regenerate code numbers according to old time step if there were changes in nodes or elements
  //
  if ((nce!=0) || (ncd!=0))
  {
    // switch on elements accroding to previous time step
    Gtt->update_elem(prev_time);
    // switch on nodes accroding to previous time step     
    Gtt->update_nodes();
    //  marking of DOFs switched on
    Gtt->update_dofs (prev_time);
    // old DOF numbering
    Ndoft = Gtt->codenum_generation(Outt);
    // destruct old system matrix
    delete Kmat;
    Kmat=NULL;
    delete Cmat;
    Cmat=NULL;
  }
}



/**
  The function updates status of all elements and nodes according 
  to their time functions.

  @param lcid - load case id
  @param istep - id of the actual time step
  @param prev_time - attained time in the previous time step

  @param lhs    - lhs vector
  @param tdlhs  - tdlhs vector

  @param tnce - the number of elements with changed status (added/removed)
  @param tncd - the number of DOFs with changed status (added/removed)

  Created by Tomas Koudelka, 31.1.2017
*/
void update_elnod_statt(long /*lcid*/, long /*istep*/, double prev_time, double *lhs, double *tdlhs, long &tnce, long &tncd)
{
  //  searching for new elements
  //  only new elements are switched on
  //  former elements are temporarily switched off
  tnce = Gtt->search_newelem (Tp->time,prev_time);
  
  //  marking of elements switched on
  Gtt->update_elem (Tp->time);
  
  //  searching for new DOFs
  //  only new DOFs are switched on
  //  former DOFs are temporarily switched off
  tncd = Gtt->search_newdofs (Tp->time,prev_time);
  
  //  marking of DOFs switched on
  Gtt->update_dofs (Tp->time);
  
  //  generation of new code numbers
  Ndoft = Gtt->codenum_generation (Outt);  
  
  if (tnce!=0 || tncd!=0){
    //  assembling of displacement vector for new code numbers
    Tt->lhs_restore (lhs,Lsrst->lhsi,tdlhs);
    
    //  cleaning of the conductivity matrix
    if (Kmat != NULL){
      delete Kmat;
      Kmat=NULL;
    }
    //  cleaning of the capacity matrix
    if (Cmat != NULL){
      delete Cmat;
      Cmat=NULL;
    }
  }
}




/**
   function solves linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   JK, 20.12.2002,
   revised 24.8.2006
*/
void linear_nonstat_solv_vform ()
{
  long i,j,n;
  double dt,end_time,alpha,*d,*p,*lhs,*tdlhs,*rhs,*f;
  
  //  number of unknowns
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  //  step number
  i=0;
  
  //  initiation of material models
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation();  
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();

  ///
  if ((Adat != NULL) && (Adat->give_step()>0)) {
    Adat->statedata_restore();
    // timectrl ??
    i = Tp->istep;
    print_initt(-1, "at");  // ??
  }
  else {
    //  inital printing of output and graphical informations
    if (Tp->hdbcont.restore_stat()){
      solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
      print_initt(-1, "at");
      //print_stept(0,i,Tp->time,NULL);
    }
    else{
      compute_req_valt (0);
      print_initt(-1, "wt");
      print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
    }
  }
  
  /// termitovo kontrola transferu state variables  DOCASNE
  if (Adat != NULL){
    if (Adat->give_step()) {
      Adat->run (2, true);
      Adat->answer = 0;
    }
  }
  
  bool breakloop = false;
  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime ();
    //  computation of backward time step
    dt = Tp->timecont.actualbacktimeincr ();
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    Tp->istep = i;
    
    //  capacity matrix    
    capacity_matrix (0);
    

   
    //  conductivity matrix
    conductivity_matrix (0);

    
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
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,f,Outt);
    
    //  nodal values computed from nodal derivatives
    //  d_{n+1} = dd + alpha dt v_{n+1}
    for (j=0;j<n;j++){

      //debug??!!
      if(lhs[j]>=0.0)
	lhs[j]=-1.0;

      lhs[j] = d[j] + alpha*dt*tdlhs[j];

      //debug??!!
      //if(lhs[j]>=0.0)
      //lhs[j]=-1.0;
    }
    
    solution_correction ();
    //  nodes - integration points interpolation
    approximation ();
    compute_req_valt (0);
    
    // pridal Madera 18.1.2008
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
    {
      if (Mesprt==1)
	fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    }
    
    /// adaptivity
    if (Tp->adaptivityflag  &&  Tp->time < end_time)
      if (i%2 == 0)
	breakloop = Adat->run (2, true);
    
  }while(Tp->time < end_time  &&  breakloop == false);
  
  delete [] p;
  delete [] d;
  delete [] f;
  
  print_closet();
}



/**
   function solves linear nonstationary transport problem
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 23.12.2002,
   revised 24.8.2006
*/
void linear_nonstat_solv_dform ()
{
  long i,j,n;
  double dt,end_time,alpha;
  double *d,*f,*p,*lhs,*tdlhs,*rhs,*fi;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //  vector of internal fluxes
  fi = new double [n];
  nullv (fi,n);


  //  nodes - integration points interpolation
  approximation ();
  actual_previous_change ();


  // vlozil JM 15.1.2008
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  //  step number
  i=0;

  //  initiation of material models
  Tm->initmaterialmodels();
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    print_initt(-1, "at");
    //print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  }

  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime ();
    //  computation of backward time step
    dt = Tp->timecont.actualbacktimeincr ();
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    
    //  capacity matrix
    capacity_matrix (0);
    
    /*
    FILE *aux;
    aux = fopen ("matice.txt","a");
    fprintf (aux,"\n\n Matice C ulohy, time %le\n",Tp->time);
    Cmat->printmat (aux);
    */
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    /*
    fprintf (aux,"\n\n Matice K ulohy, time %le\n",Tp->time);
    Kmat->printmat (aux);
    */
    
    //  predictor
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    //  C dd
    Cmat->gmxv (d,p);
    
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed nodal fluxes f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    //mechanical influence of transport problem
    trfel_right_hand_side2 (0,rhs,n);
  
    
    //  computation of the right hand side vector
    //  alpha dt f_n+1 + C dd
    for (j=0;j<n;j++){
      f[j] = rhs[j]*alpha*dt + p[j];
    }
    
    /*
    fprintf (aux,"\n\n Matice C + alpha dt K ulohy, time %le\n",Tp->time);
    Kmat->printmat (aux);
    
    fprintf (aux,"\n\n prava strana ulohy, time %le\n",Tp->time);
    for (j=0;j<n;j++){
      fprintf (aux,"% 15.12le\n",f[j]);
    }
    */
    
    //  solution of the system of algebraic equations
    //  (C + alpha dt K) d_{n+1} = alpha dt f_n+1 + C dd
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,f,Outt);
    
    //  first time derivatives of nodal values
    // v_{n+1} = 1/alpha/dt (d_{n+1}-dd)
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    
    
    /*
    fprintf (aux,"\n\n vektor reseni, time %le\n",Tp->time);
    for (j=0;j<n;j++){
      fprintf (aux,"% 15.12le\n",lhs[j]);
    }
    fclose(aux);
    */

    
    //  computation of values at integration points from nodal values
    solution_correction ();
    approximation ();
    
    // pridal Madera 18.1.2008
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    //  printing of output and graphical informations
    if (Tp->fluxcomp >= 1)
      internal_fluxes (fi,n);

    /*
    nullv (fi,n);
    nodal_energy (fi,n,dt);
    actual_previous_change ();
    double jkjk=0.0;
    for (j=0;j<n;j++){
      jkjk += fi[j]*fi[j];
    }
    fprintf (stdout,"\n jkjk   %15.10le",jkjk);      
    fprintf (Outt,"\n jkjk   %15.10le",jkjk);      
    */

    print_stept(0,i,Tp->time,rhs);// printing of output and graphical informations
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    if ((Tp->timecont.isitimptime ()==1) && (Tp->hdbcont.hdbtype))
    {
      if (Mesprt==1)
	fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    }
    
  }while(Tp->time<end_time);
  
  delete [] fi;
  delete [] f;
  delete [] p;
  delete [] d;
  
  print_closet();
}

/**
   function solves linear nonstationary transport problem
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 23.12.2002,
   revised 24.8.2006
*/
void linear_nonstat_radiation_solv_dform ()
{
  long i,j,n;
  double dt,end_time,alpha;
  double *d,*f,*p,*lhs,*tdlhs,*rhs,*fi;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //  vector of internal fluxes
  fi = new double [n];
  nullv (fi,n);


  //  nodes - integration points interpolation
  approximation ();

  // vlozil JM 15.1.2008
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  //  step number
  i=0;

  //  initiation of material models
  Tm->initmaterialmodels();
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    print_initt(-1, "at");
    //print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  }
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime ();
    //  computation of backward time step
    dt = Tp->timecont.actualbacktimeincr ();
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
        
    //  capacity matrix
    capacity_matrix (0);
    

    //  conductivity matrix
    conductivity_matrix (0);

    //fprintf (Outt,"\n\n matice vodivosti \n");
    //Kmat->printmat (Outt);
    //fprintf (Outt,"\n\n konec matice vodivosti \n\n\n");
    //fprintf (Outt,"\n\n matice kapacity \n");
    //Cmat->printmat (Outt);
    //fprintf (Outt,"\n\n konec matice kapacity \n\n\n");
    
    //  predictor
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    //  C dd
    Cmat->gmxv (d,p);
    
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed nodal fluxes f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    //mechanical influence of transport problem
    trfel_right_hand_side2 (0,rhs,n);//old

    
    //  computation of the right hand side vector
    //  alpha dt f_n+1 + C dd
    for (j=0;j<n;j++){
      f[j] = rhs[j]*alpha*dt + p[j];
    }
    
    Tt->edge_temperature ();
    Tt->heat_fluxes (f,Outt);

    //fprintf (Outt,"\n\n matice ulohy \n");
    //Kmat->printmat (Outt);
    //fprintf (Outt,"\n\n konec matice ulohy \n\n\n");
    //fprintf (Outt,"\n\n vektor prave strany \n");
    //for (j=0;j<n;j++){
    //fprintf (Outt,"\n f %le",f[j]);
    //}
    //fprintf (Outt,"\n\n konec vektoru prave strany \n\n\n");

    //  solution of the system of algebraic equations
    //  (C + alpha dt K) d_{n+1} = alpha dt f_n+1 + C dd
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,f,Outt);
    
    //  first time derivatives of nodal values
    // v_{n+1} = 1/alpha/dt (d_{n+1}-dd)
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    //  computation of values at integration points from nodal values
    solution_correction ();
    approximation ();
    
    // pridal Madera 18.1.2008
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    //  printing of output and graphical informations
    if (Tp->fluxcomp >= 1)
      internal_fluxes (fi,n);

    print_stept(0,i,Tp->time,rhs);// printing of output and graphical informations
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    if ((Tp->timecont.isitimptime ()==1) && (Tp->hdbcont.hdbtype))
    {
      if (Mesprt==1)
        fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    }
  }while(Tp->time<end_time);
  
  delete [] fi;
  delete [] f;
  delete [] p;
  delete [] d;
  
  print_closet();
}



/**
   function solves linear nonstationary transport problem
   time integration is based on subcycling
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 27. 5. 2014
*/
void linear_nonstat_solv_dform_subcycl ()
{
  long i,j,n;
  double dt,end_time,alpha;
  double *d,*f,*p,*lhs,*tdlhs,*rhs,*fi;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //  vector of internal fluxes
  fi = new double [n];
  nullv (fi,n);


  //  nodes - integration points interpolation
  approximation ();

  // vlozil JM 15.1.2008
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  //  step number
  i=0;

  //  initiation of material models
  Tm->initmaterialmodels();
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    print_initt(-1, "at");
    //print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  }
  
  //  list of elements belonging to particular subdomains is assembled
  Gtt->stop->elem_lists ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  
  //  cela uloha bude resena implicitne
  //  male podoblasti budou reseny explicitne
  //  pro male oblasti by se mela vyuzit gmatrix se systemem ukladani elem matrix
  //  ulozi se jen matice zapnutych prvku a zaroven se ukladaji jejich kodova cisla
  //  nemelo by byt nutne generovat nova kodova cisla
  do{
    //  determination of time steps on subdomains
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime ();
    //  computation of backward time step
    dt = Tp->timecont.actualbacktimeincr ();
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    
    
    // ******************************************
    //  integration on selected subdomains only
    // ******************************************
    //  explicit algorithm is used in selected subdomains
    //  matrices of the whole problem are not assembled
    //  conductivity matrices are assembled only on elements switched on
    //  they are stored in the format element_matrices=40
    //  capacity matrices are assembled only on elements switched on
    //  only diagonal entries are stored in a vector
    //  capacity matrices have to be diagonalized
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    // ***********************************************
    //  end of the integration on selected subdomain
    // ***********************************************

    
    //  capacity matrix
    capacity_matrix (0);
    
    /*
    FILE *aux;
    aux = fopen ("matice.txt","a");
    fprintf (aux,"\n\n Matice C ulohy, time %le\n",Tp->time);
    Cmat->printmat (aux);
    */
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    /*
    fprintf (aux,"\n\n Matice K ulohy, time %le\n",Tp->time);
    Kmat->printmat (aux);
    */
    
    //  predictor
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    //  C dd
    Cmat->gmxv (d,p);
    
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed nodal fluxes f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    //mechanical influence of transport problem
    trfel_right_hand_side2 (0,rhs,n);

    
    //  computation of the right hand side vector
    //  alpha dt f_n+1 + C dd
    for (j=0;j<n;j++){
      f[j] = rhs[j]*alpha*dt + p[j];
    }
    
    /*
    fprintf (aux,"\n\n Matice C + alpha dt K ulohy, time %le\n",Tp->time);
    Kmat->printmat (aux);
    
    fprintf (aux,"\n\n prava strana ulohy, time %le\n",Tp->time);
    for (j=0;j<n;j++){
      fprintf (aux,"% 15.12le\n",f[j]);
    }
    */
    
    //  solution of the system of algebraic equations
    //  (C + alpha dt K) d_{n+1} = alpha dt f_n+1 + C dd
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,f,Outt);
    
    //  first time derivatives of nodal values
    // v_{n+1} = 1/alpha/dt (d_{n+1}-dd)
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    /*
    fprintf (aux,"\n\n vektor reseni, time %le\n",Tp->time);
    for (j=0;j<n;j++){
      fprintf (aux,"% 15.12le\n",lhs[j]);
    }
    fclose(aux);
    */

    
    //  computation of values at integration points from nodal values
    solution_correction ();
    approximation ();
    
    // pridal Madera 18.1.2008
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    //  printing of output and graphical informations
    if (Tp->fluxcomp >= 1)
      internal_fluxes (fi,n);

    print_stept(0,i,Tp->time,rhs);// printing of output and graphical informations
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    if ((Tp->timecont.isitimptime ()==1) && (Tp->hdbcont.hdbtype))
    {
      if (Mesprt==1)
        fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    }
  }while(Tp->time<end_time);
  
  delete [] fi;
  delete [] f;
  delete [] p;
  delete [] d;
  
  print_closet();
}














/**
   function solves linear nonstationary transport problem
   the explicit finite difference method is used
   the forward version is used
   
   @param lcid - load case id
   
   JK, 21. 10. 2023
*/
void explicit_difference_method (long lcid)
{
  long n,i,j;
  double dt,end_time;
  double *lhs,*tdlhs,*rhs,*p,*z,*w;
  
  //  the number of degrees of freedom in the problem
  n = Ndoft;
  
  //  nodal values - displacements at next step
  //  (assigment of initial values)
  lhs = Lsrst->give_lhs (lcid);
  nullv (lhs,n);

  //  time derivatives of nodal values
  //  (assigment of initial values)
  tdlhs = Lsrst->give_tdlhs (lcid);
  nullv (tdlhs,n);
  
  //  vector of the right hand side
  rhs = Lsrst->give_rhs (2*lcid);
  nullv (rhs,n);

  //  auxiliary vector
  p = new double [n];
  w = new double [n];
  z = new double [n];


  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  //  step number
  i=0;

  //  initial nodal values
  //copyv (Lsrst->lhsi,lhs,n);
  
  //  nodes - integration points interpolation
  approximation ();
  actual_previous_change ();

  //  initiation of material models
  Tm->initmaterialmodels();
  
  // ********************************************************
  //  computation of the time derivatives of nodal values at the initial time
  // ********************************************************
  
  //  capacity matrix    
  capacity_matrix (lcid);
  
  //  conductivity matrix
  conductivity_matrix (lcid);
  
  //  right hand side
  trfel_right_hand_side (lcid,rhs,n);
  
  //  K u_0
  Kmat->gmxv (lhs,p);
  
  //  f_0 - K u_0
  cmulv (-1.0,p,n);
  addv (rhs,p,n);
  
  //  C v_0 = f_0 - K u_0
  Tp->ssle->solve_system (Gtt,Cmat,tdlhs,rhs,Outt);
  
  for (j=0;j<n;j++){
    //p[j] = lhs[j] - dt*tdlhs[j];
    p[j] = lhs[j];
  }
  
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    //solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    print_initt(-1, "at");
    //print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  }

  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime ();
    //  computation of backward time step
    dt = Tp->timecont.actualbacktimeincr ();
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    
    //  right hand side
    trfel_right_hand_side (lcid,rhs,n);
    
    //  K u_j
    Kmat->gmxv (lhs,z);
    //  C u_{j-1}
    //Cmat->gmxv (p,w);
    Cmat->gmxv (lhs,w);
    
    for (j=0;j<n;j++){
      //rhs[j]=2.0*dt*rhs[j] - 2.0*dt*z[j] + w[j];
      rhs[j]=dt*rhs[j] - dt*z[j] + w[j];
      //p[j] = lhs[j];
    }
    
    Tp->ssle->solve_system (Gtt,Cmat,lhs,rhs,Outt);    
    
    print_stept(0,i,Tp->time,rhs);// printing of output and graphical informations
    //print_stept(0,i,Tp->time,NULL);
    print_flusht();

    
  }while(Tp->time < end_time);
  
  delete [] p;
  delete [] w;
  delete [] z;

  print_closet();
}




/**
   function solves linear nonstationary transport problem
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   this is a special version for extremely large problems
   the capacity and conductivity matrices are assumed constant
   the capacity matrix is diagonalized
   the capacity matrix should be stored as a diagonal matrix
   the conductivity matrix is stored in the symmetric compressed row storage scheme
   systems of linear equations should be solved by the conjugate gradient method
   the last solution should be used as an initial approximation in the next step
   
   the coefficient alpha must not be zero (because it is the d-version)
   
   JK, 22. 10. 2023,
*/
void optim_linear_nonstat_solv_dform ()
{
  long i,j,n;
  double dt,pdt,end_time,alpha,zero;
  double *d,*f,*p,*lhs,*tdlhs,*rhs,*fi;
  
  zero=1.0e-6;
  
  //  number of unknowns / degrees of freedom
  n=Ndoft;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //  vector of internal fluxes
  fi = new double [n];
  nullv (fi,n);


  //  nodes - integration points interpolation
  approximation ();
  actual_previous_change ();
  
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  if (alpha<zero){
    print_err("alpha coefficient in the optimized trapezoidal method is zero", __FILE__, __LINE__, __func__);
    abort ();
  }
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  pdt = dt;

  //  step number
  i=0;
  
  //  initiation of material models
  Tm->initmaterialmodels();
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    print_initt(-1, "at");
    //print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  }

  //  capacity matrix
  capacity_matrix (0);
  //  conductivity matrix
  conductivity_matrix (0);
  
  //  matrix of the system of equations
  //  C + alpha.dt.K
  Kmat->scalgm (dt*alpha);
  Kmat->addgm (1.0,*Cmat);
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime ();
    //  computation of backward time step
    dt = Tp->timecont.actualbacktimeincr ();
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    
    if (fabs(pdt-dt)>zero){
      //  the time step was changed
      //  therefore the matrices have to be recalculated
      
      pdt = dt;
      
      //  capacity matrix
      capacity_matrix (0);
      
      /*
	FILE *aux;
	aux = fopen ("matice.txt","a");
	fprintf (aux,"\n\n Matice C ulohy, time %le\n",Tp->time);
	Cmat->printmat (aux);
      */
      
      //  conductivity matrix
      conductivity_matrix (0);
      
      /*
	fprintf (aux,"\n\n Matice K ulohy, time %le\n",Tp->time);
	Kmat->printmat (aux);
      */
      
      //  matrix of the system of equations
      //  C + alpha.dt.K
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
    }
    
    //  predictor
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    //  C dd
    Cmat->gmxv (d,p);
    
    
    //  prescribed nodal fluxes f_{n+1}
    trfel_right_hand_side (0,rhs,n);
    //mechanical influence of transport problem
    trfel_right_hand_side2 (0,rhs,n);
  
    
    //  computation of the right hand side vector
    //  alpha dt f_n+1 + C dd
    for (j=0;j<n;j++){
      f[j] = rhs[j]*alpha*dt + p[j];
    }
    
    /*
    fprintf (aux,"\n\n Matice C + alpha dt K ulohy, time %le\n",Tp->time);
    Kmat->printmat (aux);
    
    fprintf (aux,"\n\n prava strana ulohy, time %le\n",Tp->time);
    for (j=0;j<n;j++){
      fprintf (aux,"% 15.12le\n",f[j]);
    }
    */
    
    //  solution of the system of algebraic equations
    //  (C + alpha dt K) d_{n+1} = alpha dt f_n+1 + C dd
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,f,Outt);
    
    //  first time derivatives of nodal values
    // v_{n+1} = 1/alpha/dt (d_{n+1}-dd)
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    
    
    /*
    fprintf (aux,"\n\n vektor reseni, time %le\n",Tp->time);
    for (j=0;j<n;j++){
      fprintf (aux,"% 15.12le\n",lhs[j]);
    }
    fclose(aux);
    */

    
    //  computation of values at integration points from nodal values
    //solution_correction ();
    //approximation ();
    
    //  printing of output and graphical informations
    if (Tp->fluxcomp >= 1)
      internal_fluxes (fi,n);

    /*
    nullv (fi,n);
    nodal_energy (fi,n,dt);
    actual_previous_change ();
    double jkjk=0.0;
    for (j=0;j<n;j++){
      jkjk += fi[j]*fi[j];
    }
    fprintf (stdout,"\n jkjk   %15.10le",jkjk);      
    fprintf (Outt,"\n jkjk   %15.10le",jkjk);      
    */

    print_stept(0,i,Tp->time,rhs);// printing of output and graphical informations
    //print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    if ((Tp->timecont.isitimptime ()==1) && (Tp->hdbcont.hdbtype))
    {
      if (Mesprt==1)
	fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    }
    
  }while(Tp->time<end_time);
  
  delete [] fi;
  delete [] f;
  delete [] p;
  delete [] d;
  
  print_closet();
}
