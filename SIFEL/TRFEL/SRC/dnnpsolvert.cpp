#include "npsolvert.h"
#include "dnnpsolvert.h"
#include "nnpsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "transprint.h"
#include "saddpoint.h"
#include "backupsolt.h"
#include "elemswitcht.h"
#include <string.h>

void solve_discont_nonlin_nonstationary_problem ()
{
  nonlin_nonstat_dform ();

  //nonstat_solv_dform_comp ();
}

/**
  Function allocates and initializes global vectors used in npsolver
  It is used in one_step concept solver.

  Created by Tomas Koudelka 6.2014 according to JK
*/

void nonstat_solver_dform_init (long lcid, np_glob_vec &np_gv)
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
  np_gv.alloc_daux(n);
  

  if(Tp->tprob == nonlinear_nonstationary_problem)
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
      //  initiation of transport material models
      Tm->initmaterialmodels();
      //  initiation of transport material models at aux. int. points
      Tm->aip_initmaterialmodels();
      compute_req_valt (lcid);
      print_initt(-1, "at");
    }
    else{
      //  initiation of transport material models
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
   The function performs one step of solution of  non-linear nonstationary transport 
   problem solved by Newton-Raphson method.
   The d-form version of the generalized trapezoidal method is used -
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   Created by Tomas Koudelka 6.2014 according to JK
*/
long one_step_nonlinear_dform (long lcid,double time, double dt, long istep, long /*li*/, np_glob_vec &np_gv)
{
  long j,n,nbdof,nsad,ini,stop;
  double norfb;
  double *err,alpha,zero,*thresh;
  double *f,*d,*p,*v,*z,*lhs,*tdlhs,*rhs,*lhsb,*tdlhsb;

  //new time increment
  Tp->time = time;
  //  new step number
  Tp->istep = istep;
  Tp->jstep = -1;

  //  the number of primal unknowns / primal degrees of freedom
  //  it does not contain the number of multipliers
  n=Ndoft;
  //  the number of primal boundary/interface unknowns
  nbdof=Gtt->nbdof;
  //  the number of unknowns in reduced problem
  //  the sum of the number of primal boundary/interface unknowns and the number of Lagrange multipliers
  nsad=Gtt->nsad;

  //  computer zero
  zero=Tp->zero;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;

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

  // auxiliary vectors:
  v = np_gv.v;
  z = np_gv.z;
  //  backup of nodal values
  lhsb = np_gv.lhsb;
  //  backup of time derivatives of nodal values
  tdlhsb = np_gv.tdlhsb;

  //  parameter of the trapezoidal method
  alpha=Tp->alpha;
  //  prescribed norm of residual
  err=Tp->errarr;
  //  maximum number of iterations in one time step
  ini=Tp->nii;

  if (Mesprt==1)  fprintf (stdout,"\n\n --------------------------------------------------------------");
  if (Mesprt==1)  fprintf (stdout,"\n TRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
  if (Mesprt==1)  fprintf (stdout,"\n --------------------------------------------------------------\n");
  
  // update material properties and auxiliary values
  Tm->updateipval ();
    
  //  backup of attained nodal values and their time derivatives
  // lhsb[i] = lhs[i]
  copyv(lhs, lhsb, n);
  // tdlhsb[i] = tdlhs[i]
  copyv(tdlhs, tdlhsb, n);

  //  capacity matrix
  capacity_matrix (0);
    
  //  conductivity matrix
  conductivity_matrix (0);
    
  //  predictor d
  //  d[j] = lhs[j]+dt*(1.0-alpha)*tdlhs[j];
  addmultv(lhs, tdlhs, dt*(1.0-alpha), d, n);

  //  C.d
  Cmat->gmxv (d,p);
    
  //  matrix of the system of equations
  //  C + alpha.dt.K
  Kmat->scalgm (dt*alpha);
  Kmat->addgm (1.0,*Cmat);
    
  //  prescribed fluxes f_{n+1}
  trfel_right_hand_side (lcid,f,n);

  //  C.d + alpha dt f_{n+1}
  //  rhs[j] = f[j]*alpha*dt + p[j];
  addmultv(p, f, alpha*dt, rhs, n);

  Tt->compute_jumps (rhs);
  nullv (z,n);
  Tt->compute_jumps (z);

  Gtt->mult_localization (Kmat);
    

  Kmat->diag_check (zero,rhs);
    
  Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    
  //  computation of time derivatives of the nodal values
  //  tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
  subv(lhs, d, tdlhs, n);
  cmulv(1.0/dt/alpha, tdlhs, n);
    
  //  physically corrected solution
  solution_correction ();

  //  approximation of nodal values into ontegration points
  approximation ();
    
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  //  iteration for equilibrium
  for (j=0;j<ini;j++)
  {
    //  capacity matrix
    capacity_matrix (lcid);
      
    //  conductivity matrix
    conductivity_matrix (lcid);
      
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
      
    //  (C + alpha.dt.K).d
    Gtt->mult_localization (Kmat);
    nullv (v,n);
    Kmat->gmxv (lhs,v);
      
    //  C.d
    nullv (p,n);
    Cmat->gmxv (d,p);
      
    // rhs[k]=f[k]*alpha*dt+p[k]-v[k];
    addmultv(p, f, alpha*dt, rhs, n);
    subv(rhs, v, rhs, n);
    // z[k]+=f[k]*alpha*dt+p[k];
    addmultv(z, f, alpha*dt, n);
    addv(z, p, n);

    Tt->compute_jumps (rhs);
    
    stop = norm_computation_vec (rhs,z,err,thresh,2,1,norfb);
      
    nullv (z,n);
    Tt->compute_jumps (z);
      
    if (Mesprt==1)  fprintf (stdout,"\n inner iteration %ld ______",j);

    if (stop==1){
      break;
    }
    
    Kmat->diag_check (zero,rhs);


    Tp->ssle->solve_system (Gtt,Kmat,z,rhs,Outt);

    //	lhs[k]+=z[k];
    addv(lhs,z,n);
      
    nullv (z,n);
    Tt->compute_jumps (z);
      
    //  computation of time derivatives of the nodal values
    //  tdlhs[k]=(lhs[k]-d[k])/dt/alpha;
    subv(lhs, d, tdlhs, n);
    cmulv(1.0/dt/alpha, tdlhs, n);
      
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();
      
      
    // nulleqother ();
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
      // ***********************************
      
  }
  
  if (stop == 0)
  {
    return -1;  
  }

  return j;
}



/**
   The function solves non-linear nonstationary transport problem by Newton-Raphson method.
   The d-form version of the generalized trapezoidal method is used -
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   Created by Tomas Koudelka 06.2014 
   according to JK, 20.12.2002 and TKo mtsolver.cpp
*/
void nonstat_solv_dform_comp ()
{
  long i,li,nsts,ret,n;
  double newtime,dt,dtmin,dtmax,dtdef,end_time;
  np_glob_vec np_gv;
  long lcid = 0;

  //
  //  initialization phase
  //
  //  nodes - integration points interpolation
  approximation();  
  nonstat_solver_dform_init(lcid, np_gv);

  //  number of transport degrees of freedom
  n=Ndoft;
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

  bool breakloop = false;
  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    i++;
  
    newtime = Tp->time = Tp->timecont.newtime(dt);
    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Tp->timecont.timefun.tfunc != constant) && Tp->timecont.tct == 0)
      dt = Tp->timecont.actualbacktimeincr ();

    if(Tp->tprob == discont_nonlin_nonstat_problem){
      // non-linear algorithm
      // perform one time step with non-linear solver
      ret = one_step_nonlinear_dform(lcid, newtime, dt, i, li, np_gv);
    }
    
    // handling with time controler for nonlinear solver, only for adaptive time increment
    if(Tp->tprob == discont_nonlin_nonstat_problem && Tp->timecont.tct > 0)
    { //only for adaptive time controler
      if (ret >= 0) // equilibrium was attained
      {

	//Tm->freezing_thawing ();
	fprintf (stdout,"\n\n TADY NEMAS CO DELAT (dnnpsolvert.cpp, line 380)\n");
	abort();

        //////////////////////////////////////////////////
        //  printing of output and graphical informations
        compute_req_valt (lcid);
        print_stept(lcid,Tp->istep,Tp->time,np_gv.rhs);//print_stept(0,i,Tp->time,NULL);
        print_flusht();
  
        if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
        {
          if (Mesprt==1)
            fprintf (stdout,"\n Creating TRFEL backup file\n");
          solvert_save (np_gv.lhs,np_gv.tdlhs,np_gv.f,Tp->istep,Tp->time,dt,Tp->timecont,n);
        }

        //  backup of attained nodal values and their time derivatives
        //  lhsb[k] = np_gv.lhs[k]; //storing results from previous inner step
        copyv(np_gv.lhs, np_gv.lhsb, n);
        //  tdlhsb[k] = np_gv.tdlhs[k]; //storing results from previous inner step
        copyv(np_gv.tdlhs, np_gv.tdlhsb, n);
        //  approximation of nodal values into integration points
        approximation ();
	  
        if (ret == 0)
          nsts++; 
        else
          nsts=0;     

        dtdef = Tp->timecont.actualforwtimeincr();
        if (nsts==2)
        {
          dt*=2.0;
          nsts=0;
	      
          if (Mesprt==1)  
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
        }
        if (dt<=dtdef)//increase of time increment according to prescribed one
          dt = dtdef;
	  
        if (dt>dtmax)//maximum time increment
          dt = dtmax;
      }
      else
      {
        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
        nsts = 0;
        dt/=2.0;
        Tp->timecont.oldtime ();
        i--;
	  
        // np_gv.lhs[k] = lhsb[k]; //restoring results from previous time step
        copyv(np_gv.lhsb, np_gv.lhs, n);
        // np_gv.tdlhs[k] = tdlhsb[k]; //restoring results from previous time step
        copyv(np_gv.tdlhsb, np_gv.tdlhs, n);
        //  approximation of nodal values into ontegration points
        approximation ();
	  
        if (Mesprt==1)  
          fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
	  
        if (dt<dtmin)
        {
          if (Mesprt==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
          if (Mesprt==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
          if (Mesprt==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
          compute_req_valt (lcid);
          print_stept_forced(lcid, i, Tp->time, np_gv.rhs);
          print_flusht();
          break;
        }
      }
    }
    double ti;
    ti = total_integral (0);
    fprintf (stdout,"\n total integral %le",ti);
    
    // adaptivity
    if (Tp->adaptivityflag  &&  Tp->time < end_time)
      if (i%2 == 0)
	breakloop = Adat->run (2, true);
    
  }while(Tp->time < end_time  &&  breakloop == false);
  
  print_closet();
}



/**
   function solves nonstationary transport problem with discontinuity
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis

   JK, 19.8.2008
*/
void nonlin_nonstat_dform ()
{
  long i,j,k,n,nbdof,nsad,lcid,ani,ini,stop,nsts;
  //long *dmcn;
  double *err,dt,dtmin,dtmax,end_time,alpha,zero,norfb,*thresh;
  double *f,*d,*p,*v,*z,*lhs,*lhsi,*tdlhs,*rhs,*lhsb,*tdlhsb;
  //double *condmat,*condvect,*splhs,*sprhs;
  densemat dm;
  
  //  load case id must be equal to zero in this type of problem
  lcid=0;
  
  
  //  the number of primal unknowns / primal degrees of freedom
  //  it does not contain the number of multipliers
  n=Ndoft;
  //  the number of primal boundary/interface unknowns
  nbdof=Gtt->nbdof;
  //  the number of unknowns in reduced problem
  //  the sum of the number of primal boundary/interface unknowns and the number of Lagrange multipliers
  nsad=Gtt->nsad;
  
  //  computer zero
  zero=Tp->zero;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->give_lhsi(0);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //  auxiliary vector
  v = new double [n];
  nullv (v,n);
  //  auxiliary vector
  z = new double [n];
  nullv (z,n);
  //  backup of nodal values
  lhsb = new double [n];
  nullv (lhsb,n);
  //  backup of time derivatives of nodal values
  tdlhsb = new double [n];
  nullv (tdlhs,n);


  
  //  array containing condensed matrix
  //condmat = new double [nbdof*nbdof];
  //nullv (condmat,nbdof*nbdof);
  //  array containing condensed vector
  //condvect = new double [nbdof];
  //nullv (condvect,nbdof);
  
  //  dense matrix for resulting matrix of the system
  //dm.alloc (nsad);
  
  //  array containing solution of the saddle point problem
  //splhs = new double [nsad];
  //nullv (splhs,nsad);
  //  array containing the right hand side of the saddle point problem
  //sprhs = new double [nsad];
  //nullv (sprhs,nsad);

  //  code numbers of primal variables in the resulting system of equations
  //dmcn = new long [nbdof];
  //for (i=0;i<nbdof;i++){
  //dmcn[i]=i+1;
  //}

  
  
  //  nodes - integration points interpolation
  approximation ();

  // vlozil JM 15.1.2008
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  
  //  parameter of the trapezoidal method
  alpha=Tp->alpha;
  //  prescribed norm of residual
  err=Tp->errarr;
  //  maximum number of iterations in one time step
  ini=Tp->nii;
  //  minimum time increment
  dtmin=Tp->timecont.dtmin;
  
  //  maximum time increment
  dtmax=Tp->timecont.dtmax;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // vlozil JM 25.4.2008----*****
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  //  toto vyhodit
  //approximation (); 
  
  //  toto vyhodit
  // nulleqother ();
  //if (Tp->nvs==1 && Tp->pnvs==1)
  //actual_previous_nodval ();
  // -----*******
  //  az sem
  

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
  //  number of time step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  stop=0;

  //  toto vyhodit
  //Tm->initmaterialmodels();

  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  do{
    
    //  time update
    Tp->time=Tp->timecont.newtime (dt);
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;
    
    //  backup of attained nodal values and their time derivatives
    for (j=0;j<n;j++){
      lhsb[j]=lhs[j];
      tdlhsb[j]=tdlhs[j];
    }

    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
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

/*    
    fprintf (Outt,"\n\n kontrola prave strany");
    fprintf (Outt,"\n\n\n\n\n kontrola v case  %le",Tp->time);
    //fprintf (stdout,"\n\n kontrola prave strany");
    for (j=0;j<n;j++){
      //fprintf (stdout,"\n rhs %6ld   %le",j,rhs[j]);
      fprintf (Outt,"\n rhs-f %6ld   %20.15lf",j,f[j]);
    }
    //fprintf (stdout,"\n\n");
    fprintf (Outt,"\n\n");
*/

    //  C.d + alpha dt f_{n+1}
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
    }

 

    //Tt->compute_jumps (rhs);
    //nullv (z,n);
    //Tt->compute_jumps (z);

    //Gtt->mult_localization (Kmat);
    

    //Kmat->diag_check (zero,rhs);
    
    //if (i==1)
    //Kmat->printmat (Outt);
    
    //fprintf (Outt,"\n\n kontrola prave strany");
    //fprintf (Outt,"\n\n\n\n\n kontrola v case  %le",Tp->time);
    //fprintf (stdout,"\n\n kontrola prave strany");
    //for (j=0;j<n;j++){
      //fprintf (stdout,"\n rhs %6ld   %le",j,rhs[j]);
      //fprintf (Outt,"\n rhs %6ld   %20.15lf",j,rhs[j]);
    //}
    //fprintf (stdout,"\n\n");
    //fprintf (Outt,"\n\n");
    


      
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);
    

/*
    fprintf (Outt,"\n\n kontrola reseni");
    //fprintf (stdout,"\n\n\n\n\n kontrola reseni time %le",Tp->time);
    //fprintf (stdout,"\n\n kontrola reseni");
    for (j=0;j<n;j++){
      //fprintf (stdout,"\n lhs %6ld  %le",j,lhs[j]);
      fprintf (Outt,"\n lhs %6ld  %20.15lf",j,lhs[j]);
    }
    //fprintf (stdout,"\n\n");
    fprintf (Outt,"\n\n");
*/
    
    //Kmat->printmat (Outt);
    
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
    
    //if (Tp->nvs==1 && Tp->pnvs==1)
      //actual_previous_nodval ();
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
      Gtt->mult_localization (Kmat);
      nullv (v,n);
      Kmat->gmxv (lhs,v);
      
      //  C.d
      nullv (p,n);
      Cmat->gmxv (d,p);
      
      
      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k]-v[k];
	//	fprintf (stdout,"\n rhs   %15.10le   k je %d",rhs[k],k);
	//norres+=rhs[k]*rhs[k];
	z[k]+=f[k]*alpha*dt+p[k];
	//norrhs+=z[k]*z[k];

	//fprintf (Outt,"\n %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",f[k],p[k],v[k],z[k],rhs[k]);

	//fprintf (stdout,"\n rhs %6ld   %le  f  %15.12le   v  %15.12le  p  %15.12le  z %15.12le",k,rhs[k],f[k],v[k],p[k],z[k]);
	//fprintf (stdout,"\n rhs %6ld   %le  f  %le   v  %le  p  %le  z %le",k,rhs[k],f[k],v[k],p[k],z[k]);
      }
      //norres=sqrt(norres);
      //norrhs=sqrt(norrhs);
      
      


    


      //Tt->compute_jumps (rhs);
      
      stop = norm_computation_vec (rhs,z,err,thresh,2,1,norfb);
      
      
      //nullv (z,n);
      //Tt->compute_jumps (z);
      //fprintf (stdout,"\n pauza \n");
      
      if (stop==1){
	break;
      }
      
      if (Mesprt != 0)  fprintf (stdout,"\n iteration number %ld",j);
      
      /*
      nullv (rhs,n);
      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k];
      }
      
      Tt->compute_jumps (rhs);
      */

      //Kmat->diag_check (zero,rhs);


      Tp->ssle->solve_system (Gtt,Kmat,z,rhs,Outt);

/*      
      fprintf (Outt,"\n\n casovy krok %ld,  iterace %ld\n",i,j);
      for (k=0;k<n;k++){
	fprintf (Outt,"\n z %6ld   % 15.12le",k,z[k]);
      }
      fprintf (Outt,"\n\n konec kontroly\n");
*/      
      
      for (k=0;k<n;k++){
	lhs[k]+=z[k];
      }
      
      //nullv (z,n);
      //Tt->compute_jumps (z);
      
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
      //if (Tp->nvs==1 && Tp->pnvs==1)
      //actual_previous_nodval ();
      // ***********************************
      
    }
    
    //  actual number of performed iterations
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
      //if (Tp->nvs==1 && Tp->pnvs==1)
      //actual_previous_nodval ();
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
	print_err("time increment is less than minimum time increment",__FILE__,__LINE__,__func__);
	print_err("computation fails",__FILE__,__LINE__,__func__);
	abort ();
      }
    }else{
      
      Tm->updateipval();
      
      // other computaion
      if (Tp->othercomp==1){
	if (Tp->otherpos==1)
	  {
	    //	compute_nodeotherst_comp ();
	    compute_nodeotherst ();
	    compute_ipotherst ();
	  }
	if (Tp->otherpos==2)
	  compute_nodeotherst ();
	if (Tp->otherpos==3)//computed directly at nodes
	  compute_nodeotherst_comp ();
      }
      
      Tm->freezing_thawing ();
      

      //  time increment
      //Tp->time = Tp->timecont.newtime ();
      //dtdef = Tp->timecont.actualforwtimeincr ();
      
      if (nsts==2){
	double times, acttimes;
	times =  Tp->timecont.starttime();
	acttimes = Tp->timecont.actualtime();
	
	if ((times+10)>acttimes) {
	  dt = dt;
	}
	else{
	  dt*=2;
	  if (dt>dtmax) {
	    dt = dtmax;
	  }
	}
	//dt*=1.2;
	nsts=0;
	//	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	//fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le, dtdef = %le",dt,dtdef);
       //	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le",dt);
	//fprintf (Outt2,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt= %lf",dt);
      }
      
      //if (dt>dtdef)
      //dt=dtdef;
      
      //Tp->time = Tp->timecont.newtime (dt);
      
      
      //time+=dt;
      //Tp->time=time;
      //Tp->timecont.time=time;
      //Tp->timecont.fdt=dt;
      
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      
    }
    
    /*
    double ti;
    ti = total_integral (0);
    fprintf (stdout,"\n total integral %le",ti);

    if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
    {
      if (Mesprt==1)
        fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    }
    */

  }while(Tp->time<end_time);
  
  delete [] v;
  delete [] p;
  delete [] d;
  delete [] f;
  delete [] lhsb;
  delete [] tdlhsb;
  delete [] z;
  print_closet();
}
