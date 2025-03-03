#include "mpi.h"
#include "ppcsolver.h"
#include "pglobal.h"
#include "pglobalt.h"
#include "pglobalc.h"
#include "globalc.h"
#include "global.h"
#include "globalt.h"
#include "globmat.h"
#include "globmatt.h"
#include "globmatc.h"
#include "elemswitcht.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "mechprint.h"
#include "transprint.h"
#include "backupsol.h"
#include "backupsolt.h"
#include "seqfilesm.h"
#include "seqfilest.h"
#include "genfile.h"
#include "pnpsolvert.h"
#include "pmtsolver.h"
#include <math.h>
#include <string.h>

/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem   
   transport problems are linear or nonlinear
   mechanical analysis is always nonlinear

   4.2009 TKo
*/
void par_solve_pcouplprob ()
{
  long lcid;
  
  //  load case id must be equal to zero
  //  the problems are nonlinear and superposition method cannot be used
  lcid=0;
  
  switch (Cp->tnlinsol){
  case newtonc:{
    //par_newton_raphson_parcoupl (lcid);
    //par_newton_raphson_parcoupl_comp (lcid);
    par_newton_raphson_parcoupl_common_dt (lcid);
    break;
  }
  default:{
    print_err("unknown solver of nonlinear equation system is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function solves parallel partially coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   TKr, 19/02/2013, according to JK and TKo
*/
void par_newton_raphson_parcoupl_comp (long lcid)
{
  int mefel_convergence,trfel_convergence;
  double dt,end_time;
  long i,k,nt,ti,tli,tnsts,tret;
  double tnewtime,tdt, tdtmin,tdtmax, tdtdef, tend_time;
  double *lhsb,*tdlhsb;
  np_glob_vec np_gv;
  long tlcid = 0;//zatim??!!
  long mi, mli, mnsts, mret;
  double mnewtime,mdt, mdtmin,mdtmax, mdtdef, mend_time;
  mt_glob_vec mt_gv;
  long mlcid = 0;//zatim??!!

  //
  //  PMETR initialization phase
  //
  //  starting time
  Cp->time = Cp->timecon.starttime ();
  //  time increment
  dt = Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();  

  //this part below is not neccessary because of independent time controlers of TRFEL and MEFEL
  //  mechanical time controller is rewritten by coupled time controller
  //Mp->time = Cp->time;
  //Mp->timecon.take_values (Cp->timecon);
  //  transport time controller is rewritten by coupled time controller
  //Tp->time = Cp->time;
  //Tp->timecont.take_values (Cp->timecon);
  // end this part


  //  approximation of quatities from TRFEL nodes to MEFEL integration points
  init_trfel_mefel ();

  //
  //  PTRFEL initialization phase
  //
  //par_nonstat_solver_init(tlcid, np_gv);//tady promyslet nejde slinkovat??!!
  par_nonstat_trfel_init(tlcid, np_gv);//vytvoreni nove funkce v PARMETRu
  //  initial time increment
  tdt = Tp->timecont.initialtimeincr ();
  //  minimum time increment
  tdtmin=Tp->timecont.dtmin;
  //  maximum time increment
  tdtmax=Tp->timecont.dtmax;
  //  end time
  tend_time = Tp->timecont.endtime ();
  //  number of step
  tli = ti = np_gv.istep;
  // return status from the one_step function
  tret = 0;
  
  //
  //  PMEFEL initialization phase
  //
  //par_visco_solver_init(mlcid, mt_gv);//tady promyslet nejde slinkovat??!!
  par_visco_mefel_init(mlcid, mt_gv);//vytvoreni nove funkce v PARMETRu
  //  initial time increment
  mdt = Mp->timecon.initialtimeincr ();
  //  minimum time increment used in MEFEL inner loop
  mdtmin=Mp->timecon.dtmin;
  //  maximum time increment
  mdtmax=Mp->timecon.dtmax;
  //  end time
  mend_time = Mp->timecon.endtime ();

  //  number of step
  mli = mi = mt_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  mnsts=0;
  // return status from the one_step function
  mret = 0;
  
  // only for adaptive time increment
  if(Tp->timecont.tct > 0){
    nt=Ndoft;
    lhsb   = new double [nt];
    tdlhsb   = new double [nt];
    nullv (lhsb,nt);
    nullv (tdlhsb,nt);
  }

  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;  // number of METR step
  mefel_convergence = 1;//setup of flag for successfully performed MEFEL step to yes
  trfel_convergence = 1;//setup of flag for successfully performed TRFEL step to yes

  do{
    i++;
    if(mefel_convergence == 0 || trfel_convergence == 0){
      Tp->timecont.oldtime ();

      if (Myrank==0 && Mesprc==1)  
	fprintf (stdout,"\n\n --------------------------------------------------------------------------------------------");
      if (Myrank==0 && Mesprc==1)  
	fprintf (stdout,"\n REDUCED TIME STEP, PMETR TIME STEP = %ld : PMETR TIME = %e, pmetr time increment = %e",i,Cp->time,dt);
      if (Myrank==0 && Mesprc==1)  
	fprintf (stdout,"\n --------------------------------------------------------------------------------------------");
      i--;
    }
    else{
            
      //new time and new time increment
      // driving time is set to TRFEL
      tnewtime = Tp->time = Tp->timecont.newtime (tdt);
      Cp->time = Tp->time;

      if (Myrank==0 && Mesprc==1)  
	fprintf (stdout,"\n\n -------------------------------------------------------------------------------------");
      if (Myrank==0 && Mesprc==1)
	fprintf (stdout,"\n NEW PMETR TIME STEP = %ld: PMETR TIME = %e, pmetr time increment = %e",i,Cp->time,dt);
      if (Myrank==0 && Mesprc==1)
	fprintf (stdout,"\n -------------------------------------------------------------------------------------");
    }
    
    //PMETR is copying data
    pass_coup_data(lcid);
    
    //PTRFEL part
    trfel_convergence = 1;
    ti = i;
    
    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Tp->timecont.timefun.tfunc != constant) && Tp->timecont.tct == 0)
      tdt = Tp->timecont.actualbacktimeincr ();
    
    if(Tp->tprob == nonstationary_problem)
      // linear algorithm
      // perform one time step with linear solver
      //tret = par_one_step_linear(tlcid, tnewtime, tdt, ti, tli, np_gv);//tady promyslet nejde slinkovat??!!
      tret = par_one_step_trfel_linear(tlcid, tnewtime, tdt, ti, tli, np_gv);//vytvoreni nove funkce v PARMETRu

    if(Tp->tprob == nonlinear_nonstationary_problem)
      // non-linear algorithm
      // perform one time step with non-linear solver
      //tret = par_one_step_nonlinear(tlcid, tnewtime, tdt, ti, tli, np_gv);//tady promyslet nejde slinkovat??!!
      tret = par_one_step_trfel_nonlinear(tlcid, tnewtime, tdt, ti, tli, np_gv);//vytvoreni nove funkce v PARMETRu
    
    // handling with time controler for nonlinear solver, only for adaptive time increment
    if(Tp->tprob == nonlinear_nonstationary_problem && Tp->timecont.tct > 0){ //only for adaptive time controler
      if (tret >= 0) // equilibrium was attained
	{
	  trfel_convergence = 1;
	  
	  for (k=0;k<nt;k++){
	    lhsb[k] = np_gv.lhs[k]; //storing results from previous inner step
	    tdlhsb[k] = np_gv.tdlhs[k]; //storing results from previous inner step
	  }
	  
	  if (tret == 0)
	    tnsts++;      
	  tdtdef = Tp->timecont.actualforwtimeincr();
	  if (tnsts==2)
	    {
	      tdt*=2.0;
	      tnsts=0;
	      
	      if (Mesprt==1)  
		fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	    }
	  if (tdt<=tdtdef)//increase of time increment according to prescribed one
	    tdt = tdtdef;
	  
	  if (tdt>tdtmax)//maximum time increment
	    tdt = tdtmax;
	}
      else
	{
	  //  reduction of the time increment because
	  //  inner loop was not able to enforce equilibrium
	  trfel_convergence = 0;
	  
	  tdt/=2.0;
	  
	  if (Mesprt==1)  
	    fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
	  
	  if (tdt<tdtmin)
	    {
	      if (Mesprt==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
	      if (Mesprt==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
	      if (Mesprt==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
	      compute_req_valt (lcid);
	      print_stept(lcid, ti, Tp->time, np_gv.rhs);
	      print_flusht();
	      break;
	    }
	}
    }
    
    //PMETR cleaning trfel matrices for more memory in mefel part
    if (Cp->cleanmatrix == clean_yes){
      //  cleaning matrices for more free memmory
      //  cleaning of the conductivity matrix
      if (Kmat != NULL){
	fprintf (stderr,"\n\n Cleaning of conductivity matrix \n\n ");
	delete Kmat;
	Kmat=NULL;
      }
      //  cleaning of the capacity matrix
      if (Cmat != NULL){
	fprintf (stderr,"\n\n Cleaning of capacity matrix \n\n ");
	delete Cmat;
	Cmat=NULL;
      }
    }


    //PMEFEL part
    mefel_convergence = 1;
    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Mp->timecon.timefun.tfunc != constant) && Mp->timecon.tct == 0)
      mdt = Mp->timecon.actualbacktimeincr ();
    
    // MEFEL part is set to follow TRFEL part
    //new time increment
    if (Mp->time < Tp->time)
      mnewtime = Mp->time = Mp->timecon.newtime (mdt);
    
    if (Mp->time <= Tp->time){
      //  new mefel step number
      mi++;
      
      //copying of data
      pass_coup_data(lcid);
      
      // perform one time step
      //mret = par_one_step(mlcid, mnewtime, mdt, mi, mli, mt_gv);//tady promyslet nejde slinkovat??!!
      mret = par_one_step_mefel(mlcid, mnewtime, mdt, mi, mli, mt_gv);//vytvoreni nove funkce v PARMETRu
      
      // handling with time controler for nonlinear solver, only for adaptive time increment
      if(Mp->timecon.tct > 0){//only for adaptive time controler
	if (mret >= 0) // equilibrium was attained
	  {
	    mefel_convergence = 1;
	    
	    if (mret == 0)
	      mnsts++;      
	    mdtdef = Mp->timecon.actualforwtimeincr();
	    if (mnsts==2)
	      {
		mdt*=2.0;
		mnsts=0;
		if (Myrank==0 && Mespr==1)  
		  fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	      }
	    if (mdt>mdtdef)
	      mdt=mdtdef;
	    
	    if (mdt>mdtmax)//maximum time increment
	      mdt = mdtmax;
	  }
	else
	  {
	    //  reduction of the time increment because
	    //  inner loop was not able to enforce equilibrium
	    mefel_convergence = 0;
	    
	    mdt/=2.0;
	    Mp->timecon.oldtime ();
	    mi--;
	    
	    if (Myrank==0 && Mespr==1)  
	      fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
	    if (mdt<mdtmin)
	      {
		if (Myrank==0 && Mespr==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
		if (Myrank==0 && Mespr==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
		if (Myrank==0 && Mespr==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
		compute_req_val (mlcid);
		print_step_forced(mlcid, mi, Mp->time, mt_gv.fb);
		//print_step_forced(mlcid, mi, Mb->dlc[0].gf[0].getval(Mp->time), mt_gv.fb);
		print_flush();
		break;
	      }
	  }
      }
      
      //PMETR cleaning mefel matrices for more memory in trfel part
      if(Cp->cleanmatrix == clean_yes){
	//cleaning of stiffness matrix
	if (Smat != NULL){
	  fprintf (stderr,"\n\n Cleaning of stiffness matrix \n\n ");
	  delete Smat;
	  Smat=NULL;
	}
      }
    }
    
  }while(Cp->time<=end_time);//METR time controler is driving only start and finish of computation
  print_close ();
  print_closet();
  
  if(Tp->timecont.tct > 0){
    delete [] lhsb;
    delete [] tdlhsb;
  }
}


/**
  Function allocates and initializes global vectors used in npsolver
  It is used in one_step concept solver.

  Created by tomas Krejci according to Tomas Koudelka, 02/2013
*/

void par_nonstat_trfel_init (long lcid, np_glob_vec &np_gv)
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

  //  nodes - integration points interpolation
  approximation();  
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (np_gv.lhs, np_gv.tdlhs, np_gv.f, np_gv.istep, Tp->time, dt, Tp->timecont, n);
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (lcid);
    print_initt(-1, "at",Pcp->fnit,Pcp->feit);
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (lcid);
    print_initt(-1, "wt",Pcp->fnit,Pcp->feit);
    print_stept(lcid, np_gv.istep, Tp->time, np_gv.rhs);
  }
  print_flusht();

}

/**
  Function allocates and initializes global vectors used in pmtsolver
  It is used in one_step concept solver.

  TKr, 19/02/2013 according to Tomas Koudelka (mtsolver.cpp)
*/
void par_visco_mefel_init(long lcid, mt_glob_vec &mt_gv)
{
  long n;
  double dt;  
  //  number of mechanical degrees of freedom
  n=Ndofm;

  mt_gv.alloc(n);
  //  vector of nodal displacements
  mt_gv.r = Lsrs->give_lhs (lcid);
  //  vector of nodal prescribed forces
  mt_gv.f = Lsrs->give_rhs (lcid);

  //  initial time
  Mp->time=Mp->timecon.starttime ();

  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat())
  {
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file\n");
    solver_restore (mt_gv.r, mt_gv.fp, mt_gv.istep, Mp->time, dt, &Mp->timecon, n);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    print_init(-1, "at",Pcp->fnim,Pcp->feim);
    //print_step(lcid,i,Mp->time,f);
    print_flush();
  }
  else
  {
    mt_gv.istep=0;
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    if (Mp->eigstrains > 0)
      internal_forces(lcid, mt_gv.fp);
    // print initial step to the output files
    print_init(-1, "wt",Pcp->fnim,Pcp->feim);
    print_step(lcid, mt_gv.istep, Mp->time, mt_gv.f);
    //print_step(lcid, i, 0.0, f);
    print_flush();
  }
}


/**
   function solves parallel linear nonstationary transport problem
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   TKr, 26/02/2013 according to JK
*/
long par_one_step_trfel_linear (long lcid,double time, double dt, long istep, long li, np_glob_vec &np_gv)
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

  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n\n ------------------------------------------------------------------------");
  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n PTRFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Tp->time,dt);
  if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n ------------------------------------------------------------------------\n");
  
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

  compute_req_valt (0);  
  print_stept(0,istep,Tp->time,rhs);
  print_flusht();
  
  if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
    solvert_save (lhs,tdlhs,f,istep,Tp->time,dt,Tp->timecont,n);

  if (j == 0)
    return 0;
  
  return 1;
}


/**
   function solves parallel non-linear nonstationary transport problem  by Newton-Raphson method
   the v-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   TKr, 26/02/2013 according to JK
*/
long par_one_step_trfel_nonlinear (long lcid,double time, double dt, long istep, long li, np_glob_vec &np_gv)
{
  long j,k,n,ini;
  double alpha,*f,*d,*p,*lhs,*tdlhs,*rhs;
  double *fb,*fi;
  double norf_last;
  double zero,norf,*err,*thresh;  
  
  //new time increment
  Tp->time = time;
  //  new step number
  Tp->istep = istep;
  
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
  Kmat->scalgm (dt*alpha);
  Kmat->addgm (1.0,*Cmat);
  Kmat->copygm (*Cmat);

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

  //  physically corrected solution
  solution_correction ();
  //  nodes - integration points interpolation
  approximation ();
  
  // pridal Madera 18.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  //////////////////////////////////////////////////
  norf_last = 1.0e20;
  //inner iteration loop//////////////////
  for (j=0;j<ini;j++){
    
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
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      Kmat->copygm (*Cmat);
    }
    
    if (Tp->trestype==lrhst){
      //  correct right hand side of transport part
      trfel_right_hand_side (lcid,rhs,n);
      // Solver computes residuum from system of equations
      Kmat->gmxv (tdlhs,fi);//new fi vector
      // Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
      for (k=0;k<n;k++){
	fb[k] = rhs[k] - p[k] - fi[k];
      }
    }
    
    if (Tp->trestype==fluxest){
      // Solver computes unbalanced fluxes
      internal_fluxes (fi,n);//new fi vector
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
	fb[k]=fi[k];
      }
    }
    
    norf = Psolt->pss(fb, fb, Outt);//norf = ss (fb,fb,n);

    if (Myrank==0 && Mesprt==1)  fprintf (stdout,"\n inner iteration %ld   error %e",j,norf);
    if (norf<err[0]) break;
    
    Psolt->par_linear_solver (Gtt,Cmat,fi,fb,Outt,Mesprt);//fi is output now
    
    for (k=0;k<n;k++){
      tdlhs[k]+=fi[k];
      lhs[k]+=alpha*dt*fi[k];
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
	
	if (Myrank==0 && Mesprt==1)  
	  fprintf (stdout,"\n\n convergence control: inner iteration skiped %ld error %e\n\n",j,norf);	  
	
	return -2;//break;
      }

      norf_last = norf; //storing norf from previous inner step
    }
  }  

  if(norf >= err[0])
    return -1;
  
  //////////////////////////////////////////////////
  //  printing of output and graphical informations
  compute_req_valt (lcid);
  print_stept(lcid,istep,Tp->time,rhs);//print_stept(0,i,Tp->time,NULL);
  print_flusht();
  
  if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
    solvert_save (lhs,tdlhs,f,istep,Tp->time,dt,Tp->timecont,n);
  
  return 0;
}


/**
  Function executes one step of parallel time dependent mechanical algorithm
  (visco-elasticity, visco-plasticity, etc.)
  
  @param lcid - load case id
  @param time - actual time
  @param dt   - actual time increment
  @param istep - time step id
  @param li    - time step id of the first performed time step
  @param f - %vector of prescribed nodal forces
  @param fi - %vector of computed nodal forces
  @param fp - %vector of previous nodal forces
  @param dr - %vector of displacement increments
  @param r - %vector of total displacements
  @param lhsb - backup of %vector of nodal unknowns
  
  @retval  0 - if the equilibrium was attained with NO inner loop
  @retval  1 - if the equilibrium was attained and inner loop WAS performed
  @retval -1 - if the equilibrium was NOT attained
  
  TKr, 26/02/2013 according to JK and TKo
*/
long par_one_step_mefel (long lcid,double time, double dt, long istep, long li, mt_glob_vec &mt_gv)
{
  long j,n,ini;
  double *f,*fl,*fi,*fb,*fp,*r,*dr,*lhsb;
  double norf, err, zero;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  Mp->nlman->divc_step = new double[10];
  Mp->nlman->divc_err  = new double[10];
  
  //new time increment
  Mp->time = time;
  //  new step number
  Mp->istep = istep;
  //  indicator of strain computation
  //  it means, strains at integration points have not been computed
  Mp->strainstate=0;
  //  indicator of stress computation
  //  it means, stresses at integration points have not been computed
  Mp->stressstate=0;
  //  indicator of computation of other array
  //  it means, stresses at integration points have not been computed
  Mp->otherstate=0;
  
  //  vector of nodal displacements
  r  = mt_gv.r;
  //  vector of nodal prescribed forces
  f  = mt_gv.f;
  //  vector of increments of nodal displacements
  dr = mt_gv.dr;
  //  vector of prescribed forces from the previous step
  fp = mt_gv.fp;
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  fl = mt_gv.fl;
  //  vector of internal forces
  fi = mt_gv.fi;
  //  auxiliary force vector
  fb = mt_gv.fb;
  //  backup of the nodal displacements
  lhsb = mt_gv.lhsb;
  
  //  number of mechanical degrees of freedom
  n=Ndofm;
  
  // maximum number of steps for inner loop
  ini = Mp->nlman->niilnr;
  
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  
  //  computer zero
  zero = Mp->zero;
  
  //  backup of attained nodal values
  for (j=0;j<n;j++)
    lhsb[j]=r[j];
  
  
  if (Myrank==0 && Mespr==1)  fprintf (stdout,"\n\n ------------------------------------------------------------------------");
  if (Myrank==0 && Mespr==1)  fprintf (stdout,"\n PMEFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Mp->time,dt);
  if (Myrank==0 && Mespr==1)  fprintf (stdout,"\n ------------------------------------------------------------------------\n");
  
  //  assembling of vectors of prescribed nodal forces
  mefel_right_hand_side (lcid,f,fl);
  
  //  computation of sums of force components in particular directions
  //  it is used in special creep and consolidation models
  Mb->comp_sum (f);
  //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
  //      istep, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
  fflush(Out);
  
  //  increments of prescribed forces
  incr_internal_forces (lcid,fi);
  
  // assembling of right-hand side
  // fb[j]=f[j]-fp[j]+fi[j]
  subv(f, fp, fb, n);
  addv(fb, fi, fb, n);
  
  norf = par_gnewton_raphson_one_step_mefel(lcid, Mp->nlman, fl, r, fb, dr, fi, istep-1, j, li, ini, err);
  
  if (norf>err)
    {
      // equilibrium state was not attained
      if (Myrank==0 && Mespr==1) 
	fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);
      
      //  backup is used for actual values
      //  r[k]=lhsb[k]
      copyv(lhsb, r, n);
      
      return -1;
    }
  
  // equilibrium state was attained
  if ((Myrank==0 && Mespr==1) && (j > 0))  
    fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
  
  // actual total load vector becomes previous total load vector
  // fp[k]=f[k]
  copyv(f, fp, n);
  
  // print resulting vector of internal forces
  Mb->comp_sum (fi);
  //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
  //      istep, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
  fflush(Out);
  
  //  update of values stored at integration points
  Mm->updateipval();
  compute_req_val (lcid);
  print_step(lcid, istep, Mp->time, fl);
  //print_step(lcid, i, Mb->dlc[0].gf[0].getval(Mp->time), fl);
  print_flush();
  
  if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.save_stat()))
    {
      if (Mespr==1)
	fprintf (stdout,"\n Creating backup file\n");
      
      solver_save (r,fp,istep,Mp->time,dt,&Mp->timecon,n);
    }

  if (j == 0)
    return 0;
  
  return 1;
}


/**
  Function performs calculation of one parallel load/time step of the Newton-Raphson 
  method for the given load case. Solved equation system does not contain 
  time variable.
  
  @param lcid  - load case id
  @param nlman - pointer to structure conatining setup of the solver
  @param ra    - %vector of attained displacements
  @param fa    - attained load %vector
  @param fb    - residual %vector or load %vector increment - right-hand side %vector
  @param dr    - %vector of displacement increment - left-hand side %vector
  @param istep - time/load step id
  @param j     - inner loop step id (output parameter, set to -1 if no inner loop was performed) 
  @param li    - initial value of time/load step id
  @param ierr  - required normed error of residual %vector

  @return The function returns reached lambda parameter.

  TKr, 19/02/2013 according to Tomas Koudelka
*/
double par_gnewton_raphson_one_step_mefel(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, double *dr, double *fi, 
                                long istep, long &j, long li, long ini, double ierr)
{
  long n = Ndofm;
  double norf, norfa;
  matrix lsm_a(3,3); // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  
  //  assembling of tangent stiffness matrix
  assemble_stiffness_matrix(lcid,istep,-1,li,no);
    
  //  solution of K(r).v=F
  Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);//Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

  //  update total displacement vector
  //  ra[j]+=dr[j];
  addv(ra, dr, n);
    
  //  computation of internal forces
  internal_forces (lcid,fi);
    
  //  vector of unbalanced forces
  //  fb[j] = fa[j] - fi[j];
  subv(fa, fi, fb, n);

  // norm of attained load vector
  norfa=Psolm->pss (fa,fa,Out);//normv(fa,n);
  //  norm of vector of unbalanced forces
  norf=Psolm->pss(fb,fb,Out);//normv(fb,n);
  if (norfa != 0.0)
    norf /= norfa;
    
  if ((Myrank==0)&&(Mespr==1)) 
  fprintf (stdout,",    norm %e",norf);

  if (norf<ierr)
  {
    j = -1;
    return norf;
  }
    
  //  internal iteration loop
  fillm(0.0, lsm_a);
  fillv(0.0, lsm_r);
  for (j=0; j<ini; j++)
  {
    Mp->jstep = j;

    // re-assemble stifness matrix if it is required
    assemble_stiffness_matrix(lcid,istep,j,li,no);
    
    // solution of K(r).v=F
    Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);//Mp->ssle->solve_system(Gtm,Smat,dr,fb,Out);
    
    // ra[k]+=dr[k];
    addv(ra, dr, n);

    // computation of internal forces
    internal_forces(lcid,fi);
      
    // vector of unbalanced forces
    // fb[k]=fa[k]-fi[k]
    subv(fa, fi, fb, n);
      
    // norm of vector of unbalanced forces
    norf=Psolm->pss(fb,fb,Out);//normv(fb,n);
    if (norfa != 0.0)
      norf /= norfa;

    if ((Myrank==0) && (Mespr==1))
      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf %e",istep,j,norf);
      
    if (norf < ierr) 
      // equilibrium was attained
      return norf;

    // divergence detection with help of the least square method
    if (check_divergency(nlman, lsm_a, lsm_r, lsm_l, j, norf))
      return norf;
  }

  // equilibrium was not attained in this step
  return norf;  
}



/*************************************** This below is old ****************************************/


/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method

   @param lcid - load case id

   TKo 4.2009
*/
void par_newton_raphson_parcoupl (long lcid)
{
  switch (Tp->tprob){
  case nonstationary_problem:{
    //  linear nonstationary transport problem
    //  nonlinear mechanical problem
    par_newton_raphson_parcoupl_lin (lcid);
    break;
  }
  case nonlinear_nonstationary_problem:{
    //  nonlinear nonstationary transport problem
    //  nonlinear mechanical problem
    par_newton_raphson_parcoupl_nonlin (lcid);
    break;
  }
  default:{
    print_err("unknown METR problem type is required", __FILE__,__LINE__,__func__);
  }
  }
}



/**
   Function solves system of linear TRFEL algebraic and non-linear MEFEL equations by Newton-Raphson method
   for time-dependent problems

   @param lcid - load case id
   
   TKo 4.2009
*/
void par_newton_raphson_parcoupl_lin (long lcid)
{
  int mefel_convergence;
  long i,ii,j,k,nm,nt,ini,nsts,nii;
  double zero,dt,dtmin,dtdef,end_time,alpha,s;
  double *dr,*fb,*r,*d,*f,*fp,*fi,*fl,*p,*lhst,*tdlhst,*rhst,*lhsb,*lhstb,*tdlhstb;
  double norf, norfa, err, lsm_res;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  double alphadtn;         // alpha*dt for the actual time step
  double alphadto = 0.0;   // alpha*dt for the previous time step
  double zerodt = 1.0e-12; // alphadtn and alphadto are assumed to be identical with this tolerance

  //  maximum number of iterations in inner loop
  ini = Mp->nlman->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  //  number of mechanical degrees of freedom
  nm=Ndofm;
  //  number of transport degrees of freedom
  nt=Ndoft;

  /* *************************** */
  /*  vectors of transport part  */
  /* *************************** */

  //  vector of nodal values
  lhst=Lsrst->give_lhs (lcid);
  //  vector of time derivatives of nodal values
  tdlhst=Lsrst->give_tdlhs (lcid);
  //  vector of the right hand side
  rhst=Lsrst->give_rhs (lcid);

  //  array containing predictors
  d = new double [nt];
  nullv (d,nt);
  //  auxiliary vector
  p = new double [nt];
  nullv (p,nt);
  //  backup of the lhs vector
  lhstb = new double [nt];
  nullv (lhstb,nt);
  //  backup of the tdlhst vector
  tdlhstb = new double [nt];	  
  nullv (tdlhstb,nt);

  /* **************************** */
  /*  vectors of mechanical part  */
  /* **************************** */

  //  vector of nodal displacements
  r = Lsrs->give_lhs (lcid);
  //  vector of nodal prescribed forces
  f = Lsrs->give_rhs (lcid);

  //  vector of increments of nodal displacements
  dr = new double [nm];
  nullv (dr,nm);
  //  vector of prescribed forces from the previous step
  fp = new double [nm];
  nullv (fp,nm);
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  fl = new double [nm];
  nullv (fl,nm);
  //  vector of internal forces
  fi = new double [nm];
  nullv (fi,nm);
  //  auxiliary force vector
  fb = new double [nm];
  nullv (fb,nm);
  //  backup of the nodal displacements
  lhsb = new double [nm];
  nullv (lhsb,nm);  

  //  coefficient of trapezoidal method
  alpha=Tp->alpha;
  zero=Tp->zero;
  
  //driving time from METR
  //  starting time
  Cp->time=Cp->timecon.starttime ();
  //  time increment
  dt=Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();

  //  mechanical time controller is rewritten by coupled time controller
  Mp->time = Cp->time;
  Mp->timecon.take_values (Cp->timecon);
  //  minimum time increment used in MEFEL inner loop
  dtmin=Mp->timecon.dtmin;
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);

  //  initial temperature is approximated to the mechanical integration points
  approximation_inittemper();
  //  nodes - integration points interpolation in TRFEL part
  approximation();

  actual_previous_nodval();

  // **************************
  //  main iteration loop  ***
  // **************************  
  i=0;  // number of TRFEL step
  ii=0; // number of MEFEL step
  //  number of successful time steps in MEFEL (number of steps without inner iterations)
  nsts=0;

  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file for MEFEL\n");
    solver_restore (r,fp,ii,Mp->time,dt,&Mp->timecon,Ndofm);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    print_init(-1, "at",Pcp->fnim,Pcp->feim);
//    print_step(lcid,i,Mp->time,f);
//    print_flush();
  }
  else
  {
    // initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    print_init(-1, "wt",Pcp->fnim,Pcp->feim);
    print_step(lcid, i, Mp->time, f);
    print_flush();
  }
  if (Tp->hdbcont.restore_stat()){
    if (Mesprt==1)
      fprintf (stdout,"\n Reading of backup file for TRFEL\n");
    solvert_restore (lhst,tdlhst,f,i,Cp->time,dt,Cp->timecon,Ndoft);
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "at",Pcp->fnit,Pcp->feit);
//    print_stept(0,i,Tp->time,NULL);
//    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt",Pcp->fnit,Pcp->feit);
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }

  mefel_convergence = 1;//setup of flag for successfully performed MEFEL step to yes

  // *********************************************************
  //                     transport part
  // *********************************************************
  do{
    if(mefel_convergence == 0){
      Cp->timecon.oldtime ();
      Tp->timecont.oldtime ();
      if ((Myrank==0) && (Mesprc==1))  fprintf (stdout,"\n\n --------------------------------------------------------------------------------------------");
      if ((Myrank==0) && (Mesprc==1))  fprintf (stdout,"\n REDUCED TIME STEP, METR TIME STEP = %ld : METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
      if ((Myrank==0) && (Mesprc==1))  fprintf (stdout,"\n --------------------------------------------------------------------------------------------");
      i--;
    }

    //  time update
    //new time and new time increment
    Cp->time = Cp->timecon.newtime (dt);
    Tp->time = Tp->timecont.newtime (dt);
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);    
    i++;

    if(mefel_convergence == 1){
      if ((Myrank==0) && (Mesprc==1))  fprintf (stdout,"\n\n -------------------------------------------------------------------------------------------");
      if ((Myrank==0) && (Mesprc==1))  fprintf (stdout,"\n NEW METR TIME STEP = %ld: METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
      if ((Myrank==0) && (Mesprc==1))  fprintf (stdout,"\n -------------------------------------------------------------------------------------------");
    }
    if ((Myrank==0) && (Mesprt==1))  fprintf (stdout,"\n ------------------------------------------------------------------------------------");
    if ((Myrank==0) && (Mesprt==1))  fprintf (stdout,"\n TRFEL TIME STEP = %ld, TRFEL TIME = %e, trfel time increment = %e",i,Tp->time,dt);
    if ((Myrank==0) && (Mesprt==1))  fprintf (stdout,"\n ------------------------------------------------------------------------------------");

    Tm->updateipval ();

    alphadtn = alpha*dt;
    if (Myrank==0)
    {
      fprintf(stdout, "\n adto=%le, adtn=%le, adto/adtn-1=%le, Cmat=%p, Kmat=%p, fcsolv=%d\n", 
              alphadto, alphadtn, fabs(alphadto/alphadtn - 1.0), (void*)Cmat, (void*)Kmat, Cp->fcsolv); 
      fflush(stdout);
    }
    if ((fabs(alphadto/alphadtn - 1.0) > zerodt) || (fabs(alphadtn) < zero) || 
        (Cmat==NULL) || (Kmat == NULL) || (Cp->fcsolv != linearc))
    {
      //  conductivity matrix
      conductivity_matrix (lcid);
      
      //  capacity matrix
      capacity_matrix (lcid);
    }
    //  predictor
    for (j=0;j<nt;j++){
      p[j] = lhst[j] + (1.0-alpha)*dt*tdlhst[j];
    }
    
    //  right hand side of transport part
    trfel_right_hand_side (0,rhst,nt);
    
    //  auxiliary vector
    Kmat->gmxv (p,d);

        
    if ((fabs(alphadto/alphadtn - 1.0) > zerodt) || (fabs(alphadtn) < zero) || 
        (Cmat==NULL) || (Kmat == NULL) || (Cp->fcsolv != linearc))
    {
      //  matrix of the system of equations
      //  C + alpha.dt.K
      //Kmat->scalgm (dt*alpha);
      //Kmat->addgm (1.0,*Cmat);
      // Kmat should be saved for future use, Cmat can be modified
      Cmat->addgm(dt*alpha, *Kmat);
    }
    
    for (j=0;j<nt;j++){
      rhst[j] = rhst[j] - d[j];
      d[j]=tdlhst[j];
    }
    
    //  solution of the system of algebraic equations
    //Psolt->par_linear_solver (Gtt,Kmat,tdlhst,rhst,Outt,Mesprt);
    Psolt->par_linear_solver (Gtt,Cmat,tdlhst,rhst,Outt,Mesprt);

    //  backup of attained nodal values
    for (j=0;j<nt;j++){
      lhstb[j]=lhst[j];
      tdlhstb[j]=tdlhst[j];
    }
    
    //  new nodal values
    for (j=0;j<nt;j++){
      s=(1.0-alpha)*d[j]+alpha*tdlhst[j];
      lhst[j]+=dt*s;
    }

    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();

    actual_previous_nodval ();
    
    if (Cp->fcsolv == fullnewtonc){
      //  cleaning matrices for more free memmory
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
    // *********************************************************
    //  end of transport part and beginning of mechanical part
    // *********************************************************
    mefel_convergence = 1;
    //new time increment
    Mp->time = Mp->timecon.newtime (dt);

    if (Mp->time <= Tp->time){
      //  new step number
      ii++;
      //  indicator of strain computation
      //  it means, strains at integration points have not been computed
      Mp->strainstate=0;
      //  indicator of stress computation
      //  it means, stresses at integration points have not been computed
      Mp->stressstate=0;
      //  indicator of computation of other array
      //  it means, stresses at integration points have not been computed
      Mp->otherstate=0;

      //  approximation of temperatures from TRFEL nodes to MEFEL integration points
      approximation_temper ();
      //  approximation of humidities from TRFEL nodes to MEFEL integration points
      approximation_humid ();

      //  backup of attained nodal values
      for (j=0;j<nm;j++){
	lhsb[j]=r[j];
      }

      if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n\n -------------------------------------------------------------------------------");
      if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n MEFEL TIME STEP = %ld, MEFEL TIME = %e, mefel time increment = %e",ii,Mp->time,dt);
      if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n -------------------------------------------------------------------------------");

      //  assembling of vectors of prescribed nodal forces
      mefel_right_hand_side (lcid,f,fl);

      //  computation of sums of force components in particular directions
      //  it is used in special creep and consolidation models
      Mb->comp_sum (f);
      //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
      //i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      fflush(Out);
            
      //  computation of unbalanced forces which will be on the right hand side
      incr_internal_forces (lcid,fi);
      
      //  right hand side assembling
      for (j=0;j<nm;j++){
	fb[j]=f[j]-fp[j]+fi[j];
      }    

/*      norf = ss(f,f,nm);
      fprintf(stdout, "\n\n norm(f)=%le,", norf);
      norf = ss(fp,fp,nm);
      fprintf(stdout, " norm(fp)=%le,", norf);
      norf = ss(fi,fi,nm);
      fprintf(stdout, " norm(fi)=%le,", norf);
      norf = ss(fb,fb,nm);
      fprintf(stdout, " norm(f-fp+fi)=%le\n\n", norf);*/
      
      //  assembling of tangent stiffness matrix
      if ((Mp->nlman->stmat==tangent_stiff) || (ii == 1) || (Smat==NULL))
	stiffness_matrix (lcid);
      
      //  solution of K(r).v=F
      Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);

      //  new displacement vector
      for (j=0;j<nm;j++){
	r[j]+=dr[j];
      }

      //  computation of internal forces
      internal_forces (lcid,fi);

      // norm of total load vector      
      norfa=Psolm->pss (f,f,Out);
      // vector of unbalanced forces
      for (j=0;j<nm;j++){
	fb[j] = fl[j] - fi[j];
      }
      //  norm of vector of unbalanced forces
      norf=Psolm->pss (fb,fb,Out);
      if (norfa > 0.0)
	norf /= norfa;
      
      if ((Myrank==0)&&(Mespr==1))
        fprintf (stdout,"\n\n Norf before inner iteration loop norf = %e\n\n\n",norf);
      j = 0;
      if (Psolm->compare_on_master(norf, err) <= 0) // norf <= err
      {
        //  number of successful time steps
        nsts++;
        //  no inner iteration
        nii=1;
      }
      else
        nii=0;

      if (nii==0){
        // iteration of unbalanced forces caused by time independent models
        //  internal iteration loop
        fillm(0.0, lsm_a);
        fillv(0.0, lsm_r);
	for (j=0;j<ini;j++)
	{	    	    
          if((Cp->fcsolv == fullnewtonc) || 
             ((Mp->nlman->stmat==tangent_stiff) && (j%5 == 0)))
          {
            //assembling of stiffness matrix for fullnewton
            stiffness_matrix (lcid);
          }
          //  solution of K(r).v=F
          Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);

          // new displacement vector
          for (k=0;k<nm;k++)
            r[k]+=dr[k];

          //  computation of internal forces
          internal_forces (lcid,fi);
	    
          //  vector of unbalanced forces
          for (k=0;k<nm;k++)
            fb[k]= fl[k] - fi[k];
	    
          //  norm of vector of unbalanced forces
          norf=Psolm->pss (fb,fb,Out);
          if (norfa > 0.0)
            norf /= norfa;

//          if ((Myrank==0)&&(Mespr==1)&&(j%10==0)) // print each 10-th step
          if ((Myrank==0)&&(Mespr==1)) 
            fprintf (stdout,"\n Inner loop number j = %ld,   norf = %e \n",j,norf);
          if (Psolm->compare_on_master(norf, err) <= 0) // norf <= err, divergence detection
          {
            // convergence attained
            // number of successful time steps
            nsts++;
            //  stop inner iteration
            nii=1;
            break;
          }
          else
            nii=0;

          // divergence detection with help of least square method
          if (j > 1)
          {
            lsm_res = lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero,1);
            if (Psolm->compare_on_master(lsm_res, 0.0) > 0) // lsm_res > 0.0, divergence detection
            {
              if ((Myrank==0) && (Mespr==1))
                fprintf (stdout,"\n\n Divergence of inner iteation loop is detected. Inner loop is stopped.\n\n");
              //  stop inner iteration
              nii=0;
              break;
            }
          }
          else
            lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero,0);
        }
      }

      if(Cp->fcsolv == fullnewtonc){
	//cleaning of stiffness matrix
	if (Smat != NULL){
	  delete Smat;
	  Smat=NULL;
	}
      }

      if (nii==0)
      {
        // equilibrium state was not attained
        if ((Myrank==0) && (Mespr==1))
          fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);
        mefel_convergence = 0;

        //  backup is used for actual MEFEL values
        for (j=0;j<nm;j++)
          r[j]=lhsb[j];

        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
        dt/=2.0;
        Mp->timecon.oldtime ();
        ii--;

        if ((Myrank==0) && (Mespr==1))
          fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
        if (Psolm->compare_on_master(dt, dtmin) < 0){ // dt<dtmin
          if ((Myrank==0) && (Mespr==1)) fprintf (stderr,"\n\n time increment is less than minimum time increment");
          if ((Myrank==0) && (Mespr==1)) fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
          compute_req_val (lcid);
          print_step_forced(lcid, i, Mp->time, f);
          print_flush();
          break;
        }
      }
      else
      {
        // equilibrium state was attained
        if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
        mefel_convergence = 1;

        //  actual total load vector becomes previous total load vector
        for (j=0;j<nm;j++)
          fp[j]=f[j];

        dtdef = Mp->timecon.actualforwtimeincr ();
        if (nsts==2){
          dt*=2.0;
          nsts=0;
          if ((Myrank==0) && (Mespr==1)) fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
        }
        if (dt>dtdef)
          dt=dtdef;

        // print resulting vector of internal forces
        Mb->comp_sum (fi);
        //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
	//      i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
        fflush(Out);

        //  update of values stored at integration points
        Mm->updateipval();
        compute_req_val (lcid);
        print_step(lcid, ii, Mp->time, f);
        print_flush();

        if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.hdbtype))
        {
          if (Mespr==1)
            fprintf (stdout,"\n Creating backup file for MEFEL\n");
          solver_save (r,f,ii,Mp->time,dt,&Mp->timecon,Ndofm);
        }
      }
    }

    if(mefel_convergence == 0){
      //  backup is used for actual TRFEL values
      for (j=0;j<nt;j++){
        lhst[j]=lhstb[j];
        tdlhst[j]=tdlhstb[j];
      }
    }
    else{
      alphadto = alphadtn;
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      if ((Tp->timecont.isitimptime()==1) && Tp->hdbcont.save_stat())
      {
        if (Mesprt==1)
          fprintf (stdout,"\n Creating backup file for TRFEL\n");
        solvert_save (lhst,tdlhst,rhst,i,Cp->time,dt,Cp->timecon,Ndoft);
      }
    }
  } while(Cp->time<=end_time);
  
  print_close ();
  print_closet ();
  
  delete [] fi;
  delete [] fb;
  delete [] fp;
  delete [] dr;
  delete [] fl;
  delete [] lhsb;

  delete [] p;
  delete [] d;

  delete [] lhstb;
  delete [] tdlhstb;
}



/**
   Function solves system of non-linear TRFEL algebraic and non-linear MEFEL equations by Newton-Raphson method
   for time-dependent problems.

   @param lcid - load case id
   
   TKr 5.4.2007
*/
void par_newton_raphson_parcoupl_nonlin (long lcid)
{
  int mefel_convergence;
  long i,ii,j,k,nm,nt,ini,init,nsts,nii;
  double zero,dt,dtmin,dtdef,end_time,alpha;
  double *dr,*fb,*r,*d,*f,*fp,*fi,*fl,*p,*lhst,*tdlhst,*rhst,*lhsb,*lhstb,*tdlhstb;
  double *fbt,*fit,*lhst_last;
  //double s;
  double norf_last, lsm_res;
  double norf, norfa,err,errt;  
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);

  //  maximum number of iterations in inner loop
  ini = Mp->nlman->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  //  maximum number of iterations in inner loop
  init = Tp->nii;
  //  required norm of vector of unbalanced forces
  errt = Tp->errarr[0];//docasne, nutne doplnit??!!
  //  number of mechanical degrees of freedom
  nm=Ndofm;
  //  number of transport degrees of freedom
  nt=Ndoft;

  /* *************************** */
  /*  vectors of transport part  */
  /* *************************** */
  //  left hand side
  lhst=Lsrst->give_lhs (lcid);
  //  vector of time derivatives of nodal values
  tdlhst=Lsrst->give_tdlhs (lcid);
  //  right hand side
  rhst=Lsrst->give_rhs (lcid);
  
  //  array containing predictors
  d = new double [nt];
  nullv (d,nt);
  //  auxiliary vector
  p = new double [nt];
  nullv (p,nt);
  fbt = new double [nt];
  nullv (fbt,nt);
  fit = new double [nt];
  nullv (fit,nt);
  lhst_last = new double [nt];
  
  //  backup of the lhst vector
  lhstb = new double [nt];
  nullv (lhstb,nt);
  //  backup of the tdlhst vector
  tdlhstb = new double [nt];	  
  nullv (tdlhstb,nt);

  /* **************************** */
  /*  vectors of mechanical part  */
  /* **************************** */

  //  vector of nodal displacements
  r = Lsrs->give_lhs (lcid);
  //  vector of nodal prescribed forces
  f = Lsrs->give_rhs (lcid);
  
  //  vector of increments of nodal displacements
  dr = new double [nm];
  nullv (dr,nm);
  //  vector of prescribed forces from the previous step
  fp = new double [nm];
  nullv (fp,nm);
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  fl = new double [nm];
  nullv (fl,nm);
  //  vector of internal forces
  fi = new double [nm];
  nullv (fi,nm);
  //  auxiliary force vector
  fb = new double [nm];
  nullv (fb,nm);
  //  backup of the nodal displacements
  lhsb = new double [nm];
  nullv (lhsb,nm);

  //  coefficient of trapezoidal method
  alpha=Tp->alpha;
  zero=Tp->zero;

  //driving time from METR
  //  starting time
  Cp->time=Cp->timecon.starttime ();
  //  time increment
  dt=Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();
  
  //  mechanical time controller is rewritten by coupled time controller
  Mp->time = Cp->time;
  Mp->timecon.take_values (Cp->timecon);
  //  minimum time increment used in MEFEL inner loop
  dtmin=Mp->timecon.dtmin;
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);

  //  initial temperature is approximated to the mechanical integration points
  approximation_inittemper ();
  
  //  nodes - integration points interpolation in TRFEL part
  approximation ();
  
  actual_previous_nodval ();

  // **************************
  //  main iteration loop   ***
  // **************************
  i=0;  // number of TRFEL step
  ii=0; // number of MEFEL step
  //  number of successful time steps in MEFEL (number of steps without inner iterations)
  nsts=0;
  
  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file for MEFEL\n");
    solver_restore (r,fp,ii,Mp->time,dt,&Mp->timecon,Ndofm);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    print_init(-1, "at",Pcp->fnim,Pcp->feim);
//    print_step(lcid,i,Mp->time,f);
//    print_flush();
  }
  else
  {
    // initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    print_init(-1, "wt",Pcp->fnim,Pcp->feim);
    print_step(lcid, i, Mp->time, f);
    print_flush();
  }
  if (Tp->hdbcont.restore_stat()){
    if (Mesprt==1)
      fprintf (stdout,"\n Reading of backup file for TRFEL\n");
    solvert_restore (lhst,tdlhst,f,i,Cp->time,dt,Cp->timecon,Ndoft);
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "at",Pcp->fnit,Pcp->feit);
    //    print_stept(0,i,Tp->time,NULL);
    //    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt",Pcp->fnit,Pcp->feit);
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }

  mefel_convergence = 1;//setup of flag for successfully performed MEFEL step to yes

  // *********************************************************
  //                     transport part
  // *********************************************************
  do{
    if(mefel_convergence == 0){
      Cp->timecon.oldtime ();
      Tp->timecont.oldtime ();
      if (Mesprc==1)  fprintf (stdout,"\n\n --------------------------------------------------------------------------------------------");
      if (Mesprc==1)  fprintf (stdout,"\n REDUCED TIME STEP, METR TIME STEP = %ld : METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
      if (Mesprc==1)  fprintf (stdout,"\n --------------------------------------------------------------------------------------------");
      i--;
    }

    //  time update
    //new time and new time increment
    Cp->time = Cp->timecon.newtime (dt);
    Tp->time = Tp->timecont.newtime (dt);
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);    
    i++;

    if(mefel_convergence == 1){
      if (Mesprc==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------------------");
      if (Mesprc==1)  fprintf (stdout,"\n NEW METR TIME STEP = %ld: METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
      if (Mesprc==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------------------");
    }
    if (Mesprt==1)  fprintf (stdout,"\n ------------------------------------------------------------------------------------");
    if (Mesprt==1)  fprintf (stdout,"\n TRFEL TIME STEP = %ld, TRFEL TIME = %e, trfel time increment = %e",i,Tp->time,dt);
    if (Mesprt==1)  fprintf (stdout,"\n ------------------------------------------------------------------------------------");
   
    Tm->updateipval ();

    //  conductivity matrix
    conductivity_matrix (lcid);
    
    //  capacity matrix
    capacity_matrix (lcid);
    
    //  predictor
    for (j=0;j<nt;j++){
      p[j] = lhst[j] + (1.0-alpha)*dt*tdlhst[j];
    }
    
    //  right hand side of transport part
    trfel_right_hand_side (0,rhst,nt);
    
    //  auxiliary vector
    Kmat->gmxv (p,d);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    for (j=0;j<nt;j++){  // !!!!!!!!! lisi se od linearni verze
      //verze1:
      rhst[j] = rhst[j] - d[j];
      //verze2:
      //d[j]=tdlhst[j];
    }
    
    //  solution of the system of algebraic equations
    Psolt->par_linear_solver (Gtt,Kmat,tdlhst,rhst,Outt,Mesprt);

    //  backup of attained nodal values
    for (j=0;j<nt;j++){
      lhstb[j]=lhst[j];
      tdlhstb[j]=tdlhst[j];
    }

    //  new nodal values
    for (j=0;j<nt;j++){  // !!!!!!!!! lisi se od linearni verze
      //verze1:
      //s=(1.0-alpha)*d[j]+alpha*tdlhst[j];
      //lhst[j]+=dt*s;
      //verze2:
      lhst[j] = p[j] + alpha*dt*tdlhst[j];
    }
    
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();
    
    norf_last = 1.0e20;    
    //inner iteration loop for trfel
    for (j=0;j<init;j++){      
      //  physically corrected solution
      solution_correction ();
      //  approximation of nodal values into ontegration points
      approximation ();
      
      // full newton-raphson
      if (Tp->trsolv == fullnewtont){

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
	//  right hand side of transport part
	trfel_right_hand_side (0,rhst,nt);
	// Solver computes residuum from system of equations
	Kmat->gmxv (tdlhst,fit);
	// Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
	for (k=0;k<nt;k++){
	  fbt[k] = rhst[k] - d[k] - fit[k];
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
      
      //  norm of vector of unbalanced fluxes
      norf=Psolt->pss (fbt,fbt, Out);
      
      if ((Myrank==0)&&(Mesprt==1))  
        fprintf (stdout,"\n TRFEL inner loop number j = %ld,   norf = %e",j,norf);
      
      if (Psolt->compare_on_master(norf, errt) <= 0) // norf <= err
        break;
      
      Psolt->par_linear_solver (Gtt,Kmat,d,fbt,Outt,Mesprt);
      
      for (k=0;k<nt;k++){
	tdlhst[k]+=d[k];
	lhst[k]+=alpha*dt*d[k];
      }
      
      //  physically corrected solution
      solution_correction ();
      //  approximation of nodal values into ontegration points
      approximation ();
      
      //  convergence control for newton-raphson
      //  condition of decreesing of error
      if (Tp->convergcontrolt==yes){
	if (Psolt->compare_on_master(norf, norf_last)>0){ // norf > norf_last
	  for (k=0;k<nt;k++){
	    lhst[k] = lhst_last[k]; //storing results from previous inner step
	  }
	  
	  //  physically corrected solution
	  solution_correction ();
	  //  approximation of nodal values into ontegration points
	  approximation ();	  

          if ((Myrank==0)&&(Mesprt==1))  fprintf (stdout,"\n\nTRFEL convergence control: inner iteration skipped %ld error %e\n\n",j,norf);
	  
	  break;
	}
	for (k=0;k<nt;k++){
	  lhst_last[k] = lhst[k]; //storing results from previous inner step
	}	
	norf_last = norf; //storing norf from previous inner step
      }
    }

    if (Cp->fcsolv == fullnewtonc){
      //cleaning matrices for more free memmory
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

    // *********************************************************
    //  end of transport part and beginning of mechanical part
    // *********************************************************
    mefel_convergence = 1;
    //new time increment
    Mp->time = Mp->timecon.newtime (dt);

    if (Mp->time <= Tp->time){
      //  new step number
      ii++;
      //  indicator of strain computation
      //  it means, strains at integration points have not been computed
      Mp->strainstate=0;
      //  indicator of stress computation
      //  it means, stresses at integration points have not been computed
      Mp->stressstate=0;
      //  indicator of computation of other array
      //  it means, stresses at integration points have not been computed
      Mp->otherstate=0;

      //  approximation of temperatures from TRFEL nodes to MEFEL integration points
      approximation_temper ();
      //  approximation of humidities from TRFEL nodes to MEFEL integration points
      approximation_humid ();
            
      //  backup of attained nodal values
      for (j=0;j<nm;j++){
	lhsb[j]=r[j];
      }
      
      if (Mespr==1)  fprintf (stdout,"\n\n -----------------------------------------------------------------------------------");
      if (Mespr==1)  fprintf (stdout,"\n MEFEL TIME STEP = %ld, MEFEL TIME = %e, mefel time increment = %e",ii,Mp->time,dt);
      if (Mespr==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------");
            
      //  assembling of the right hand side (prescribed load)
      mefel_right_hand_side (lcid,f,fl);
      
      //  computation of sums of force components in particular directions
      //  it is used in special creep and consolidation models
      Mb->comp_sum (f);
      fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
              i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      fflush(Out);
      
      //  computation of unbalanced forces which will be on the right hand side
      incr_internal_forces (lcid,fi);
      
      //  right hand side assembling
      for (j=0;j<nm;j++){
	fb[j]=f[j]-fp[j]+fi[j];
      }    
      
/*    norf = ss(f,f,nm);
      fprintf(stdout, "\n\n norm(f)=%le,", norf);
      norf = ss(fp,fp,nm);
      fprintf(stdout, " norm(fp)=%le,", norf);
      norf = ss(fi,fi,nm);
      fprintf(stdout, " norm(fi)=%le,", norf);
      norf = ss(fb,fb,nm);
      fprintf(stdout, " norm(f-fp+fi)=%le\n\n", norf);*/

      //  assembling of tangent stiffness matrix
      if ((Mp->nlman->stmat==tangent_stiff) || (ii == 1) || (Smat==NULL))
	stiffness_matrix (lcid);
      
      //  solution of K(r).v=F
      Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);

      //  new displacement vector
      for (j=0;j<nm;j++){
	r[j]+=dr[j];
      }
      
      //  computation of internal forces
      internal_forces (lcid,fi);

      // norm of total load vector      
      norfa=Psolm->pss (f,f,Out);
      // vector of unbalanced forces
      for (j=0;j<nm;j++){
	fb[j] = fl[j] - fi[j];
      }
      //  norm of vector of unbalanced forces
      norf=Psolm->pss (fb,fb,Out);
      if (norfa > 0.0)
	norf /= norfa;
      
      if ((Myrank==0)&&(Mespr==1))  fprintf (stdout,"\n\n Norf before inner iteration loop norf = %e\n\n\n",norf);
      j = 0;
      if (Psolm->compare_on_master(norf, err) <= 0) // norf <= err
      {
        //  number of successful time steps
        nsts++;
        //  no inner iteration
        nii=1;
      }
      else
        nii=0;

      if (nii==0){
        //  iteration of unbalanced forces caused by time independent models
        //  internal iteration loop
        fillm(0.0, lsm_a);
        fillv(0.0, lsm_r);
	for (j=0;j<ini;j++)
        {
          if((Cp->fcsolv == fullnewtonc) || 
             ((Mp->nlman->stmat==tangent_stiff) && (j%5 == 0)))
          {
            //assembling of stiffness matrix for fullnewton
            stiffness_matrix (lcid);
          }
          //  solution of K(r).v=F
          Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);
	    
          // new displacement vector
          for (k=0;k<nm;k++)
            r[k]+=dr[k];

          //  computation of internal forces
          internal_forces (lcid,fi);
	    
          //  vector of unbalanced forces
          for (k=0;k<nm;k++)
            fb[k]= fl[k] - fi[k];
	    
          //  norm of vector of unbalanced forces
          norf=Psolm->pss (fb,fb,Out);
          if (norfa > 0.0)
            norf /= norfa;
	    
//          if ((Mespr==1)&&(Mespr==1)&&(j%10==0)) // print each 10-th step
          if ((Myrank==0)&&(Mespr==1)) 
            fprintf (stdout,"\n Inner loop j =  %ld     norf = %e\n",j,norf);

          if (Psolm->compare_on_master(norf, err) <= 0) // norf <= err, divergence detection
          {
            // convergence attained
            // number of successful time steps
            nsts++;
            //  stop inner iteration
            nii=1;
            break;
          }
          else
            nii=0;	

          // divergence detection with help of least square method
          if (j > 1)
          {
            lsm_res = lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero,1);
            if (Psolm->compare_on_master(lsm_res, 0.0) > 0) // lsm_res > 0.0, divergence detection
            {
              if (Myrank==0)
                fprintf (stdout,"\n\n Divergence of inner iteation loop is detected. Inner loop is stopped.\n\n");
              //  stop inner iteration
              nii=0;
              break;
            }
          }
          else
            lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero,0);
        }
      }

      if(Cp->fcsolv == fullnewtonc){
	//cleaning of stiffness matrix
	if (Smat != NULL){
	  delete Smat;
	  Smat=NULL;
	}
      }
      
      if (nii==0)
      {
        // equilibrium state was not attained
        if ((Myrank==0) && (Mespr==1))
          fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);

        mefel_convergence = 0;

        //  backup is used for actual MEFEL values
        for (j=0;j<nm;j++)
          r[j]=lhsb[j];

        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
        dt/=2.0;
        Mp->timecon.oldtime ();
        ii--;

        if ((Myrank==0) && (Mespr==1))
          fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
        if (Psolm->compare_on_master(dt, dtmin) < 0){ // dt<dtmin
          if ((Myrank==0) && (Mespr==1)) fprintf (stderr,"\n\n time increment is less than minimum time increment");
          if ((Myrank==0) && (Mespr==1)) fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
          compute_req_val (lcid);
          print_step_forced(lcid, i, Mp->time, f);
          print_flush();
          break;
        }
      }
      else
      {
        // equilibrium state was attained
        if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
        mefel_convergence = 1;

        //  actual total load vector becomes previous total load vector
        for (j=0;j<nm;j++)
          fp[j]=f[j];

        dtdef = Mp->timecon.actualforwtimeincr ();
        if (nsts==2){
          dt*=2.0;
          nsts=0;
          if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
        }
        if (dt>dtdef)
          dt=dtdef;

        // print resulting vector of internal forces
        Mb->comp_sum (fi);
        fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
                i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
        fflush(Out);

        //  update of values stored at integration points
        Mm->updateipval();
        compute_req_val (lcid);
        print_step(lcid, ii, Mp->time, f);
        print_flush();

        if ((Mp->timecon.isitimptime ()==1) && (Mp->hdbcont.hdbtype))
        {
          if (Mespr==1)
            fprintf (stdout,"\n Creating backup file for MEFEL\n");
          
          solver_save (r,f,ii,Mp->time,dt,&Mp->timecon,Ndofm);
        }
      }
    }
    if(mefel_convergence == 0){
      //  backup is used for actual TRFEL values
      for (j=0;j<nt;j++){
	lhst[j]=lhstb[j];
	tdlhst[j]=tdlhstb[j];
      }
    }
    else{
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat()){
	if (Mesprt==1)
	  fprintf (stdout,"\n Creating backup file for TRFEL\n");
        solvert_save (lhst,tdlhst,rhst,i,Cp->time,dt,Cp->timecon,Ndoft);
      }
    }
  }while(Cp->time<=end_time);
  
  print_close ();
  print_closet ();
  
  delete [] fi;
  delete [] fb;
  delete [] fp;
  delete [] dr;
  delete [] fl;
  delete [] lhsb;

  delete [] p;
  delete [] d;

  delete [] fit;
  delete [] fbt;
  delete [] lhst_last;

  delete [] lhstb;
  delete [] tdlhstb;
}




/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   TKo, 27/08/2014
*/
void par_newton_raphson_parcoupl_common_dt (long lcid)
{
  long mefel_convergence,trfel_convergence,stop;
  double dt,end_time;
  long i,ti,tli,tnsts, btnsts, tret;
  double tnewtime, tdt;
  np_glob_vec np_gv;
  long tlcid = 0;//zatim??!!
  long mi, mli, mnsts, bmnsts, mret;
  double mnewtime, mdt, mdtr;
  mt_glob_vec mt_gv;
  long mlcid = 0;//zatim??!!
  long rest_calcm, rest_calct;
  double prev_timem, prev_timet;
  vector gv;
  vector grhs;

  //
  //  METR initialization phase
  //
  //  starting time
  Cp->time = Cp->timecon.starttime ();
  //  time increment
  dt = Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();  

  Mp->time = Cp->time;
  Tp->time = Cp->time;
  //this part below is not neccessary because of independent time controlers of TRFEL and MEFEL
  //  mechanical time controller is rewritten by coupled time controller
  //Mp->timecon.take_values (Cp->timecon);
  //  transport time controller is rewritten by coupled time controller
  //Tp->timecont.take_values (Cp->timecon);
  // end this part


  //  initial approximation of quantities from TRFEL nodes to MEFEL integration points
  initapproximation();
  init_trfel_mefel ();
  approximation();

  //
  //  TRFEL initialization phase
  //
  par_nonstat_solver_init(tlcid, rest_calct, np_gv);
  //  initial time increment
  tdt = Tp->timecont.initialtimeincr ();
  //  number of step
  tli = ti = np_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  tnsts = btnsts = 0;
  // return status from the one_step function
  tret = 0;
  // store initial value for attained time in the previous time step
  // (prev_time is exploited in the growing transport problems)
  prev_timet = Tp->time;
  
  if (Cp->dpt == pass_by_aux_ip){
    aip_approximation(Mm->tnip, MTipmap);
    actualize_aip_nonmechq(Tm->tnip, TMipmap);
  }

  //
  //  MEFEL initialization phase
  //
  par_visco_solver_init(mlcid, rest_calcm, mt_gv);

  //  initial time increment
  mdt = Mp->timecon.initialtimeincr ();
  //  number of step
  mli = mi = mt_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  mnsts = bmnsts = 0;
  // return status from the one_step function
  mret = 0;
  // store initial value for attained time in the previous time step
  // (prev_time is exploited in the growing mechanical problems)
  prev_timem = Mp->time;
  
  // ***************************
  //  main iteration loop  ****
  // ***************************

  // number of METR step
  i=0;  
  //setup of flag for successfully performed MEFEL step to yes
  mefel_convergence = 1;
  //setup of flag for successfully performed TRFEL step to yes
  trfel_convergence = 1;
  //  auxiliary indicator of convergence
  stop=0;
  
  //saving values for the first iteration
  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    Tt->lhs_save (np_gv.lhs,Lsrst->lhsi,np_gv.tdlhs);
  if (Mp->tprob == growing_mech_structure)  
    Mt->save_nodval(lcid);

  //  METR is copying data between TRFEL and MEFEL in both directions
  pass_coup_data(lcid);

  do{
    i++;

    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Mp->timecon.timefun.tfunc != constant) && Mp->timecon.tct == 0)
      mdt = Mp->timecon.actualbacktimeincr ();

    //  TRFEL time
    tnewtime = Tp->time = Tp->timecont.newtime (tdt);
    //  MEFEL time
    mnewtime = Mp->time = Mp->timecon.newtime (mdt);
    Mp->dtime = mdt;


    if (Tp->time > Mp->time)
      Cp->time = Mp->time;
    else
      Cp->time = Tp->time;
    
    
    //  time synchronization
    if (Tp->time > Mp->time){
      Tp->timecont.oldtime ();
      tnewtime = Tp->time = Tp->timecont.newtime (mdt);
      tdt = mdt; // actual TRFEL dt for linear solvers (due to recovery of decreased time step in MEFEL)
      if (btnsts)
        tnsts = btnsts;
    }
    
    // (prev_time is exploited in the growing transport problems)
    if (((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin)) &&
        (prev_timet == Cp->timecon.starttime()))
      prev_timet = Tp->time;
    
    //  METR is copying data
    //  copy_data();
    //  copying of mechanical data to TRFEL
    //  mefel_trfel();
    //mefel_trfel_by_nodes_comp ();  //debug??

    //  TRFEL solution
    trfel_convergence = 0;

    //for testing:
    //fprintf (Outt,"\n\n --------------------------------------------------------------");
    //fprintf (Outt,"\n TRFEL Time step = %ld,  Time %e,  Time increment = %e",i,Tp->time,tdt);
    //fprintf (Outt,"\n --------------------------------------------------------------\n");

    ti++;
      
    if((Tp->tprob == nonstationary_problem) || (Tp->tprob == growing_np_problem)){
      // linear algorithm
      // perform one time step with linear solver
      tret = par_one_step_linear(tlcid, tnewtime, tdt, ti, tli, np_gv);
      trfel_convergence = 1;
    }

    if(Tp->tprob == discont_nonstat_problem){
      print_err("the solver has not yet been implemented", __FILE__, __LINE__, __func__);
      // perform one time step with linear solver
      // tret = one_step_linear(tlcid, tnewtime, tdt, ti, tli, np_gv);
      // trfel_convergence = 1;
    }

    if(Tp->tprob == discont_nonlin_nonstat_problem){
      // perform one time step with nonlinear solver for discontinuities dform
      print_err("the solver has not yet been implemented", __FILE__, __LINE__, __func__);
      //if (Myrank == 0){
      //  // global vectors of coarse problem for used computation of selective norm of residual vector
      //  reallocv();
      //  reallocv();
      //}
      //tret = one_step_nonlinear_dform(tlcid, tnewtime, tdt, ti, tli, np_gv);
      trfel_convergence = 1;
    }
      
    if((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin)){
      // non-linear algorithm
      // perform one time step with non-linear solver
      if (Myrank == 0){
        // global vectors of coarse problem for used computation of selective norm of residual vector
        reallocv(Psolt->schcom->ndofcp, gv);
        reallocv(Psolt->schcom->ndofcp, grhs);
      }
      tret = par_one_step_nonlinear(tlcid, tnewtime, tdt, prev_timet, rest_calct, ti, tli, np_gv, gv.a, grhs.a);
      trfel_convergence = 1;
    }
      
    //  handling with time controler for nonlinear solver
    if ((((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin)) || 
         (Tp->tprob == discont_nonlin_nonstat_problem))  && Tp->timecont.tct > 0)
    {
      if (tret >= 0){
        // equilibrium has been attained
	
        //  approximation of nodal values into ontegration points
        approximation ();
        if (Cp->dpt == pass_by_aux_ip)
          aip_approximation(Mm->tnip, MTipmap);
        Tm->updateipval();
	//Tm->freezing_thawing ();

        trfel_convergence = 1;
	  
        if ((tret == 0) || (tret < Tp->nii/4))
          tnsts++;
        else
          tnsts = btnsts = 0; 
     
        if (tnsts>1){
          if (tdt < Tp->timecont.dtmax){
            tdt*=2.0;
            if (tdt > Tp->timecont.dtmax)
              tdt = Tp->timecont.dtmax;
            tnsts=0;
            btnsts=1;
	    
            if ((Mesprt==1) && (Myrank == 0))  
              fprintf (stdout,"\n\n TRFEL time increment may be enlarged because no inner loop was neccessary in previous steps");
          }
        }
      }
      else{
        // equilibrium has not been attained
        tnsts = btnsts = 0;

        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
        trfel_convergence = 0;
	  
        tdt/=2.0;
        ti--;
	  
        if(Tp->timecont.tct > 0) // for adaptive time controler
        {
          //restoring results from previous time step
          copyv(np_gv.lhsb, np_gv.lhs, Ndoft);
          copyv(np_gv.tdlhsb, np_gv.tdlhs, Ndoft);
        }
        //  approximation of nodal values into integration points
        approximation ();
        if (Cp->dpt == pass_by_aux_ip)
          aip_approximation(Mm->tnip, MTipmap);
	  
        if ((Mesprt==1) && (Myrank == 0)) 
          fprintf (stdout,"\n\n TRFEL time increment is reduced to dt=%le because the inner loop was not able to enforce equilibrium", tdt);
        
        if (tdt<Tp->timecont.dtmin){
          if ((Mesprt==1) && (Myrank == 0))  fprintf (stderr," TRFEL time increment is less than the minimum time increment");
          if ((Mesprt==1) && (Myrank == 0))  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
          if ((Mesprt==1) && (Myrank == 0))  fprintf (stderr,"\n FORCED output of results from this step is performed\n");

          compute_req_valt (lcid);
          print_stept_forced(lcid, ti, Tp->time, np_gv.rhs);
          print_flusht();
	  
          compute_req_val (mlcid);
          print_step_forced(mlcid, mi, Mp->time, mt_gv.fb);
          print_flush();

          stop=1;
          break;
        }
	  
        Tp->timecont.oldtime ();
        Mp->timecon.oldtime ();

        // change status of elements and nodes according to the previous time step
        if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
        {
          update_elemnod_prev_timet(prev_timet, np_gv.ncd, np_gv.nce);
        }
	
        continue;
      }
    }

    
    if (stop==1){
      //  TRFEL requires smaller time step than the minimum time step
      //  the code will be terminated
      break;
    }
    
    //  time synchronization
    if (Mp->time > Tp->time){
      Mp->timecon.oldtime ();
      mdt=Tp->time - Mp->timecon.time;
      mnewtime = Mp->timecon.newtime (mdt);
      Mp->time = mnewtime;
      Mp->dtime = mdt;
      if (bmnsts == 1)
        mnsts = bmnsts;
    }
    
    //  METR cleans TRFEL matrices for more memory in MEFEL part
    if (Cp->cleanmatrix == clean_yes){
      //  cleaning matrices for more free memmory
      //  cleaning of the conductivity matrix
      if (Kmat != NULL){
	fprintf (stdout,"\n\n Cleaning of conductivity matrix \n\n ");
	delete Kmat;
	Kmat=NULL;
      }
      //  cleaning of the capacity matrix
      if (Cmat != NULL){
	fprintf (stdout,"\n\n Cleaning of capacity matrix \n\n ");
	delete Cmat;
	Cmat=NULL;
      }
    }
    //print_stept(lcid,Tp->istep,Tp->time,np_gv.rhs);//print_stept(0,i,Tp->time,NULL);//tady??!!
    //  MEFEL solution
    mefel_convergence = 0;
    //  new mefel step number
    mi++;
      
    
    //fprintf (Out,"\n\n -------------------------------------------------------------------------------------------");//debug??!!
    //fprintf (Out,"\n NEW MEFEL TIME STEP = %ld: MEFEL TIME = %e, mefel time increment = %e",i,Mp->time,mdt);
    //fprintf (Out,"\n -------------------------------------------------------------------------------------------\n");
    //fflush(Out);

    //  copying of data
    //  copying of transport data to MEFEL
    trfel_mefel ();
      
    if ((Mp->tprob == growing_mech_structure) && (prev_timem == Cp->timecon.starttime()))
      prev_timem = Mp->time;
    //  perform one time step
    mret = par_one_step (mlcid, mnewtime, mdt, mdtr, prev_timem, rest_calcm, mi, mli, mt_gv);
      
    if (mret >= 0){
      //  equilibrium has been attained
	
      if ((mret == 0) || (mret < Mp->nlman->niilnr/4))
        mnsts++;      
      else
        mnsts = bmnsts = 0;
	
      if (mnsts > 1){
        if (mdt < Mp->timecon.dtmax){
          mdt*=2.0;
          if (mdt > Mp->timecon.dtmax)
            mdt = Mp->timecon.dtmax;
          mnsts=0;
          if ((Mespr==1) && (Myrank == 0))  
            fprintf (stdout,"\n\n MEFEL time increment is enlarged because no inner loop was neccessary in previous steps\n");
        }
      }
      //surface_fluxes (Outt); //printing in TRFEL part moved into outdriver
    }
    else{
      // equilibrium has not been attained
      mnsts = bmnsts = 0;
      
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      mefel_convergence = 0;
      
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      if (mdtr < 1.0)
        mdt *= mdtr;
      else
        mdt/=2.0;
      
      Mp->timecon.oldtime ();
      Tp->timecont.oldtime ();
      mi--;
      ti--;	
      if(Tp->timecont.tct > 0) // for adaptive time controler
      {
        //restoring results from previous time step
        copyv(np_gv.lhsb, np_gv.lhs, Ndoft);
        copyv(np_gv.tdlhsb, np_gv.tdlhs, Ndoft);
      }
      //  approximation of nodal values into integration points
      approximation ();
      if (Cp->dpt == pass_by_aux_ip)
        aip_approximation(Mm->tnip, MTipmap);

      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();

      if ((Mespr==1) && (Myrank == 0))
        fprintf (stdout,"\n\n MEFEL time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (mdt<Mp->timecon.dtmin){
        if ((Mespr==1) && (Myrank == 0))  fprintf (stderr, " MEFEL time increment is less than the minimum time increment");
        if ((Mespr==1) && (Myrank == 0))  fprintf (stderr, "\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
        if ((Mespr==1) && (Myrank == 0))  fprintf (stderr, "\n FORCED output of results from this step is performed\n");
       
 
        compute_req_valt (lcid);
        print_stept_forced(lcid, ti, Tp->time, np_gv.rhs);
        print_flusht();

        compute_req_val (mlcid);
        print_step_forced(mlcid, mi, Mp->time, mt_gv.fb);
        print_flush();
	
        stop=1;
        break;
      }
      // change status of elements and nodes according to the previous time step
      if (Mp->tprob == growing_mech_structure)
      {
        update_elemnod_prev_time(prev_timem, mt_gv.ncd, mt_gv.nce);
      }
      continue;
    }
    
    if (stop==1){
      //  MEFEL requires smaller time step than the minimum time step
      //  the code will be terminated
      break;
    }
    
    //    
    // Equilibrium was attained both in TRFEL and MEFEL
    //

    // For nonlinear TRFEL problems
    //
    // backup nodal values, lshsb = lhs
    copyv(np_gv.lhs, np_gv.lhsb, Ndoft);
    // backup time derivatives of nodal values, tdlhsb = tdlhs
    copyv(np_gv.tdlhs, np_gv.tdlhsb, Ndoft);
    
    //  printing of TRFEL output and graphical informations
    compute_req_valt (lcid);
    print_stept(lcid,Tp->istep,Tp->time,np_gv.rhs);//print_stept(0,i,Tp->time,NULL);
    print_flusht();
    // perform TRFEL backup
    if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
    {
      if ((Mesprt==1) && (Myrank == 0))
        fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (np_gv.lhs,np_gv.tdlhs,np_gv.f,Tp->istep,Tp->time,tdt,Tp->timecont,Ndoft);
    }
    
    // for TRFEL nonlinear problems, save vectors of unknowns and unknown increments
    if (((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin) || 
         (Tp->tprob == discont_nonlin_nonstat_problem))  && Tp->timecont.tct > 0)
    {
      // storage of results from previous inner step
      copyv(np_gv.lhs, np_gv.lhsb, Ndoft);
      copyv(np_gv.tdlhs, np_gv.tdlhsb, Ndoft);
    }

    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
    {
      // Following two commands must be performed here because in case that
      // elements or nodes are changed in the next step and the step would not converge,
      // the lhs and tdlhs could not be recovered correctly
	
      // backup of nodal values because code numbers may be changed
      //debug??!!Tt->lhs_save (np_gv.lhs,Lsrst->lhsi,np_gv.tdlhs);
      //  backup of nodal forces because code numbers may be changed//from mefel??!!
      //Mt->save_nodforce(mt_gv.f);//from mefel???!!!
	
      // new elements can be definitely changed to old ones
      Gtt->update_auxinf();
      prev_timet = Tp->time;
    }

    // for MEFEL nonlinear problems
    //  Backup of r is used in case that the Newton-Raphson procedure does not converge.
    //  lhsb[j]=r[j]
    copyv(mt_gv.r, mt_gv.lhsb, Ndofm);
    // actual total load vector becomes previous total load vector
    // fp[k]=f[k]
    copyv(mt_gv.f, mt_gv.fp, Ndofm);

    if (Mp->tprob == growing_mech_structure)
    {
      // Following two commands must be performed here because in case that
      // elements or nodes are changed in the next step and the step would not converge,
      // the nodal values and nodal forces could not be recovered correctly
      // nodval due to changes in Mt->nodedispl and nodforce due to changes in vector f

      // backup of nodal values because code numbers may be changed due to
      Mt->save_nodval(lcid);
      //  backup of nodal forces because code numbers may be changed
      Mt->save_nodforce(mt_gv.f);

      // new elements can be definitely changed to old ones
      Gtm->update_auxinf();
      prev_timem = Mp->time;
    }

    //  update of values stored at integration points
    Mm->updateipval();
    
    mefel_trfel(lcid);

    //  printing of MEFEL output and graphical informations
    compute_req_val (lcid);
    print_step(lcid, Mp->istep, Mp->time, mt_gv.fl);
    print_flush();
	
    // perform MEFEL backup
    if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.save_stat()))
    {
      if ((Mespr==1) && (Myrank == 0))
        fprintf (stdout,"\n Creating backup file\n");
      solver_save (mt_gv.r,mt_gv.fp,Mp->istep,Mp->time,dt,&Mp->timecon,Ndofm);
    }

    //  METR cleans MEFEL matrices for more memory in TRFEL part
    if(Cp->cleanmatrix == clean_yes){
      //  cleaning of stiffness matrix
      if (Smat != NULL){
	fprintf (stdout,"\n\n Cleaning of stiffness matrix \n\n ");
	delete Smat;
	Smat=NULL;
      }
    }
    
    if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
      prev_timet = tnewtime;

  }while (Cp->time < end_time); 
  
  print_close ();
  print_closet();
}
