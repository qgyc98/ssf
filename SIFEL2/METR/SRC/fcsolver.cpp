#include "fcsolver.h"
#include "globalc.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "lhsrhs.h"
#include "globalt.h"
#include "globmat.h"
#include "globmatt.h"
#include "globmatc.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "elemswitch.h"
#include "dloadcase.h"
#include "transprint.h"
#include "mechprint.h"
#include "backupsol.h"
#include "backupsolt.h"

#include "intpoints.h"
#include "intpointst.h"


/**
   function solves fully coupled thermo-hydro-mechanical problem
   
   JK
*/
void solve_fcouplprob ()
{
  //------------------------------------
  //myslim, ze toto tady byt nemusi TKr 12/9/2008
  double *rhs;
  
  Lsrs->lhs=Lsrsc->lhs;
  Lsrs->rhs=Lsrsc->rhs;
  Lsrst->lhs=Lsrsc->lhs;
  Lsrst->rhs=Lsrsc->rhs;
  
  rhs = Lsrsc->rhs;
  
  //  vectors of loading
  // zde musi byt dlc a sesteveni pro pocatecni stav
  //Mb->dlc[0].assemble (0,rhs,NULL,Ndofc,Cp->timecon.starttime ());
  //Mb->lc[1].assemble (0,rhs+Ndofc);
  //-------------------------------------


  //  solver
  nonlinear_solver_coupl (0);
  
  Lsrs->lhs=NULL;
  Lsrs->rhs=NULL;
  Lsrst->lhs=NULL;
  Lsrst->rhs=NULL;
}


void nonlinear_solver_coupl (long lcid)
{
  switch (Cp->tnlinsol){
    case newtonc:{
      //newton_raphson_coupl (lcid);
      newton_raphson_coupl_vform (lcid);
      //newton_raphson_coupl_dform (lcid);
      break;
    }
    case newtondform:{
      //newton_raphson_coupl (lcid);
      //newton_raphson_coupl_vform (lcid);
      newton_raphson_coupl_dform (lcid);
      break;
    }
    default:
      print_err("unknown solver type (tnlinsol=%d) of nonlinear equation system is required", __FILE__, __LINE__, __func__, int(Cp->tnlinsol));
  }
}


/**
   function solves fully coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   created 22/7/2020, 
   revised by TKr 12/9/2008, 
   renamed, modified by TKo 4.2019
*/
void newton_raphson_coupl_tkr (long lcid)
{
  long i,j,n,ini,nsts;
  double zero, cnewtime, dt, dtmin, end_time, alpha, err;
  // allocated vectors
  double *d, *p, *fi, *dfi, *fl, *rhst, *fb, *tdlhsi, *auxtdlhs;
  double *tdlhsb, *lhsb;
  double *fp;
  // pointers to global vectors
  double *lhs, *lhsi, *tdlhs, *rhs;
  // vectors of selective norms
  vector vnorfa(1+Tp->ntm), vnorfb(1+Tp->ntm), tol(1+Tp->ntm), aterr(1+Tp->ntm);

  //  coefficient of trapezoidal method
  alpha=Cp->alpha;
  zero=Cp->zero;
  //  maximum number of iterations in inner loop
  ini=Cp->niilnr;
  //  required norm of vector of unbalanced forces and fluxes
  err=Cp->errnr;

  // prepare vector of required norms of unbalanced flux/stress resultants
  tol[0] = Mp->nlman->errnr;
  for(i=0; i<Tp->ntm; i++)
    tol[i+1] = Tp->errarr[i];

  //  total number of degrees of freedom of problem
  //  it is sum of DOFs from MEFEL and TRFEL
  n=Ndofc;
  
  //  pointers assigment to local arrays
  //  arrays are allocated in METR in Lsrsc

  //  vector of nodal values
  lhs = Lsrsc->give_lhs (lcid);
  //  vector of time derivatives of nodal values
  tdlhs = Lsrsc->give_tdlhs (lcid);
  //  vector of initial values
  lhsi = Lsrsc->lhsi;
  //  vector of the right hand side
  rhs = Lsrsc->give_rhs (lcid);
  
  //  initial values of arryas
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);

  //-----------------------------------
  //this below should be changed??!!
  //  pointers assigment to MEFEL part
  //  arrays are allocated in METR
  Lsrs->lhs = lhs;
  Lsrs->rhs = rhs;
  
  //  pointers assigment to TRFEL part
  //  arrays are allocated in METR
  Lsrst->lhs = lhs;
  Lsrst->tdlhs = tdlhs;
  Lsrst->rhs = rhs;
  //-----------------------------------

  //  vector of unknown predictor d_n + dt*(1-alpha)*v_n
  d = new double [n];
  nullv (d,n);
  //  vector of rhs predictor K_n*(u_n+dt(1-alpha)*v_n = K_n*d
  p = new double [n];
  nullv (p,n);
  //  vector of corrector for TRFEL part, i.e. flux increment
  dfi = new double [n];
  nullv (dfi,n);
  //  residuum vector
  fb = new double [n];
  nullv (fb,n);
  //  residuum vector
  fp = new double [n];
  nullv (fp,n);
  //  vector of internal forces for MEFEL part or flux resultant vector
  fi = new double [n];
  nullv (fi,n);
  //  vector of prescribed force loads from MEFEL, it does not contain forces caused by prescribed displacements
  fl = new double [n];
  nullv (fl,n);
  // vector for trfel right hand side
  rhst = new double [n];
  nullv(rhst, n);
  // vector of tdlhs increments
  tdlhsi = new double [n];
  nullv(tdlhsi, n);

  // vector of auxiliary tdlhs 
  auxtdlhs = new double [n];
  nullv(auxtdlhs, n);

  if(Cp->timecon.tct > 0){
    //  backup of the lhs vector
    lhsb   = new double [n];
    nullv (lhsb,n);
    //  backup of the tdlhst vector
    tdlhsb   = new double [n];
    nullv (tdlhsb,n);
  }
  else
  {
    lhsb   = NULL;
    tdlhsb = NULL;
  }
 
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
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);
  //  minimum time increment used in METR inner loop
  dtmin=Cp->timecon.dtmin;
  //  number of successful time steps (number of steps without inner iterations)
  nsts = 0;

  //  initial approximation of quantities from TRFEL nodes to MEFEL integration points
  //trfel_mefel_by_nodes_comp(); // je to potreba k Cmu->initmaterialmodels ?

  //  initial approximation of quantities from TRFEL nodes to MEFEL integration points
  initapproximation();

  //  initiation of coupled material models
  Cmu->initmaterialmodels();  
  //Cml->initmaterialmodels();  

  /**************************/
  /*  main iteration loop  **/
  /**************************/
  i=0;
  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    init_trfel_mefel();
    Mm->initmaterialmodels(lcid, true);
    Mm->aip_initmaterialmodels(lcid, true);
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file for MEFEL\n");
    solver_restore (lhs,rhs,i,Cp->time,dt,&Cp->timecon,n);
    Mm->updateother();
    Mp->timecon.take_values (Cp->timecon);
    print_init(-1, "at");
    print_flush();
  }
  else
  {
    // initiation of mechanical material models
    //  approximation of temperatures and humidities from TRFEL nodes to MEFEL integration points
    init_trfel_mefel ();
    //init_trfel_mefel_by_nodes_comp ();
    Mm->initmaterialmodels(lcid, false);
    Mm->aip_initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    if (Mp->eigstrains > 0)
      internal_forces(lcid, rhs);
    print_init(-1, "wt");
    print_step(lcid, i, Cp->time, rhs);
    print_flush();
  }
  if (Tp->hdbcont.restore_stat()){
    Tm->initmaterialmodels();
    Tm->aip_initmaterialmodels();
    if (Mesprt==1)
      fprintf (stdout,"\n Reading of backup file for TRFEL\n");
    solvert_restore (lhs,tdlhs,rhs,i,Cp->time,dt,Cp->timecon,n);
    Tp->timecont.take_values (Cp->timecon);
    compute_req_valt (0);
    print_initt(-1, "at");
    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    Tm->aip_initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Cp->time,NULL);
    print_flusht();
  }
  
  //  nodes - integration points interpolation
  approximationc ();  
  //  METR is copying data between TRFEL and MEFEL in both directions
  pass_coup_data(lcid);

  do{
    
    i++;
    
    Mp->istep = i;
    Mp->jstep = -1;
    Tp->istep = i;
    Tp->jstep = -1;
    //  new time and new time increment
    cnewtime = Cp->time = Cp->timecon.newtime (dt);
    //dt = Cp->timecon.actualbacktimeincr ();
    Mp->time = Cp->time;
    Mp->timecon.take_values (Cp->timecon);
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);
    
    if (Mesprc==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------------------");
    if (Mesprc==1)  fprintf (stdout,"\n NEW METR TIME STEP = %ld: METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
    if (Mesprc==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------------------");
    
    //  stiffness and conductivity matrices
    zero_order_matrix (lcid);

    //D0mat->printmat(Outc);//debug??!!

    //  capacity matrix
    first_order_matrix (lcid);

    //D1mat->printmat(Outc);//debug??!!
    //  previous solution
    for (j=0;j<n;j++){
      d[j] = lhs[j];
      fp[j] = rhs[j];
    }
    
    //  assembling of vectors of prescribed nodal forces and fluxes
    metr_right_hand_side (lcid,rhs,fl);
    
    //C.1/dt
    D1mat->scalgm (1.0/dt);

    //  matrix
    //  C.1/dt - (1-alpha).K
    D1mat->addgm (alpha-1.0,*D0mat);

    //  auxiliary vector  (C/dt - K(1-alpha)d)
    D1mat->gmxv (d,p);
    
    //  capacity matrix once again
    first_order_matrix (lcid);
    //  matrix of the system of equations
    //C.1/dt
    D1mat->scalgm (1.0/dt);
    //  C.1/dt + alpha.K
    D1mat->addgm (alpha,*D0mat);

    //  right hand side assembling
    for (j=0;j<n;j++){
      fb[j] = (1.0-alpha)*fp[j] + alpha*rhs[j] + p[j];
    }
    
    //  solution of the system of algebraic equations
    Cp->ssle->solve_system (Gtu,D1mat,tdlhs,fb,Outc);
    
    //  new left hand side vector (unknowns) vector
    for (j=0;j<n;j++){
      lhs[j] = tdlhs[j];
    }
    
    //  update of values stored at integration points before results printing
    updateval();
    // actualize passed quantities according to attained values
    pass_coup_data(lcid);
    //  printing of MEFEL output and graphical informations
    compute_req_val (0);
    print_step(0, i, Cp->time, fi);
    print_flush();
    //  printing of TRFEL output and graphical informations
    compute_req_valt (0);
    print_stept(0, i, Cp->time, fi);
    print_flusht();
    
    if ((Mp->timecon.isitimptime()==1) && Mp->hdbcont.save_stat()){
      if (Mespr==1)
	fprintf (stdout,"\n Creating backup file for MEFEL\n");
      solver_save (lhs,rhs,i,Cp->time,dt,&Cp->timecon,n);
    }
    if ((Tp->timecont.isitimptime()==1) && Tp->hdbcont.save_stat()){
      if (Mesprt==1)
	fprintf (stdout,"\n Creating TRFEL backup file\n");
      solvert_save (lhs,tdlhs,rhs,i,Cp->time,dt,Cp->timecon,n);
    }
    
  }while(Cp->time<end_time);
  
  print_close ();
  print_closet ();
  
  delete [] fl;
  delete [] fi;
  delete [] fb;
  delete [] fp;
  delete [] p;
  delete [] rhst;
  delete [] d;
  delete [] dfi;
  delete [] tdlhsi;
  delete [] lhsb;
  delete [] tdlhsb;
}



/**
   function solves fully coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   cerated 24/4/2007, 
   revised by TKr 12/9/2008, 
   renamed, modified by TKo 4.2019
*/
void newton_raphson_coupl_vform (long lcid)
{
  long i,j,k,n,ini,nsts, kk, ii;
  long totiter = 0;
  double zero, cnewtime, dt, dtmin, end_time, alpha, err, dtr;
  // allocated vectors
  double *d, *p, *pm, *fp, *fi, *dfi, *fl, *rhst, *fb, *tdlhsi, *auxtdlhs;
  double *tdlhsb, *lhsb;
  // pointers to global vectors
  double *lhs, *lhsi, *tdlhs, *rhs;
  // vector of DOF block id
  ivector dofbid(Ndofc);
  // vectors of selective norms
  vector vnorfa(1+Tp->ntm), vnorfb(1+Tp->ntm), tol(1+Tp->ntm), aterr(1+Tp->ntm);
  long attol_flag;

  // step id where results will be printed out at each NR iteration
  // set to -1 for no detailed output
  long det_out_i = -1;

  //  coefficient of trapezoidal method
  alpha=Cp->alpha;
  zero=Cp->zero;
  //  maximum number of iterations in inner loop
  ini=Cp->niilnr;
  //  required norm of vector of unbalanced forces and fluxes
  err=Cp->errnr;

  // prepare vector of required norms of unbalanced flux/stress resultants
  assemble_dof_block_id(dofbid);
  tol[0] = Mp->nlman->errnr;
  for(i=0; i<Tp->ntm; i++)
    tol[i+1] = Tp->errarr[i];

  //  total number of degrees of freedom of problem
  //  it is sum of DOFs from MEFEL and TRFEL
  n=Ndofc;
  
  //  pointers assigment to local arrays
  //  arrays are allocated in METR in Lsrsc

  //  vector of nodal values
  lhs = Lsrsc->give_lhs (lcid);
  //  vector of time derivatives of nodal values
  tdlhs = Lsrsc->give_tdlhs (lcid);
  //  vector of initial values
  lhsi = Lsrsc->lhsi;
  //  vector of the right hand side
  rhs = Lsrsc->give_rhs (lcid);
  
  //  initial values of arryas
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);

  //-----------------------------------
  //this below should be changed??!!
  //  pointers assigment to MEFEL part
  //  arrays are allocated in METR
  Lsrs->lhs = lhs;
  Lsrs->rhs = rhs;
  
  //  pointers assigment to TRFEL part
  //  arrays are allocated in METR
  Lsrst->lhs = lhs;
  Lsrst->tdlhs = tdlhs;
  Lsrst->rhs = rhs;
  //-----------------------------------

  //  vector of unknown predictor d_n + dt*(1-alpha)*v_n
  d = new double [n];
  nullv (d,n);
  //  vector of rhs predictor K_n*(u_n+dt(1-alpha)*v_n = K_n*d
  p = new double [n];
  nullv (p,n);
  //  vector of corrector for TRFEL part, i.e. flux increment
  dfi = new double [n];
  nullv (dfi,n);
  //  residuum vector
  fb = new double [n];
  nullv (fb,n);
  //  vector of internal forces for MEFEL part or flux resultant vector
  fi = new double [n];
  nullv (fi,n);
  //  vector of prescribed force loads from MEFEL, it does not contain forces caused by prescribed displacements
  fl = new double [n];
  nullv (fl,n);
  // vector for trfel right hand side
  rhst = new double [n];
  nullv(rhst, n);
  // vector of tdlhs increments
  tdlhsi = new double [n];
  nullv(tdlhsi, n);
  // vector of rhs from previous time step
  fp = new double [n];
  nullv(fp, n);
  // vector of mechanical part for predictor
  pm = new double [n];
  nullv(pm, n);

  // vector of auxiliary tdlhs 
  auxtdlhs = new double [n];
  nullv(auxtdlhs, n);

  if(Cp->timecon.tct > 0){
    //  backup of the lhs vector
    lhsb   = new double [n];
    nullv (lhsb,n);
    //  backup of the tdlhst vector
    tdlhsb   = new double [n];
    nullv (tdlhsb,n);
  }
  else
  {
    lhsb   = NULL;
    tdlhsb = NULL;
  }
 
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
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);
  //  minimum time increment used in METR inner loop
  dtmin=Cp->timecon.dtmin;
  //  number of successful time steps (number of steps without inner iterations)
  nsts = 0;

  //  initial approximation - nodes - integration points interpolation at TRFEL and METR
  approximationc ();
  if (Cp->dpt == pass_by_aux_ip){
    // initial approximation of quantities from TRFEL nodes to TRFEL INT. points
    initapproximation();
    // pass initial values computed in TRFEL INT. points to the MEFEL array aip_nonmechq (i.e. initial values of nonmechq for AUX. points)
    actualize_aip_nonmechq(Tm->tnip, TMipmap);
  }

  //  initiation of coupled material models
  Cmu->initmaterialmodels();  
  //Cml->initmaterialmodels();  

  /**************************/
  /*  main iteration loop  **/
  /**************************/
  i=0;
  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    // Generally, init_trfel_mefel() computes initial TRFEL values at MEFEL int. points and stores them in MEFEL nonmechq array.
    // In the case of aux. point concept:
    // passes initial values from TRFEL AUX. points to nonmechq array in MEFEL (associated with the regular MEFEL INT. points)
    // initial values of TRFEL AUX. points are calculated in the body of init_trfel_mefel() by aip_initapproximation
    init_trfel_mefel();
    Mm->initmaterialmodels(lcid, true);
    Mm->aip_initmaterialmodels(lcid, true);
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file for MEFEL\n");
    solver_restore (lhs,rhs,i,Cp->time,dt,&Cp->timecon,n);
    Mm->updateother();
    Mp->timecon.take_values (Cp->timecon);
    print_init(-1, "at");
    print_flush();
  }
  else
  {
    // initiation of mechanical material models
    // Generally, init_trfel_mefel() computes initial TRFEL values at MEFEL int. points and stores them in MEFEL nonmechq array.
    // In the case of aux. point concept:
    // passes initial values from TRFEL AUX. points to nonmechq array in MEFEL (associated with the regular MEFEL INT. points)
    // initial values of TRFEL AUX. points are calculated in the body of init_trfel_mefel() by aip_initapproximation
    init_trfel_mefel ();
    Mm->initmaterialmodels(lcid, false);
    Mm->aip_initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    if ((Mp->eigstrains > 0) && (Mp->eigstrains < 4))
      internal_forces(lcid, rhs);
    print_init(-1, "wt");
    print_step(lcid, i, Cp->time, rhs);
    print_flush();
  }
  if (Tp->hdbcont.restore_stat()){
    Tm->initmaterialmodels();
    Tm->aip_initmaterialmodels();
    if (Mesprt==1)
      fprintf (stdout,"\n Reading of backup file for TRFEL\n");
    solvert_restore (lhs,tdlhs,rhs,i,Cp->time,dt,Cp->timecon,n);
    Tp->timecont.take_values (Cp->timecon);
    compute_req_valt (0);
    print_initt(-1, "at");
    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    Tm->aip_initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Cp->time,NULL);
    print_flusht();
  }
  
  //  nodes - integration points interpolation
  approximationc ();  
  //  METR is copying data between TRFEL and MEFEL in both directions
  pass_coup_data(lcid);

  do{
    
    i++;
    
    Mp->istep = i;
    Mp->jstep = -1;
    Tp->istep = i;
    Tp->jstep = -1;
    //  new time and new time increment
    cnewtime = Cp->time = Cp->timecon.newtime (dt);
    //dt = Cp->timecon.actualbacktimeincr ();
    Mp->time = Cp->time;
    Mp->timecon.take_values (Cp->timecon);
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);
    
    if (Mesprc==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------------------");
    if (Mesprc==1)  fprintf (stdout,"\n NEW METR TIME STEP = %ld: METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
    if (Mesprc==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------------------");
    
    //  stiffness and conductivity matrices
    zero_order_matrix (lcid);

    //D0mat->printmat(Outc);//debug??!!
    
    //  capacity matrix
    first_order_matrix (lcid);

    //D1mat->printmat(Outc);//debug??!!
    
    //  predictor   d+(1-alpha)*dt*v
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //fprintf(Outc,"\nVector tdlhs at Time = %e\n",Tp->time);//debug??!!
    //for (j=0;j<n;j++){
    //fprintf(Outc,"tdlhs[%ld] = %e\n",j,tdlhs[j]);//debug??!!
    //}
    //fflush(Outt);

    //  assembling of vectors of prescribed nodal forces and fluxes
    metr_right_hand_side (lcid,rhs,fl);
    
    //fflush(Outc);//debug??!!
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    D0mat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    D1mat->addgm (dt*alpha,*D0mat);

    //D1mat->printmat(Outc);//debug??!!

    //  increments of prescribed forces
    //incr_internal_gforces (lcid,fi);

    // compute MEFEL contributions to predictor vector in the form of internal forces
    nullv(fi, n);
    // set predictor lhs for MEFEL
    //Lsrs->lhs = d;
    // compute MEFEL contributions for predictor
    //internal_forces(lcid, fi);
    // restor original lhs
    //Lsrs->lhs = lhs;
    

    //  auxiliary vector  pm = K (1-alpha)*dt*v
    D0mat->gmxv (tdlhs,pm);
    cmulv((1-alpha)*dt,pm,n);

    // set mechanical components of TRFEL-style predictor
    for (k=0; k<Mt->nn; k++){
      for (kk=0; kk<Mt->give_ndofn(k); kk++){
	ii = Mt->give_dof(k, kk);
	if (ii > 0)
	  // p[ii-1] = 0.0;
	  //p[ii-1] = fi[ii-1];
	  p[ii-1] = fp[ii-1] + pm[ii-1];
      }
    }

    //  right hand side assembling
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j]-p[j];
    }
    
    //  solution of the system of algebraic equations
    Cp->ssle->solve_system (Gtu,D1mat,tdlhs,fb,Outc);
    
    //  new left hand side vector (unknowns) vector
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }

    //  inner iteration loop
    for (j=0;j<ini;j++){
      
      Mp->jstep = j;
      Tp->jstep = j;

      solution_correction ();
      //updating of int. points material models and approximations
      approximationc ();
      //copying of data
      pass_coup_data(lcid);
      
      if (Cp->fcsolv == fullnewtonc){
	//  stiffness and conductivity matrices are computed at each inner step
	zero_order_matrix (lcid);
	
	//  capacity matrix
	first_order_matrix (lcid);
	
	//  auxiliary vector  K (d+(1-alpha)*dt*v)
	D0mat->gmxv (d,p);
	
	//  matrix of the system of equations
	//  C + alpha.dt.K
	D1mat->addgm (dt*alpha,*D0mat);
      }
      
      //
      // rebuild right hand side vector rhs separately from METR+TRFEL part and MEFEL part
      // because it was rewritten in the procedure of solution of equation system 
      //
      copyv(fl, rhs, n);                    // store fl MEFEL componets to rhs, TRFEL components are zero in the fl vector
      trfel_right_hand_side(lcid, rhst, n); // rebuild TRFEL right hand side in separate vector because it zeroes sent argument
      //mechanical influence of transport problem
      trfel_right_hand_side2 (lcid,rhst,n);
      addv(rhs, rhst, rhs, n);              // add the TRFEL contributions to the rhs
      right_hand_side (lcid,rhs,n);         // METR part is added (without zeroing)


      //  vector of unbalanced fluxes and forces not tested yet!!
      if (Cp->restype == fluxesc){
	// Solver computes unbalanced fluxes and forces
	//  computation of internal forces in METR, MEFEL, TRFEL
	internal_gforces (lcid,fi);
	// check whether material model requires step length decrement
	dtr = Mm->dstep_red_mat();
	if (dtr < 1.0){
          Mp->jstep = Tp->jstep = j = -1;
          attol_flag = 0;
          break;
        }

	for (k=0;k<n;k++)
	  fb[k] = rhs[k] - fi[k];
      }
      
       // residuum vector computed from system of equation in TRFEEL and METR; from unbalanced forces in MEFEL
      if (Cp->restype == lrhsc){
	/////////////////////////////////////////////////
	// new version:

        // TRFEL part residuum from the system equations
	// Res. = F      - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v     
	//      = rhs    -    p                 - dfi
        // MEFEL part residuum can be calculated from the internal forces only
        // Res  = Fa - Fi = rhs - fi
        
        // compute residuum from system of equations for TRFEL components
        //
        // compute corrector vector, i.e. increment of flux resultants
	D1mat->gmxv (tdlhs,dfi);        
	// computation of internal forces in MEFEL
        nullv(fi, n);
	internal_forces(lcid, fi);
	// check whether material model requires step length decrement
	dtr = Mm->dstep_red_mat();
	if (dtr < 1.0){
          Mp->jstep = Tp->jstep = j = -1;
          attol_flag = 0;
          aterr[0] = 2.0*tol[0];
          break;
        }
        // set mechanical components of TRFEL-style corrector to zero
        // substitute MEFEL components of predictor vector by components of internal forces for attained displacements
	for (k=0; k<Mt->nn; k++){
          for (kk=0; kk<Mt->give_ndofn(k); kk++){
            ii = Mt->give_dof(k, kk);
            if (ii > 0){
              dfi[ii-1] = 0.0;            
              p[ii-1]   = fi[ii-1];  
              rhs[ii-1] = fl[ii-1]; // change rhs component to component due to mechanical force load
            }
          }
        }
        

        // compute resulting residual vector 
	for (k=0;k<n;k++)
	  fb[k] = rhs[k] - p[k] - dfi[k];
          /*
	/////////////////////////////////////////////////
	  //old version
	  //  computation of internal forces
	  internal_gforces (lcid,fii);
	  // check whether material model requires step length decrement
	  dtr = Mm->dstep_red_mat();
	  if (dtr < 1.0)
	  {
	  Mp->jstep = Tp->jstep = j = -1;
	  attol_flag = 0;
	  break;
	  }
	  //correct right hand side
	  metr_right_hand_side(lcid,rhs,fl);
	  // Solver computes residuum from system of equations
	  D1mat->gmxv (tdlhs,fp);//new fp vector
	  // Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v     
	  for (k=0;k<n;k++){
	  fb[k] = rhs[k]  - p[k] - fp[k];
	  }
        */
	///////////////////////////////////////////////// 
      }

      // compute norms of the attained load vector, MEFEL components
      compute_vector_sel_norms(rhs, dofbid, vnorfa);
      if (vnorfa[0] == 0.0){ 
        // in the case of zero load vector fa in mechanical problem,
        // try to compute norm of attained reactions
        compute_reactions(lcid);
        vnorfa[0] += Mt->compute_react_norm();
      }
      // compute norms of the residual vector 
      compute_vector_sel_norms(fb, dofbid, vnorfb);
      attol_flag = check_tolerance(vnorfa, vnorfb, tol, aterr);      
      if (Mesprc==1){      
	fprintf (stdout,"\n Time increment i = %ld;  Inner loop j =  %ld", i, j);
        for (k=0; k<aterr.n; k++){
          fprintf(stdout, "; rnorres_%ld = %e", k+1, aterr(k));
          //          fprintf(stdout, "; norf_%ld=%le; norres_%ld = %e", k+1, vnorfa(k), k+1, vnorfb(k));
        }
        fprintf(stdout, "\n");
      }

      totiter++;
      if (attol_flag) break;

      if (i == det_out_i){
        //  printing of MEFEL output and graphical informations
        compute_req_val (0);
        print_step_forced(0, i, Cp->time + dtmin/ini*j, fb);
        print_flush();
        //  printing of TRFEL output and graphical informations
        compute_req_valt (0);
        print_stept_forced(0, i, Cp->time + dtmin/ini*j, fb);
        print_flusht();
      }
      
      Cp->ssle->solve_system (Gtu,D1mat,tdlhsi,fb,Outc);
      
      for (k=0;k<n;k++){
	tdlhs[k] += tdlhsi[k];
	lhs[k]   += alpha*dt*tdlhsi[k];
      }
    }
    
    if (attol_flag==0){
      // equilibrium state was not attained
      nsts = 0;
      if (Mesprc==1){ 
        fprintf (stdout,"\n\n Iteration process failed at %ld-th iteration:\n", j+1);
        fprintf(stdout,"  MEFEL toler = %e, rerr = %e\n", tol[0], aterr[0]);//new
        for (k=1; k<aterr.n; k++)
          fprintf(stdout,"  TRFEL toler_m%ld = %e, rerr_m%ld = %e\n", k, tol[k], k, aterr[k]);//new
      }
	
      //  backup is used for actual METR values
      if(Cp->timecon.tct > 0){
        for (j=0;j<n;j++){
          lhs[j]=lhsb[j];
          tdlhs[j]=tdlhsb[j];
        }
      }
      approximationc();
	
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      Cp->timecon.oldtime ();
      i--;
	
      if (Mesprc==1)
        fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
        if (Mesprc==1){
          fprintf (stderr,"\n\n time increment is less than minimum time increment");
          fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
        }
        compute_req_val (0);
        print_step_forced(0, i, Cp->time, fb);
        print_flush();
        compute_req_valt (0);
        print_stept_forced(0, i, Cp->time, fb);
        print_flusht();
        break;
      }
      // update strains, stresses and state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      internal_forces(lcid,fi);

      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      // replace by above internal_forces call
      //Mm->updateother();

      pass_coup_data(lcid);
    }
    else{
      // equilibrium state was attained
      if (Mesprc==1)  
        fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,normv(aterr));
      nsts++;
      if(Cp->timecon.tct > 0){        
        subv(lhs, lhsb, auxtdlhs, n);
        cmulv(1.0/dt, auxtdlhs, n);
        compute_vector_sel_norms(auxtdlhs, dofbid, vnorfa);
        compute_vector_sel_norms(tdlhs, dofbid, vnorfb);
        //fprintf(Outc, "istep=%ld, tdlhs[0]=%.10le, auxtdlhs[0]=%.10le, tdlhs[1]=%.10le, auxtdlhs[1]=%.10le", i+1, vnorfb[0], vnorfa[0], vnorfb[1], vnorfa[1]);//debug??!!
        //fflush(Outc);
        //copyv(auxtdlhs, tdlhs, n);
        subv(lhs, lhsb, auxtdlhs, n);
        compute_vector_sel_norms(auxtdlhs, dofbid, vnorfa);
        //fprintf(Outc, ", dlhs[0]=%.10le, dlhs[1]=%.10le\n", vnorfa[0], vnorfa[1]);//debug??!!
        
        // backup of attained vector of unknowns
        copyv(lhs, lhsb, n);
        //backup of attained rate of vector of unknowns
        copyv(tdlhs, tdlhsb, n);
      }
      // backup right hand side of the mechanical part and store it as the attained load vector in the last equilibrium step
      nullv(fp, Ndofc); // this must be done here, because mefel_right_hand_side zeroes only Ndofm components of the vectors
      mefel_right_hand_side(lcid, fp, NULL);
	
      if (nsts==2){
        if (dt < Cp->timecon.dtmax){
          dt*=2.0;
          if (dt > Cp->timecon.dtmax)
            dt = Cp->timecon.dtmax;
          nsts=0;
          if (Mespr==1)  
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 2 steps");
        }
      }
      
      //  update of values stored at integration points before results printing
      updateval();
      // actualize passed quantities according to attained values
      pass_coup_data(lcid);
      //  printing of MEFEL output and graphical informations
      compute_req_val (0);
      print_step(0, i, Cp->time, fi);
      print_flush();
      //  printing of TRFEL output and graphical informations
      compute_req_valt (0);
      print_stept(0, i, Cp->time, fi);
      print_flusht();
      
      if ((Mp->timecon.isitimptime()==1) && Mp->hdbcont.save_stat()){
        if (Mespr==1)
          fprintf (stdout,"\n Creating backup file for MEFEL\n");
        solver_save (lhs,rhs,i,Cp->time,dt,&Cp->timecon,n);
      }
      if ((Tp->timecont.isitimptime()==1) && Tp->hdbcont.save_stat()){
        if (Mesprt==1)
          fprintf (stdout,"\n Creating TRFEL backup file\n");
        solvert_save (lhs,tdlhs,rhs,i,Cp->time,dt,Cp->timecon,n);
      }
    }
  }while(Cp->time<end_time);

  fprintf(stdout, "\nTotal number of performed iterations: %ld\n", totiter);
  print_close ();
  print_closet ();
  
  delete [] fl;
  delete [] fi;
  delete [] fb;
  delete [] fp;
  delete [] p;
  delete [] pm;
  delete [] rhst;
  delete [] d;
  delete [] dfi;
  delete [] tdlhsi;
  delete [] lhsb;
  delete [] tdlhsb;
}


/**
   function solves fully coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method, incremental formulation (dform) is being used 
   rather than velocity formulation (vform)
   
   @param lcid - load case id
   
   24/4/2007, revised by TKr 12/9/2008
*/
void newton_raphson_coupl_dform (long lcid)
{
  long i,j,k,n,ini,nsts, kk, ii;
  double zero, cnewtime, dt, dtmin, end_time, alpha, err, dtr;
  // allocated vectors
  double *d, *p, *fi, *fl, *rhst, *fb, *dlhsi, *auxtdlhs;
  double *tdlhsb, *lhsb;
  // pointers to global vectors
  double *lhs, *lhsi, *tdlhs, *rhs;
  // vector of DOF block id
  ivector dofbid(Ndofc);
  // vectors of selective norms
  vector vnorfa(1+Tp->ntm), vnorfb(1+Tp->ntm), tol(1+Tp->ntm), aterr(1+Tp->ntm);
  long attol_flag;

  //  coefficient of trapezoidal method
  alpha=Cp->alpha;
  zero=Cp->zero;
  //  maximum number of iterations in inner loop
  ini=Cp->niilnr;
  //  required norm of vector of unbalanced forces and fluxes
  err=Cp->errnr;

  // prepare vector of required norms of unbalanced flux/stress resultants
  assemble_dof_block_id(dofbid);
  tol[0] = Mp->nlman->errnr;
  for(i=0; i<Tp->ntm; i++)
    tol[i+1] = Tp->errarr[i];

  //  total number of degrees of freedom of problem
  //  it is sum of DOFs from MEFEL and TRFEL
  n=Ndofc;
  
  //  pointers assigment to local arrays
  //  arrays are allocated in METR in Lsrsc

  //  vector of nodal values
  lhs = Lsrsc->give_lhs (lcid);
  //  vector of time derivatives of nodal values
  tdlhs = Lsrsc->give_tdlhs (lcid);
  //  vector of initial values
  lhsi = Lsrsc->lhsi;
  //  vector of the right hand side
  rhs = Lsrsc->give_rhs (lcid);
  
  //  initial values of arryas
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);

  //-----------------------------------
  //this below should be changed??!!
  //  pointers assigment to MEFEL part
  //  arrays are allocated in METR
  Lsrs->lhs = lhs;
  Lsrs->rhs = rhs;
  
  //  pointers assigment to TRFEL part
  //  arrays are allocated in METR
  Lsrst->lhs = lhs;
  Lsrst->tdlhs = tdlhs;
  Lsrst->rhs = rhs;
  //-----------------------------------

  //  vector of unknown predictor d_n + dt*(1-alpha)*v_n
  d = new double [n];
  nullv (d,n);
  //  vector of rhs predictor K_n*(u_n+dt(1-alpha)*v_n = K_n*d
  p = new double [n];
  nullv (p,n);
  //  residuum vector
  fb = new double [n];
  nullv (fb,n);
  //  vector of internal forces for MEFEL part or flux resultant vector
  fi = new double [n];
  nullv (fi,n);
  //  vector of prescribed force loads from MEFEL, it does not contain forces caused by prescribed displacements
  fl = new double [n];
  nullv (fl,n);
  // vector for trfel right hand side
  rhst = new double [n];
  nullv(rhst, n);
  // vector of lhs increments
  dlhsi = new double [n];
  nullv(dlhsi, n);

  // vector of auxiliary tdlhs 
  auxtdlhs = new double [n];
  nullv(auxtdlhs, n);

  if(Cp->timecon.tct > 0){
    //  backup of the lhs vector
    lhsb   = new double [n];
    nullv (lhsb,n);
    //  backup of the tdlhst vector
    tdlhsb   = new double [n];
    nullv (tdlhsb,n);
  }
  else
  {
    lhsb   = NULL;
    tdlhsb = NULL;
  }
 
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
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);
  //  minimum time increment used in METR inner loop
  dtmin=Cp->timecon.dtmin;
  //  number of successful time steps (number of steps without inner iterations)
  nsts = 0;

  //  initial approximation of quantities from TRFEL nodes to MEFEL integration points
  //trfel_mefel_by_nodes_comp(); // je to potreba k Cmu->initmaterialmodels ?

  //  initiation of coupled material models
  Cmu->initmaterialmodels();  
  //Cml->initmaterialmodels();  

  /**************************/
  /*  main iteration loop  **/
  /**************************/
  i=0;
  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    init_trfel_mefel();
    Mm->initmaterialmodels(lcid, true);
    Mm->aip_initmaterialmodels(lcid, true);
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file for MEFEL\n");
    solver_restore (lhs,rhs,i,Cp->time,dt,&Cp->timecon,n);
    Mm->updateother();
    Mp->timecon.take_values (Cp->timecon);
    print_init(-1, "at");
    print_flush();
  }
  else
  {
    // initiation of mechanical material models
    //  approximation of temperatures and humidities from TRFEL nodes to MEFEL integration points
    //init_trfel_mefel ();
    init_trfel_mefel_by_nodes_comp ();
    Mm->initmaterialmodels(lcid, false);
    Mm->aip_initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    if (Mp->eigstrains > 0)
      internal_forces(lcid, rhs);
    print_init(-1, "wt");
    print_step(lcid, i, Cp->time, rhs);
    print_flush();
  }
  if (Tp->hdbcont.restore_stat()){
    Tm->initmaterialmodels();
    Tm->aip_initmaterialmodels();
    if (Mesprt==1)
      fprintf (stdout,"\n Reading of backup file for TRFEL\n");
    solvert_restore (lhs,tdlhs,rhs,i,Cp->time,dt,Cp->timecon,n);
    Tp->timecont.take_values (Cp->timecon);
    compute_req_valt (0);
    print_initt(-1, "at");
    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    Tm->aip_initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt");
    print_stept(0,i,Cp->time,NULL);
    print_flusht();
  }
  
  //  nodes - integration points interpolation
  approximationc ();  
  //  METR is copying data between TRFEL and MEFEL in both directions
  pass_coup_data(lcid);

  do{
    
    i++;
    
    Mp->istep = i;
    Mp->jstep = -1;
    Tp->istep = i;
    Tp->jstep = -1;
    //  new time and new time increment
    cnewtime = Cp->time = Cp->timecon.newtime (dt);
    //dt = Cp->timecon.actualbacktimeincr ();
    Mp->time = Cp->time;
    Mp->timecon.take_values (Cp->timecon);
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);
    
    if (Mesprc==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------------------");
    if (Mesprc==1)  fprintf (stdout,"\n NEW METR TIME STEP = %ld: METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
    if (Mesprc==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------------------");
    
    //  predictor   d+(1-alpha)*dt*v
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  conductivity matrix
    zero_order_matrix_fc ();
    
    //  capacity and stiffness matrices
    first_order_matrix_fc ();

    //  assembling of vectors of prescribed nodal forces and fluxes
    metr_right_hand_side (lcid,rhs,fl);

    if (Cp->time > 276000.0)
      Mp->time = Cp->time;
      
    
    //fflush(Outc); //debug??!!
    
    //  auxiliary vector  C*(r+(1-alpha)*dt*v) = C*d = p
    D1mat->gmxv (d,p);
          // replace by above internal_forces call

    //  matrix of the system of equations
    //  C + alpha*dt*K
    D0mat->scalgm (dt*alpha);
    D0mat->addgm (1.0,*D1mat);
    
    //  increments of prescribed forces
    //incr_internal_gforces (lcid,fi);

    // compute MEFEL contributions to predictor vector in the form of internal forces
    nullv(fi, n);
    // set predictor lhs for MEFEL
    Lsrs->lhs = d;
    // compute MEFEL contributions for predictor
    internal_forces(lcid, fi);
    // restor original lhs
    Lsrs->lhs = lhs;
    
    //  right hand side assembling of the TRFEL predictor
    //  alpha*dt*F + C*(r+(1-alpha)*dt*v) = alpha*dt*F + C*d = alpha*dt*F + p
    addmultv(p, rhs, alpha*dt, fb, n);

    // set mechanical components of TRFEL predictor
    /*    for (k=0; k<Mt->nn; k++){
      for (kk=0; kk<Mt->give_ndofn(k); kk++){
        ii = Mt->give_dof(k, kk);
        if (ii > 0)
          // p[ii-1] = 0.0;
          fb[ii-1] = (rhs[ii-1]-fi[ii-1]);
      }
      }*/
    
    //  solution of the system of algebraic equations
    Cp->ssle->solve_system (Gtu,D0mat,lhs,fb,Outc);
    
    //  computation of time derivatives of the nodal values
    //  tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    subv(lhs, d, tdlhs, n);
    cmulv(1.0/dt/alpha, tdlhs, n);

    //  inner iteration loop
    for (j=0;j<ini;j++){
      
      Mp->jstep = j;
      Tp->jstep = j;

      solution_correction ();
      //updating of int. points material models and approximations
      approximationc ();
      //copying of data
      pass_coup_data(lcid);
      
      if ((Cp->fcsolv == fullnewtonc) || (Cp->restype == lrhsc)){
	//  stiffness and conductivity matrices are computed at each inner step
	zero_order_matrix_fc ();
	
	//  capacity matrix
	first_order_matrix_fc ();

        /*	
	//  auxiliary vector  C*(r+(1-alpha)*dt*v) = C*d = p
	D1mat->gmxv (d,p);
	
	//  matrix of the system of equations
	//  C + alpha.dt.K
        D0mat->scalgm (dt*alpha);
        D0mat->addgm (1.0,*D1mat);
        */
      }
      
      //
      // rebuild right hand side vector rhs separately from METR+TRFEL part and MEFEL part
      // because it was rewritten in the procedure of solution of equation system 
      //
      copyv(fl, rhs, n);                    // store fl MEFEL componets to rhs, TRFEL components are zero in the fl vector
      trfel_right_hand_side(lcid, rhst, n); // rebuild TRFEL right hand side in separate vector because it zeroes its argument
      addv(rhs, rhst, rhs, n);              // add the TRFEL contributions to the rhs
      right_hand_side (lcid,rhs,n);         // METR part is added (without zeroing)
 


      //  vector of unbalanced fluxes and forces not tested yet!!
      if (Cp->restype == fluxesc){
	// Solver computes unbalanced fluxes and forces
	//  computation of internal forces in METR, MEFEL, TRFEL
	internal_gforces (lcid,fi);
	// check whether material model requires step length decrement
	dtr = Mm->dstep_red_mat();
	if (dtr < 1.0){
          Mp->jstep = Tp->jstep = j = -1;
          attol_flag = 0;
          break;
        }

	for (k=0;k<n;k++)
	  fb[k] = rhs[k] - fi[k];
      }
      
       // residuum vector computed from system of equation in TRFEEL and METR; from unbalanced forces in MEFEL
      if (Cp->restype == lrhsc){
	/////////////////////////////////////////////////
	// new version:

        // TRFEL part residuum from the system equations
	// Res. = F      - C*v_{n+1} + K*r_{n+1}
        // MEFEL part residuum can be calculated from the internal forces only
        // Res  = Fa - Fi = rhs - fi
        
        // compute residuum from system of equations for TRFEL components
        //
	D0mat->gmxv(lhs, fi);
        D1mat->gmxv(tdlhs, p);
        addv(fi, p, fi, n);
        // compute resulting residual vector fro TRFEL part
	for (k=0;k<n;k++)
	  fb[k] = rhs[k] - fi[k];
	// computation of internal forces in MEFEL
        nullv(fi, n);
	internal_forces(lcid, fi);
	// check whether material model requires step length decrement
	dtr = Mm->dstep_red_mat();
	if (dtr < 1.0){
          Mp->jstep = Tp->jstep = j = -1;
          attol_flag = 0;
          break;
        }
        // set mechanical components of TRFEL-style corrector to zero
        // substitute MEFEL components of predictor vector by components of internal forces for attained displacements
	for (k=0; k<Mt->nn; k++){
          for (kk=0; kk<Mt->give_ndofn(k); kk++){
            ii = Mt->give_dof(k, kk);
            if (ii > 0){
              fb[ii-1] = fl[ii-1] - fi[ii-1];
            }
          }
        }        
      }

      // compute norms of the attained load vector, MEFEL components
      compute_vector_sel_norms(rhs, dofbid, vnorfa);
      if (vnorfa[0] == 0.0){ 
        // in the case of zero load vector fa in mechanical problem,
        // try to compute norm of attained reactions
        compute_reactions(lcid);
        vnorfa[0] += Mt->compute_react_norm();
      }
      // compute norms of the residual vector 
      compute_vector_sel_norms(fb, dofbid, vnorfb);
      attol_flag = check_tolerance(vnorfa, vnorfb, tol, aterr);
      
      if (Mesprc==1){      
	fprintf (stdout,"\n Time increment i = %ld;  Inner loop j =  %ld", i, j);
        for (k=0; k<aterr.n; k++){
          fprintf(stdout, "; rnorres_%ld = %e", k+1, aterr(k));
          fprintf(stdout, "; norf_%ld=%le; norres_%ld = %e", k+1, vnorfa(k), k+1, vnorfb(k));
        }
        fprintf(stdout, "\n");
      }

      if (attol_flag) break;

      // solve new increments of lhs in the form
      // (C_{n+1,j} + alpha*dt*K_{n+1,j})*dlhs = alpha*dt*r_{n+1,j}
      // r_{n+1,j} is the residual vector
      D0mat->scalgm (dt*alpha);
      D0mat->addgm (1.0,*D1mat);

      Cp->ssle->solve_system (Gtu,D0mat, dlhsi,fb,Outc);
      
      // add new increments to lhs, calculate new time derivative of lhs vector
      for (k=0;k<n;k++){
	lhs[k] += dlhsi[k];
        tdlhs[k] = (lhs[k]-d[k])/dt/alpha;
      }
    }
    
    if (attol_flag==0){
      // equilibrium state was not attained
      nsts = 0;
      if (Mesprc==1) 
        fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,normv(aterr));//new
	
      //  backup is used for actual METR values
      if(Cp->timecon.tct > 0){
        for (j=0;j<n;j++){
          lhs[j]=lhsb[j];
          tdlhs[j]=tdlhsb[j];
        }
      }
      approximationc();
      pass_coup_data(lcid);
	
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      Cp->timecon.oldtime ();
      i--;
	
      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();

      if (Mesprc==1)
        fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
        if (Mesprc==1){
          fprintf (stderr,"\n\n time increment is less than minimum time increment");
          fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
        }
        compute_req_val (0);
        print_step_forced(0, i, Cp->time, fb);
        print_flush();
        compute_req_valt (0);
        print_stept_forced(0, i, Cp->time, fb);
        print_flusht();
        break;
      }
    }
    else{
      // equilibrium state was attained
      if (Mesprc==1)  
        fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,normv(aterr));
      nsts++;
      if(Cp->timecon.tct > 0){        
        subv(lhs, lhsb, auxtdlhs, n);
        cmulv(1.0/dt, auxtdlhs, n);
        compute_vector_sel_norms(auxtdlhs, dofbid, vnorfa);
        compute_vector_sel_norms(tdlhs, dofbid, vnorfb);
        //fprintf(Outc, "istep=%ld, tdlhs[0]=%.10le, auxtdlhs[0]=%.10le, tdlhs[1]=%.10le, auxtdlhs[1]=%.10le", i+1, vnorfb[0], vnorfa[0], vnorfb[1], vnorfa[1]);//debug??!!
        //fflush(Outc);//debug??!!
        //copyv(auxtdlhs, tdlhs, n);
        subv(lhs, lhsb, auxtdlhs, n);
        compute_vector_sel_norms(auxtdlhs, dofbid, vnorfa);
        //fprintf(Outc, ", dlhs[0]=%.10le, dlhs[1]=%.10le\n", vnorfa[0], vnorfa[1]);//debug??!!
        
        // backup of attained vector of unknowns
        copyv(lhs, lhsb, n);
        //backup of attained rate of vector of unknowns
        copyv(tdlhs, tdlhsb, n);
        //  actual total load vector becomes previous total load vector
        //copyv(rhs, fp, n);//toto zkotrolovat??!!
      }
	
      if (nsts==2){
        if (dt < Cp->timecon.dtmax){
          dt*=2.0;
          if (dt > Cp->timecon.dtmax)
            dt = Cp->timecon.dtmax;
          nsts=0;
          if (Mespr==1)  
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 2 steps");
        }
      }
      
      //  update of values stored at integration points before results printing
      updateval();
      // actualize passed quantities according to attained values
      pass_coup_data(lcid);
      //  printing of MEFEL output and graphical informations
      compute_req_val (0);
      print_step(0, i, Cp->time, fi);
      print_flush();
      //  printing of TRFEL output and graphical informations
      compute_req_valt (0);
      print_stept(0, i, Cp->time, fi);
      print_flusht();
      
      if ((Mp->timecon.isitimptime()==1) && Mp->hdbcont.save_stat()){
        if (Mespr==1)
          fprintf (stdout,"\n Creating backup file for MEFEL\n");
        solver_save (lhs,rhs,i,Cp->time,dt,&Cp->timecon,n);
      }
      if ((Tp->timecont.isitimptime()==1) && Tp->hdbcont.save_stat()){
        if (Mesprt==1)
          fprintf (stdout,"\n Creating TRFEL backup file\n");
        solvert_save (lhs,tdlhs,rhs,i,Cp->time,dt,Cp->timecon,n);
      }
    }
  }while(Cp->time<end_time);
  
  print_close ();
  print_closet ();
  
  delete [] fl;
  delete [] fi;
  delete [] fb;
  delete [] p;
  delete [] rhst;
  delete [] d;
  delete [] dlhsi;
  delete [] lhsb;
  delete [] tdlhsb;
}


/**
   function solves system of nonlinear algebraic
   equations by Newton-Raphson method
   
   @param lcid - load case id
   
   5.4.2003
*/
void newton_raphson_coupl (long lcid)
{
  long i,j,k,n,ini;
  double zero,dt,end_time,alpha,norfb,err;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;
  
  //  coefficient of trapezoidal method
  alpha=Cp->alpha;
  zero=Cp->zero;
  //  maximum number of iterations in inner loop
  ini=Cp->niilnr;
  //  required norm of vector of unbalanced forces and fluxes
  err=Cp->errnr;

  //  total number of degrees of freedom of problem
  //  it is sum of DOFs from MEFEL and TRFEL
  n=Ndofc;
  
  //  pointers assigment to local arrays
  //  arrays are allocated in METR in Lsrsc

  //  vector of nodal values
  lhs = Lsrsc->give_lhs (lcid);
  //  vector of time derivatives of nodal values
  tdlhs = Lsrsc->give_tdlhs (lcid);
  //  vector of initial values
  lhsi = Lsrsc->lhsi;
  //  vector of the right hand side
  rhs = Lsrsc->give_rhs (lcid);
  
  //  initial values of arryas
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);

  //-----------------------------------
  //this below should be changed??!!
  //  pointers assigment to MEFEL part
  //  arrays are allocated in METR
  Lsrs->lhs = lhs;
  Lsrs->rhs = rhs;
  
  //  pointers assigment to TRFEL part
  //  arrays are allocated in METR
  Lsrst->lhs = lhs;
  Lsrst->tdlhs = tdlhs;
  Lsrst->rhs = rhs;
  //-----------------------------------

  //  array containing predictors
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //....
  arhs = new double [n];
  nullv (arhs,n);
  //....
  fb = new double [n];
  nullv (fb,n);
  //....
  fi = new double [n];
  nullv (fi,n);
  
  //  starting time
  Cp->time=Cp->timecon.starttime ();
  //  time increment
  dt=Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();
  
  //  nodes - integration points interpolation
  approximationc ();  
  
  //  initiation of coupled material models
  //Cm->initmaterialmodels();  
  //  initiation of mechanical material models
  //Mm->initmaterialmodels();
  //  initiation of transport material models
  //Tm->initmaterialmodels();

  //computation of required values in MEFEL
  compute_req_valt (0);

  /**************************/
  /*  main iteration loop  **/
  /**************************/
  i=0;
  print_init(-1, "wt");
  print_step(0,i,Cp->time,NULL);
  print_initt(-1, "wt");
  print_stept(0,i,Cp->time,NULL);
  do{
    
    if (Mesprc==1)  fprintf (stdout,"\n time %e",Cp->time);
    
    //  stiffness and conductivity matrices
    zero_order_matrix (lcid);

    //  capacity matrix
    first_order_matrix (lcid);
    


    //fprintf (Outc,"\n\n\n MATICE VODIVOSTI (TUHOSTI) \n\n\n");
    //D0mat->printmat (Outc);

    //fprintf (Outc,"\n\n\n MATICE KAPACITY \n\n\n");
    //D1mat->printmat (Outc);

    //  predictor   d+(1-alpha)*dt*v
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    metr_right_hand_side (lcid,rhs,p);

    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    D0mat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    D0mat->scalgm (dt*alpha);
    D0mat->addgm (1.0,*D1mat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    //fprintf (Outc,"\n\n\n\n\n kontrola prave strany");
    //for (ijk=0;ijk<n;ijk++){
    //fprintf (Outc,"\n %ld %e %e %e",ijk,fb[ijk],rhs[ijk],p[ijk]);
    //}

    //  solution of the system of algebraic equations
    //D0mat->solve_system (Gtt,tdlhs,fb);
    Cp->ssle->solve_system (Gtt,D0mat,tdlhs,fb,Outt);

    //fprintf (Outc,"\n\n\n\n\n kontrola reseni");
    //for (ijk=0;ijk<n;ijk++){
    //fprintf (Outc,"\n %ld %e",ijk,tdlhs[ijk]);
    //}

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    //fprintf (Outc,"\n\n\n\n\n kontrola obnovenych promennych");
    //for (ijk=0;ijk<n;ijk++){
    //fprintf (Outc,"\n %ld %e",ijk,lhs[ijk]);
    //}
    
    solution_correction ();
    approximationc ();
    
    //fprintf (stdout,"\n iterace cislo %ld",i);

    internal_gforces (lcid,fi);

    //fprintf (Outc,"\n\n kontrola vnitrnich toku\n");
    //for (j=0;j<n;j++){
    //fprintf (Outc,"\n %ld  %30.20f",j,fi[j]);
    //}
    
    //  vector of unbalanced fluxes
    if (Cp->restype == fluxesc){
      for (k=0;k<n;k++){
	fb[k]=fi[k];
      }
    }
    if (Cp->restype == lrhsc){
      D0mat->gmxv (tdlhs,fi);
      for (k=0;k<n;k++){
	fb[k] = rhs[k] - p[k] - fi[k];
      }
    }

    norfb = ss (fb,fb,n);
    
    if (Mesprc==1)  fprintf (stdout,"\n%e %e",norfb,err);
    
    if (norfb<err){
      Cp->time = Cp->timecon.newtime ();
      dt = Cp->timecon.actualbacktimeincr ();
      i++;
      
      Mp->time = Cp->time;
      Mp->timecon.take_values (Cp->timecon);
      //Mp->timecon.seconds_days ();
      
      Tp->time = Cp->time;
      Tp->timecont.take_values (Cp->timecon);
      
      print_step(0,i,Cp->time,NULL);//zatim??!!
      print_flush();
      print_stept(0,i,Cp->time,NULL);//zatim??!!
      print_flusht();
      continue;
    }
    
    //  inner iteration loop
    for (j=0;j<ini;j++){
      
      if (Cp->fcsolv == fullnewtonc){
	//  stiffness and conductivity matrices
	zero_order_matrix (lcid);
	
	//  capacity matrix
	first_order_matrix (lcid);
	
	//  auxiliary vector  K (d+(1-alpha)*dt*v)
	D0mat->gmxv (d,p);
	
	//  matrix of the system of equations
	//  C + alpha.dt.K
	D0mat->scalgm (dt*alpha);
	D0mat->addgm (1.0,*D1mat);
      }
      
      //D0mat->solve_system (Gtt,p,fb);
      Cp->ssle->solve_system (Gtt,D0mat,p,fb,Outt);

      for (k=0;k<n;k++){
        tdlhs[k]+=p[k];
        lhs[k]+=alpha*dt*p[k];
      }
      
      solution_correction ();
      approximationc ();
      internal_gforces (lcid,fi);
      
      //fprintf (Outc,"\n\n kontrola vnitrnich toku ve vnitrim cyklu\n");
      //for (jj=0;jj<n;jj++){
      //fprintf (Outc,"\n %ld  %30.20f",jj,fi[jj]);
      //}
      
      //  vector of unbalanced fluxes
      if (Cp->restype == fluxesc){
	for (k=0;k<n;k++){
	  fb[k]=fi[k];
	}
      }
      if (Cp->restype == lrhsc){
	D0mat->gmxv (tdlhs,fi);
	for (k=0;k<n;k++){
	  fb[k] = rhs[k] - p[k] - fi[k];
	}
      }
      
      norfb = ss (fb,fb,n);
      
      if (Mesprc==1)  fprintf (stdout,"\n iterace %ld   chyba %e",j,norfb);
      
      if (norfb<err){
        break;
      }
    }
    


    //  new time and new time increment
    Cp->time = Cp->timecon.newtime ();
    dt = Cp->timecon.actualbacktimeincr ();
    i++;
    
    Mp->time = Cp->time;
    Mp->timecon.take_values (Cp->timecon);
    //Mp->timecon.seconds_days ();
    
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);
    
    print_step(0,i,Cp->time,NULL);//zatim??!!
    print_flush();
    print_stept(0,i,Cp->time,NULL);//zatim??!!
    print_flusht();
    
    //fprintf (Outc,"\n\n\n posun konce tyce d %f",lhs[Gtm->gnodes[4].cn[0]]);//debug??!!

  }while(Cp->time<end_time);
  
  
  delete [] fi;
  delete [] fb;
  delete [] p;
  delete [] arhs;
  delete [] d;
}




/**
  The function assembles vector of dof block identfiers. These identifiers 
  denote which problem the given dof is connected with.
  dofbid[i] = 0 means that the i-th dof represents mechanical unknown
  dofbid[i] > 0 means that the i-th dof represents unknown

  @param dofbid[in,out] - vector of dof block identfiers, dimension must be Ndofc

  Created by Tomas Koudelka, 02.2022
*/
void assemble_dof_block_id(ivector &dofbid)
{
  long i, j, ii;

  fillv(-1, dofbid);
  
  for(i=0; i<Mt->nn; i++){
    for(j=0; j<Mt->give_ndofn(i); j++){
      ii = Mt->give_dof(i, j);
      if (ii > 0)
        dofbid(ii-1) = 0;
    }
  }
    
  for(i=0; i<Tt->nn; i++){
    for(j=0; j<Tt->give_ndofn(i); j++){
      ii = Tt->give_dof(i, j);
      if (ii > 0)
        dofbid(ii-1) = j+1;
    }
  }

  for (i=0; i<Ndofc; i++){
    if (dofbid(i) < 0)
      print_err("block identfier of DOF number %ld cannot be recognized.", __FILE__, __LINE__, __func__, i+1);
  }    
}



/**
  The function computes selective vector norms for particular blocks of global vector c.

  @param c [in]      - array of the global %vector components, dimension is Ndofc
  @param dofbid[in]  - vector/array of block identfiers of particular DOFs, its dimension is Ndofc, 
                       dofbid[i] = 0 means that i*th dof represents a mechanical unknown
                       dofbid[i] > 0 means that i-th dof represents an unknown of the dofbid[i]-th transported media
  @param norms [out] - %vector of the resulting norms, its dimension is Tp->ntm + 1
  
  Created by Tomas Koudelka, 24.5.2018
*/
void compute_vector_sel_norms(const double *c, const ivector &dofbid, vector &norms)
{
  long i;

  if (norms.n == 1+Tp->ntm)
  {
    nullv(norms);  

    for (i=0; i<Ndofc; i++){
      norms(dofbid(i)) += c[i]*c[i];
    }
    for(i=0; i<1+Tp->ntm; i++)
      norms(i) = sqrt(norms(i));
  }
  else
  {
    print_err("wrong number of components (norms.n=%ld) of the vector for norm storage\n"
              " norms.n must be 1+Tp->ntm=%ld", __FILE__, __LINE__, __func__, norms.n, 1+Tp->ntm);
    abort();
  }
}



/**
  The function checks required tolerances of selective norms of the residual %vector normalized 
  to selective righte hand side %vector norms.
  
  @param vnorfa[in] - %vector of selective norms of attained right hand side %vector
  @param vnorfb[in] - %vector of selective norms of the residual %vector
  @param tol[in]    - %vector of required tolerances for particular selective norms of the 
                      residual %vector relative to vnorfa components
  @param aterr[out] - %vector of selective norms of the residual %vector relative to vnorfa components

  @retval 0 - required tolerance was NOT attained for all selective norms
  @retval 1 - required tolerance was attained for all selective norms

  Created by Tomas Koudelka, 24.5.2018
*/
long check_tolerance(const vector &vnorfa, const vector &vnorfb, const vector &tol, vector &aterr)
{
  long i, ret = 0;

  for (i=0; i<vnorfa.n; i++)
  {
    if (vnorfa[i] > 0.0)
    {
      aterr[i] = vnorfb[i]/vnorfa[i];
      if (aterr[i] < tol[i])
        ret++;
    }
    else
    {
      aterr[i] = vnorfb[i];
      if (aterr[i] < tol[i])
        ret++;
    }
  }
  if (ret == tol.n) // all norms of the residual vector are in tolerance
    return 1;

  return 0; // required tolerance was not attained for some components
}



void vector_assemb (double *c,double *m,double *t)
{
  long i,j;
  
  for (i=0;i<Gtm->ndof;i++){
    j=Gtm->cngtopcorr[i];
    c[j]=m[i];
  }
  for (i=0;i<Gtt->ndof;i++){
    j=Gtt->cngtopcorr[i];
    c[j]=t[i];
  }

}

void vector_decomp (double *c,double *m,double *t)
{
  long i,j;
  
  for (i=0;i<Gtm->ndof;i++){
    j=Gtm->cngtopcorr[i];
    m[i]=c[j];
  }
  for (i=0;i<Gtt->ndof;i++){
    j=Gtt->cngtopcorr[i];
    t[i]=c[j];
  }

}



void print_log()
{
  long i, j, k;

  fprintf(Out, "\nMech. int. points:\n");
  for(i=0; i<Mm->tnip; i++)
  {
    fprintf(Out, "%3ld: strain: ", i);
    for(j=0; j<Mm->ip[i].ncompstr; j++)
      fprintf(Out, "%le ", Mm->ip[i].strain[j]);
    fprintf(Out, "\n%3s: stress: ", "");
    for(j=0; j<Mm->ip[i].ncompstr; j++)
      fprintf(Out, "%le ", Mm->ip[i].stress[j]);
    fprintf(Out, "\n%3s: other: ", "");
    for(j=0; j<Mm->ip[i].ncompother; j++)
      fprintf(Out, "%le ", Mm->ip[i].other[j]);
    fprintf(Out, "\n%3s: eqother: ", "");
    for(j=0; j<Mm->ip[i].ncompother; j++)
      fprintf(Out, "%le ", Mm->ip[i].eqother[j]);
    fprintf(Out, "\n");
  }

  fprintf(Out, "\nNonmechq array:\n");
  for(i=0; i<Mm->tnip; i++)
  {
    for(j=0; j<Mm->nnmq; j++)
      fprintf(Out, "%le ", Mm->nonmechq[j][i]);
    fprintf(Out, "\n");
  }

  fprintf(Out, "\nTrans. int. points:\n");
  for(i=0; i<Tm->tnip; i++)
  {
    fprintf(Out, "%3ld: av: ", i);
    for(j=0; j<Tp->ntm; j++)
      fprintf(Out, "%le ", Tm->ip[i].av[j]);
    fprintf(Out, "%3s: pv: ", "");
    for(j=0; j<Tp->ntm; j++)
      fprintf(Out, "%le ", Tm->ip[i].pv[j]);
    fprintf(Out, "\n%3s: grad: ", "");
    for(j=0; j<Tp->ntm; j++)
    {
      for(k=0; k<Tm->ip[i].ncompgrad; k++)
        fprintf(Out, "%le ", Tm->ip[i].grad[j][k]);
      if (j < Tp->ntm-1)
        fprintf(Out, "\n%11s", "");
    }
    fprintf(Out, "\n%3s: flux: ", "");
    for(j=0; j<Tp->ntm; j++)
    {
      for(k=0; k<Tm->ip[i].ncompgrad; k++)
        fprintf(Out, "%le ", Tm->ip[i].fluxes[j][k]);
      if (j < Tp->ntm-1)
        fprintf(Out, "\n%11s", "");
    }
    fprintf(Out, "\n%3s: other: ", "");
    for(j=0; j<Tm->ip[i].ncompother; j++)
      fprintf(Out, "%le ", Tm->ip[i].other[j]);
    fprintf(Out, "\n%3s: eqother: ", "");
    for(j=0; j<Tm->ip[i].ncompeqother; j++)
      fprintf(Out, "%le ", Tm->ip[i].eqother[j]);
    fprintf(Out, "\n");
  }

  fprintf(Out, "\nNontransq array:\n");
  for(i=0; i<Tm->tnip; i++)
  {
    for(j=0; j<Tm->nntq; j++)
      fprintf(Out, "%le ", Tm->nontransq[j][i]);
    fprintf(Out, "\n");
  }
}
