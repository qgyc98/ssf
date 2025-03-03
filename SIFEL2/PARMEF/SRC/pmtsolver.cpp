#include "mpi.h"
#include "pnewtonraph.h"
#include "pmtsolver.h"
#include "pglobal.h"
#include "seqfilesm.h"
#include "genfile.h"
#include <string.h>


void par_solve_timemech_prob ()
{
  long lcid;
  
  //  load case must be equal to zero, no superposition can be used
  //  in this type of analysis
  lcid=0;
  
  //  solver of the problem
  //par_solve_timemech_prob (lcid);
  par_solve_timemech_prob2 (lcid);
  
}

/**
  Function solves parallel time dependent problem (viscoplasticity).
  It uses one_step concept instead of standard one.
  
  @param lcid - load case id
   
  @return The function does not return anything.
   
  TKr, 19/02/2013 according to Tomas Koudelka (mtsolver.cpp)
  Actualized by TKo, 09.2020
*/
void par_solve_timemech_prob2 (long lcid)
{
  long i, li, nsts, ret, rest_calc;
  double dt, dtr, dtmin, dtmax, end_time, prev_time, time;
  mt_glob_vec mt_gv;

  //
  //  initialization phase
  //
  par_visco_solver_init(lcid, rest_calc, mt_gv);

  // store initial value for attained time in the previous time step
  // (prev_time is exploited in the growing mechanical problems)
  prev_time = Mp->time;

  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  //  minimum time increment
  dtmin=Mp->timecon.dtmin;
  //  maximum time increment
  dtmax=Mp->timecon.dtmax;
  //  end time
  end_time = Mp->timecon.endtime ();

  //  number of step
  li = i = mt_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;
  // return status from the one_step function
  ret = 0;

  // just for evaluation of hypoplasticity model ???!!!  
  Neval = 0.0;

  if (Mp->tprob == growing_mech_structure)  
    //saving values for the first iteration
    Mt->save_nodval(lcid);

  // ***************************
  //   main iteration loop  ****
  // ***************************
  do
  {
    //  new step number
    i++;
        
    //  computation of backward time step - only for adaptive time controler or nonconstant time increment
    if((Mp->timecon.timefun.tfunc != constant) && Mp->timecon.tct == 0)
      dt = Mp->timecon.actualbacktimeincr ();

    // perform one time step
    time=Mp->timecon.newtime(dt);
    Mp->dtime = dt;
    if ((Mp->tprob == growing_mech_structure) && (prev_time == Mp->timecon.starttime()))
      prev_time = time;
    
    // perform one time step
    ret = par_one_step(lcid, time, dt, dtr, prev_time, rest_calc, i, li, mt_gv);
    
    if (ret >= 0){ // equilibrium was attained
    
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
        prev_time = Mp->time;
      }

      //  update of values stored at integration points
      Mm->updateipval();
      compute_req_val (lcid);
      print_step(lcid, Mp->istep, Mp->time, mt_gv.fl);
      print_flush();
	
      if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.save_stat()))
      {
        if (Myrank==0 && Mespr==1)
          fprintf (stdout,"\n Creating MEFEL backup file\n");
        solver_save (mt_gv.r,mt_gv.fp,Mp->istep,Mp->time,dt,&Mp->timecon,Ndofm);
      }
      if (ret == 0)
        nsts++;      
      if (nsts==2)
      {
        if (dt < dtmax)
        {
          dt*=2.0;
          if (dt>dtmax)//maximum time increment
            dt = dtmax;
          nsts=0;
          if (Myrank==0 && Mespr==1)
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 2 steps");
        }
      }
    }
    else{
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      nsts = 0;
      
      if (dtr < 1.0)
        dt *= dtr;
      else
        dt/=2.0;
      Mp->timecon.oldtime ();
      i--;
      
      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();
      
      if (Myrank==0 && Mespr==1)
        fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
        if (Myrank==0 && Mespr==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
        if (Myrank==0 && Mespr==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
        if (Myrank==0 && Mespr==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
        //  update of values stored at integration points
        Mm->updateipval();
        compute_req_val (lcid);
        print_step_forced(lcid, i, Mp->time, mt_gv.fb);
        //print_step_forced(lcid, i, Mb->dlc[0].gf[0].getval(Mp->time), fb);
        print_flush();
        break;
      }
      // change status of elements and nodes according to the previous time step
      if (Mp->tprob == growing_mech_structure)
      {
        update_elemnod_prev_time(prev_time, mt_gv.ncd, mt_gv.nce);
      }
    }
  }while(Mp->time<end_time);
  
  print_close ();
}

/**
  Function allocates and initializes global vectors used in pmtsolver
  It is used in one_step concept solver.

  TKr, 19/02/2013 according to Tomas Koudelka (mtsolver.cpp)
  Actualized by TKo, 09.2020
*/
void par_visco_solver_init(long lcid, long &rest_calc, mt_glob_vec &mt_gv)
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
  Mp->nlman->divc_step = new double[10];
  Mp->nlman->divc_err  = new double[10];

  rest_calc = 0;

  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat())
  {
    // initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    // initiation of material models on auxiliary integration points
    Mm->aip_initmaterialmodels(lcid, true);
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file\n");
    solver_restore (mt_gv.r, mt_gv.fp, mt_gv.istep, Mp->time, dt, &Mp->timecon, n);
    Mm->updateother();

    //  Backup of r is used in case that the Newton-Raphson procedure does not converge.
    //  lhsb[j]=r[j]
    copyv(mt_gv.r, mt_gv.lhsb, n);

    if (Mp->tprob == growing_mech_structure)    
      update_elnod_stat_after_hdbrest(lcid);

    // update actual material index for all elements
    Mm->update_actual_mat_id(lcid, 0, true);
    // set indicator of calculation started from restored backup data
    rest_calc = 1;
    print_init(-1, "at",Pmp->fni,Pmp->fei);
  }
  else
  {
    mt_gv.istep=0;
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    // initiation of material models on auxiliary integration points
    Mm->aip_initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    // compute stresses from eigenstrains or initial irreversible strains
    if ((Mp->eigstrains > 0) && (Mp->eigstrains < 4))
      internal_forces(lcid, mt_gv.fp);
    // print initial step to the output files
    print_init(-1, "wt",Pmp->fni,Pmp->fei);
    print_step(lcid, mt_gv.istep, Mp->time, mt_gv.f);
    //print_step(lcid, i, 0.0, f);
    print_flush();
  }
}

/**
  Function executes one step of parallel time dependent mechanical algorithm
  (visco-elasticity, visco-plasticity, etc.)
  
  @param lcid - load case id
  @param time - actual time
  @param dt   - actual time increment
  @param dtr   - the minimum ratio of required time step size to the actual one required by material models (output)
  @param prev_time - attained time in the previous time step (needed only for growing structures)
  @param rest_calc - inidicator whether step follows restorage from the backup immediately (=1) or no (=0)
  @param istep - time step id
  @param li    - time step id of the first performed time step
  @param mt_gv - structure with pointers to %vectors of righthand and lefthand side
  
  @retval  0 - if the equilibrium was attained with NO inner loop
  @retval  1 - if the equilibrium was attained and inner loop WAS performed
  @retval -1 - if the equilibrium was NOT attained
  
  TKr, 19/02/2013 according to JK and TKo
  Actualized by TKo, 09.2020
*/
long par_one_step (long lcid,double time, double dt, double &dtr, double prev_time, long rest_calc, long istep, long li, mt_glob_vec &mt_gv)
{
  long j,n,ini,um, mncd, mnce, mnae;
  double *f,*fl,*fi,*fb,*fp,*r,*dr,*lhsb, *flp;
  long *ifn;
  double norf, err, zero;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  answertype fusm;
  
  //new time increment
  Mp->time = time;
  //  new step number
  Mp->istep = istep;
  Mp->jstep = -1;
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
  //  vector of prescribed forces due to force load from the previous step
  flp = mt_gv.flp;
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  fl = mt_gv.fl;
  //  vector of internal forces
  fi = mt_gv.fi;
  //  auxiliary force vector
  fb = mt_gv.fb;
  //  backup of the nodal displacements
  lhsb = mt_gv.lhsb;
  // indicators of nodes on interfaces between old and new parts of the structure
  ifn = mt_gv.ifn;
  
  //  number of mechanical degrees of freedom
  n=Ndofm;
  
  // update actual material index for all old elements and call initval for those whose indeces have changed
  if (rest_calc == 0)
  {
    um = Mm->update_actual_mat_id(lcid, 1, false);
    if (um)   fusm = yes;
    else      fusm = no;
  }  
  else  // do not update if the calculation was restored from HD
    fusm = no;

  // update statuses of all nodes, DOFS and elements for growing structures,
  // calculation of initial displacements of new parts of the structure
  if (Mp->tprob == growing_mech_structure)
  {    
    update_elnod_stat(lcid, istep, prev_time, ifn, r, fb, fp, mnce, mncd, mnae);
    mt_gv.ncd = mncd;
    mt_gv.nce = mnce;
    n = Ndofm;
  }

  // maximum number of steps for inner loop
  ini = Mp->nlman->niilnr;
  
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  
  //  computer zero
  zero = Mp->zero;
  
  //  backup of attained nodal values
  for (j=0;j<n;j++)
    lhsb[j]=r[j];
  
  
  if (Myrank==0 && Mespr==1)  fprintf (stdout,"\n\n --------------------------------------------------------------");
  if (Myrank==0 && Mespr==1)  fprintf (stdout,"\n PMEFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Mp->time,dt);
  if (Myrank==0 && Mespr==1)  fprintf (stdout,"\n --------------------------------------------------------------\n");
  
  //  assembling of vectors of prescribed nodal forces
  mefel_right_hand_side (lcid,f,fl);
  
  //  computation of sums of force components in particular directions
  //  it is used in special creep and consolidation models
  Mb->comp_sum (f);
  //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
  //      istep, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
  fflush(Out);
  
  // assemble force vector due to removed elements and apply it gradually if required
  if (Mp->tprob == growing_mech_structure){
    forces_due2removed_elem(lcid, mnce-mnae, istep, prev_time, fi, fp, flp, fb, r);
  }

  //  increments of prescribed forces
  incr_internal_forces (lcid,fi);
  
  // assembling of right-hand side
  // fb[j]=f[j]-fp[j]+fi[j]
  subv(f, fp, fb, n);
  addv(fb, fi, fb, n);

  // perform one load step with the help of Newton-Raphson procedure
  j = 0;
  dtr = 1.0;
  norf = par_gnewton_raphson_one_step(lcid, Mp->nlman, fl, r, fb, dr, fi, dtr, istep-1, j, li, ini, err);
  
  if (norf>err){
    // equilibrium state was not attained
    if ((Myrank == 0) && (Mespr==1)) fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);
    
    //  backup is used for actual values
    //  r[k]=lhsb[k]
    copyv(lhsb, r, n);
    return -1;
  }
  
  // equilibrium state was attained
    if ((Myrank == 0) && (Mespr==1) && (j > 0))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
  
  
  // print resulting vector of internal forces
  Mb->comp_sum (fi);
  //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
  //      istep, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
  //fflush(Out);
  
  
  if (j < 0)
    return 0;
  
  //return j+1;
  if (Mp->jstep > 1)
    return 1;
  else
    return 0;
}


/**
   function solves time dependent mechanical problem in parallel
   inertial forces are negligible
   
   JK, 24.7.2008
*/
void par_solve_timemech_prob (long lcid)
{
  long i,j,k,n,ini,nsts,nii;
  double end_time,zero,err,dt,dtmin,dtdef,norfa,norf;
  double *f,*fi,*fp,*fb,*fl,*dr,*r,*lhsb;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  double lsm_res;
  
  //  initialization phase
  
  //  maximum number of iterations in inner loop
  ini = Mp->nlman->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  //  number of rows of the matrix
  n = Ndofm;
  //  computer zero
  zero = Mp->zero;
  
  
  //  vector of nodal displacements
  r = Lsrs->give_lhs (lcid);
  //  vector of nodal prescribed forces
  f = Lsrs->give_rhs (lcid);

  // allocation of vectors  

  //  vector of increments of nodal displacements
  dr   = new double [n];
  nullv (dr,n);
  //  vector of prescribed forces from the previous step
  fp   = new double [n];
  nullv (fp,n);
  //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  fl   = new double [n];
  nullv (fl,n);
  //  vector of internal forces
  fi   = new double [n];
  nullv (fi,n);
  //  auxiliary force vector
  fb   = new double [n];
  nullv (fb,n);
  //  backup of the nodal displacements
  lhsb = new double [n];
  nullv (lhsb,n);
    
  //  initial time
  Mp->time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  //  minimum time increment
  dtmin=Mp->timecon.dtmin;
  
  //  number of step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;


  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file\n");
    solver_restore (r,fp,i,Mp->time,dt,&Mp->timecon,n);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    print_init(-1, "at",Pmp->fni,Pmp->fei);
    print_step(lcid, i, Mp->time, f);
    print_flush();
  }
  else
  {
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    print_init(-1, "wt",Pmp->fni,Pmp->fei);
    print_step(lcid, i, Mp->time, f);
    print_flush();
  }
  //print_close();


  // ***************************
  //   main iteration loop  ****
  // ***************************
  do{
    //new time increment
    Mp->time = Mp->timecon.newtime (dt);
    //  new step number
    i++;
    //  indicator of strain computation
    //  it means, strains at integration points have not been computed
    Mp->strainstate=0;
    //  indicator of stress computation
    //  it means, stresses at integration points have not been computed
    Mp->stressstate=0;
    //  indicator of computation of other array
    //  it means, stresses at integration points have not been computed
    Mp->otherstate=0;
    
    //  backup of attained nodal values
    for (j=0;j<n;j++){
      lhsb[j]=r[j];
    }
    
    if ((Myrank==0)&&(Mespr==1))  fprintf (stdout,"\n\n -------------------------------------------------------------");
    if ((Myrank==0)&&(Mespr==1))  fprintf (stdout,"\n Time step = %ld,  Time %e,  Time increment = %e",i,Mp->time,dt);
    if ((Myrank==0)&&(Mespr==1))  fprintf (stdout,"\n -------------------------------------------------------------\n");

    //  assembling of vectors of prescribed nodal forces
    mefel_right_hand_side (lcid,f,fl);

    //  computation of sums of force components in particular directions
    //  it is used in special creep and consolidation models
    Mb->comp_sum (f);
    //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
    //i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
    fflush(Out);
    
    //  increments of prescribed forces
    incr_internal_forces (lcid,fi);
    
    //  right hand side assembling
    for (j=0;j<n;j++){
      fb[j]=f[j]-fp[j]+fi[j];
    }    

    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == 1) || (Smat==NULL))
      stiffness_matrix (lcid);

    //  solution of K(r).v=F
    Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);

    //  new displacement vector
    for (j=0;j<n;j++){
      r[j]+=dr[j];
    }
    
    //  computation of internal forces
    internal_forces (lcid,fi);
        
    // norm of total load vector      
    norfa=Psolm->pss (f,f,Out);       
    // vector of unbalanced forces
    for (j=0;j<n;j++){
      fb[j] = fl[j] - fi[j];
    }
    //  norm of vector of unbalanced forces
    norf=Psolm->pss(fb,fb,Out);

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
        //  assembling of tangent stiffness matrix
        if ((Mp->nlman->stmat==tangent_stiff) && (j%5 == 0))
          stiffness_matrix (lcid);
	//  solution of K(r).v=F
	Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);
	
        //  new displacement vector
	for (k=0;k<n;k++)
	  r[k]+=dr[k];
	
	//  computation of internal forces
	internal_forces (lcid,fi);
	
	//  vector of unbalanced forces
	for (k=0;k<n;k++)
	  fb[k]=fl[k]-fi[k];
	
	//  norm of vector of unbalanced forces
	norf=Psolm->pss(fb,fb,Out);	
        if (norfa > 0.0)
	  norf /= norfa;

//	if ((Myrank==0) && (Mespr==1) && (j%10==0)) // print each 10-th step
	if ((Myrank==0) && (Mespr==1))
	  fprintf (stdout,"\n Time increment %lf   inner loop j %ld     norf=%e",Mp->time,j,norf);
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
    
    if (nii==0)
    {
      // equilibrium state was not attained
      if ((Myrank==0) && (Mespr==1))
        fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);

      //  backup is used for actual values
      for (j=0;j<n;j++){
	r[j]=lhsb[j];
      }

      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      Mp->timecon.oldtime ();
      i--;
      
      if ((Myrank==0) && (Mespr==1))
        fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (Psolm->compare_on_master(dt, dtmin) < 0){ // dt<dtmin
        fflush(stdout);
	if ((Myrank==0) && (Mespr==1)) fprintf (stderr,"\n\n time increment is less than minimum time increment");
	if ((Myrank==0) && (Mespr==1)) fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
        compute_req_val (lcid);
        print_step_forced(lcid, i, Mp->time, f);
        print_flush();
        break;
      }      
    }
    else{
      // equilibrium state was attained
      if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);

      //  actual total load vector becomes previous total load vector
      for (j=0;j<n;j++){
        fp[j]=f[j];
      }

      dtdef = Mp->timecon.actualforwtimeincr ();
      if (nsts==2){
	dt*=2.0;
	nsts=0;
	if ((Myrank==0) && (Mespr==1))
          fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
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
      print_step(lcid, i, Mp->time, f);
      print_flush();

      if ((Mp->timecon.isitimptime ()==1) && (Mp->hdbcont.hdbtype))
      {
        if ((Myrank==0) && (Mespr==1))
	  fprintf (stdout,"\n Creating backup file\n");
	solver_save (r,fp,i,Mp->time,dt,&Mp->timecon,n);      
      }
    }

  }while(Mp->time<=end_time);
    
  print_close ();
  

  delete [] fi;  
  delete [] fb;  
  delete [] fp;  
  delete [] dr; 
  delete [] fl;  
  delete [] lhsb;
  

}
