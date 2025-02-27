#include "mpi.h"
#include "pglobalt.h"
#include "pcpcsolver.h"
#include "pglobalc.h"
#include "pcpsolver.h"
#include "pglobal.h"
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
#include <math.h>
#include <string.h>

/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   for growing constructions,  transport problems are linear or nonlinear
   mechanical analysis is always nonlinear

   4.2009 TKo
*/
void par_solve_gpcouplprob ()
{
  long lcid;
  
  //  load case id must be equal to zero
  //  the problems are nonlinear and superposition method cannot be used
  lcid=0;
  
  switch (Cp->tnlinsol){
  case newtonc:{
    par_newton_raphson_gparcoupl (lcid);
    break;
  }
  default:{
    print_err("unknown solver of nonlinear equation system is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   for growing constructions by the Newton-Raphson method

   @param lcid - load case id

   TKo 4.2009
*/
void par_newton_raphson_gparcoupl (long lcid)
{
  switch (Tp->tprob){
  case growing_np_problem:{
    //  linear nonstationary transport problem
    //  nonlinear mechanical problem
    par_newton_raphson_gparcoupl_lin (lcid);
    break;
  }
  case growing_np_problem_nonlin:{
    //  nonlinear nonstationary transport problem
    //  nonlinear mechanical problem
    par_newton_raphson_gparcoupl_nonlin (lcid);
    break;
  }
  default:{
    print_err("unknown TRFEL problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   Function solves system of linear TRFEL algebraic and non-linear MEFEL equations by Newton-Raphson method
   for time-dependent problems for growing structure

   @param lcid - load case id
   
   TKo 4.2009
*/
void par_newton_raphson_gparcoupl_lin (long lcid)
{
  int mefel_convergence;
  long i,ii,j,k,nm,nt,ini,nsts,mncd,tncd,mnce,tnce,nii,mnae,nifn;
  long *ifn;
  double zero,dt,dtmin,dtdef,end_time,alpha,s,prev_timem,prev_timet,btime;
  double *dr,*fb,*r,*d,*f,*fp,*fi,*fl,*p,*lhst,*tdlhst,*rhst,*lhsb,*lhstb,*tdlhstb;
  double norf, norfa, err, lsm_res;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);

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
  // array with marked nodes on interface between new and old part of structure
  ifn = new long[Mt->nn];
  memset(ifn, 0, sizeof(*ifn)*Mt->nn);

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
  //  backup of time
  prev_timem = Mp->time;
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);
  //  backup of time
  prev_timet = Tp->time;

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
    Mt->rhs_save (fp);
    //  approximation of temperatures from TRFEL nodes to MEFEL integration points
    approximation_inittemper ();
    //  approximation of humidities from TRFEL nodes to MEFEL integration points
    //approximation_inithumid ();
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    print_init(-1, "at",Pcp->fnim,Pcp->feim);
//    print_step(lcid,i,Mp->time,f);
//    print_flush();
  }
  else
  {
    //  approximation of temperatures from TRFEL nodes to MEFEL integration points
    approximation_inittemper ();
    //  approximation of humidities from TRFEL nodes to MEFEL integration points
    //approximation_inithumid ();
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
    Tm->initmaterialmodels ();
    compute_req_valt (0);
    print_initt(-1, "at",Pcp->fnit,Pcp->feit);
//    print_stept(0,i,Tp->time,NULL);
//    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels ();
    compute_req_valt (0);
    print_initt(-1, "wt",Pcp->fnit,Pcp->feit);
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }
  //print_close();

  mefel_convergence = 1;//setup of flag for successfully performed MEFEL step to yes

  //saving values for the first iteration
  Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
  //saving values for the first iteration
  Mt->save_nodval(lcid);
  //Mt->lhs_save (r,1);

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
    Tp->timecont = Cp->timecon;//this should be corrected in future
    i++;

    if(mefel_convergence == 1){
      if (Mesprc==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------------------");
      if (Mesprc==1)  fprintf (stdout,"\n NEW METR TIME STEP = %ld: METR TIME = %e, metr time increment = %e",i,Cp->time,dt);
      if (Mesprc==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------------------");
    }
    if (Mesprt==1)  fprintf (stdout,"\n ------------------------------------------------------------------------------------");
    if (Mesprt==1)  fprintf (stdout,"\n TRFEL TIME STEP = %ld, TRFEL TIME = %e, trfel time increment = %e",i,Tp->time,dt);
    if (Mesprt==1)  fprintf (stdout,"\n ------------------------------------------------------------------------------------");
    
    //  searching for new elements
    //  only new elements are switched on
    //  former elements are temporarily switched off
    tnce = Gtt->search_newelem (Tp->time,prev_timet);
    
    //  assignment of initial nodal values
    //  it works for problem with growing number of elements
    //  it does not work generally for problems with increasing and decreasing number of elements
    if (tnce>0){
      Tt->initial_nodval ();
    }

    //  searching for new DOFs
    //  only new DOFs are switched on
    //  former DOFs are temporarily switched off
    tncd = Gtt->search_newdofs (Tp->time,prev_timet);
 
    //  marking of elements switched on
    Gtt->update_elem (Tp->time);
    
    //  marking of nodes switched on
    Gtt->update_nodes ();
    
    //  marking of DOFs switched on
    Gtt->update_dofs (Tp->time);
    
    //  generation of new code numbers
    // Ndoft = Gtt->codenum_generation (Outt);
    nt=Ndoft;
            
    //  cleaning of the conductivity matrix in case of changes
    if (tncd>0){
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

    if (tnce!=0 || tncd!=0){
      //  assembling of displacement vector for new code numbers
      Tt->lhs_restore (lhst,Lsrst->lhsi,tdlhst);
    }

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
    
    for (j=0;j<nt;j++){
      rhst[j] = rhst[j] - d[j];
      d[j]=tdlhst[j];
    }
    
    //  solution of the system of algebraic equations
    Psolt->par_linear_solver (Gtt,Kmat,tdlhst,rhst,Outt,Mesprt);

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

    Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);

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

      //  backup of attained nodal values
      for (j=0;j<nm;j++){
	lhsb[j]=r[j];
      }

      if (Mespr==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------");
      if (Mespr==1)  fprintf (stdout,"\n MEFEL TIME STEP = %ld, MEFEL TIME = %e, mefel time increment = %e",ii,Mp->time,dt);
      if (Mespr==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------");

      //  searching for changed elements
      //  status of all elements is updated
      mnce = Gtm->search_changed_elem(Mp->time, mnae);    
      //  searching for new DOFs and
      //  status of nodes is updated
      mncd = Gtm->search_changed_dofs(Mp->time,prev_timem);
      
      nifn = 0;
      ifn = NULL;
      if (mnce!=0){
	// create list of nodes that belong to interfaces
	// between old and new elements
	nifn = Gtm->search_iface_nodes(ifn);
      }  
      
      
      if ((mnce!=0) || (mncd!=0))
	{
	  // nodal initial displacements investigation in case of 
	  // new supports, zero nodal values of switched off nodes
	  Mt->save_node_inidispl(lcid, Mp->time, prev_timem);        
	  
	  //
	  // CALCULATION OF INITIAL DISPLACEMENTS FOR NEW PARTS OF STRUCTURE
	  //
	  
	  if (nifn) // if there are some interfaces between old and new part of structure
	    {
	      // switch on only new elements
	      Gtm->switch_new_elem();
	      
	      // switch on nodes on the new part of structure
	      Gtm->update_nodes();
	      
	      // prepare active nodes (including nodes on interface)
	      // for generation of code numbers
	      Gtm->update_active_dofs(Mp->time);
	      
	      // new DOF numbering
	      Ndofm = Gtm->codenum_generation(Out);
	      
	      // restore nodal displacements for interface nodes 
	      // (zero displacements will be in the rest of nodes)
	      Mt->restore_nodval(r, Ndofm);
	      
	      // calculate stresses due to displacements on the interface nodes
	      stress_initdispl(lcid);
	      
	      // prepare active nodes (including nodes on interface)
	      // for generation of code numbers
	      Gtm->update_active_dofs(Mp->time);
	      // EXCLUDE INTERFACE NODES
	      Gtm->clear_intf_dofs(ifn);
	      
	      // new DOF numbering
	      Ndofm = Gtm->codenum_generation(Out);
	      // use integration point stresses that have been already calculated
	      // to establish nodal forces due to prescribed displacements
	      Mp->strcomp = 0;
	      internal_forces(lcid,fb);
	      cmulv(-1.0, fb, Ndofm);
	      Mp->strcomp = 1;
	      
	      
	      // assamble new stiffness matrix, reset left hand side vector
	      delete Smat;
	      Smat=NULL;
	      stiffness_matrix (lcid);
	      nullv(r, Ndofm);
	      
	      // solution of system of equations
	      //Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);
	      //  solution of K(r).v=F
	      Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);
	      
	      // update vector of total displacements
	      addv(r, dr, nm);
	      
	      // store initial displacements on the new part of structure
	      // in order to prescribe stress free state on required new elements
	      Mt->save_elem_inidispl(lcid, ifn);
	      
	      // switch off interface nodes in order to prevent 
	      // the values in the interface nodes to be rewritten
	      Gtm->remove_nodes(ifn);
	      
	      // store attained nodal values due to prescribed values on interfaces
	      // on the new part of structure
	      Mt->save_nodval(lcid);
	    }
	  
	  // update true status of elements
	  Gtm->update_elem(Mp->time);
	  
	  // switch on nodes on the new part of structure
	  Gtm->update_nodes();
	  
	  //  marking of DOFs switched on
	  Gtm->update_dofs (Mp->time);
	  
	  // new DOF numbering
	  nm = Ndofm = Gtm->codenum_generation(Out);
	  
	  // destruct old system matrix
	  delete Smat;
	  Smat=NULL;
	  
	  // restore vector of attained unknowns
	  Mt->restore_nodval(r,nm);
	  
	  //
	  // Reassemble vector fp of attained load vector from the previous time step
	  //
	  
	  // switch on elements accroding to previous time step
	  Gtm->update_elem(prev_timem);
	  // switch on nodes accroding to previous time step     
	  Gtm->update_nodes();
	  // backup actual time
	  btime = Mp->time;
	  // set previous time in order to assemble load vector fp
	  Mp->time = prev_timem;
	  // DOF numbers remains the same as for the actual time step
	  mefel_right_hand_side(lcid, fp, fl);
	  Mp->time = btime;
	  // switch on elements accroding to ACTUAL time step
	  Gtm->update_elem(Mp->time);
	  // switch on nodes accroding to ACTUAL time step     
	  Gtm->update_nodes();
	  
	  // update reaction status indicators at nodes and elements
	  Mt->comreac();
	  // update status indicators of prescribed displacements on elements
	  Mt->elemprescdisp();
	  // prescribe stress free state on the integration points of new parts
	  Mt->clean_ip_new_elem();
	}
      
      //  approximation of temperatures from TRFEL nodes to MEFEL integration points
      approximation_temper ();
      //  approximation of humidities from TRFEL nodes to MEFEL integration points
      approximation_humid ();

      //  assembling of vectors of prescribed nodal forces
      mefel_right_hand_side (lcid,f,fl);
      
      //  computation of sums of force components in particular directions
      //  it is used in special creep and consolidation models
      Mb->comp_sum (f);
      Mb->comp_sum_react();
      Mb->comp_sum_pdreact();
      fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
              i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      fflush(Out);
      
      if ((mnce - mnae) > 0){
	// there were some removed elements =>
	// calculate nodal forces due to removed elements
	
	// use integration point stresses that have been already calculated
	Mp->strcomp = 0;
	internal_forces(lcid,fi);
	Mp->strcomp = 1;
	
	//  fb[j] = fp[j] - fi[j];
	subv(fp, fi, fb, nm);
      }
      else
	nullv(fb, nm);
      
      //  computation of unbalanced forces which will be on the right hand side
      incr_internal_forces (lcid,fi);
      
      //  right hand side assembling
      for (j=0;j<nm;j++){
	fb[j]=f[j]-fp[j]+fi[j];
      }    
      
      //  assembling of tangent stiffness matrix
      if ((Mp->nlman->stmat==tangent_stiff) || (ii == 1) || (Smat==NULL))
	stiffness_matrix (lcid);
      
      //  solution of K(r).v=F
      Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);
      
      //  update total displacement vector  r[j]+=dr[j]
      addv(r, dr, nm);
      
      //  computation of internal forces
      internal_forces (lcid,fi);

      // norm of total load vector      
      norfa=Psolm->pss (f,f,Out);
      
      // vector of unbalanced forces fb[j] = fl[j] - fi[j]
      subv(fl, fi, fb, nm);      
      
      //  norm of vector of unbalanced forces
      norf=Psolm->pss (fb,fb,Out);
      if (norfa > 0.0)
	norf /= norfa;
      
      if ((Myrank==0)&&(Mespr==1))
        fprintf (stdout,"\n\n Norf before inner iteration loop, norf = %e\n\n\n",norf);
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
	    
	  //  vector of unbalanced forces fb[k] = fl[k] - fi[k]
	  subv(fl, fi, fb, nm);
	  
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
	
	//  backup of nodal values and time because code numbers may be changed
	Mt->save_nodval(lcid);
	
	// new elements can be definitely changed to old ones
	Gtm->update_auxinf();
	prev_timem = Mp->time;
	
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
        fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
                i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
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
  delete [] ifn;

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
void par_newton_raphson_gparcoupl_nonlin (long lcid)
{
  int mefel_convergence;
  long i,ii,j,k,nm,nt,ini,init,nsts,mncd,tncd,mnce,tnce,nii,mnae,nifn;
  long *ifn;
  double zero,dt,dtmin,dtdef,end_time,alpha,prev_timem,prev_timet,btime;
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
  errt = Tp->errarr[0];//docasne, nutne upravit??!!
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
  // array with marked nodes on interface between new and old part of structure
  ifn = new long[Mt->nn];
  memset(ifn, 0, sizeof(*ifn)*Mt->nn);

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
  //  backup of time
  prev_timem = Mp->time;
  //  transport time controller is rewritten by coupled time controller
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);
  //  backup of time
  prev_timet = Tp->time;

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
    Mt->rhs_save (fp);
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
    Tm->initmaterialmodels ();
    compute_req_valt (0);
    print_initt(-1, "at",Pcp->fnit,Pcp->feit);
    //    print_stept(0,i,Tp->time,NULL);
    //    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels ();
    compute_req_valt (0);
    print_initt(-1, "wt",Pcp->fnit,Pcp->feit);
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }
  //print_close();

  mefel_convergence = 1;//setup of flag for successfully performed MEFEL step to yes

  //saving values for the first iteration
  Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
  //saving values for the first iteration
  Mt->save_nodval(lcid);
  //Mt->lhs_save (r,1);

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
    
    //  searching for new elements
    //  only new elements are switched on
    //  former elements are temporarily switched off
    tnce = Gtt->search_newelem (Tp->time,prev_timet);
    
    //  assignment of initial nodal values
    //  it works for problem with growing number of elements
    //  it does not work generally for problems with increasing and decreasing number of elements
    if (tnce>0){
      Tt->initial_nodval ();
    }

    //  searching for new DOFs
    //  only new DOFs are switched on
    //  former DOFs are temporarily switched off
    tncd = Gtt->search_newdofs (Tp->time,prev_timet);
 
    //  marking of elements switched on
    Gtt->update_elem (Tp->time);
    
    //  marking of nodes switched on
    Gtt->update_nodes ();
    
    //  marking of DOFs switched on
    Gtt->update_dofs (Tp->time);
    
    //  generation of new code numbers
    // Ndoft = Gtt->codenum_generation (Outt);
    nt=Ndoft;
            
    //  cleaning of the conductivity matrix in case of changes
    if (tncd>0){
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

    if (tnce!=0 || tncd!=0){
      //  assembling of displacement vector for new code numbers
      Tt->lhs_restore (lhst,Lsrst->lhsi,tdlhst);
    }

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

    Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);

    //  approximation of nodal values into ontegration points
    approximation ();
    
    norf_last = 1.0e20;    
    //inner iteration loop for trfel
    for (j=0;j<init;j++){      
      //  physically corrected solution
      solution_correction ();
      
      Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
      
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
      
      Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
      
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

	  Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);

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
      
      //  backup of attained nodal values
      for (j=0;j<nm;j++){
	lhsb[j]=r[j];
      }
      
      if (Mespr==1)  fprintf (stdout,"\n\n -----------------------------------------------------------------------------------");
      if (Mespr==1)  fprintf (stdout,"\n MEFEL TIME STEP = %ld, MEFEL TIME = %e, mefel time increment = %e",ii,Mp->time,dt);
      if (Mespr==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------");
      
      //  searching for changed elements
      //  status of all elements is updated
      mnce = Gtm->search_changed_elem(Mp->time, mnae);    
      //  searching for new DOFs and
      //  status of nodes is updated
      mncd = Gtm->search_changed_dofs(Mp->time,prev_timem);
      
      nifn = 0;
      if (mnce!=0){
	// create list of nodes that belong to interfaces
	// between old and new elements
	nifn = Gtm->search_iface_nodes(ifn);
      }  
      
      if ((mnce!=0) || (mncd!=0))
	{
	  // nodal initial displacements investigation in case of 
	  // new supports, zero nodal values of switched off nodes
	  Mt->save_node_inidispl(lcid, Mp->time, prev_timem);        
	  
	  //
	  // CALCULATION OF INITIAL DISPLACEMENTS FOR NEW PARTS OF STRUCTURE
	  //
	  
	  if (nifn) // if there are some interfaces between old and new part of structure
	    {
	      // switch on only new elements
	      Gtm->switch_new_elem();
	      
	      // switch on nodes on the new part of structure
	      Gtm->update_nodes();
	      
	      // prepare active nodes (including nodes on interface)
	      // for generation of code numbers
	      Gtm->update_active_dofs(Mp->time);
	      
	      // new DOF numbering
	      Ndofm = Gtm->codenum_generation(Out);
	      
	      // restore nodal displacements for interface nodes 
	      // (zero displacements will be in the rest of nodes)
	      Mt->restore_nodval(r, Ndofm);
	      
	      // calculate stresses due to displacements on the interface nodes
	      stress_initdispl(lcid);
	      
	      // prepare active nodes (including nodes on interface)
	      // for generation of code numbers
	      Gtm->update_active_dofs(Mp->time);
	      // EXCLUDE INTERFACE NODES
	      Gtm->clear_intf_dofs(ifn);
	      
	      // new DOF numbering
	      Ndofm = Gtm->codenum_generation(Out);
	      // use integration point stresses that have been already calculated
	      // to establish nodal forces due to prescribed displacements
	      Mp->strcomp = 0;
	      internal_forces(lcid,fb);
	      cmulv(-1.0, fb, Ndofm);
	      Mp->strcomp = 1;
	      
	      
	      // assamble new stiffness matrix, reset left hand side vector
	      delete Smat;
	      Smat=NULL;
	      stiffness_matrix (lcid);
	      nullv(r, Ndofm);
	      
	      //  solution of K(r).v=F
	      Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);
	      
	      // update vector of total displacements
	      addv(r, dr, nm);
	      
	      // store initial displacements on the new part of structure
	      // in order to prescribe stress free state on required new elements
	      Mt->save_elem_inidispl(lcid, ifn);
	      
	      // switch off interface nodes in order to prevent 
	      // the values in the interface nodes to be rewritten
	      Gtm->remove_nodes(ifn);
	      
	      // store attained nodal values due to prescribed values on interfaces
	      // on the new part of structure
	      Mt->save_nodval(lcid);
	    }
	  
	  // update true status of elements
	  Gtm->update_elem(Mp->time);
	  
	  // switch on nodes on the new part of structure
	  Gtm->update_nodes();
	  
	  //  marking of DOFs switched on
	  Gtm->update_dofs (Mp->time);
	  
	  // new DOF numbering
	  nm = Ndofm = Gtm->codenum_generation(Out);
	  
	  // destruct old system matrix
	  delete Smat;
	  Smat=NULL;
	  
	  // restore vector of attained unknowns
	  Mt->restore_nodval(r,nm);
	  
	  //
	  // Reassemble vector fp of attained load vector from the previous time step
	  //
	  
	  // switch on elements accroding to previous time step
	  Gtm->update_elem(prev_timem);
	  // switch on nodes accroding to previous time step     
	  Gtm->update_nodes();
	  // backup actual time
	  btime = Mp->time;
	  // set previous time in order to assemble load vector fp
	  Mp->time = prev_timem;
	  // DOF numbers remains the same as for the actual time step
	  mefel_right_hand_side(lcid, fp, fl);
	  Mp->time = btime;
	  // switch on elements accroding to ACTUAL time step
	  Gtm->update_elem(Mp->time);
	  // switch on nodes accroding to ACTUAL time step     
	  Gtm->update_nodes();
	  
	  // update reaction status indicators at nodes and elements
	  Mt->comreac();
	  // update status indicators of prescribed displacements on elements
	  Mt->elemprescdisp();
	  // prescribe stress free state on the integration points of new parts
	  Mt->clean_ip_new_elem();
	}
      
      //  approximation of temperatures from TRFEL nodes to MEFEL integration points
      approximation_temper ();
      //  approximation of humidities from TRFEL nodes to MEFEL integration points
      approximation_humid ();
      
      //  assembling of the right hand side (prescribed load)
      mefel_right_hand_side (lcid,f,fl);
      
      //  computation of sums of force components in particular directions
      //  it is used in special creep and consolidation models
      Mb->comp_sum (f);
      Mb->comp_sum_react();
      Mb->comp_sum_pdreact();
      fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
              i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      fflush(Out);
      
      if ((mnce - mnae) > 0) 
	{
	  // there were some removed elements =>
	  // calculate nodal forces due to removed elements
	  
	  // use integration point stresses that have been already calculated
	  Mp->strcomp = 0;
	  internal_forces(lcid,fi);
	  Mp->strcomp = 1;
	  
	  //  fb[j] = fp[j] - fi[j];
	  subv(fp, fi, fb, nm);
	}
      else
	nullv(fb, nm);
      
      
      //  computation of unbalanced forces which will be on the right hand side
      incr_internal_forces (lcid,fi);
      
      //  right hand side assembling
      for (j=0;j<nm;j++){
        fb[j]=f[j]-fp[j]+fi[j];
      }    
      
      //  assembling of tangent stiffness matrix
      if ((Mp->nlman->stmat==tangent_stiff) || (ii == 1) || (Smat==NULL))
        stiffness_matrix (lcid);
      
      //  solution of K(r).v=F
      Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);

      //  update total displacement vector  r[j]+=dr[j]
      addv(r, dr, nm);
      
      //  computation of internal forces
      internal_forces (lcid,fi);
      
      // norm of total load vector      
      norfa=Psolm->pss (f,f,Out);
      
      // vector of unbalanced forces fb[j] = fl[j] - fi[j]
      subv(fl, fi, fb, nm);
      
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
        // iteration of unbalanced forces caused by time independent models
        // internal iteration loop
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
            
	  //  vector of unbalanced forces fb[k] = fl[k] - fi[k]
	  subv(fl, fi, fb, nm);
	  
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
              if ((Myrank==0)&&(Mespr==1))
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

	//  backup of nodal values and time because code numbers may be changed
	Mt->save_nodval(lcid);
        prev_timem = Mp->time;

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
  delete [] ifn;

  delete [] p;
  delete [] d;

  delete [] fit;
  delete [] fbt;
  delete [] lhst_last;

  delete [] lhstb;
  delete [] tdlhstb;
}
