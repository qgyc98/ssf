#include "mpi.h"
#include <string.h>
#include "pcpsolver.h"
#include "pglobal.h"
#include "seqfilesm.h"
#include "genfile.h"
#include <string.h>


void par_solve_prob_constr_phases ()
{
  long i=0;
  for (i=0; i<Lsrs->nlc; i++)
    par_solve_prob_constr_phases(i);
}


/**
   Function solves non-linear MEFEL equations by Newton-Raphson method for time-dependent problems
   for growing structure
   
   @param lcid - load case id
   
   15.6.2007, TKr
*/
void par_solve_prob_constr_phases (long lcid)
{
  long i,j,k,n,mncd,mnce,ini,nsts,nii,mnae;
  long nifn;
  long *ifn;
  double dt,dtmin,prev_timem,end_time,btime;
  double *dr,*fb,*r,*f,*fp,*fl,*fi,*lhsb;
  double norf, norfa, err, dtdef,lsm_res;  
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  
  // array with marked nodes on interface between new and old part of structure
  ifn = new long[Mt->nn];
  memset(ifn, 0, sizeof(*ifn)*Mt->nn);
  //  initialization phase

  //  maximum number of iterations in inner loop
  ini = Mp->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->errnr;
  //  number of mechanical degrees of freedom
  n=Ndofm;
  //  minimum time increment
  dtmin=Mp->timecon.dtmin;

  //  vector of nodal displacements
  r = Lsrs->give_lhs (lcid);
  //  vector of nodal prescribed forces
  f = Lsrs->give_rhs (lcid);

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
  Mp->time = Mp->timecon.starttime ();
  //  backup of time
  prev_timem = Mp->time;
  //  end time
  end_time = Mp->timecon.endtime ();
  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  //  minimum time increment
  dtmin=Mp->timecon.dtmin;

  // **************************
  //  main iteration loop  ***
  // **************************

  //  number of step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file\n");
    solver_restore (r,fp,i,Mp->time,dt,&Mp->timecon,n);
    Mt->rhs_save (fp);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    print_init(-1, "at",Pmp->fni,Pmp->fei);
    //print_step(lcid,i,Mp->time,f);
    //print_flush();
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

  //saving values for the first iteration
  //  Mt->lhs_save (r,1);

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
            
    //  searching for changed elements
    //  status of all elements is updated
    mnce = Gtm->search_changed_elem(Mp->time, mnae);    
    //  searching for new DOFs and
    //  status of nodes is updated
    mncd = Gtm->search_changed_dofs(Mp->time,prev_timem);

    nifn = 0;
    if (mnce!=0){
      //  attained displacements have to be saved on new elements
      //  otherwise stresses will not be correct
      //      Mt->initial_displ ();
    }
      
    //  searching for new DOFs
    //  only new DOFs are switched on
    //  former DOFs are temporarily switched off
    mncd = Gtm->search_newdofs (Mp->time,prev_timem);
      
    //  marking of elements switched on
    Gtm->update_elem (Mp->time);
      
    //  marking of nodes switched on
    Gtm->update_nodes ();
      
    //  marking of DOFs switched on
    Gtm->update_dofs (Mp->time);
      
    //  definition of DOFs
    Ndofm = Gtm->codenum_generation (Out);
    n=Ndofm;
      
    //  update of nodes and DOFs
    if (mncd>0){
      if (Smat != NULL){

        delete Smat;
        Smat=NULL;
        stiffness_matrix (lcid);
        nullv(r, Ndofm);
        
        // solution of system of equations
        //Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);
	//  solution of K(r).v=F
	Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);
        
        // update vector of total displacements
        addv(r, dr, n);
      
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
      n = Ndofm = Gtm->codenum_generation(Out);

      // destruct old system matrix
      delete Smat;
      Smat=NULL;

      // restore vector of attained unknowns
      Mt->restore_nodval(r,n);
 
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
    
    //  assembling of the right hand side (prescribed loading)
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
      subv(fp, fi, fb, n);
    }
    else
      nullv(fb, n);

    //  computation of unbalanced forces which will be on the right hand side
    incr_internal_forces (lcid,fi);

    //  right hand side assembling
    for (j=0;j<n;j++){
      fb[j]+=f[j]-fp[j]+fi[j];
    }
    
    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == 1) || (Smat==NULL))
      stiffness_matrix (lcid);
      
    //  solution of K(r).v=F
    Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);

    //  update total displacement vector  r[j]+=dr[j]
    addv(r, dr, n);
    
    //  computation of internal forces
    internal_forces (lcid,fi);

    // norm of total load vector      
    norfa=Psolm->pss (f,f,Out);    

    // vector of unbalanced forces fb[j] = fl[j] - fi[j]
    subv(fl, fi, fb, n);

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
        Mp->jstep = j;
        //assembling of stiffness matrix for fullnewton
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
        subv(fl, fi, fb, n);
	
        //  norm of vector of unbalanced forces
	norf=Psolm->pss(fb,fb,Out);	
        if (norfa > 0.0)
          norf /= norfa;
	    
//	if ((Myrank==0) && (Mespr==1) && (j%10==0)) // print each 10-th step
	if ((Myrank==0) && (Mespr==1))
  	  fprintf (stdout,"\n Inner loop number j = %ld,     norf = %e \n",j,norf);

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

    if (nii==0)
    {
      // equilibrium state was not attained
      if ((Myrank==0) && (Mespr==1))
        fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);

      //  backup is used for actual values
      for (j=0;j<n;j++)
        r[j]=lhsb[j];
      	    	    
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
    else
    {
      // equilibrium state was attained
      if ((Myrank==0) && (Mespr==1))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);

      //  actual total load vector becomes previous total load vector
      for (j=0;j<n;j++)
        fp[j]=f[j];
            
      //  backup of nodal values and time because code numbers may be changed
      //      Mt->lhs_save (r,1);
      Mt->rhs_save (f);

      prev_timem = Mp->time;
	    
      dtdef = Mp->timecon.actualforwtimeincr ();
      if (nsts==2){
        dt*=2.0;
        nsts=0;
	if ((Myrank==0) && (Mespr==1))
          fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
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
      
      Mb->comp_sum_react();
      Mb->comp_sum_pdreact();
      
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
