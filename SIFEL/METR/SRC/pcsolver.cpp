#include "pcsolver.h"
#include "globalc.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "outdriverm.h"
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
#include "mtsolver.h"
#include "npsolvert.h"
#include "dnnpsolvert.h"

/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem   
   transport problems are linear or nonlinear
   mechanical analysis is always nonlinear
*/
void solve_pcouplprob ()
{

  switch (Cp->tnlinsol){
  case newtonc:{
    //newton_raphson_parcoupl (0);
    //newton_raphson_parcoupl_comp (0);
    newton_raphson_parcoupl_common_dt (0);
    break;
  }
  default:{
    print_err("unknown solver of nonlinear equation system is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   TKo, 27/08/2014
*/
void newton_raphson_parcoupl_common_dt (long lcid)
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
  if(Tp->tprob == discont_nonlin_nonstat_problem)
    nonstat_solver_dform_init(tlcid, np_gv);
  else
    nonstat_solver_init(tlcid, rest_calct, np_gv);
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
  visco_solver_init(mlcid, rest_calcm, mt_gv);

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

    //  TRFEL time
    tnewtime = Tp->time = Tp->timecont.newtime (tdt);
    //  MEFEL time
    mnewtime = Mp->time = Mp->timecon.newtime (mdt);
 

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

    //for debug testing:
    //fprintf (Outt,"\n\n --------------------------------------------------------------");
    //fprintf (Outt,"\n TRFEL Time step = %ld,  Time %e,  Time increment = %e",i,Tp->time,tdt);
    //fprintf (Outt,"\n --------------------------------------------------------------\n");

    ti++;
      
    if((Tp->tprob == nonstationary_problem) || (Tp->tprob == growing_np_problem)){
      // linear algorithm
      // perform one time step with linear solver
      tret = one_step_linear(tlcid, tnewtime, tdt, prev_timet, rest_calct, ti, tli, np_gv);
      trfel_convergence = 1;
    }

    if(Tp->tprob == discont_nonstat_problem){
      print_err("the solver has not yet been implemented", __FILE__, __LINE__, __func__);
      // perform one time step with linear solver
      // tret = one_step_linear(tlcid, tnewtime, tdt, ti, tli, np_gv);
      // trfel_convergence = 1;
    }

    if(Tp->tprob == discont_nonlin_nonstat_problem){
      // perform one time step with nonlinear solver
      tret = one_step_nonlinear_dform(tlcid, tnewtime, tdt, ti, tli, np_gv);
      trfel_convergence = 1;
    }
      
    if((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin)){
      // non-linear algorithm
      // perform one time step with non-linear solver
      tret = one_step_nonlinear(tlcid, tnewtime, tdt, prev_timet, rest_calct, ti, tli, np_gv);
      trfel_convergence = 1;
    }
      
    //zde debug??!!
    //if(Tp->time >= 1.029600e+06){
    //compute_req_valt (lcid);
    //print_stept_forced(lcid, ti, Tp->time, np_gv.rhs);
    //print_flusht();
    //
    //compute_req_val (mlcid);
    //print_step_forced(mlcid, mi, Mp->time, mt_gv.fb);
    //print_flush();
    //exit(0);
    //}

    //  handling with time controler for nonlinear solver
    if ((((Tp->tprob == nonlinear_nonstationary_problem) || (Tp->tprob == growing_np_problem_nonlin)) || 
         (Tp->tprob == discont_nonlin_nonstat_problem))  && Tp->timecont.tct > 0)
    {
      if (tret >= 0){
        // equilibrium has been attained

        //  approximation of nodal values into integration points
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
	    
            if (Mesprt==1)  
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
	  
        if (Mesprt==1)  
          fprintf (stdout,"\n\n TRFEL time increment is reduced to dt=%le because the inner loop was not able to enforce equilibrium", tdt);
        
        if (tdt<Tp->timecont.dtmin){
          if (Mesprt==1)  fprintf (stderr," TRFEL time increment is less than the minimum time increment");
          if (Mesprt==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
          if (Mesprt==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");

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
    mret = one_step (mlcid, mnewtime, mdt, mdtr, prev_timem, rest_calcm, mi, mli, mt_gv);
      
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
          if (Mespr==1)  
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

      if (Mespr==1)  
        fprintf (stdout,"\n\n MEFEL time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (mdt<Mp->timecon.dtmin){
        if (Mespr==1)  fprintf (stderr, " MEFEL time increment is less than the minimum time increment");
        if (Mespr==1)  fprintf (stderr, "\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
        if (Mespr==1)  fprintf (stderr, "\n FORCED output of results from this step is performed\n");
       
 
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
      if (Mesprt==1)
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
      if (Mespr==1)
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



/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   TKr, 06/11/2012
   this function is necessary to correct???!!
*/
void newton_raphson_parcoupl_comp (long lcid)
{
  long mefel_convergence,trfel_convergence,stop;
  double dt,end_time;
  long i,k,nt,ti,tli,tnsts, btnsts,tret;
  double tnewtime,tdt, atdt, tdtmin;
  double *lhsb,*tdlhsb;
  np_glob_vec np_gv;
  long tlcid = 0;//zatim??!!
  long mi, mli, mnsts, bmnsts, mret;
  double mnewtime, mdt, mdtr, mdtmin;
  mt_glob_vec mt_gv;
  long mlcid = 0;//zatim??!!
  long rest_calc,rest_calct;
  double prev_timem,prev_timet;

  //
  //  METR initialization phase
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


  //  initial approximation of quatities from TRFEL nodes to MEFEL integration points
  init_trfel_mefel ();
  approximation();
  pass_coup_data (lcid);

  //
  //  TRFEL initialization phase
  //
  if(Tp->tprob == discont_nonlin_nonstat_problem)
    nonstat_solver_dform_init(tlcid,np_gv);
  else
    nonstat_solver_init(tlcid,rest_calct,np_gv);
  //  initial time increment
  atdt = tdt = Tp->timecont.initialtimeincr ();
  //  minimum time increment
  tdtmin=Tp->timecont.dtmin;
  //  maximum time increment
  //tdtmax=Tp->timecont.dtmax;
  //  end time
  //tend_time = Tp->timecont.endtime ();
  //  number of step
  tli = ti = np_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  tnsts = btnsts = 0;
  // return status from the one_step function
  tret = 0;
  
  //
  //  MEFEL initialization phase
  //
  visco_solver_init(mlcid, rest_calc, mt_gv);
  //  initial time increment
  mdt = Mp->timecon.initialtimeincr ();
  //  minimum time increment used in MEFEL inner loop
  mdtmin=Mp->timecon.dtmin;
  //  maximum time increment
  //mdtmax=Mp->timecon.dtmax;
  //  end time
  //mend_time = Mp->timecon.endtime ();
  //  number of step
  mli = mi = mt_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  mnsts = bmnsts = 0;
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
  else
  {
    lhsb   = NULL;
    tdlhsb = NULL;
  }
  
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
  
  if (Mp->tprob == growing_mech_structure)  
    //saving values for the first iteration
    Mt->save_nodval(lcid);

  do{
    i++;
    
    //  TRFEL time
    tnewtime = Tp->time = Tp->timecont.newtime (tdt);
    //  MEFEL time
    mnewtime = Mp->time = Mp->timecon.newtime (mdt);
    
    if (Tp->time > Mp->time)
      Cp->time = Mp->time;
    else
      Cp->time = Tp->time;
    
    
    //  time synchronization
    if (Tp->time > Mp->time){
      Tp->timecont.oldtime ();
      tnewtime = Tp->time = Tp->timecont.newtime (mdt);
      atdt = mdt; // actual TRFEL dt for linear solvers (due to recovery of decreased time step in MEFEL)
      if (btnsts == 1)
        tnsts = btnsts;
    }
    else
    {
      atdt = tdt; // actual TRFEL dt for linear solvers (due to recovery of decreased time step in MEFEL)
    }
    
    // store initial value for attained time in the previous time step
    // (prev_time is exploited in the growing mechanical problems)
    prev_timem = Mp->time;
    prev_timet = Tp->time;

    //  METR is copying data
    pass_coup_data(lcid);

    //  TRFEL solution
    trfel_convergence = 0;
    do{
      
      ti++;

      if ((Tp->tprob == growing_np_problem) && (prev_timet == Tp->timecont.starttime()))
        prev_timet = Tp->time;
      if ((Tp->tprob == growing_np_problem_nonlin) && (prev_timet == Tp->timecont.starttime()))
        prev_timet = Tp->time;      

      if(Tp->tprob == nonstationary_problem){
	// linear algorithm
	// perform one time step with linear solver
	//tret = one_step_linear(tlcid, tnewtime, atdt, ti, tli, np_gv);
	trfel_convergence = 1;
      }

      if(Tp->tprob == discont_nonstat_problem){
        print_err("the solver has not yet been implemented", __FILE__, __LINE__, __func__);
	// perform one time step with linear solver
	// tret = one_step_linear(tlcid, tnewtime, atdt, ti, tli, np_gv);
	// trfel_convergence = 1;
      }

      if(Tp->tprob == discont_nonlin_nonstat_problem){
	// perform one time step with linear solver
	//tret = one_step_nonlinear_dform(tlcid, tnewtime, atdt, rest_calct, ti, tli, np_gv);
	trfel_convergence = 1;
      }
      
      if(Tp->tprob == nonlinear_nonstationary_problem) 
	// non-linear algorithm
	// perform one time step with non-linear solver
	//tret = one_step_nonlinear(tlcid, tnewtime, atdt, rest_calct, ti, tli, np_gv);
      
      //  handling with time controler for nonlinear solver
      if (((Tp->tprob == nonlinear_nonstationary_problem) || 
           (Tp->tprob == discont_nonlin_nonstat_problem))  && Tp->timecont.tct > 0)
      {
	if (tret >= 0){
	  // equilibrium has been attained
	  
	  trfel_convergence = 1;
	  
          if(Tp->timecont.tct > 0){
            for (k=0;k<Ndoft;k++){
              lhsb[k] = np_gv.lhs[k]; //storing results from previous inner step
              tdlhsb[k] = np_gv.tdlhs[k]; //storing results from previous inner step
            }
          }
	  
	  if (tret == 0)
	    tnsts++;
          else
            tnsts = btnsts = 0;      

	  if (tnsts==2){
	    atdt = tdt*=2.0;
	    tnsts=0;
	    
	    if (Mesprt==1)  
	      fprintf (stdout,"\n\n TRFEL time increment is enlarged because no inner loop was neccessary in previous 2 steps");
	  }
	}
	else{
	  // equilibrium has not been attained
           tnsts = btnsts = 0;      

	  //  reduction of the time increment because
	  //  inner loop was not able to enforce equilibrium
	  trfel_convergence = 0;
	  
	  atdt = tdt/=2.0;
	  ti--;
	  
	  if (Mesprt==1)  
	    fprintf (stdout,"\n\n TRFEL time increment is reduced to dt=%le because the inner loop was not able to enforce equilibrium", tdt);
	  
	  if (tdt<tdtmin){
	    if (Mesprt==1)
	      print_err(" TRFEL time increment is less than the minimum time increment", __FILE__, __LINE__, __func__);
	    
	    compute_req_valt (lcid);
	    print_stept(lcid, ti, Tp->time, np_gv.rhs);
	    print_flusht();
	    
	    stop=1;
	    break;
	  }
	 
	  // change status of elements and nodes according to the previous time step
	  if ((Tp->tprob == growing_np_problem) || (Tp->tprob == growing_np_problem_nonlin))
	  {
	    update_elemnod_prev_timet(prev_timet, np_gv.ncd, np_gv.nce);//????????????????tady?????!!!!!
	  }
 
	  Tp->timecont.oldtime ();
 	  //  TRFEL time
	  tnewtime = Tp->time = Tp->timecont.newtime (tdt);
	}
      }
    }while (trfel_convergence == 0);
    
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
    
    
    //  MEFEL solution
    mefel_convergence = 0;
    do{
      
      //  new mefel step number
      mi++;
      
      //  copying of data
      pass_coup_data(lcid);
      
      if ((Mp->tprob == growing_mech_structure) && (prev_timem == Mp->timecon.starttime()))
        prev_timem = Mp->time;
      //  perform one time step
      mret = one_step (mlcid, mnewtime, mdt, mdtr, prev_timem, rest_calc, mi, mli, mt_gv);
      
      if (mret >= 0){
	//  equilibrium has been attained
	
	if (Mp->time >= Tp->time){
	  mefel_convergence = 1;
	}
	
	if ((mret == 0) || (mret < Mp->nlman->niilnr/4))
	  mnsts++;      
        else
          mnsts = bmnsts = 0;
	
	if (mnsts==2){
	  mdt*=2.0;
	  mnsts=0;
          bmnsts = 1;
	  if (Mespr==1)  
	    fprintf (stdout,"\n\n MEFEL time increment is enlarged because no inner loop was neccessary in previous 2 steps");
	}
	
	if (mefel_convergence == 0){
	  //  MEFEL time
	  if (Tp->time-Mp->time < mdt)
	    mdt = Tp->time-Mp->time;
	  mnewtime = Mp->time = Mp->timecon.newtime (mdt);
	}
	
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
	mi--;
	
    // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
    Mm->updateother();

	if (Mespr==1)  
	  fprintf (stdout,"\n\n MEFEL time increment is reduced because the inner loop was not able to enforce equilibrium");
	if (mdt<mdtmin){
	  if (Mespr==1)
	    print_err(" MEFEL time increment is less than the minimum time increment", __FILE__, __LINE__, __func__);
	  
	  compute_req_val (mlcid);
	  print_step_forced(mlcid, mi, Mp->time, mt_gv.fb);
	  //print_step_forced(mlcid, mi, Mb->dlc[0].gf[0].getval(Mp->time), mt_gv.fb);
	  print_flush();
	  
	  stop=1;
	  break;
	}
        // change status of elements and nodes according to the previous time step
        if (Mp->tprob == growing_mech_structure)
        {
          update_elemnod_prev_time(prev_timem, mt_gv.ncd, mt_gv.nce);
        }
	
	//  MEFEL time
	mnewtime = Mp->time = Mp->timecon.newtime (mdt);
	
      }
    }while(mefel_convergence == 0);
    
    
    if (stop==1){
      //  MEFEL requires smaller time step than the minimum time step
      //  the code will be terminated
      break;
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
  }while (Cp->time < end_time); 
  
  print_close ();
  print_closet();
  
  delete [] lhsb;
  delete [] tdlhsb;

}



/**
   function solves partially coupled thermo-hydro-mechanical time-dependent problem
   by the Newton-Raphson method
   
   @param lcid - load case id
   
   11.6.2008
*/
void newton_raphson_parcoupl (long lcid)
{
  switch (Tp->tprob){
  case nonstationary_problem:{
    //  linear nonstationary transport problem
    //  nonlinear mechanical problem
    newton_raphson_parcoupl_lin (lcid);
    break;
  }
  case nonlinear_nonstationary_problem:{
    //  nonlinear nonstationary transport problem
    //  nonlinear mechanical problem
    newton_raphson_parcoupl_nonlin (lcid);
    break;
  }
  default:{
    print_err("unknown TRFEL problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   Function solves system of linear TRFEL algebraic and non-linear MEFEL equations by Newton-Raphson method
   for time-dependent problems

   @param lcid - load case id
   
   15.7.2005 JK, corrected by TKr 5.4.2007
   12/06/2012 TKr updated
*/
void newton_raphson_parcoupl_lin (long lcid)
{
  int mefel_convergence;
  long i,ii,j,k,nm,nt,ini,nsts;
  double zero,dt,dtmin,dtdef,end_time,alpha,s;
  double *dr,*fb,*r,*d,*f,*fp,*fi,*fl,*p,*lhst,*tdlhst,*rhst,*lhsb,*lhstb,*tdlhstb;
  double norf, norfa, err;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  Mp->nlman->divc_step = new double[10];
  Mp->nlman->divc_err  = new double[10];

  //  maximum number of iterations in inner loop
  ini = Mp->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->errnr;
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
    //  initiation of mechanical material models
    //  approximation of temperatures and humidities from TRFEL nodes to MEFEL integration points
    init_trfel_mefel ();
    Mm->initmaterialmodels(lcid, true);
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file for MEFEL\n");
    solver_restore (r,fp,ii,Mp->time,dt,&Mp->timecon,Ndofm);
    Mm->updateother();
    print_init(-1, "at");
    //    print_step(lcid,i,Mp->time,f);
    //    print_flush();
  }
  else{
    // initiation of mechanical material models
    //  approximation of temperatures and humidities from TRFEL nodes to MEFEL integration points
    init_trfel_mefel ();
    Mm->initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    if (Mp->eigstrains > 0)
      internal_forces(lcid, fp);
    print_init(-1, "wt");
    print_step(lcid, i, Mp->time, f);
    print_flush();
  }
  if (Tp->hdbcont.restore_stat()){
    //  initiation of transport material models
    Tm->initmaterialmodels();
    if (Mesprt==1)
      fprintf (stdout,"\n Reading of backup file for TRFEL\n");
    solvert_restore (lhst,tdlhst,f,i,Cp->time,dt,Cp->timecon,Ndoft);
    compute_req_valt (0);
    print_initt(-1, "at");
    //    print_stept(0,i,Tp->time,NULL);
    //    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt");
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

    //copying of data
    pass_coup_data(lcid);

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
    Tp->ssle->solve_system (Gtt,Kmat,tdlhst,rhst,Outt);

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

      //copying of data
      pass_coup_data(lcid);

      //  backup of attained nodal values
      for (j=0;j<nm;j++){
	lhsb[j]=r[j];
      }

      if (Mespr==1)  fprintf (stdout,"\n\n -------------------------------------------------------------------------------");
      if (Mespr==1)  fprintf (stdout,"\n MEFEL TIME STEP = %ld, MEFEL TIME = %e, mefel time increment = %e",ii,Mp->time,dt);
      if (Mespr==1)  fprintf (stdout,"\n -------------------------------------------------------------------------------");

      //  assembling of vectors of prescribed nodal forces
      mefel_right_hand_side (lcid,f,fl);

      //  computation of sums of force components in particular directions
      //  it is used in special creep and consolidation models
      Mb->comp_sum (f);
      //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
      //              i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      //fflush(Out);
            
      //  computation of unbalanced forces which will be on the right hand side
      incr_internal_forces (lcid,fi);
      
      //  right hand side assembling
      for (j=0;j<nm;j++){
	fb[j]=f[j]-fp[j]+fi[j];
      }    

      /* norf = ss(f,f,nm);
	 fprintf(stdout, "\n\n norm(f)=%le,", norf);
	 norf = ss(fp,fp,nm);
	 fprintf(stdout, " norm(fp)=%le,", norf);
	 norf = ss(fi,fi,nm);
	 fprintf(stdout, " norm(fi)=%le,", norf);
	 norf = ss(fb,fb,nm);
	 fprintf(stdout, " norm(f-fp+fi)=%le\n\n", norf);
      */

      //  assembling of tangent stiffness matrix
      if ((Mp->nlman->stmat==tangent_stiff) || (ii == 1) || (Smat==NULL))
	stiffness_matrix (lcid);
      
      //  solution of K(r).v=F
      Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

      //  new displacement vector
      for (j=0;j<nm;j++){
	r[j]+=dr[j];
      }
      
      //  computation of internal forces
      internal_forces (lcid,fi);
      
      // norm of total load vector      
      norfa=ss(f,f,nm);
      // vector of unbalanced forces
      for (j=0;j<nm;j++){
	fb[j] = fl[j] - fi[j];
      }
      //  norm of vector of unbalanced forces
      norf=ss(fb,fb,nm);
      if (norfa > 0.0)
	norf /= norfa;
      
      if (Mespr==1)  fprintf (stdout,"\n\n Norf before inner iteration loop, norf = %e\n\n\n",norf);	    
      j = 0;
      if (norf<err){
        //  number of successful time steps
        nsts++;
      }
      else
	{
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
	      Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);
	      
	      // new displacement vector
	      for (k=0;k<nm;k++)
		r[k]+=dr[k];
	      
	      //  computation of internal forces
	      internal_forces (lcid,fi);
	      
	      //  vector of unbalanced forces
	      for (k=0;k<nm;k++)
		fb[k]= fl[k] - fi[k];
	      
	      //  norm of vector of unbalanced forces
	      norf=ss(fb,fb,nm);
	      if (norfa > 0.0)
		norf /= norfa;
	      
	      //          if ((Mespr==1)&&(j%10==0)) // print each 10-th step
	      if (Mespr==1)
		fprintf (stdout,"\n Inner loop number j = %ld,   norf = %e \n",j,norf);
	      
	      if (norf<err)  break; // convergence attained
	      
	      // divergence detection with help of least square method
	      if (j > 200)
		{
		  if (lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), Mp->nlman->divc_step[j%10], 
			       Mp->nlman->divc_err[j%10], norf, Mp->zero,1) > 0.0)
		    {
		      if (Mespr==1)
			fprintf (stdout,"\n\n Divergence of inner iteation loop is detected. Inner loop is stopped.\n\n");
		      break;
		    }
		  Mp->nlman->divc_step[j%10] = double(j+1);
		  Mp->nlman->divc_err[j%10] = err;
		}
	      else
		{
		  lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero,0);
		  Mp->nlman->divc_step[j%10] = double(j+1);
		  Mp->nlman->divc_err[j%10] = err;
		}
	    }
	}
      
      if(Cp->fcsolv == fullnewtonc){
	//cleaning of stiffness matrix
	if (Smat != NULL){
	  delete Smat;
	  Smat=NULL;
	}
      }

      if (j==ini || norf>err)
	{
	  // equilibrium state was not attained
	  if (Mespr==1) fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);
	  
	  mefel_convergence = 0;
	  
	  //  backup is used for actual MEFEL values
	  for (j=0;j<nm;j++)
	    r[j]=lhsb[j];
	  
	  if (Mespr==1)  fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
	  if (dt<dtmin){
	    if (Mespr==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
	    if (Mespr==1)  fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
	    compute_req_val (lcid);
	    print_step_forced(lcid, i, Mp->time, f);
	    print_flush();
	    break;
	  }
	}
      else
	{
	  // equilibrium state was attained
	  if (Mespr==1)  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
	  mefel_convergence = 1;
	  
	  //  actual total load vector becomes previous total load vector
	  for (j=0;j<nm;j++)
	    fp[j]=f[j];
	  
	  dtdef = Mp->timecon.actualforwtimeincr ();
	  if (nsts==2){
	    dt*=2.0;
	    nsts=0;
	    if (Mespr==1)  fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	  }
	  if (dt>dtdef)
	    dt=dtdef;
	  
	  // print resulting vector of internal forces
	  Mb->comp_sum (fi);
	  //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
	  //i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
	  //fflush(Out);
	  
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
      if ((Tp->timecont.isitimptime()==1) && (Tp->hdbcont.hdbtype))
	{
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
  
  delete [] lhstb;
  delete [] tdlhstb;
}



/**
   Function solves system of non-linear TRFEL algebraic and non-linear MEFEL equations by Newton-Raphson method
   for time-dependent problems.

   @param lcid - load case id
   
   TKr 5.4.2007, pozn.: 31/08/2012 TKr stare aktualizovat??!!!
*/
void newton_raphson_parcoupl_nonlin (long lcid)
{
  int mefel_convergence;
  long i,ii,j,k,nm,nt,ini,init,nsts;
  double zero,dt,dtmin,dtdef,end_time,alpha;
  double *dr,*fb,*r,*d,*f,*fp,*fi,*fl,*p,*lhst,*tdlhst,*rhst,*lhsb,*lhstb,*tdlhstb;
  double *fbt,*fit,*lhst_last;
  //double s;
  double norf_last;
  double norf, norfa,err,errt;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);

  //  maximum number of iterations in inner loop
  ini = Mp->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->errnr;
  //  maximum number of iterations in inner loop
  init = Tp->nii;
  //  required norm of vector of unbalanced forces
  errt = Tp->err;
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
    print_init(-1, "at");
//    print_step(lcid,i,Mp->time,f);
//    print_flush();
  }
  else
  {
    // initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    internal_forces(lcid, fp);
    print_init(-1, "wt");
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
    print_initt(-1, "at");
//    print_stept(0,i,Tp->time,NULL);
//    print_flusht();
  }
  else{
    //  initiation of transport material models
    Tm->initmaterialmodels();
    compute_req_valt (0);
    print_initt(-1, "wt");
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
    Tp->ssle->solve_system (Gtt,Kmat,tdlhst,rhst,Outt);

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
      
      norf = ss (fbt,fbt,nt);
      
      if (Mesprt==1)  fprintf (stdout,"\n TRFEL inner loop number j = %ld,   norf = %e",j,norf);
      
      if (norf<errt) break;
      
      Tp->ssle->solve_system (Gtt,Kmat,d,fbt,Outt);
      
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
	if (norf > norf_last){
	  for (k=0;k<nt;k++){
	    lhst[k] = lhst_last[k]; //storing results from previous inner step
	  }
	  
	  //  physically corrected solution
	  solution_correction ();
	  //  approximation of nodal values into ontegration points
	  approximation ();	  

          if (Mesprt==1)  fprintf (stdout,"\n\nTRFEL convergence control: inner iteration skipped %ld error %e\n\n",j,norf);
	  
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
      //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
      //            i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      //fflush(Out);
      
      //  computation of unbalanced forces which will be on the right hand side
      incr_internal_forces (lcid,fi);
      
      //  right hand side assembling
      for (j=0;j<nm;j++){
	fb[j]=f[j]-fp[j]+fi[j];
      }    
      
      norf = ss(f,f,nm);
      fprintf(stdout, "\n\n norm(f)=%le,", norf);
      norf = ss(fp,fp,nm);
      fprintf(stdout, " norm(fp)=%le,", norf);
      norf = ss(fi,fi,nm);
      fprintf(stdout, " norm(fi)=%le,", norf);
      norf = ss(fb,fb,nm);
      fprintf(stdout, " norm(f-fp+fi)=%le\n\n", norf);
      
      //  assembling of tangent stiffness matrix
      if ((Mp->nlman->stmat==tangent_stiff) || (ii == 1) || (Smat==NULL))
	stiffness_matrix (lcid);
      
      //  solution of K(r).v=F
      Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

      //  new displacement vector
      for (j=0;j<nm;j++){
	r[j]+=dr[j];
      }
      
      //  computation of internal forces
      internal_forces (lcid,fi);

      // norm of total load vector      
      norfa=ss(f,f,nm);
      // vector of unbalanced forces
      for (j=0;j<nm;j++){
	fb[j] = fl[j] - fi[j];
      }
      //  norm of vector of unbalanced forces
      norf=ss(fb,fb,nm);
      if (norfa > 0.0)
	norf /= norfa;
      
      if (Mespr==1)  fprintf (stdout,"\n\n Norf before inner iteration loop norf = %e\n\n\n",norf);
      j = 0;
      if (norf<err){
        nsts++;
      }
      else
      {
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
          Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);
	    
          // new displacement vector
          for (k=0;k<nm;k++)
            r[k]+=dr[k];

          //  computation of internal forces
          internal_forces (lcid,fi);
	    
          //  vector of unbalanced forces
          for (k=0;k<nm;k++)
            fb[k]= fl[k] - fi[k];
	    
          //  norm of vector of unbalanced forces
          norf=ss(fb,fb,nm);
          if (norfa > 0.0)
            norf /= norfa;
	    
//          if ((Mespr==1)&&(j%10==0)) // print each 10-th step
          if (Mespr==1)
            fprintf (stdout,"\n Inner loop j =  %ld     norf = %e\n",j,norf);

          if (norf<err)  break; // convergence attained

          // divergence detection with help of least square method
          if (j > 1)
          {
            if (lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero, 1) > 0.0)
            {
              if (Mespr==1)
                fprintf (stdout,"\n\n Divergence of inner iteation loop is detected. Inner loop is stopped.\n\n");
              break;
            }
          }
          else
            lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero, 0);
        }
      }

      if(Cp->fcsolv == fullnewtonc){
	//cleaning of stiffness matrix
	if (Smat != NULL){
	  delete Smat;
	  Smat=NULL;
	}
      }
      
      if (j==ini || norf>err)
      {
        // equilibrium state was not attained
        if (Mespr==1) fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);

        mefel_convergence = 0;

        //  backup is used for actual MEFEL values
        for (j=0;j<nm;j++)
          r[j]=lhsb[j];

        //  reduction of the time increment because
        //  inner loop was not able to enforce equilibrium
        dt/=2.0;
        Mp->timecon.oldtime ();
        ii--;

        // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
        Mm->updateother();

        if (Mespr==1)  fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
        if (dt<dtmin){
          if (Mespr==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
          if (Mespr==1)  fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
          compute_req_val (lcid);
          print_step_forced(lcid, i, Mp->time, f);
          print_flush();
          break;
        }
      }
      else
      {
        // equilibrium state was attained
        if (Mespr==1)  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
        mefel_convergence = 1;

        //  actual total load vector becomes previous total load vector
        for (j=0;j<nm;j++)
          fp[j]=f[j];

        dtdef = Mp->timecon.actualforwtimeincr ();
        if (nsts==2){
          dt*=2.0;
          nsts=0;
          if (Mespr==1)  fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
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










/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */



/**
   function solves system of nonlinear algebraic
   equations by Newton-Raphson method
   
   @param lcid - load case id
   
   19.5.2004 JK, corrected by TKr 5.4.2007
*/
void newton_raphson_parcoupl_nonlin_old (long lcid)
{
  long i,ii,j,k,nm,nt,ini;
  double zero,dt,end_time,alpha,norfb,err;
  double *r,*dr,*d,*f,*fb,*fp,*fim,*p,*lhs,*tdlhs,*rhs;
  
  char fname[1020];
  FILE *aux;

  filename_decomposition (Outdm->outfn,Mp->path,Mp->filename,Mp->suffix);
  sprintf (fname,"%s%s.bac",Mp->path, Mp->filename);

  aux = fopen (fname,"w");

  //  number of mechanical degrees of freedom
  nm=Ndofm;
  //  number of transport degrees of freedom
  nt=Ndoft;
  
  
  r = new double [nm];
  dr = new double [nm];

  /* *************************** */
  /*  vectors of transport part  */
  /* *************************** */
  //  array containing predictors
  d = new double [nt];
  nullv (d,nt);
  //  auxiliary vector
  p = new double [nt];
  nullv (p,nt);
  //  right hand side
  fb = new double [nt];
  nullv (fb,nt);

  /* **************************** */
  /*  vectors of mechanical part  */
  /* **************************** */
  //  vector of nodal displacements
  r = new double [nm];
  nullv (r,nm);
  //  vector of prescribed forces
  f = new double [nm];
  nullv (f,nm);
  //  vector of prescribed forces from the previous step
  fp = new double [nm];
  nullv (fp,nm);
  //  vector of internal forces
  fim = new double [nm];
  nullv (fim,nm);
  
  
  alpha=Tp->alpha;
  zero=Tp->zero;
  err=Tp->err;
  ini=Tp->nii;
  
  
  //  starting time
  Cp->time=Cp->timecon.starttime ();
  //  time increment
  dt=Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();
  
  //  nodes - integration points interpolation
  approximation ();
  
  // **************************
  //  main iteration loop  ***
  // **************************
    
  Mp->time = Cp->time;
  Mp->timecon.take_values (Mp->timecon);
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);

  i=0;
  ii=0;
  print_init(-1, "wt");
  print_step(lcid,ii,Mp->time,f);
  print_flush();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  print_flusht();
  

  approximation_temper ();
  approximation_humid ();
  
  //  initiation of mechanical material models
  Mm->initmaterialmodels(lcid, false);
  
  // *********************************************************
  //                     transport part
  // *********************************************************
  do{
    
    if (Mesprc==1)  fprintf (stdout,"\n\n -----------------------------------------------------------------------------");
    if (Mesprc==1)  fprintf (stdout,"\n metr time %e,    trfel time %e,   mefel time %e",Cp->time,Tp->time,Mp->time);
    if (Mesprc==1)  fprintf (stdout,"\n ---------------------------------------------------------------------------\n");

    if (Mesprc==1)  fprintf (stdout,"\n\n provadi se TRFEL,  trfel time %e",Tp->time);

    lhs=Lsrst->lhs;
    tdlhs=Lsrst->tdlhs;
    rhs=Lsrst->rhs;
    
    //  conductivity matrix
    conductivity_matrix (lcid);
    
    //  capacity matrix
    capacity_matrix (lcid);
    
    //  predictor   d+(1-alpha)*dt*v
    for (j=0;j<nt;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side of transport part
    trfel_right_hand_side (0,rhs,nt);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<nt;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    //  solution of the system of algebraic equations
    //Kmat->solve_system (Gtt,tdlhs,fb);
    Cp->ssle->solve_system (Gtt,Kmat,tdlhs,fb,Outt);

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<nt;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    

    solution_correction ();
    //  approximation of nodal values into ontegration points
    approximation ();
    
    //  computation of actual fluxes
    //  fb contains unbalanced fluxes
    internal_fluxes (fb,nt);
    
    //  norm of unbalanced fluxes
    norfb = ss (fb,fb,nt);
    
    if (Mesprc==1)  fprintf (stdout,"\n%e %e",norfb,err);    
    
    if (norfb>err){
      
      //  inner iteration loop
      for (j=0;j<ini;j++){
	
	//  correction, the right hand side is created by unbalanced fluxes
	//Kmat->solve_system (Gtt,p,fb);
	Cp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

	for (k=0;k<nt;k++){
	  tdlhs[k]+=p[k];
	  lhs[k]+=alpha*dt*p[k];
	}

	solution_correction ();
	//  approximation of nodal values into ontegration points
	approximation ();
	
	//  computation of actual fluxes
	//  fb contains unbalanced fluxes
	internal_fluxes (fb,nt);
	
	//  norm of unbalanced fluxes
	norfb = ss (fb,fb,nt);
	
	if (Mesprc==1)  fprintf (stdout,"\n iterace %ld   chyba %e",j,norfb);
	
	if (norfb<err){
	  break;
	}
      }
    }
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    
    // *********************************************************
    //  end of transport part and beginning of mechanical part
    // *********************************************************

    if (Mp->time <= Tp->time){

      if (Mesprc==1)  fprintf (stdout,"\n provadi se MEFEL,  mefel time %e",Mp->time);
    
    approximation_temper ();
    approximation_humid ();
    
    
    Mm->updateipval();
    
    //  mechanical part of the problem
    lhs=Lsrs->lhs;
    rhs=Lsrs->rhs;
    
    //  assembling of the right hand side (prescribed loading)
    mefel_right_hand_side (lcid,f);
    
    //  computation of sums of force components in particular directions
    //  it is used in special creep and consolidation models
    Mb->comp_sum (f);
    
    //Mp->phase=1;
    //Mm->est=eigstrain;
    //  computation of unbalanced forces which will be on the right hand side
    internal_forces (lcid,fim);
    
    
    for (j=0;j<nm;j++){
      rhs[j]=f[j]-fp[j]+fim[j];
      fp[j]=f[j];
    }    
    
    //  assembling of tangent stiffness matrix
    stiffness_matrix (lcid);
    
    //  solution of K(r).v=F
    //Smat->solve_system (Gtm,r,rhs);
    Cp->ssle->solve_system (Gtm,Smat,r,rhs,Out);

    //  new displacement vector
    for (j=0;j<nm;j++){
      lhs[j]+=r[j];
    }
    
    //Mp->phase=2;
    //  update of computed values
    internal_forces (lcid,fim);

    compute_req_val (lcid);

    Mp->time = Mp->timecon.newtime ();
    Mp->timecon.take_values (Mp->timecon);

    if (Mp->timecon.isitimptime ()==1){
      solver_save (lhs,fp,ii,Mp->time,dt,&Mp->timecon,nm);
    }
    
    print_step(lcid, ii, Mp->time, f);
    print_flush();
    ii++;
    }

    //  new time and new time increment
    Cp->time = Cp->timecon.newtime ();
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);
    
    dt = Cp->timecon.actualbacktimeincr ();
    i++;
    
  }while(Cp->time<=end_time);
  
  print_close ();
  print_closet ();
  
  delete [] fim;
  delete [] fb;
  delete [] fp;
  delete [] f;
  delete [] p;
  delete [] d;
  delete [] dr;
  delete [] r;

}




/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */
/* -------------------------------------------------------------- */


/**
   function solves system of nonlinear algebraic
   equations by Newton-Raphson method
   
   @param lcid - load case id
   
   19.5.2004 - modified by TKr 31.10.2005
*/
void newton_raphson_parcoupl_nonlin_new (long lcid)
{
  long i,j,k,nm,nt,ini;
  double zero,dt,end_time,alpha,norfb,err;
  double *r,*dr,*d,*f,*fb,*fp,*fim,*fi,*p,*lhs,*tdlhs,*rhs;
  
  char fname[1020];
  FILE *aux;

  filename_decomposition (Outdm->outfn,Mp->path,Mp->filename,Mp->suffix);
  sprintf (fname,"%s%s.bac",Mp->path, Mp->filename);

  aux = fopen (fname,"w");

  //  number of mechanical degrees of freedom
  nm=Ndofm;
  //  number of transport degrees of freedom
  nt=Ndoft;
  
  
  r = new double [nm];
  dr = new double [nm];

  /* *************************** */
  /*  vectors of transport part  */
  /* *************************** */
  //  array containing predictors
  d = new double [nt];
  nullv (d,nt);
  //  auxiliary vector
  p = new double [nt];
  nullv (p,nt);
  //  right hand side
  fb = new double [nt];
  nullv (fb,nt);
  //  vector of internal fluxes
  fi = new double [nt];
  nullv (fi,nt);
  
  /* **************************** */
  /*  vectors of mechanical part  */
  /* **************************** */
  //  vector of nodal displacements
  r = new double [nm];
  nullv (r,nm);
  //  vector of prescribed forces
  f = new double [nm];
  nullv (f,nm);
  //  vector of prescribed forces from the previous step
  fp = new double [nm];
  nullv (fp,nm);
  //  vector if internal forces
  fim = new double [nm];
  nullv (fim,nm);
  
  
  alpha=Tp->alpha;
  zero=Tp->zero;
  err=Tp->err;
  ini=Tp->nii;
  
  
  //  starting time
  Cp->time=Cp->timecon.starttime ();
  //  time increment
  dt=Cp->timecon.initialtimeincr ();
  //  end time
  end_time = Cp->timecon.endtime ();
  
  //  nodes - integration points interpolation
  approximation ();
  
  // **************************
  //  main iteration loop  ***
  // **************************
    
  Mp->time = Cp->time;
  Mp->timecon.take_values (Mp->timecon);
  Tp->time = Cp->time;
  Tp->timecont.take_values (Cp->timecon);

  i=0;
  print_init(-1, "wt");
  print_step(lcid,i,Mp->time,f);
  print_flush();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  print_flusht();
  

  approximation_temper ();
  approximation_humid ();
  
  //  initiation of mechanical material models
  Mm->initmaterialmodels(lcid, false);
  
  // *********************************************************
  //                     transport part
  // *********************************************************
  do{
    
    if (Mesprc==1)  fprintf (stdout,"\n\n -----------------------------------------------------------------------------");
    if (Mesprc==1)  fprintf (stdout,"\n metr time %e,    trfel time %e,   mefel time %e",Cp->time,Tp->time,Mp->time);
    if (Mesprc==1)  fprintf (stdout,"\n ---------------------------------------------------------------------------\n");

    if (Mesprc==1)  fprintf (stdout,"\n\n provadi se TRFEL,  trfel time %e",Tp->time);

    lhs=Lsrst->lhs;
    tdlhs=Lsrst->tdlhs;
    rhs=Lsrst->rhs;
    
    //  conductivity matrix
    conductivity_matrix (lcid);
    
    //  capacity matrix
    capacity_matrix (lcid);
    
    //  predictor   d+(1-alpha)*dt*v
    for (j=0;j<nt;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side of transport part
    trfel_right_hand_side (0,rhs,nt);
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<nt;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    //  solution of the system of algebraic equations
    //Kmat->solve_system (Gtt,tdlhs,fb);
    Cp->ssle->solve_system (Gtt,Kmat,tdlhs,fb,Outt);

    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<nt;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }

        
    //inner iteration loop
    for (j=0;j<ini;j++){
      
      solution_correction ();
      approximation ();
      
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
      
      
      if (Tp->trestype==lrhst){
	Kmat->gmxv (tdlhs,fi);
	//  Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
	for (k=0;k<nt;k++){
	  fb[k] = rhs[k] - p[k] - fi[k];
	}
      }
      if (Tp->trestype==fluxest){
	internal_fluxes (fi,nt);
	//  vector of unbalanced fluxes
	for (k=0;k<nt;k++){
	  fb[k]=fi[k];
	}
      }
      
      norfb = ss (fb,fb,nt);
      
      if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,norfb);

      if (norfb<err){
	break;
      }

      //Kmat->solve_system (Gtt,p,fb);
      Cp->ssle->solve_system (Gtt,Kmat,p,fb,Outt);

      for (k=0;k<nt;k++){
	tdlhs[k]+=p[k];
	lhs[k]+=alpha*dt*p[k];
      }
    }
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
    
    // *********************************************************
    //  end of transport part and beginning of mechanical part
    // *********************************************************

    if (Mp->time <= Tp->time){

      if (Mesprc==1)  fprintf (stdout,"\n provadi se MEFEL,  mefel time %e",Mp->time);
    
    approximation_temper ();
    approximation_humid ();
    
    
    Mm->updateipval();
    
    //  mechanical part of the problem
    lhs=Lsrs->lhs;
    rhs=Lsrs->rhs;
    
    //  assembling of the right hand side (prescribed loading)
    mefel_right_hand_side (lcid,f);
    
    //  computation of sums of force components in particular directions
    //  it is used in special creep and consolidation models
    Mb->comp_sum (f);
    
    //Mp->phase=1;
    //Mm->est=eigstrain;
    //  computation of unbalanced forces which will be on the right hand side
    internal_forces (lcid,fim);
    
    
    for (j=0;j<nm;j++){
      rhs[j]=f[j]-fp[j]+fim[j];
      fp[j]=f[j];
    }    
    
    //  assembling of tangent stiffness matrix
    stiffness_matrix (lcid);
    
    //  solution of K(r).v=F
    //Smat->solve_system (Gtm,r,rhs);
    Cp->ssle->solve_system (Gtm,Smat,r,rhs,Out);

    //  new displacement vector
    for (j=0;j<nm;j++){
      lhs[j]+=r[j];
    }
    
    //Mp->phase=2;
    //  update of computed values
    internal_forces (lcid,fim);

    compute_req_val (lcid);

    Mp->time = Mp->timecon.newtime ();
    Mp->timecon.take_values (Mp->timecon);

    if (Mp->timecon.isitimptime ()==1){
      solver_save (lhs,fp,i,Mp->time,dt,&Mp->timecon,nm);
    }
    
    print_step(lcid, i, Mp->time, f);
    print_flush();
        
    }

    //  new time and new time increment
    Cp->time = Cp->timecon.newtime ();
    Tp->time = Cp->time;
    Tp->timecont.take_values (Cp->timecon);
    
    dt = Cp->timecon.actualbacktimeincr ();
    i++;
    
  }while(Cp->time<=end_time);
  
  print_close ();
  print_closet ();
  
  delete [] fim;
  delete [] fi;
  delete [] fb;
  delete [] fp;
  delete [] f;
  delete [] p;
  delete [] d;
  delete [] dr;
  delete [] r;

}
