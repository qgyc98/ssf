#include <string.h>
#include "cpsolver.h"
#include "mtsolver.h"
#include "newtonraph.h"
#include "nssolver.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "lhsrhs.h"
#include "mechbclc.h"
#include "globmat.h"
#include "mechprint.h"
#include "matrix.h"
#include "vector.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "dloadcase.h"
#include "node.h"
#include "element.h"
#include "slipsurf.h"



/**
  The function calls solver of time dependent problems of growing structures for each load case.

  @return The function does not return anything.

  Created by TKr,
*/
void solve_prob_constr_phases ()
{
  long i=0;
  for (i=0; i<Lsrs->nlc; i++)
    visco_solver2 (i);
  //solve_prob_constr_phases(i);

  if (Mm->mcoul)
  {
    fprintf(stdout, "\n Mohr-Coulomb material detected -> starting detection of slip surfaces\n");
    long *plast_ip = new long[Mm->tnip];
    long nplz = detect_plastic_zones(plast_ip, 1.0e-3);
    fprintf(stdout, "\n There are detected %ld of plastic zones\n", nplz);
    for(i=0; i<nplz; i++)
    {
      long err;
      double a, b, minx, maxx;
      err = approx_slip_surf(i+1, plast_ip, a, b, minx, maxx);
      if (err)
        fprintf(stdout, "  - Approximation of the %ld slip surface could not be calculated\n", i+1);
      else
        fprintf(stdout, "  - Approximation of the %ld slip surface - a=%le, b=%le, minx=%le, maxx=%le\n", i+1, a, b, minx, maxx);
    }
    delete [] plast_ip;
  }
  
}



/**
  Function solves non-linear MEFEL equations by Newton-Raphson method for time-dependent problems
  for growing structure
   
  @param lcid - load case id

  @return The function does not return anything.
   
  Created by TKr, 15.6.2007
  Modified by TKo,
*/
void solve_prob_constr_phases (long lcid)
{
  long i,j,n,ini,nsts,mncd,mnce,mnae,li;
  long nifn, um;
  answertype fusm;
  long rest_calc=0;  
  long *ifn;
  double dt,dtmin,dtdef,end_time,prev_timem;
  double *f,*fl,*fi,*fb,*fp,*flp,*r,*dr,*lhsb;
  double norf, norfa, err, zero;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  nonlinman cpnlman;
  
  // array with marked nodes on interface between new and old part of structure
  ifn = new long[Mt->nn];
  memset(ifn, 0, sizeof(*ifn)*Mt->nn);
  Mp->nlman->divc_step = new double[10];
  Mp->nlman->divc_err  = new double[10];
  //  initialization phase

  //  maximum number of iterations in inner loop
  ini = Mp->nlman->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  //  number of mechanical degrees of freedom
  // at the begining, it represents all possible dofs at nodes 
  // because all elements and all nodes are switched on
  n=Ndofm;
  //  computer zero
  zero = Mp->zero;

  // zero norm of unbalanced force vector because initial state is in equilibrium
  // and norf is used for proper handling of output graphics GiD files
  norf = 0.0;

  //  vector of nodal displacements
  r = Lsrs->give_lhs (lcid);
  //  vector of nodal prescribed forces
  f = Lsrs->give_rhs (lcid);

  // allocation of vectors  

  //  vector of increments of nodal displacements
  dr   = new double [n];
  nullv (dr,n);
  //  vector of prescribed forces from the previous step including temperature forces, etc.
  fp   = new double [n];
  nullv (fp,n);
  //  vector of prescribed forces from the previous step due to prescribed force load
  flp   = new double [n];
  nullv (flp,n);
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
  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  minimum time increment
  dtmin=Mp->timecon.dtmin;
  //  backup of time
  prev_timem = Mp->time;

  //  number of step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  //  initial printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    if (Mespr==1)
      fprintf (stdout,"\n Reading of backup file\n");
    solver_restore (r,fp,i,Mp->time,dt,&Mp->timecon,n);
    li = i;
    //  backup of time
    prev_timem = Mp->time;

    //  Backup of r is used in case that the Newton-Raphson procedure does not converge.
    //  lhsb[j]=r[j]
    copyv(r, lhsb, n);

    //  marking of elements switched on
    Gtm->update_elem (Mp->time);
    //  at Mp->time=prev_timem must be auxinf=leso
    Gtm->update_auxinf();
    //  marking of nodes switched on
    Gtm->update_nodes ();
    //  marking of DOFs switched on
    Gtm->update_dofs (Mp->time);
    // new DOF numbering
    n = Ndofm = Gtm->codenum_generation(Out);
    // save nodal values 
    Mt->save_nodval(lcid);
    // update reaction status indicators at nodes and elements
    Mt->comreac();
    // update status indicators of prescribed displacements on elements
    Mt->elemprescdisp();
    // update actual material index for all elements
    Mm->update_actual_mat_id(lcid, 0, true);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    
    rest_calc = 1;
    print_init(-1, "at");
    //print_step(lcid,i,Mp->time,f);
    //print_flush();
  }
  else
  {
    // compute stresses from eigenstrains or initial irreversible strains
    //internal_forces(lcid, fp);
    li = 0;
    print_init(-1, "wt");
    print_step(lcid, i, Mp->time, f);
    print_flush();
  }
  //print_close();

  //saving values for the first iteration
  Mt->save_nodval(lcid);

  // ***************************
  //   main iteration loop  ****
  // ***************************
  do{
    //new time increment
    Mp->time = Mp->timecon.newtime (dt);
    if (prev_timem == Mp->timecon.starttime())
      prev_timem = Mp->time;
    //  new step number
    i++;
    Mp->istep = i;
    //  indicator of strain computation
    //  it means, strains at integration points have not been computed
    Mp->strainstate=0;
    //  indicator of stress computation
    //  it means, stresses at integration points have not been computed
    Mp->stressstate=0;
    //  indicator of computation of other array
    //  it means, stresses at integration points have not been computed
    Mp->otherstate=0;

    if (Mespr==1)  fprintf (stdout,"\n\n -------------------------------------------------------------");
    if (Mespr==1)  fprintf (stdout,"\n Time step = %ld,  Time %e,  Time increment = %e",i,Mp->time,dt);
    if (Mespr==1)  fprintf (stdout,"\n -------------------------------------------------------------\n");

    // update actual material index for all old elements and call initval for those whose indeces have changed
    if (rest_calc == 0)
    {
      um = Mm->update_actual_mat_id(lcid, 1, false);
      if (um)   fusm = yes;
      else      fusm = no;
    }  
    else  // do not update if the calculation was restored from HD
    {
      rest_calc = 0;
      fusm = no;
    }

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

      // actualize displacement vctors of axis definition nodes for active axes      
      if (Mp->rot_inidispl == yes)
        Mb->actualize_displ_def_node(lcid, Mp->time);

      // switch on only new elements
      Gtm->switch_new_elem();
        
      // switch on nodes on the new part of structure
      Gtm->update_nodes();
        
      // prepare active nodes (including nodes on interface)
      // for generation of code numbers
      Gtm->update_active_dofs(Mp->time);
        
      // new DOF numbering including interface nodes
      Ndofm = Gtm->codenum_generation(Out);
        
      // restore nodal displacements for interface nodes 
      // (zero displacements will be in the rest of nodes)
      Mt->restore_nodval(r, Ndofm);

      // update actual material index for all new elements
      // but do not initialize new active material models
      Mm->update_actual_mat_id(lcid, 0, false);

      //  initiation of mechanical material models
      Mm->initmaterialmodels(lcid, false);
      
          
      //
      // CALCULATION OF INITIAL DISPLACEMENTS FOR NEW PARTS OF STRUCTURE
      //

      if (nifn) // if there are some interfaces between old and new part of structure
      {
        if (Mp->comp_inidispl) // initial displacements will be calculated for all new elements
        {
          // store prescribed initial displacements into the vector of unknowns directly
          // in order to determine the forces due to these prescribed values
          Mb->store_inidispl(ifn, r);
          if (Mp->rot_inidispl == yes)
            Mb->store_rotinidispl(ifn, Mp->time, r);

          // calculate stresses due to displacements on the interface nodes
          stress_initdispl(lcid);
          
          // prepare active nodes (including nodes on interface)
          // for generation of code numbers
          Gtm->update_active_dofs(Mp->time);
          // EXCLUDE INTERFACE NODES
          Gtm->clear_intf_dofs(ifn);
          // EXCLUDE DOFs at nodes with prescribed initial displacements 
          Mb->clear_inidispl(ifn);
          // EXCLUDE DOFs at nodes with prescribed initial displacements due to rotation about given axis
          if (Mp->rot_inidispl == yes)
            Mb->clear_rotinidispl(ifn, Mp->time);
          
          // new DOF numbering
          Ndofm = Gtm->codenum_generation(Out);

          //
          // Apply prescribed initial displacements from the initial conditions at new active nodes 
          //
          // Detect the number of prescribed initial displacements
          Mb->dlc[lcid].npid = Mb->num_dofs_inidispl(ifn);
          if (Mp->rot_inidispl == yes)
            // Detect the number of prescribed initial displacements by rotations at the current time
            Mb->dlc[lcid].nrpid = Mb->num_dofs_rotinidispl(ifn, Mp->time);
          // Total number of prescribed initial displacements
          Mb->dlc[lcid].tnpid = Mb->dlc[lcid].npid + Mb->dlc[lcid].nrpid;
          // Detect the maximum of prescribed displacements in the given dloadcase
          Mb->dlc[lcid].compute_max_npd();
          Mb->dlc[lcid].realloc_pid_array();

          if (Mp->rot_inidispl == yes)
            Mb->apply_inidispl(lcid, ifn);
          // Apply prescribed initial displacements due to rotation about given axis
          if (Mp->rot_inidispl == yes)
            Mb->apply_rotinidispl(lcid, ifn, Mp->time);
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
          fprintf(Out, "\n\n\nAfter generation of code numbers for new part at time: %le\n", Mp->time);
          Gtm->codenum_print(Out);
          Mp->ssle->solve_system (Gtm,Smat,r,fb,Out);
          
          // update vector of total displacements
          //addv(r, dr, n);
          
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
        else 
        {
          // initial displacements will be taken only form nodes on the interface 
          // between new and old parts, all new nodes have zero displacements or 
          // prescribed values
              
          // store initial displacements on the new part of structure
          // in order to prescribe stress free state on required new elements
          Mt->save_elem_inidispl(lcid, ifn);
        }
      }

      // switch on elements accroding to ACTUAL time step
      Gtm->update_elem(Mp->time);

      // switch on nodes accroding to ACTUAL time step     
      Gtm->update_nodes();

      //  marking of DOFs switched on
      Gtm->update_dofs (Mp->time);
      
      // new DOF numbering
      n = Ndofm = Gtm->codenum_generation(Out);

      // destruct old system matrix
      delete Smat;
      Smat=NULL;

      // restore vector of attained unknowns from the previous time step + eventually calculated 
      // displacements on the new parts
      Mt->restore_nodval(r,n);
 
      //
      // Restore vector fp of attained load vector from the previous time step
      //
      if (i > 1)
        Mt->restore_nodforce(fp, n);

      // update reaction status indicators at nodes and elements
      Mt->comreac();
      // update status indicators of prescribed displacements on elements
      Mt->elemprescdisp();
      // prescribe stress free state on the integration points of new parts
      Mt->clean_ip_new_elem();

      /* prepared for groups of results
        if (norf<err) // equilibrium was attained in the previous time step
        {
          // close groups of results and mesh from previous time step
          print_close();
          // export new mesh and add new group of results 
          print_init(-1, "at");
        }
      */
    }

    //  assembling of vectors of prescribed nodal forces
    mefel_right_hand_side (lcid,f,fl);

    //  computation of sums of force components in particular directions
    //  it is used in special creep and consolidation models

    Mb->comp_sum (f);
    Mb->comp_sum_react();
    Mb->comp_sum_pdreact();
    if (Mt->give_ndofn(0) > 2)
      fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
              i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
    else
      fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le", 
              i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1]);
    fflush(Out);

    if ((mnce - mnae) > 0) 
    {
      // There were some removed elements =>
      // calculate nodal forces due to removed elements:
      // f^{re}_i - f^{re}_lp      

      // Stored total load vector fp from the previous time step contains unwanted contributions from the removed elements
      // and therefore it must be reassembled again on actual active elements only.
      double aux = Mp->time;
      Mp->time = prev_timem;           
      mefel_right_hand_side (lcid,fp,flp);

      // switch on removed elements in the actual time step
      // the previous time may be left in Mp->time because the gtopology procedure does not depend on it
      Gtm->switch_removed_elem();

      // Calculate contribution f^{re}_lp which contains contributions due to force load from removed elements.
      // The result is stored in vector fb while the vector fi contains the total load vector 
      // (including eigenstrains and eigenstresses contributions, temperatures, etc.). The vector fi is not used 
      // in further computations and therefor will be rewritten in the next statements.
      Mp->time = prev_timem;           
      mefel_right_hand_side (lcid,fi,fb);  // fb = f^{re}_lp

      // Calculate contribution f^{re}_i with the help of integration point stresses that have been already calculated
      // in the previous time step.
      Mp->strcomp = 0;
      internal_forces(lcid,fi);
      Mp->strcomp = 1;

      // fb = f^{re}_i - f^{re}_{lp}
      subv(fi, fb, fb, n);

      // restore the actual time and update status of active elements
      Mp->time = aux;
      Gtm->update_elem(Mp->time);
      if (Mp->cpsmooth == yes)
      {
        // load vector which is equivalent to vector of internal forces on the remaining part of element domain
        // fc = f_{lp} + f^{re}_{lp} - f^{re}_i
        subv(flp, fb, n);
        cpnlman.tnlinsol = newtong;
        // maximum number of steps in the equilibrium iteration in the case of contruction phase change
        cpnlman.ninr = long(1.0/Mp->cpminincrnr)+1;
        // maximum number of iterations in inner loop in the equilibrium iteration in the case of contruction phase change
        cpnlman.niilnr = Mp->nlman->niilnr;
        //  required norm of vector of unbalanced forces in the equilibrium iteration in the case of contruction phase change
        cpnlman.errnr = Mp->nlman->errnr;
        // default increment in the equilibrium iteration in the case of contruction phase change
        cpnlman.incrnr = Mp->cpincrnr;
        // minimum increment in the equilibrium iteration in the case of contruction phase change
        cpnlman.minincrnr = Mp->cpminincrnr;
        // maximum increment in the equilibrium iteration in the case of contruction phase change
        cpnlman.maxincrnr = 1.0;
        // there will be defined required value of the load coefficient to be attained
        cpnlman.check_lambdar = on;
        // required value of the load coefficient to be attained
        cpnlman.lambdar = 1.0;
        // load coefficient starts form zero
        Mp->lambda = 0.0;
        // fi vector will be used for the storage of attained load vector fa in the Newton-Raphson (NR) procedure temporarily
        nullv(fi, n);
        /* The vector flp represents load vector on the remaining part of the domain from the previous time step which becomes 
           constant load vector fc in NR. The vector fb represents load increment of traction due to element removal on the 
           interface nodes between remaining and removed domain parts and thus it becomes proportional load vector flp in NR.
           No change of internal forces due to eigenstrain/stresses or temperatures will be assumed in course of smooth element 
           removal and therefor fc = flc and fp = flp in NR. */
        if (Mespr==1) fprintf(stdout, "\n\n -----------------------------------------------------------------------\n");
        if (Mespr==1) fprintf(stdout, " Smoothed element removal has been started at time %le\n\n", Mp->time);
        aux = gnewton_raphson2(lcid, &cpnlman, r, fi, flp, fb, flp, fb, 0, 0.0, no);
        // if the required load coefficient value lambdar = 1.0 has not been attained => error, break the computation
        if (aux < 1.0)
        {
          print_err("Equilibrium cannot be established after smoothed element removal in time step %le.\n"
                    " Prescribe lower minimum load increment cpminincrnr (actual cpminincrnr=%le)\n"
                    " or remove less elements.", __FILE__, __LINE__, __func__, Mp->time, Mp->cpminincrnr);
          if (Mespr==1)  fprintf (stderr,"\n Computation failed (file %s, line %d)\n",__FILE__,__LINE__);
          if (Mespr==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
          compute_req_val (lcid);
          print_step_forced(lcid, i, Mp->time+0.1, fb);
          print_flush();
          break;
        }
        else
        {
          if (Mespr==1) fprintf(stdout, "\n\n Smoothed element removal has been finished successfully\n");
          if (Mespr==1) fprintf(stdout, " -----------------------------------------------------------------------\n");
          nullv(fb, n);
        }
      }
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
    assemble_stiffness_matrix(lcid, i, -1, li, fusm);

    //  solution of K(r).v=F
    Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

    //  update total displacement vector  r[j]+=dr[j]
    addv(r, dr, n);
        
    //  computation of internal forces
    internal_forces (lcid,fi);

    // norm of load vector      
    norfa=normv(fl,n);

    // vector of unbalanced forces fb[j] = fl[j] - fi[j]
    subv(fl, fi, fb, n);

    /*
    if (Mp->time > 1.0)
    {
      compute_req_val (lcid);
      print_step_forced(lcid, i, Mp->time+(j+1)*0.001, fb);
      print_flush();
      }*/

    //  norm of vector of unbalanced forces
    norf=normv(fb,n);
    if (norfa > 0.0)
      norf /= norfa;

    if (Mespr==1)  fprintf (stdout,"\n\n Norf before inner iteration loop norf = %e\n\n\n",norf);
    j = -1;
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
        Mp->jstep = j;

        // re-assemble stifness matrix if it is required
        assemble_stiffness_matrix(lcid,i,j+1,li,no);

        //  solution of K(r).v=F
        Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

        //  new displacement vector r[k]+=dr[k]
        addv(r, dr, n);

        //  computation of internal forces
        internal_forces (lcid,fi);

        //  vector of unbalanced forces fb[k] = fl[k] - fi[k]
        subv(fl, fi, fb, n);
        //  norm of vector of unbalanced forces
        norf=normv(fb,n);
        if (norfa > 0.0)
          norf /= norfa;

//        if ((Mespr==1)&&(j%10==0)) // print each 10-th step
	if (Mespr==1)
  	  fprintf (stdout,"\n Inner loop number j = %ld,     norf = %e \n",j,norf);

        if (norf<err)  break; // convergence attained

        // divergence detection with help of least square method
        /*
        if (j > 10)
        {
          if (lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), Mp->nlman->divc_step[j%10], 
                       Mp->nlman->divc_err[j%10], norf, Mp->zero,1) > 0.0)
	  {
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
        */
      }
    }

    if (j==ini || norf>err)
    {
      // equilibrium state was not attained
      if (Mespr==1) fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);


      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      if (Mespr==1)  fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin)
      {
        if (Mespr==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
        if (Mespr==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
        if (Mespr==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
        compute_req_val (lcid);
        print_step_forced(lcid, i, Mp->time+0.1, fb);
        print_flush();
        break;
      }

      Mp->timecon.oldtime ();
      i--;
        
      //    
      // regenerate code numbers according to old time step if there were changes in nodes or elements
      //
      if ((mnce!=0) || (mncd!=0))
      {
        // switch on elements accroding to previous time step
        Gtm->update_elem(prev_timem);
        // switch on nodes accroding to previous time step     
        Gtm->update_nodes();
        //  marking of DOFs switched on
        Gtm->update_dofs (prev_timem);
        // old DOF numbering
        n = Ndofm = Gtm->codenum_generation(Out);
        // destruct old system matrix
        delete Smat;
        Smat=NULL;
      }

      //  backup is used for attained values
      for (j=0;j<n;j++)
        r[j]=lhsb[j];

    }
    else
    {
      // equilibrium state was attained
      if ((Mespr==1) && (j >= 0))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);

      //  Backup of r is used in case that the Newton-Raphson procedure does not converge.
      //  lhsb[j]=r[j]
      copyv(r, lhsb, n);

      //  actual total load vector becomes previous total load vector
      //  fp[j]=f[j]
      copyv(f, fp, n);

      // Following two commands must be performed here because in case that
      // elements or nodes are changed in the next step and the step would not converge,
      // the nodal values and nodal forces could not be recovered correctly
      // nodval due to changes in Mt->nodedispl and nodforce due to changes in vector f

      // backup of nodal values because code numbers may be changed due to
      Mt->save_nodval(lcid);
      //  backup of nodal forces because code numbers may be changed
      Mt->save_nodforce(f);

      // new elements can be definitely changed to old ones
      Gtm->update_auxinf();
      prev_timem = Mp->time;

      dtdef = Mp->timecon.actualforwtimeincr ();
      if (nsts==2)
      {
        dt*=2.0;
        nsts=0;
        if (Mespr==1)  fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
      }
      if (dt>dtdef)
        dt=dtdef;

      // print resulting vector of internal forces      
      //Mb->comp_sum (fi);
      //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
      //      i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      //fflush(Out);

      //  update of values stored at integration points
      Mm->updateipval();
      compute_req_val (lcid);
      print_step(lcid, i, Mp->time, f);
      print_flush();
	
      Mb->comp_sum_react();
      Mb->comp_sum_pdreact();
      if (Mt->give_ndofn(0) > 2)
      {
        fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], rx=%le, ry=%le, rz=%le", 
                i, Mp->time, Mp->time/86400, Mb->reactsumcomp[0], Mb->reactsumcomp[1], Mb->reactsumcomp[2]);
        fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], pdx=%le, pdy=%le, pdz=%le", 
                i, Mp->time, Mp->time/86400, Mb->pd_reactsumcomp[0], Mb->pd_reactsumcomp[1], Mb->pd_reactsumcomp[2]);
      }
      else
      {
        fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], rx=%le, ry=%le", 
                i, Mp->time, Mp->time/86400, Mb->reactsumcomp[0], Mb->reactsumcomp[1]);
        fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], pdx=%le, pdy=%le", 
                i, Mp->time, Mp->time/86400, Mb->pd_reactsumcomp[0], Mb->pd_reactsumcomp[1]);
      }
      fflush(Out);

      if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.hdbtype))
      {
        if (Mespr==1)
	  fprintf (stdout,"\n Creating backup file\n");
        solver_save (r,fp,i,Mp->time,dt,&Mp->timecon,n);
      }
    }
  }while(Mp->time<end_time);
  
  print_close ();
  
  delete [] fi;  
  delete [] fb;  
  delete [] fp;  
  delete [] dr; 
  delete [] fl;  
  delete [] flp;
  delete [] lhsb;
  delete [] ifn;
}
