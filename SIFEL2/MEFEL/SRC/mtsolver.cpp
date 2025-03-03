#include "mtsolver.h"
#include "backupsol.h"
#include "mtglvec.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "mechprint.h"
#include "matrix.h"
#include "vector.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "dloadcase.h"
#include "elemswitch.h"
#include "nssolver.h"
#include "newtonraph.h"
#include "mathem.h"
// for evaluation of hypoplasticity ???!!!
#include "intpoints.h"
//#include "hypoplunsatexptherm.h"
#include "hypoplunsatexptherm2.h"
//#include "hypoplunsatexptherm2.h"
#include "slipsurf.h"
#include <string.h>
#include <string>

using std::string;
//using std::to_string;


#ifndef FNAMELEN
 #define FNAMELEN 1001
#endif



/**
  The function calls solver of time dependent problems (viscoplasticity)
  for each load case.

  @return The function does not return anything.

  Created by JK,
*/
void solve_time_dep_prob ()
{
  long lcid;
  
  //  load case must be equal to zero, no superposition can be used
  //  in this type of analysis
  lcid=0;
  
  //  solver of the problem
  visco_solver2 (lcid);
}



/**
  Function solves time dependent problem (viscoplasticity)
   
  @param lcid - load case id
   
  @return The function does not return anything.
   
  Created by JK, 27.10.2001
  Modified by Tomas Koudelka
*/
void visco_solver (long lcid)
{
  long i,j,k,n,ini,nsts;
  double dt,dtmin,dtdef,end_time;
  double *f,*fl,*fi,*fb,*fp,*r,*dr,*lhsb;
  double norf, norfa, err, zero;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  Mp->nlman->divc_step = new double[10];
  Mp->nlman->divc_err  = new double[10];
  
  //  initialization phase

  //  maximum number of iterations in inner loop
  ini = Mp->nlman->niilnr;
  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;
  //  number of mechanical degrees of freedom
  n=Ndofm;
  //  computer zero
  zero = Mp->zero;

  //  vector of nodal displacements
  r = Lsrs->give_lhs (lcid);
  nullv (r,n);
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
  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  Mp->dtime = dt;
  //  end time
  end_time = Mp->timecon.endtime ();
  //  minimum time increment
  dtmin=Mp->timecon.dtmin;

  //  number of step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  //  inital printing of output and graphical informations
  if (Mp->hdbcont.restore_stat()){
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, true);
    if (Mespr==1)
      fprintf (stdout,"\n Reading of MEFEL backup file\n");
    solver_restore (r,fp,i,Mp->time,dt,&Mp->timecon,n);
    Mm->updateother();
    print_init(-1, "at");
    //print_step(lcid,i,Mp->time,f);
    //print_flush();
  }
  else
  {
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    print_init(-1, "wt");
    print_step(lcid, i, Mp->time, f, fb);
    //print_step(lcid, i, 0.0, f);
    print_flush();
  }
  //print_close();

  // ***************************
  //   main iteration loop  ****
  // ***************************
  do{
    //new time increment
    Mp->time = Mp->timecon.newtime (dt);
    Mp->dtime = dt;
    //  new step number
    i++;
    Mp->istep = i;
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

    //  backup of attained nodal values
    for (j=0;j<n;j++){
      lhsb[j]=r[j];
    }

    if (Mespr==1)  fprintf (stdout,"\n\n -------------------------------------------------------------");
    if (Mespr==1)  fprintf (stdout,"\n MEFEL Time step = %ld,  Time %e,  Time increment = %e",i,Mp->time,dt);
    if (Mespr==1)  fprintf (stdout,"\n -------------------------------------------------------------\n");
    
    //  assembling of vectors of prescribed nodal forces
    mefel_right_hand_side (lcid,f,fl);
    
    //  computation of sums of force components in particular directions
    //  it is used in special creep and consolidation models
    Mb->comp_sum (f);
    //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
    //      i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
    //fflush(Out);

    //  increments of prescribed forces
    incr_internal_forces (lcid,fi);

    //  right hand side assembling
    for (j=0;j<n;j++){
      fb[j]=f[j]-fp[j]+fi[j];
    }
    
/*    norf = ss(f,f,nm);
    fprintf(stdout, "\n\n norm(f)=%le,", norf);
    norf = ss(fp,fp,nm);
    fprintf(stdout, " norm(fp)=%le,", norf);
    norf = ss(fi,fi,nm);
    fprintf(stdout, " norm(fi)=%le,", norf);
    norf = ss(fb,fb,nm);
    fprintf(stdout, " norm(f-fp+fi)=%le\n\n", norf); */

    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == 1) || (Smat==NULL))
      stiffness_matrix (lcid);

    //  solution of K(r).v=F
    Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

    //  new displacement vector
    for (j=0;j<n;j++){
      r[j]+=dr[j];
    }

    //  computation of internal forces
    internal_forces (lcid,fi);

    // norm of load vector      
    norfa=ss(fl,fl,n);
    // vector of unbalanced forces
    for (j=0;j<n;j++){
      fb[j] = fl[j] - fi[j];
    }

    //  norm of vector of unbalanced forces
    norf=ss(fb,fb,n);
    if (norfa > 0.0)
      norf /= norfa;

    if (Mespr==1)  fprintf (stdout,"\n\n Norf before inner iteration loop norf = %e\n\n\n",norf);
    j = 0;
    if (norf<err){
      nsts++;      
    }
    else{
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
        Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

        //  new displacement vector
        for (k=0;k<n;k++)
	  r[k]+=dr[k];

        //  computation of internal forces
        internal_forces (lcid,fi);

        //  vector of unbalanced forces
        for (k=0;k<n;k++)
          fb[k]=fl[k] - fi[k];

        //  norm of vector of unbalanced forces
        norf=ss(fb,fb,n);
        if (norfa > 0.0)
          norf /= norfa;

//        if ((Mespr==1)&&(j%10==0)) // print each 10-th step
	if (Mespr==1)
  	  fprintf (stdout,"\n Inner loop number j = %ld,     norf = %e \n",j,norf);

        if (norf<err)  break; // convergence attained

        // divergence detection with help of least square method
        if (j > 10)
        {
          if (lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), Mp->nlman->divc_step[j%10], 
                       Mp->nlman->divc_err[j%10], norf, zero,1) > 0.0)
	  {
	    fprintf (stdout,"\n\n Divergence of inner iteation loop is detected. Inner loop is stopped.\n\n");
	    break;
	  }
          Mp->nlman->divc_step[j%10] = double(j+1);
          Mp->nlman->divc_err[j%10] = err;
	}
        else
        {
          lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, zero,0);
          Mp->nlman->divc_step[j%10] = double(j+1);
          Mp->nlman->divc_err[j%10] = err;
        }
      }
    }

    if (j==ini || norf>err)
    {
      // equilibrium state was not attained
      if (Mespr==1) fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);

      nsts = 0;

      //  backup is used for actual values
      for (j=0;j<n;j++)
        r[j]=lhsb[j];

      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      Mp->timecon.oldtime ();
      i--;

      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();

      if (Mespr==1)  fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin)
      {
        if (Mespr==1)  fprintf (stderr,"\n\n time increment is less than minimum time increment");
        if (Mespr==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
        if (Mespr==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
        //  update of values stored at integration points due to print of unequilibriated eqother values
        Mm->updateipval();
        compute_req_val (lcid);
	print_step_forced(lcid, i, Mp->time, fl, fb);
        //print_step_forced(lcid, i, Mb->dlc[0].gf[0].getval(Mp->time), fb);
        print_flush();
        break;
      }
    }
    else
    {
      // equilibrium state was attained
      if (Mespr==1)  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);

      //  actual total load vector becomes previous total load vector
      for (j=0;j<n;j++)
        fp[j]=f[j];

      dtdef = Mp->timecon.actualforwtimeincr ();
      if (nsts==2)
      {
        dt*=2.0;
        nsts=0;
        if (Mespr==1)  fprintf (stdout,"\n\n time increment is enlarged to %le because no inner loop was neccessary in previous 3 steps", dt);
      }
//      if (dt>dtdef)
//        dt=dtdef;

      // print resulting vector of internal forces
      Mb->comp_sum (fi);
      //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fix=%le, fiy=%le, fiz=%le", 
      //      i, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
      //fflush(Out);

      //  update of values stored at integration points
      Mm->updateipval();
      compute_req_val (lcid);
      print_step(lcid, i, Mp->time, fi, fb);
      //print_step(lcid, i, Mb->dlc[0].gf[0].getval(Mp->time), fl);
      print_flush();
	
      if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.save_stat()))
      {
        if (Mespr==1)
	  fprintf (stdout,"\n Creating MEFEL backup file\n");
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



/**
  Function solves time dependent problem (viscoplasticity).
  It uses one_step concept instead of standard one.
  
  @param lcid - load case id
   
  @return The function does not return anything.
  
  Created by Tomas Koudelka, 11.2011, corrected by TKr 27/03/2013
*/
void visco_solver2 (long lcid)
{
  long i, li, nsts, ret, rest_calc;
  double time, dt, dtmin, dtmax, dtr, end_time, prev_time;
  mt_glob_vec mt_gv;
  
  //
  //  initialization phase
  //
  visco_solver_init(lcid, rest_calc, mt_gv);
  
  // store initial value for attained time in the previous time step
  // (prev_time is exploited in the growing mechanical problems)
  prev_time = Mp->time;

  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  //  end time
  end_time = Mp->timecon.endtime ();
  
  //  number of step
  li = i = mt_gv.istep;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;
  
  /*
  // test MuPiF
  vector tt(Mt->nn);
  fillv(10.0, tt);
  intpointval(tt.a, temperature, 1.0);
  fillv(0.0, tt);
  intpointval(tt.a, initial_temperature, 1.0);
  
  ret = one_step(lcid, 1.0, dt, i, mt_gv);
  */

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
      
    // perform one time step
    time=Mp->timecon.newtime(dt);
    if (time<0.0)
      break;

    //  minimum time increment
    dtmin=Mp->timecon.dtmin;
    //  maximum time increment
    dtmax=Mp->timecon.dtmax;
    Mp->dtime = dt;
    if ((Mp->tprob == growing_mech_structure) && (prev_time == Mp->timecon.starttime()))
      prev_time = time;

    ret = one_step(lcid, time, dt, dtr, prev_time, rest_calc, i, li, mt_gv);
    rest_calc = 0;
      
    if (ret >= 0) // equilibrium was attained
    {
 
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
      print_step(lcid, Mp->istep, Mp->time, mt_gv.fl, mt_gv.fb);
      print_flush();
	
      if ((Mp->timecon.isitimptime()==1) && (Mp->hdbcont.save_stat()))
      {
        if (Mespr==1)
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
          if (Mespr==1)  
            fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 2 steps");
        }
      }
    }
    else
    {
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      nsts = 0;

      if (dt > dtmin){
        // time step can be reduced
        if (dtr < 1.0)
          dt *= dtr;  // use reduction factor proposed by material models if they proposed ever reduction
        else
          dt/=2.0;   // use general reduction factor

        // check the minumum length of time step
        if (dt < dtmin)
          dt = dtmin;

        Mp->timecon.oldtime ();
        i--;
      
        // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
        Mm->updateother();

        if (Mespr==1) {
          if (dt == dtmin)
            fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium\n");
          else
            fprintf (stdout,"\n\n time increment is reduced because the inner loop was not able to enforce equilibrium");
        }
        // change status of elements and nodes according to the previous time step
        if (Mp->tprob == growing_mech_structure){
          update_elemnod_prev_time(prev_time, mt_gv.ncd, mt_gv.nce);
        }
      }
      else{
        // time step should be reduced but it was attained the minimum time step length
        if (Mespr==1)  fprintf (stderr,"\n\n reduced time increment is less than minimum time increment");
        if (Mespr==1)  fprintf (stderr,"\n computation failed (file %s, line %d)\n",__FILE__,__LINE__);
        if (Mespr==1)  fprintf (stderr,"\n FORCED output of results from this step is performed\n");
        //  update of values stored at integration points
        Mm->updateipval();
        compute_req_val (lcid);
        print_step_forced(lcid, i, Mp->time, mt_gv.fl, mt_gv.fb);
        print_flush();
        break;
      }
    }
  }while(Mp->time<end_time);

  print_close ();
  /*
  // experimental detection of slip surface
  if (Mm->mcoul)
  {
    fprintf(stdout, "\n Mohr-Coulomb material detected -> starting detection of slip surfaces\n");
    long *plast_ip = new long[Mm->tnip];
    long nplz = detect_plastic_zones(plast_ip, 5.0e-2);
    fprintf(stdout, "\n There are detected %ld of plastic zones\n", nplz);
    for(i=0; i<nplz; i++)
    {
      long err, npip = 0;
      double a, b, minx, maxx;
      fprintf(Out, "\nPlastic zone: %ld\nElements:\n", i+1);
      for (long j=0; j<Mm->tnip; j++){
        if (plast_ip[j] == i+1){
          npip++;
          fprintf(Out, "%ld\n", Mm->elip[j]+1);
        }
      }
      err = approx_slip_surf(i+1, plast_ip, a, b, minx, maxx);
      if (err)
        fprintf(stdout, "  - Approximation of the %ld slip surface could not be calculated\n", i+1);
      else
        fprintf(stdout, "  - Approximation of the %ld slip surface - a=%le, b=%le, minx=%le, maxx=%le, npip=%ld\n", i+1, a, b, minx, maxx, npip);
    }
    delete [] plast_ip;
    }*/
}



/**
  Function allocates and initializes global vectors used in mtsolver
  It is used in one_step concept solver.

  @param lcid      - load case id
  @param rest_calc - indicator of calculation restorage from backup 
                     (0=no restorage was required, 1=restorage was performed)
  @param mt_gv     - structure with global vectors used in the solver

  Created by Tomas Koudelka, 11.2011
*/
void visco_solver_init(long lcid, long &rest_calc, mt_glob_vec &mt_gv)
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
      fprintf (stdout,"\n Reading of MEFEL backup file\n");
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
    // reopen output files
    print_init(-1, "at");
  }
  else
  {
    mt_gv.istep=0;
    // initialize actual chain of time switch material models
    Mm->update_actual_mat_id(lcid, 0, false);
    // Switch on only elements active at the initial time,
    // it is needed for the proper initialization of creep models
    if (Mp->tprob == growing_mech_structure)
      Gtm->update_elem(Mp->time);
    //  initiation of mechanical material models
    Mm->initmaterialmodels(lcid, false);
    // initiation of material models on auxiliary integration points
    Mm->aip_initmaterialmodels(lcid, false);
    // compute stresses from eigenstrains or initial irreversible strains
    if ((Mp->eigstrains > 0) && (Mp->eigstrains < 4))
      internal_forces(lcid, mt_gv.fp);
    // switch all elements and nodes on for the export of whole FE mesh
    if (Mp->tprob == growing_mech_structure)
      Gtm->lneso_init();
    // print initial step to the output files
    print_init(-1, "wt");
    if (Mp->tprob == growing_mech_structure)
      Gtm->update_elem(Mp->time);
    print_step(lcid, mt_gv.istep, Mp->time, mt_gv.f, mt_gv.fb);
    //print_step(lcid, i, 0.0, f);
    print_flush();
  }
}



/**
  Function executes one step of time dependent mechanical algorithm
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
  @retval  >0 - if the equilibrium was attained and inner loop WAS performed
  @retval -1 - if the equilibrium was NOT attained
   
  Created by JK, 19.9.2004
  Modified by TKo, 11.2011
*/
long one_step (long lcid,double time, double dt, double &dtr, double prev_time, long rest_calc, long istep, long li, mt_glob_vec &mt_gv)
{
  long j,n,um, mncd, mnce, mnae;
  double *f,*fl,*fi,*fb,*fp,*r,*dr,*lhsb, *flp;
  long *ifn;
  double norf, err;
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

  // update statuses of all nodes, DOFs and elements for growing structures,
  // calculation of initial displacements of new parts of the structure
  if (Mp->tprob == growing_mech_structure)
  {    
    update_elnod_stat(lcid, istep, prev_time, ifn, r, fb, fp, mnce, mncd, mnae);
    mt_gv.ncd = mncd;
    mt_gv.nce = mnce;
    n = Ndofm;
  }

  //  required norm of vector of unbalanced forces
  err = Mp->nlman->errnr;

  if (Mespr==1)  fprintf (stdout,"\n\n ------------------------------------------------------------------------");
  if (Mespr==1)  fprintf (stdout,"\n MEFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Mp->time,dt);
  if (Mespr==1)  fprintf (stdout,"\n ------------------------------------------------------------------------\n");

  //if (Mespr==1)  fprintf (Out,"\n\n ------------------------------------------------------------------------");
  //if (Mespr==1)  fprintf (Out,"\n MEFEL Time step = %ld,  Time %e,  Time increment = %e",istep,Mp->time,dt);
  //if (Mespr==1)  fprintf (Out,"\n ------------------------------------------------------------------------\n");

  //  assembling of vectors of prescribed nodal forces
  mefel_right_hand_side (lcid,f,fl);
  
    
  //  computation of sums of force components in particular directions
  //  it is used in special creep and consolidation models
  Mb->comp_sum (f);
  //fprintf(Out, "\n# i=%ld, time=%le[s], time=%le[d], fx=%le, fy=%le, fz=%le", 
  //      istep, Mp->time, Mp->time/86400, Mb->sumcomp[0], Mb->sumcomp[1], Mb->sumcomp[2]);
  //fflush(Out);


  // assemble force vector due to removed elements and apply it gradually if required
  if (Mp->tprob == growing_mech_structure)
  {    
    forces_due2removed_elem(lcid, mnce-mnae, istep, prev_time, fi, fp, flp, fb, r);
  }

  //  increments of prescribed internal forces
  incr_internal_forces (lcid,fi);

  // assembling of right hand side
  // fb[j]=f[j]-fp[j]+fi[j]
  subv(f, fp, fb, n);
  addv(fb, fi, fb, n);

  // perform one load step with the help of Newton-Raphson procedure
  j = 0;  

  norf = gnewton_raphson_one_step(lcid, Mp->nlman, fl, r, fb, dr, fi, dtr, istep-1, j, li, fusm);

  // just for evaluation of hypoplasticity model ???!!!  
  /* if (Mm->hypoplustherm) //commented by TKr 11/12/2017
     {
     long i, ncompstr, nstatev = Mm->hypoplustherm[0].nstatev;
     for(i=0; i<Mm->tnip; i++)
     {
     ncompstr = Mm->ip[i].ncompstr;
     Neval += Mm->ip[i].eqother[2*ncompstr+3+nstatev+5];
     }
     }
  */

  if (norf>err)
  {
    // equilibrium state was not attained
    if (Mespr==1) fprintf (stdout,"\n\n Iteration process does not converge req. err = %e, reached err = %e", err,norf);
    
    //  backup is used for actual values
    //  r[k]=lhsb[k]
    copyv(lhsb, r, n);
    
    return -1;
  }

  // equilibrium state was attained
  if ((Mespr==1) && (j > 0))  fprintf (stdout,"\n\n Last inner iteration step number j = %ld,     norf = %e \n\n",j,norf);
  
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
  The function updates status of all elements and nodes according 
  to their time functions after restorage from the backup file.
  The function is being used in the growing mechanical problems.

  @param lcid - load case id

  Created by Tomas Koudelka 31.1.2017
*/
void update_elnod_stat_after_hdbrest(long lcid)
{
  //  marking of elements switched on
  Gtm->update_elem(Mp->time);
  //  at Mp->time=prev_timem must be auxinf=leso
  Gtm->update_auxinf();
  //  marking of nodes switched on
  Gtm->update_nodes();
  //  marking of DOFs switched on
  Gtm->update_dofs(Mp->time);
  // new DOF numbering
  Ndofm = Gtm->codenum_generation(Out);
  // save nodal values 
  Mt->save_nodval(lcid);
  // update reaction status indicators at nodes and elements
  Mt->comreac();
  // update status indicators of prescribed displacements on elements
  Mt->elemprescdisp();
}



/**
  Regenerate DOF numbers according to old time step if there were changes
  in nodes or elements.

  @param prev_time - attained time from previous time step
  @param ncd - the number of changed DOFs in the actual time step
  @param nce - the number of changed elements in the actual time step

  Created by Tomas Koudelka, 31.1.2017
*/
void update_elemnod_prev_time(double prev_time, long ncd, long nce)
{
  //    
  // regenerate code numbers according to old time step if there were changes in nodes or elements
  //
  if ((nce!=0) || (ncd!=0))
  {
    // switch on elements accroding to previous time step
    Gtm->update_elem(prev_time);
    // switch on nodes accroding to previous time step     
    Gtm->update_nodes();
    //  marking of DOFs switched on
    Gtm->update_dofs (prev_time);
    // old DOF numbering
    Ndofm = Gtm->codenum_generation(Out);
    // destruct old system matrix
    delete Smat;
    Smat=NULL;
  }
}



/**
  The function updates status of all elements and nodes according 
  to their time functions.

  @param lcid - load case id
  @param istep - id of the actual time step
  @param prev_time - attained time in the previous time step
  @param ifn - array of interface node indicators, ifn[i] = 1 means that i-th node 
               is on the interface between new and current part of the structure.
  @param r   - %vector of attained nodal displacements
  @param fb  - auxiliary load %vector used in the determination of initial displacements
  @param fp  - load %vector from previous time step
  @param mnce - the number of elements with changed status (added/removed)
  @param mncd - the number of DOFs with changed status (added/removed)
  @param mnae - the number of added elements

  Created by Tomas Koudelka, 31.1.2017
*/
void update_elnod_stat(long lcid, long istep, double prev_time, long *ifn, double *r, double *fb, double *fp,
                       long &mnce, long &mncd, long &mnae)
{
  long nifn, n;

  //  searching for changed elements
  //  status of all elements is updated
  mnce = Gtm->search_changed_elem(Mp->time, mnae);    
  //  searching for new DOFs and
  //  status of nodes is updated
  mncd = Gtm->search_changed_dofs(Mp->time, prev_time);

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
    Mt->save_node_inidispl(lcid, Mp->time, prev_time);        

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
    } // end of nifn condition

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
    if (istep > 1)
      Mt->restore_nodforce(fp, n);

    // update reaction status indicators at nodes and elements
    Mt->comreac();
    // update status indicators of prescribed displacements on elements
    Mt->elemprescdisp();
    // prescribe stress free state on the integration points of new parts
    Mt->clean_ip_new_elem();

/*
    string fname = "noddofnum"+to_string(Mp->istep)+".log";
    FILE *log = fopen(fname.c_str(), "wt");
    long tnan = 0;
    for (long i=0; i<Mt->nn; i++){
      if (Gtm->lnso[i] == 1){
        fprintf (log, "%6ld:", i+1);
        long ndofn = Gtm->give_ndofn(i);
        long n = Mt->give_ndofn(i);
        ivector cn;
        Gtm->give_gnode_code_numbers(i, cn);
        for(long j = 0; j<cn.n; j++){
          fprintf(log, " %6ld", cn[j]);
          if ((ndofn < 0) && (j%n == (n-1)))
            fprintf(log, " :");        
        }
        fprintf(log, "\n");
        tnan++;
      }
    }
    fprintf(log, "Total number of active DOFs: %ld\n", Ndofm);
    fprintf(log, "Total number of active nodes: %ld\n", tnan);
    fprintf(log, "List of active nodes:\n");
    for (long i=0; i<Mt->nn; i++){
      if (Gtm->lnso[i] == 1){
        fprintf(log, "%6ld\n", i+1);
      }
    }
    fclose(log);
*/
  
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
}


/**
  The function assembles %vector of internal forces f^{re}_i - f^{re}_lp due to removed elements. 
  Additionally if it is required, it applies the given %vector in several substeps in order to avoid the convergence problem.

  @param lcid - load case id
  @param nrel - the number of removed elements  = mnce - mnae
  @param istep - id of the actual time step
  @param prev_time - attained time from the previous time step
  @param fi   - auxiliary %vector of internal forces
  @param fp   - load %vector from previous time step (including eigenstresses, prescribed displacements, ...)
  @param flp  - load %vector from previous time step due to force load only
  @param fb   - resulting %vector of f^{re}_i - f^{re}_lp (output)
  @param r    - displacement %vector

  @retval 0 - on success
  @retval 1 - smooth application of f^{re}_i - f^{re}_lp in substeps failed

  Created by Tomas Koudelka, 31.1.2017
*/
long forces_due2removed_elem(long lcid, long nrel, long istep, double prev_time, 
                             double *fi, double *fp, double *flp, double *fb, double *r)
{
  nonlinman cpnlman;
  long n = Ndofm;

  if (nrel > 0)
  {
    // There were some removed elements =>
    // calculate nodal forces due to removed elements:
    // f^{re}_i - f^{re}_lp      

    // Stored total load vector fp from the previous time step contains unwanted contributions from the removed elements
    // and therefore it must be reassembled again on actual active elements only.
    double aux = Mp->time;
    Mp->time = prev_time;           
    mefel_right_hand_side (lcid,fp,flp);

    // switch on removed elements in the actual time step
    // the previous time may be left in Mp->time because the gtopology procedure does not depend on it
    Gtm->switch_removed_elem();

    // Calculate contribution f^{re}_lp which contains contributions due to force load from removed elements.
    // The result is stored in vector fb while the vector fi contains the total load vector 
    // (including eigenstrains and eigenstresses contributions, temperatures, etc.). The vector fi is not used 
    // in further computations and therefor will be rewritten in the next statements.
    Mp->time = prev_time;           
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
        print_step_forced(lcid, istep, Mp->time+0.1, fb);
        print_flush();
        return 1;
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

  return 0;
}
