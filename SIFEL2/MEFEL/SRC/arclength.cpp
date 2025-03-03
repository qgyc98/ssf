#include "arclength.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "nssolver.h"
#include "mechprint.h"
#include "mathem.h"
#include "global.h"
#include "globmat.h"
#include "loadcase.h"
#include "newtonraph.h"

#include "mtsolver.h"

#include <stdlib.h>



/**
  The function solves system of nonlinear algebraic equations
  by the arc-length method. This method was modified to reach the given value lambdar 
  of the lambda parameter. lambdar value is reached with required error errl. 
  Only one right hand side %vector is supported with respect to nonlinearity and absence of superposition. 
  Proportional %vector grows until the nial number of steps is reached.

  @param lcid    - load case id
  @param nlman   - pointer to structure conatining setup of the solver
  @param ra      - %vector of attained displacements
  @param fa      - attained load %vector
  @param fc      - constant part of load %vector (it is not influenced by lambda)
  @param fp      - proportional part of load %vector (it is increased by load coefficient lambda up to the nlman->lambdar)
  @param flc -  constant component of load %vector due to forces
  @param flp -  proportional component of load %vector due to forces
  @param li      - initial value of step id (default is 0)
  @param ilambda - initial value of load coefficient (default is 0)
  @param outres  - flag for performing of output of results (if yes -> print_step procedure is called for the each step)

  @return The function returns reached lambda parameter.

  Created by JK,  16.8.2001
  Rewritten by TKo, JK 08.2011
*/
double garclength(long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp, double *flc, double *flp, 
                  long li, double ilambda, answertype outres)
{
  long i,j,n,ni,ini,stop,stopj,modif;
  nonlinsolvertype nlst; // type of nonlinear solver
  double dl;       // actual value of arc length increment
  double dlmin;    // minimum value of arc length increment
  double dlmax;    // maximum value of arc length increment
  double psi;      // proportional coefficient between load vector and displacement vector
  double psi0;     // initial value of proportional coefficient psi
  double lambda;   // actual value of load coefficient
  double blambda;  // backupped value of lambda
  double dlambda;  // load coefficient increment from the inner loop
  double ddlambda; // cummulative value of load coefficent increments ddlambda = \sum(dlambda)
  double ddlambda0;// initial/reference value of ddlambda
  double norf;     // norm of residual vector
  double norfp;    // norm of proportional load vector
  double norv;     // norm of displacement increments
  double norfa;    // norm of attained load vector
  double ierr;     // relative error of residual vector
  double dtr;      // time step size reduction coeffcient required by material models
  double tmp;      // auxiliary variable
  double *r;   // displacements from previous step (backup copy)
  double *ddr; // vector of displacement increments
  double *u;   // vector of corrections for computed displacements v
  double *v;   // vector of computed displacements for actual arclength step
  double *f;   // auxiliary vector of actual load (call of solve_system rewrite it during solution)
  double *fi;  // vector of internal forces
  double totl    = 0.0; // the total(cumulative) length of arc
  long sm_change = 0;   // flag for changed content of the stiffness matrix
  long modif_psi; // flag for automatic modification psi in dependence on the ddlambda
  flagsw check_rv_fp = nlman->check_lambdar;  // flag for checking of attainment of required value for load coefficient lambdar 
  flagsw check_max_tot_al = nlman->check_tot_al; // flag fro checking maximu value of the length of arc
  resnormt normt;     // type of residual vector norm
  matrix lsm_a(3,3);         // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  double nom, denom; // auxiliary variables for the dissipation increment control
  double tau = nlman->tau_ini;  // actual value of required dissipation increment
  
  

  // type of nonlinear solver
  nlst = nlman->tnlinsol;
  //  number of rows of the matrix
  n = Ndofm;
  // initial value of load coefficient
  lambda = ilambda;
  //  maximum number of increments
  ni = nlman->nial;
  //  maximum number of iterations in one increment
  ini = nlman->niilal;
  //  required error in inner loop
  ierr = nlman->erral;
  //  length of arc
  dl = nlman->dlal;
  //  maximum length of arc
  dlmax = nlman->dlmaxal;
  //  minimum length of arc
  dlmin = nlman->dlminal;
  //  displacement-loading driving switch
  psi = nlman->psial;
  //type of residual vector norm
  normt = nlman->rnormtnr;  
  
  //  allocation of auxiliary arrays
  r   = new double[n];
  ddr = new double[n];
  u   = new double[n];
  v   = new double[n];
  f   = new double[n];
  fi  = new double[n];
  
  //  Delta lambda - cumulative increment in one load increment
  //  at the beginning of computation is always zero
  dlambda=0.0;

  // in case that no restorage from the backup was performed (lambda == 0.0)
  if (lambda == 0.0)
  {
    sm_change = assemble_stiffness_matrix (lcid,li,-1,li,no);
    //  backup of the fc because the right hand side will be destroyed during the solution of the system of equations
    copyv(fc, f, n);
    //  solution of K(r).v=F
    Mp->ssle->solve_system(Gtm,Smat,ra,f,Out);
  }

  //  norm of proportionality vector due to forces
  // norfp = normv(flp, n);
  norfp = normv(fp, n);
  if (norfp == 0.0){
    print_err("Norm of proportional load vector is zero, solution of quadratic equation for lambda would failed\n", __FILE__, __LINE__, __func__);
    abort();
  }
  modif=0;
  psi0 = psi;
  if (psi<0.0){
    //  parameter psi will be modified during the computation
    modif_psi=1;
  }else{
    //  parameter psi will be constant during the computation
    modif_psi=0;
  }


  // ***************************
  //  main iteration loop   ****
  // ***************************
  for (i=li;i<ni;i++) 
  { // loop i

    //  indicator of strain computation
    //  it means, strains at integration points have not been computed
    Mp->strainstate=0;
    //  indicator of stress computation
    //  it means, stresses at integration points have not been computed
    Mp->stressstate=0;
    //  indicator of computation of other array
    //  it means, stresses at integration points have not been computed
    Mp->otherstate=0;

    Mp->istep = i;
    Mp->jstep = -1;
    Mp->lambda = lambda;
    Mp->dlambda = 0.0;

    //  backup of left hand side vector rb[j]=ra[j];
    copyv(ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;

    if (Mespr==1)
    {
      fprintf (stdout,"\n increment %ld,   lambda=%e   dl=%e  psi=%e",i,lambda,dl,psi);
//      fprintf (Out,"\n increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
    }
    
    // assembling of the stiffness matrix
    sm_change = assemble_stiffness_matrix (lcid,i,-1,li,no);

    //  backup of the fp because the right hand side will be destroyed during the solution of the system of equations
    copyv(fp, f, n);

    //  solution of K(r).v=F
    Mp->ssle->solve_system(Gtm,Smat,v,f,Out);

    //  square of generalized norm of displacement increments
    norv = displincr(lcid, i, nlman, v, v, n);

    if (nlst == dissip_incr){
      nom = scprd(fp, ra, n) - 2.0*nlman->tau;
      denom = scprd(fp, v, n);
      if (denom == 0.0)
        dlambda = 0.0; // proceed with alternative setup for arclength
      else
        dlambda = nom/denom;
      if (dlambda < nlman->tau_lim) // dissipated energy increment is too small => switch to arc-length method with alternative setup
        dlambda = 0.0;
    }
    if ((nlst != dissip_incr) || (dlambda == 0.0)){
      if ((i == 0) && (psi < 0.0) && (norfp > Mp->zero))
        psi0 = psi = sqrt(norv/norfp/norfp);

      //  compute new dlambda increment
      tmp = norv+psi*psi*norfp*norfp;
      //  compute new dlambda increment
      dlambda = dl/sqrt(tmp);
    }
 
    // check required value of load coefficient lambdar
    if ((check_rv_fp == on) && (dlambda+lambda > nlman->lambdar))
    {
      check_rv_fp = off;
      dlambda = nlman->lambdar - lambda;      
      dl = dlambda*sqrt(tmp);
    }
    // checking of maximum value of the total (cumulative) arc length
    if ((check_max_tot_al == on) && (dl+totl > nlman->max_tot_al))
    {
      check_max_tot_al = off;
      dl = nlman->max_tot_al - totl;      
    }

    //  ddr[j]=dlambda*v[j];
    cmulv(dlambda, v, ddr, n);

    //  ra[j]+=ddr[j];
    addv(ra, ddr, n);    

/*
    // update attained load vector
    // fa[j] = fc + (lambda+dlambda)*fp[j] is calculated
    addmultv(fc, fp, lambda+dlambda, fa, n);
*/    
    // update attained load vector
    // fa[j] = flc + (lambda+dlambda)*flp[j] is calculated
    addmultv(flc, flp, lambda+dlambda, fa, n);

    // initialize value of ddlambda before inner iteration loop
    ddlambda=dlambda;
    
    if (Mespr == 1){
      fprintf(Out, "Step %ld: ddlambda = %e, norv = %e, norfp = %e\n", i, ddlambda, norv, norfp);
      fflush(Out);
    }
    
    Mp->dlambda = ddlambda;
    Mp->lambda += ddlambda; // due to prescribed displacements
    //  computation of internal forces
    internal_forces (lcid,fi);
    Mp->lambda -= ddlambda; // restore state
    fflush(Out);

    // compute residual vector f[j] = fa[j]-fi[j]
    subv(fa, fi, f, n);    

    // compute norm of the residual vector
    norf = compute_res_norm_nr(lcid, normt, f, fa, n, norfa);

    // check time step size requirements which originate in material models, 
    // after check, the dtr will contain the minimum ratio of required time step size to the actual one 
    // ( dtr \in (0.0; 1.0) ) collected from all integration points or dtr = 1.0 for no change 
    // in time step size
    dtr = Mm->dstep_red_mat();
    if (dtr < 1.0) // reduction of time step size was required by some material model
    {
      norf = 2.0*ierr;
      if (dl == dlmin) // dl cannot be decreased
        break;
      //  modification of the arc length
      dl*=dtr;
      if (dl<dlmin)
        dl=dlmin;  

      if (Mespr==1)
      {
        if (dl == dlmin)
        {
          fprintf (stdout,"\n arc length was decreased due to mat. models dl=dlmin=%e",dl);
          // fprintf (Out,"\n arc length was decreased dl=dlmin=%e",dl);
        }
        else
        {
          fprintf (stdout,"\n arc length was decreased due to mat. models dl=%e",dl);
          // fprintf (Out,"\n arc length was decreased dl=%e",dl);
        }
      }
      //  restoring of left hand side vector
      copyv (r, ra, n);
      //  restoring of lambda parameter
      lambda=blambda;
      continue;
    }

    if (norf<ierr) 
    {
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      if (Mespr==1) 	  
        fprintf(stdout,"\n increment %ld   norf=%e  ierr=%e  ddlambda=%e  dlambda=%e",i,norf,ierr, ddlambda, dlambda);

      if (nlst == dissip_incr)
        totl += tau;
      else
        totl += dl;

      modif++;

      if (modif>1) 
      {
        //  arc length modification
        dl*=2.0;
        if (dl>dlmax)  
          dl=dlmax;
        modif=0;
        if (Mespr==1)
          fprintf(stdout,"\n arc-length is modified (dl increase) dl=%e",dl);
      }
      lambda+=dlambda;

      //  this part modifies parameter psi
      if (modif_psi)
        adapt_psi_coeff(i, ddlambda, ddlambda0, psi0, psi);
      
      Mm->updateipval();
      if (outres == yes)
      {
        compute_req_val(lcid);
        print_step(lcid, i+1, lambda, fa, f);
//        print_step(lcid, (i+1)*100, lambda, fa, f);
        print_flush();
      }

      // required value of lambda was attained for the load vector -> end of arclength
      if ((nlman->check_lambdar == on) && (check_rv_fp == off))
        break;

      // required total(cumulative) arc length has been attained => end of arclength
      if ((nlman->check_tot_al == on) && (check_max_tot_al==off)) 
        break;

      continue;
    }

    /*
    else{
      if (outres == yes){
        compute_req_val(lcid);
        print_step(lcid, (i+1)*100, lambda, fa, f);
        print_flush();
      }
    }*/

    // ***************************************
    // *       inner iteration loop          *
    // ***************************************
    if (Mespr==1) 	  
      fprintf (stdout,"\n inner loop: 1 norf=%e  ierr=%e  ddlambda=%e  dlambda=%e",norf,ierr,ddlambda, dlambda);
    stop=0;
    stopj=0;
    fillm(0.0, lsm_a);
    fillv(0.0, lsm_r);
    for (j=0;j<ini;j++)
    { //loop j
      //  copy of actual step to problem description
      //  it is used in connection with printing of output values
      Mp->jstep = j;

      //  eventual update of stiffness matrix
      Mespr = 0; // switch off the message printing for assembling of stiffness %matrix
      sm_change = assemble_stiffness_matrix (lcid,i,j,li,no);
      Mespr = 1; // switch on the message printing for assembling of stiffness %matrix
      
      // solve displacement increments due to vector of unbalanced forces
      Mp->ssle->solve_system(Gtm,Smat,u,f,Out);
      if (sm_change)
      {
        //  backup of the fp because the right hand side will be destroyed during the solution of the system of equations
        copyv(fp, f, n);
        Mp->ssle->solve_system(Gtm,Smat,v,f,Out);
      }
      
      // backup vector of displacement increments from the previous step
      // f[j] = ddr[j]
      copyv(ddr, f, n);

      // add computed correction of displacement increments due to unbalanced forces
      // ddr[j] = ddr[j]+u[j]
      addv(ddr, u, n);

      //  determination of increment of load parameter
      determine_dlambda (ddr, u, v, f, n, ddlambda,psi,norfp,dl,dlambda,stopj, nlman, lcid, i);
      if (stopj)
        break;

      // check required value of load coefficient lambdar
      if ((nlman->check_lambdar == on) && (check_rv_fp == off) && 
          (dlambda+ddlambda+lambda > nlman->lambdar))
        dlambda = nlman->lambdar - ddlambda - lambda;      
      
      // actualize vector of displacements increments
      // ddr[j] += dlambda*v[j]
      addmultv(ddr, v, dlambda, n);
        
      // actualize vector of attained displacements ra
      // ra[j] += u[j] + dlambda*v[j]
      addv(ra, u, n);
      addmultv(ra, v, dlambda, n);
/*
      // update attained load vector
      // fa[k]+=dlambda*fp[k] is calculated
      addmultv(fa, fp, dlambda, n);*/

      // update attained load vector
      // fa[k]+=dlambda*flp[k] is calculated
      addmultv(fa, flp, dlambda, n);
      //update_attained_load_vector(fa, fp, n, lambda+ddlambda, dlambda, nlman);

      ddlambda+=dlambda;
	
      Mp->dlambda = ddlambda;
      Mp->lambda += ddlambda;
      //  computation of internal forces
      internal_forces (lcid,fi);
      Mp->lambda -= ddlambda;

      // check time step size requirements which originate in material models, 
      // after check, the dtr will contain the minimum ratio of required time step size to the actual one 
      // ( dtr \in (0.0; 1.0) ) collected from all integration points or dtr = 1.0 for no change 
      // in time step size
      dtr = Mm->dstep_red_mat();
      if (dtr < 1.0) // reduction of time step size was required by some material model
      {
        stopj = 1;
        norf = 2.0*ierr;
        break;
      }

      //  computation of residual forces f[k] = fa[k] - fi[k]
      subv(fa, fi, f, n);	
      // compute norm of the residual vector
      norf = compute_res_norm_nr(lcid, normt, f, fa, n, norfa);

      if (Mespr==1)
      {
        fprintf (stdout,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
        // fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e   lambda=%le,  ddlambda=%le",i,j,norf, lambda, ddlambda);
      }

      // equilibrium was attained
      if (norf<ierr)
      {
        if (j < 0.25*ini) // only few steps were necessary => set flag for step length increasing
          modif++;
        else  // inner loop was necessary with too much steps => no step length increasing
          modif = 0;

        if (modif>1){
          //  arc length modification
          dl*=2.0;
          if (dl>dlmax)  
            dl=dlmax;
          modif=0;
          if (Mespr==1)
            fprintf(stdout,"\n arc-length is modified (dl increase) dl=%e",dl);
        }
        lambda+=ddlambda;
        Mm->updateipval();
        compute_req_val (lcid);
        print_step(lcid, i+1, lambda, fa, f);
        //print_step(lcid, (i+1)*100+j+1, lambda, fa, f);
        print_flush();
        totl += dl; 
        stop=1;
        break; 
      }
      /*
      else{
        print_step(lcid, (i+1)*100+j+1, lambda, fa, f);
        print_flush();
      }*/

      // divergence detection with help of least square method
      stopj = check_divergency(nlman, lsm_a, lsm_r, lsm_l, j, norf);
      if (stopj)
        break;

    } // end of j

    //  there is no modification of the length of arc because at least one iteration was performed
    modif=0;
      
    // check of the total(cumulative) attained length of arc
    if ((nlman->check_tot_al == on) && (check_max_tot_al==off))
      break;

    if (stop==1) // the equilibrium was attained
    {  
      //  this part modifies parameter psi
      if (modif_psi)
        adapt_psi_coeff(i, ddlambda, ddlambda0, psi0, psi);
            
      if ((nlman->check_lambdar == on) && (check_rv_fp == off) && (fabs(lambda-nlman->lambdar) > nlman->errl))
        check_rv_fp = on;
      
    }

    if ((stop==0) || (stopj)) // stop  == 0 the equilibrium was not attained in the prescribed number of iterations
                              // stopj != 0 the dlambda could not be determined or divergence detected or 
                              //            material models required time step reduction
    {
      if (dl == dlmin) // dl cannot be decreased
        break;
      //  modification of the arc length
      if (dtr < 1.0) // reduction of time step size was required by some material model
        dl *= dtr;
      else
        dl/=2.0;

      if (dl<dlmin)
        dl=dlmin;  

      if (Mespr==1)
      {        
        if (dl == dlmin)
        {
          if (dtr < 1.0)
            fprintf (stdout,"\n arc length was decreased due to mat. models dl=dlmin=%e",dl);
          else
            fprintf (stdout,"\n arc length was decreased dl=dlmin=%e",dl);
          // fprintf (Out,"\n arc length was decreased dl=dlmin=%e",dl);
        }
        else
        {
          if (dtr < 1.0)
            fprintf (stdout,"\n arc length was decreased due to mat. models dl=%e",dl);
          else
            fprintf (stdout,"\n arc length was decreased dl=%e",dl);
          // fprintf (Out,"\n arc length was decreased dl=%e",dl);
        }
      }
      //  restoring of left hand side vector
      copyv (r, ra, n);
      //  restoring of lambda parameter
      lambda=blambda;
      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();
    }

    // ------------------------------------
    //  finish of main iteration loop  ----
    // ------------------------------------
  }

  fprintf(stdout, "\n The total attained lambda =% le\n", lambda);
  
  if ((i < ni-1) && (norf >= ierr))
  {    
    if (Mespr==1)
    {
      print_err("the step length was decreased to the minimum required value\n"
                " but the equilibrium could not be attained.\n"
                " FORCED output of the attained results was performed in this step.", 
                __FILE__, __LINE__, __func__);
    }
    compute_req_val (lcid);
    print_step_forced(lcid, i+1, lambda, fa, f);
    //print_step_forced(lcid, (i+1)*100+j+1, lambda, fa, f);
    print_flush();
  }

  if (Mp->hdbcont.save_stat())
  {
    //  backup of the attained state is required
    if (Mespr==1)
      fprintf (stdout,"\n Creating backup file\n");
    solver_save (ra, fa, i, lambda, dl, NULL, n);
  }  
  

  delete [] r;		    
  delete [] fi;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;

  return lambda;
}



/**
  The function solves system of nonlinear algebraic equations
  by the arc-length method. This method was modified to reach the given value lambdar 
  of the lambda parameter. lambdar value is reached with required error errl. 
  Only one right hand side %vector is supported with respect to nonlinearity and absence of superposition. 
  Proportional %vector grows until the nial number of steps is reached.

  @param lcid    - load case id
  @param nlman   - pointer to structure conatining setup of the solver
  @param ra      - %vector of attained displacements
  @param fa      - attained load %vector
  @param fc      - constant part of load %vector (it is not influenced by lambda)
  @param fp      - proportional part of load %vector (it is increased by load coefficient lambda up to the nlman->lambdar)
  @param li      - initial value of step id (default is 0)
  @param ilambda - initial value of load coefficient (default is 0)
  @param outres  - flag for performing of output of results (if yes -> print_step procedure is called for the each step)

  @return The function returns reached lambda parameter.

  Created by JK,  16.8.2001
  Rewritten by TKo, JK 08.2011
*/
/*
double garclength2(long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp, 
                   long li, double ilambda, answertype outres)
{
  long i,j,n,ni,ini,stop,stopj,modif;
  double dl;       // actual value of arc length increment
  double dlmin;    // minimum value of arc length increment
  double dlmax;    // maximum value of arc length increment
  double psi0;     // initial value of proportional coefficient psi
  double blambda;  // backupped value of lambda
  double zero;     // computer zero
  double ierr;     // relative error of residual vector
  double *r;   // displacements from previous step (backup copy)
  double *ddr; // vector of displacement increments
  double *u;   // vector of corrections for computed displacements v
  double *v;   // vector of computed displacements for actual arclength step
  double *f;   // auxiliary vector of actual load (call of solve_system rewrite it during solution)
  double *fi;  // vector of internal forces
  double totl    = 0.0; // the total(cumulative) length of arc
  long sm_change = 0;   // flag for changed content of the stiffness matrix
  long modif_psi; // flag for automatic modification psi in dependence on the ddlambda
  matrix lsm_a(3,3);         // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  flagsw check_rv_fp = nlman->check_lambdar;  // flag for checking of attainment of required value for load coefficient lambdar 
  flagsw check_max_tot_al = nlman->check_tot_al; // flag fro checking maximu value of the length of arc
  
  
  
  //  number of rows of the matrix
  n = Ndofm;
  //  computer zero
  zero = Mp->zero;
  // initial value of load coefficient
  lambda = ilambda;
  //  maximum number of increments
  ni = nlman->nial;
  //  displacement-loading driving switch
  psi = nlman->psial;
  
  //  allocation of auxiliary arrays
  r   = new double[n];
  ddr = new double[n];
  u   = new double[n];
  v   = new double[n];
  f   = new double[n];
  fi  = new double[n];
  

  
  // assemble initial load vector, fa[j]=fc[j]+lambda*fp[j]
  // in case that no restorage from the backup was performed (lambda == 0.0)
  if (lambda == 0.0)
    assemble_attained_load_vector(fa, fc, fp, n, lambda, nlman);
  //  norm of proportionality vector
  norfp = normv(fp, n);
  //  norfp = loadincr(lcid, nlman, lambda, flp, n);
  modif=0;
  psi0 = psi;
  if (psi<0.0){
    //  parameter psi will be modified during the computation
    modif_psi=1;
  }else{
    //  parameter psi will be constant during the computation
    modif_psi=0;
  }


  // ***************************
  //  main iteration loop   ****
  // ***************************
  for (i=li;i<ni;i++) 
  { // loop i

    //  indicator of strain computation
    //  it means, strains at integration points have not been computed
    Mp->strainstate=0;
    //  indicator of stress computation
    //  it means, stresses at integration points have not been computed
    Mp->stressstate=0;
    //  indicator of computation of other array
    //  it means, stresses at integration points have not been computed
    Mp->otherstate=0;

    Mp->istep = i;
    Mp->jstep = -1;

    //  backup of left hand side vector rb[j]=ra[j];
    copyv(ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;

    if (Mespr==1)
    {
      fprintf (stdout,"\n increment %ld,   lambda=%e   dl=%e  psi=%e",i,lambda,dl,psi);
//      fprintf (Out,"\n increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
    }

    //  backup of the fp because the right hand side will be destroyed during the solution of the system of equations
    copyv(fp, f, n);
    

    //  there is no modification of the length of arc because at least one iteration was performed
    modif=0;
      
    norf = graclength_one_step()
    if (norf<ierr)       //  no inner iteration loop was required  ***

    { 
      totl += dl;
    
      if (j < 0) // no inner loop was required
        modif++;

      if (modif>1) // two successive steps did not require inner loop
      {
        //  length of arc can be increased
        dl*=2.0;
        if (dl>dlmax)  
          dl=dlmax;
        modif=0;
        if (Mespr==1)
          fprintf(stdout,"\n arc-length is modified (dl increase) dl=%e",dl);
      }
      lambda+=dlambda;

      //  this part modifies parameter psi
      if (modif_psi)
        adapt_psi_coeff(i, ddlambda, ddlambda0, psi0, psi);
    
      Mm->updateipval();
      if (outres == yes)
      {
        compute_req_val(lcid);
        print_step(lcid, i+1, lambda, fa);
        print_flush();
      }

      // required value of lambda was attained for the load vector -> end of arclength
      if ((nlman->check_lambdar == on) && (check_rv_fp == off))
        break;

      // required total(cumulative) arc length has been attained => end of arclength
      if ((nlman->check_tot_al == on) && (check_max_tot_al==off)) 
        break;
    }


    if (stop==1) // the equilibrium was attained
    {  
      //  this part modifies parameter psi
      if (modif_psi)
        adapt_psi_coeff(i, ddlambda, ddlambda0, psi0, psi);
            
      if ((nlman->check_lambdar == on) && (check_rv_fp == off) && (fabs(lambda-nlman->lambdar) > nlman->errl))
        check_rv_fp = on;
    }

    if ((stop==0) || (stopj)) // stop  == 0 the equilibrium was not attained in the prescribed number of iterations
                              // stopj != 0 the dlambda could not be determined or divergence detected
    {
      if (dl == dlmin) // dl cannot be decreased
        break;
      //  modification of the arc length
      dl/=2.0;
      if (dl<dlmin)
        dl=dlmin;  

      if (Mespr==1)
      {
        if (dl == dlmin)
        {
          fprintf (stdout,"\n arc length was decreased dl=dlmin=%e",dl);
          // fprintf (Out,"\n arc length was decreased dl=dlmin=%e",dl);
        }
        else
        {
          fprintf (stdout,"\n arc length was decreased dl=%e",dl);
          // fprintf (Out,"\n arc length was decreased dl=%e",dl);
        }
      }
      //  restoring of left hand side vector
      copyv (r, ra, n);
      //  restoring of lambda parameter
      lambda=blambda;
      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();
    }

    // ------------------------------------
    //  finish of main iteration loop  ----
    // ------------------------------------
  }
  
  if ((i < ni-1) && (norf >= ierr))
  {    
    if (Mespr==1)
    {
      print_err("the step length was decreased to the minimum required value\n"
                " but the equilibrium could not be attained.\n"
                " FORCED output of the attained results was performed in this step.", 
                __FILE__, __LINE__, __func__);
    }
    compute_req_val (lcid);
    print_step_forced(lcid, i+1, lambda, f);
    print_flush();
  }

  if (Mp->hdbcont.save_stat())
  {
    //  backup of the attained state is required
    if (Mespr==1)
      fprintf (stdout,"\n Creating backup file\n");
    solver_save (ra, fp, i, lambda, dl, NULL, n);
  }  
  

  delete [] r;		    
  delete [] fi;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;

  return lambda;
}
*/


/**
  Function performs calculation of one load/time step of the Newton-Raphson 
  method for the given load case. Solved equation system does not contain 
  time variable.
  
  @param lcid  - load case id
  @param nlman - pointer to structure conatining setup of the solver
  @param fa    - attained load %vector
  @param ra    - %vector of attained displacements
  @param fb    - residual %vector or load %vector increment - righthand side %vector
  @param dr    - %vector of displacement increment - lefthand side %vector
  @param fb    - %vector of internal forces
  @param f     - auxiliary %vector of actual load (call of solve_system rewrite it during solution)
  @param istep - time/load step id
  @param j     - inner loop step id (output parameter, set to -1 if no inner loop was performed) 
  @param li    - initial value of time/load step id

  @return The function returns reached lambda parameter.

  Created by Tomas Koudelka, 11.2012
*/
/*
double garclength_one_step(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, double *dr, double *fi, 
                           long istep, long &j, long li)
{
  double *r;   // displacements from previous step (backup copy)
  double *ddr; // vector of displacement increments
  double *u;   // vector of corrections for computed displacements v
  double *fi;  // vector of internal forces
  double *ra;  // %vector of attained displacements
  double *fa;  // attained load %vector
  double *fc;  // constant part of load %vector (it is not influenced by lambda)
  double *fp;  // proportional part of load %vector (it is increased by load coefficient lambda up to the nlman->lambdar)

  double  norfp;            // norm of proportional load vector
  double  psi;              // proportional coefficient between load vector and displacement vector
  long    modif_psi;        // flag for automatic modification psi in dependence on the ddlambda
  double  lambda;           // actual value of load coefficient
  flagsw &check_rv_fp;      // flag for checking of attainment of required value for load coefficient lambdar 
  flagsw &check_max_tot_al; // flag fro checking maximu value of the length of arc

  double *v;   // vector of computed displacements for actual arclength step
  double *f;   // auxiliary vector of actual load (call of solve_system rewrite it during solution)

  double norv;     // norm of displacement increments
  double norf;     // norm of residual vector
  double norfa;    // norm of attained load vector
  double dlambda;  // load coefficient increment from the inner loop
  double ddlambda; // cummulative value of load coefficent increments ddlambda = \sum(dlambda)
  double ddlambda0;// initial/reference value of ddlambda
  double tmp;      // auxiliary variable

  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of iterations in one increment
  ini = nlman->niilal;
  //  required error in inner loop
  ierr = nlman->erral;
  //  maximum number of iterations in one increment
  ini = nlman->niilal;
  //  required error in inner loop
  ierr = nlman->erral;
  //  length of arc
  dl = nlman->dlal;
  //  maximum length of arc
  dlmax = nlman->dlmaxal;
  //  minimum length of arc
  dlmin = nlman->dlminal;

  //  Delta lambda - cumulative increment in one load increment
  //  at the beginning of computation is always zero
  dlambda=0.0;

  // assembling of the stiffness matrix
  sm_change = assemble_stiffness_matrix (lcid, istep, -1, li, no);

  //  solution of K(r).v=F
  Mp->ssle->solve_system(Gtm, Smat, v, f, Out);

  //  square of generalized norm of displacement increments
  norv = displincr(lcid, istep, nlman, v, v, n);

  if ((istep == 0) && (psi < 0.0) && (norfp > Mp->zero))
    psi0 = psi = sqrt(norv/norfp/norfp);

  //  compute new dlambda increment
  tmp = norv+psi*psi*norfp*norfp;
  //  compute new dlambda increment
  dlambda = dl/sqrt(tmp);
 
  // check required value of load coefficient lambdar
  if ((check_rv_fp == on) && (dlambda+lambda > nlman->lambdar))
  {
    check_rv_fp = off;
    dlambda = nlman->lambdar - lambda;      
    dl = dlambda*sqrt(tmp);
  }
  // checking of maximum value of the total (cumulative) arc length
  if ((check_max_tot_al == on) && (dl+totl > nlman->max_tot_al))
  {
    check_max_tot_al = off;
    dl = nlman->max_tot_al - totl;      
  }

  if (Mespr == 1)
    fprintf(Out, "Step %ld: ddlambda = %e, norv = %e, norfp = %e\n", i, ddlambda, norv, norfp);
    
  //  ddr[j]=dlambda*v[j];
  copymultv(v, ddr, dlambda, n);

  //  ra[j]+=ddr[j];
  addv(ra, ddr, n);    

  // update attained load vector
  // fa[j] = fc + (lambda+dlambda)*fp[j] is calculated
  addmultv(fc, fp, lambda+dlambda, fa, n);
    

  // initialize value of ddlambda before inner iteration loop
  ddlambda=dlambda;
    
  //  computation of internal forces
  internal_forces (lcid, fi);

  // f[j] = fa[j]-fi[j]
  subv(fa, fi, f, n);    

  // computation of norms
  norf  = normv(f, n);
  norfa = normv(fa, n);

  // proportional norm is used in case of nonzero load vector
  // otherwise absolute norm is used for vector of unbalanced forces
  if (fabs(norfa) > Mp->zero)
    norf /= norfa;

  if (norf<ierr) 
  {
    // ******************************************
    //  no inner iteration loop is required  ***
    // ******************************************
    if (Mespr==1) 	  
      fprintf(stdout,"\n increment %ld   norf=%e  ierr=%e  ddlambda=%e  dlambda=%e",i,norf,ierr, ddlambda, dlambda);

    j = -1;
    return norf;
  }    

  // ***************************************
  // *       inner iteration loop          *
  // ***************************************
  if (Mespr==1) 	  
    fprintf (stdout,"\n inner loop: 1 norf=%e  ierr=%e  ddlambda=%e  dlambda=%e",norf,ierr,ddlambda, dlambda);
  stop =0;
  stopj=0;
  fillm(0.0, lsm_a);
  fillv(0.0, lsm_r);
  for (j=0; j<ini; j++)
  { //loop j
    //  copy of actual step to problem description
    //  it is used in connection with printing of output values
    Mp->jstep = j;

    //  eventual update of stiffness matrix
    Mespr = 0; // switch off the message printing for assembling of stiffness %matrix
    sm_change = assemble_stiffness_matrix (lcid,i,j,li,no);
    Mespr = 1; // switch on the message printing for assembling of stiffness %matrix
      
    // solve displacement increments due to vector of unbalanced forces
    Mp->ssle->solve_system(Gtm, Smat, u, f, Out);
    if (sm_change)
    {
      //  backup of the fp because the right hand side will be destroyed during the solution of the system of equations
      copyv(fp, f, n);
      Mp->ssle->solve_system(Gtm, Smat, v, f, Out);
    }
      
    // backup vector of displacement increments from the previous step
    // f[j] = ddr[j]
    copyv(ddr, f, n);

    // add computed correction of displacement increments due to unbalanced forces
    // ddr[j] = ddr[j]+u[j]
    addv(ddr, u, n);

    //  determination of increment of load parameter
    determine_dlambda (ddr, u, v, f, n, ddlambda,psi,norfp,dl,dlambda,stopj, nlman,lcid,istep);
    if (stopj)
      break;

    // check required value of load coefficient lambdar
    if ((nlman->check_lambdar == on) && (check_rv_fp == off) && 
        (dlambda+ddlambda+lambda > nlman->lambdar))
      dlambda = nlman->lambdar - ddlambda - lambda;      
      
    // actualize vector of displacements increments
    // ddr[j] += dlambda*v[j]
    addmultv(ddr, v, dlambda, n);
        
    // actualize vector of attained displacements ra
    // ra[j] += u[j] + dlambda*v[j]
    addv(ra, u, n);
    addmultv(ra, v, dlambda, n);
    // update attained load vector
    // fa[k]+=dlambda*fp[k] is calculated
    addmultv(fa, fp, dlambda, n);
    //update_attained_load_vector(fa, fp, n, lambda+ddlambda, dlambda, nlman);
    
    ddlambda+=dlambda;
	

    //  computation of internal forces
    internal_forces (lcid,fi);
    //  computation of residual forces f[k] = fa[k] - fi[k]
    subv(fa, fi, f, n);	
    //  norm of the vector of residuals
    norf=normv(f,n);
    //  norm of the vector of attained forces
    norfa = normv(fa, n);
    if (norfa > zero)
    {
      //  if the attained forces are nonzero, the norm of residuals is scaled by the norm of forces
      //  this scaling is reasonable because nondimensional quantity is obtained
      norf /= norfa;
    }

    if (Mespr==1)
    {
      fprintf (stdout,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
      // fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
    }

    // equilibrium was attained
    if (norf<ierr)
    { 
      lambda+=ddlambda;
      Mm->updateipval();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, fa);
      print_flush();
      totl += dl; 
      stop=1;
      break; 
    }

    // divergence detection with help of least square method
    stopj = check_divergency(nlman, lsm_a, lsm_r, lsm_l, j, norf);
    if (stopj)
      break;
    
  } // end of j
}
*/


/**
  The function initializes dof numbers for the particular
  types of displacement control. It is used in the arclength method.

  @return The function does not return anything, it performs initialization of 
          Mp->nlman.
  
  Created by JK,
*/
void seldofinit ()
{
  long i,j,k,l,ndofn;
  double x1,x2,y1,y2,z1,z2,length;
  
  switch (Mp->nlman->displnorm){
  case alldofs:{  break; }
  case seldofs:
  case seldofscoord:{
    for (i=0;i<Mp->nlman->nsdofal;i++){
      j=Mp->nlman->seldofal[i];
      Mp->nlman->seldofal[i]=Mt->give_dof (Mp->nlman->selnodal[i],j);
    }
    break;
  }
  case selmstr:
    break;
  case selecnodes:{
    Mp->nlman->nsdofal=0;
    for (i=0;i<Mp->nlman->nsnal;i++){
      Mp->nlman->nsdofal+=Mt->give_ndofn (Mp->nlman->selnodal[i]);
    }
    Mp->nlman->seldofal = new long [Mp->nlman->nsdofal];
    k=0;
    for (i=0;i<Mp->nlman->nsnal;i++){
      ndofn=Mt->give_ndofn (Mp->nlman->selnodal[i]);
      for (j=0;j<ndofn;j++){
	l=Mt->give_dof (Mp->nlman->selnodal[i],j);
	if (l!=0){
	  Mp->nlman->seldofal[k]=l;
	  k++;
	}
      }
    }
    Mp->nlman->nsdofal=k;
    break;
  }
  case nodesdistincr:{
    Mp->nlman->nsdofal=0;
    for (i=0;i<Mp->nlman->nsnal;i++){
      Mp->nlman->nsdofal+=Mt->give_ndofn (Mp->nlman->selnodal[i]);
    }
    Mp->nlman->seldofal = new long [Mp->nlman->nsdofal];
    k=0;
    for (i=0;i<Mp->nlman->nsnal;i++){
      ndofn=Mt->give_ndofn (Mp->nlman->selnodal[i]);
      for (j=0;j<ndofn;j++){
	Mp->nlman->seldofal[k]=Mt->give_dof (Mp->nlman->selnodal[i],j);
	k++;
      }
    }
    x1=Gtm->gnodes[Mp->nlman->selnodal[0]].x;    x2=Gtm->gnodes[Mp->nlman->selnodal[1]].x;
    y1=Gtm->gnodes[Mp->nlman->selnodal[0]].y;    y2=Gtm->gnodes[Mp->nlman->selnodal[1]].y;
    z1=Gtm->gnodes[Mp->nlman->selnodal[0]].z;    z2=Gtm->gnodes[Mp->nlman->selnodal[1]].z;
    length=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
    Mp->nlman->nxal=(x2-x1)/length;  Mp->nlman->nyal=(y2-y1)/length;  Mp->nlman->nzal=(z2-z1)/length;
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown displacement norm is required in function seldofinit (%s, line %d)",__FILE__,__LINE__);
  }
  }
}



/**
  The function computes generalized norm of load %vector.
  It is used in the arclength method.

  @param lcid - load case id
  @param fp - array with components of proportional load %vector due to force load
  @param n - the number of components of fp %vector
  @param lambda - attained value of the load coefficient
  @param nlman - structure with setup of handling of lambda and lambda tresholds

  @return The function returns required norm.

  Created by JK,
*/
double loadincr (long lcid, nonlinman *nlman, double lambda, double *fp, long n)
{
  long i;
  double norm;
  double *loadincr = new double[n];
  double *mstre = NULL;
  strastre *ssa = NULL;

  norm = 0.0;
  
  nullv(loadincr, n);
  if ((nlman->check_lambdar == off) || (lambda > nlman->lambdar))
    addv(loadincr, fp, n);

  switch (Mp->nlman->displnorm)
  {
    case alldofs:
    case nodesdistincr:
      norm = ss (loadincr,loadincr,n);
      break;
    case seldofs:
    case seldofscoord:
    case selecnodes:
      for (i=0;i<Mp->nlman->nsdofal;i++)
        norm+=loadincr[Mp->nlman->seldofal[i]-1]*loadincr[Mp->nlman->seldofal[i]-1];
      break;
    case selmstr:
      if (Mb->give_num_mstress_comp(lcid)){
        mstre = Mb->lc[lcid].mstress;        
        ssa   = Mb->lc[lcid].give_mstrastre();
      }
      else
        break;
      for(i=0; i<nlman->nmstrcomp; i++){
        if (nlman->mstrastre[i] == stress){
          if (nlman->mstrid[i] < Mt->max_ncompstr){
            if (ssa[nlman->mstrid[i]] == stress)
              norm += sqr(mstre[nlman->mstrid[i]]);
            else{
              print_err("invalid macro-value type of %ld-th component is required in arc-length computation.\nIt must be %d but it is %d\n",
                        __FILE__, __LINE__, __func__, nlman->mstrid[i]+1, ssa[nlman->mstrid[i]], nlman->mstrastre[i]);
              abort();
            }
          }
          else{
            print_err("invalid macro-stress component id %ld is required in arc-length computation.\nIt must be in [1;%ld]\n",
                      __FILE__, __LINE__, __func__, nlman->mstrid[i]+1, Mt->max_ncompstr);
            abort();
          }
        }
      }
      break;
    default:
      print_err("unknown type of load increment %d is required", __FILE__, __LINE__, __func__, Mp->nlman->displnorm);
      abort();
  }
  norm = sqrt(norm);
  delete [] loadincr;

  return norm;
}



/**
  The function computes square of generalized norm of %vector of 
  displacemnt increments. The function returns either square of %vector norm computed
  from the selected componets or it returns increment of distance between 
  selected nodes. It is used in the arclength method.

  @param lcid[in] - load case id
  @param istep[in] - load increment id
  @param nlman[in] - pointer to the nonlinear manager structure with the setup of displacement increment calculation
  @param dr1[in] - array with components of the first %vector of displacement increments
  @param dr2[in] - array with components of the second %vector of displacement increments
  @param n[in] - number of components of displincr array

  @return The function returns required norm.

  Created by JK,
*/
double displincr (long lcid, long istep, nonlinman *nlman, double *dr1, double *dr2, long n)
{
  long i;
  double norm=0.0,u1,u2,v1,v2,w1,w2,aux;
  double *mstra = NULL;
  strastre *ssa = NULL;
  
  switch (Mp->nlman->displnorm){
  case alldofs:{
    norm = scprd(dr1,dr2,n);
    break;
  }
  case seldofs:
  case seldofscoord:
  case selecnodes:{
    norm=0.0;
    for (i=0;i<Mp->nlman->nsdofal;i++){
      if (Mp->nlman->seldofal[i]>0)
        norm+=dr1[Mp->nlman->seldofal[i]-1]*dr2[Mp->nlman->seldofal[i]-1];
      else
      {
        aux=prdisplincr(i, lcid, istep);
        norm += aux*aux;
      }
    }
    break;
  }
  case selmstr:
    norm = 0.0;
    if (Mb->give_num_mstrain_comp(lcid)){
      mstra = Mb->lc[lcid].mstrain;
      ssa   = Mb->lc[lcid].give_mstrastre();
    }
    else
      break;
    for (i=0; i<nlman->nmstrcomp; i++){
      if (nlman->mstrastre[i] == strain){
        if ((nlman->mstrid[i] < Mt->max_ncompstr)){
          if (ssa[nlman->mstrid[i]] == strain)
            norm += sqr(mstra[nlman->mstrid[i]]);
          else{
            print_err("invalid macro-value type of %ld-th component is required in arc-length computation.\nIt must be %d but it is %d\n",
                      __FILE__, __LINE__, __func__, nlman->mstrid[i]+1, ssa[nlman->mstrid[i]], nlman->mstrastre[i]);
            abort();
          }
        }
        else{
          print_err("invalid macro-strain component id %ld is required in arc-length computation.\nIt must be in [1;%ld]\n",
                    __FILE__, __LINE__, __func__, nlman->mstrid[i]+1, Mt->max_ncompstr);
          abort();
        }
      }
    }
    break;
  case nodesdistincr:{
    if (Mp->nlman->probdimal==2){
      u1=0.0;  u2=0.0;  v1=0.0;  v2=0.0;
      if (Mp->nlman->seldofal[0]>0)  u1=dr1[Mp->nlman->seldofal[0]-1];
      else u1 = prdisplincr(0, lcid, istep);

      if (Mp->nlman->seldofal[1]>0)  v1=dr1[Mp->nlman->seldofal[1]-1];
      else v1 = prdisplincr(1, lcid, istep);

      if (Mp->nlman->seldofal[2]>0)  u2=dr2[Mp->nlman->seldofal[2]-1];
      else u2 = prdisplincr(2, lcid, istep);

      if (Mp->nlman->seldofal[3]>0)  v2=dr2[Mp->nlman->seldofal[3]-1];
      else v2 = prdisplincr(3, lcid, istep);
      
      norm=(u2-u1)*Mp->nlman->nxal + (v2-v1)*Mp->nlman->nyal;
      norm *= norm;
    }
    if (Mp->nlman->probdimal==3){
      u1=0.0;  u2=0.0;  v1=0.0;  v2=0.0;  w1=0.0;  w2=0.0;
      if (Mp->nlman->seldofal[0]>0)  u1=dr1[Mp->nlman->seldofal[0]-1];
      else u1 = prdisplincr(0, lcid, istep);

      if (Mp->nlman->seldofal[1]>0)  v1=dr1[Mp->nlman->seldofal[1]-1];
      else v1 = prdisplincr(1, lcid, istep);

      if (Mp->nlman->seldofal[2]>0)  w1=dr1[Mp->nlman->seldofal[2]-1];
      else w1 = prdisplincr(2, lcid, istep);

      if (Mp->nlman->seldofal[3]>0)  u2=dr2[Mp->nlman->seldofal[3]-1];
      else u2 = prdisplincr(3, lcid, istep);

      if (Mp->nlman->seldofal[4]>0)  v2=dr2[Mp->nlman->seldofal[4]-1];
      else v2 = prdisplincr(4, lcid, istep);

      if (Mp->nlman->seldofal[5]>0)  w2=dr2[Mp->nlman->seldofal[5]-1];
      else w2 = prdisplincr(5, lcid, istep);
      
      norm=(u2-u1)*Mp->nlman->nxal + (v2-v1)*Mp->nlman->nyal + (w2-w1)*Mp->nlman->nzal;
      norm *= norm;
    }
    break;
  }
  default:
    print_err("unknown norm of displacement increment is required", __FILE__, __LINE__, __func__);
  }
  return norm;
}



/**
  The function returns in the case of prescribed displacement in the i-th selected dof from nonlinman.

  @param i[in]     - id of selected dof
  @param lcid[in]  - load case id
  @param istep[in] - load step id

  @return The function returns contribution of to the length of arc from prescribed displacement
          of i-th selected dof.

  Created by Tomas Koudelka, 5.2014
*/
double prdisplincr(long i, long lcid, long istep)
{
  double ret = 0.0;

  if (Mp->nlman->seldofal[i] < 0)
  {
    if (Mb->lc[lcid].pd)
      ret = Mb->lc[lcid].pd[0-Mp->nlman->seldofal[i]-1];
    
    if ((Mb->lc[lcid+1].pd) && (istep == 0))
      ret += Mb->lc[lcid+1].pd[0-Mp->nlman->seldofal[i]-1];
  }
  return ret;
}



/**
  Function computes coefficients of the quadratic equation used in the arclength method.
  The coefficients are computed with respect to type of load/displacement control.

  @param[in] lcid     - load case id
  @param[in] istep    - load step id
  @param[in] nlman    - pointer to the structure with the setup of displacement increment calculation
  @param[in] ddr      - increment of load coefficient lambda for the given step
  @param[in] v        - %vector of displacement increments
  @param[in] n        - total number of DOFs (i.e. size of ddr and v arrays)
  @param[in] ddlambda - correction of increment of load coefficient lambda
  @param[in] psi      - proportional weight coeffient between load and displacement increment influence on lambda increment determination (the switching between load/displacement control)
  @param[in] norfp    - norm of proportional load %vector
  @param[in] dl       - length of arc
  @param[out] a0      - zero order coefficient
  @param[out] a1      - first order coefficient
  @param[out] a2      - second order coefficient

  @return The function returns evaluated coefficient in the parameters a0, a1 and a2.
 
  Created by JK,
*/
void quadeqcoeff (long lcid, long istep, nonlinman *nlman, double *ddr, double *v, long n, double ddlambda, double psi, double norfp, double dl,
		  double &a0, double &a1, double &a2)
{
  double norddr, norv, norddrv;
  
  switch (Mp->nlman->displnorm){
  case alldofs:
  case seldofs:
  case seldofscoord:
  case selmstr:
  case selecnodes:{
    a0 = displincr(lcid, istep, nlman, ddr, ddr, n);
    a2 = displincr(lcid, istep, nlman, v, v, n);
    norddrv = displincr(lcid, istep, nlman, ddr, v, n); 
    a1 = 2.0*norddrv;

    a0 += ddlambda*ddlambda*psi*psi*norfp*norfp - dl*dl;
    a1 += 2.0*ddlambda*psi*psi*norfp*norfp;
    a2 += psi*psi*norfp*norfp;
    break;
  }
  case nodesdistincr:{
    norddr = displincr(lcid, istep, nlman, ddr, ddr, n);
    norv   = displincr(lcid, istep, nlman, v, v, n);
    norddrv = sqrt(norddr)*sqrt(norv); // it is possible to calculate it in this way because ddr and v are scalars in this case
    a0 = norddr + ddlambda*ddlambda*norfp*norfp*psi*psi-dl*dl;
    a1 = 2.0*(norddrv + ddlambda*norfp*norfp*psi*psi);
    a2 = norv + norfp*norfp*psi*psi;
    break;
  }
  default:
    print_err("unknown norm of displacement increment is required", __FILE__, __LINE__, __func__);
  }
  // quadeqcoef_log(); // debugging log
}



/**
   The function determines increment of lambda (delta lambda).
   
   @param ddr     - Delta r + u 
   @param u       - u = K^{-1}*(fa-fi)
   @param v       - v = K^{-1}*fp
   @param ddrprev - Delta r (not updated by the %vector u)
   @param n       - number of unknowns
   @param ddlambda - Delta lambda
   @param psi   - proportional coefficient between load and displacement vectors
   @param norfp - the norm of proportional %vector
   @param dl    - length of the arc
   @param dlambda - delta lambda (output of this function)
   @param stop - indicator of end of iteration (output of this function)
   @param nlman - pointer to structure with arclength setup (it is used for full arclength method)
   @param lcid  - load case id
   @param istep - load step id
   
   @return The function returns delta lambda in the parameter dlambda and
           it sets the stop parameter to indicate the end of the arclength due
           to nonsolvable constrained conditions.

   Created by JK, 30.8.2010
   Modified
*/
void determine_dlambda (double *ddr, double *u, double *v, double *ddrprev, long n, 
                        double ddlambda, double psi, double norfp, double dl, 
                        double &dlambda, long &stop, nonlinman *nlman, long lcid, long istep)
{
  long k,numr;
  double l1,l2,aux,a0,a1,a2,ss1,ss2,ss3,ss4,ss5,nom,denom;
  
  switch (Mp->nlman->dlam)
  {
    case nodetermination:{
      print_err("some reasonable type of lambda determination has to be selected",__FILE__,__LINE__,__func__);
      break;
    }
    case minvalue:
    case maxvalue:
    case minangle:{
      //  coefficient of quadratic equation
      quadeqcoeff (lcid, istep, nlman, ddr, v, n, ddlambda, psi, norfp, dl, a0, a1, a2);
      //  solution of quadratic equation
      numr = solv_polynom_2(a2, a1, a0, l1, l2);
      switch (numr){
        case -1:{
          print_err("infinite number of solution of constrained condition\n"
                    "(all coefficients of quadratic equation are equal to zero", __FILE__, __LINE__, __func__);
          abort();
          break;
        }
        case 0:{
          print_err("nonsolvable constrained condition in function arclength", __FILE__, __LINE__, __func__);
          stop=1;
          break;
        }
        case 1:{
          dlambda = l1;
          break;
        }
        default:{}
      }
      break;
    }
    case linearizedmeth:{
      break;
    }
    case fullmethod:
      // unbalanced length of arc
      // a0 = ddrprev*ddrprev + ddlambda^2 * psi^2 * fp * fp - dl^2
      a0 = displincr(lcid, istep, nlman, ddrprev, ddrprev, n);
      a0 += ddlambda*ddlambda*psi*psi*norfp*norfp - dl*dl;
      break;
      /*    case dissip_incr:
            a0 = */
    default:{
      print_err("unknown type of detemination of lambda parameter is required",__FILE__,__LINE__,__func__);
    }
  }
  
  
  switch (Mp->nlman->dlam){
    case nodetermination:{
      print_err("some reasonable type of lambda determination has to be selected",__FILE__,__LINE__,__func__);
      break;
    }
    case minvalue:{
      if (fabs(l1)>fabs(l2))  
        dlambda=l2;
      else                    
        dlambda=l1;
      break;
    }
    case maxvalue:{
      if (fabs(l1)>fabs(l2))  
        dlambda=l1;
      else 
        dlambda=l2;
      break;
    }
    case minangle:{
      ss1=0.0;  ss2=0.0;
      ss3=0.0;  ss4=0.0;  ss5=0.0;
      for (k=0;k<n;k++){
        ss1+=(ddr[k]+l1*v[k])*ddrprev[k];
        ss2+=(ddr[k]+l2*v[k])*ddrprev[k];
        ss3+=(ddr[k]+l1*v[k])*(ddr[k]+l1*v[k]);
        ss4+=(ddr[k]+l2*v[k])*(ddr[k]+l2*v[k]);
        ss5+=ddrprev[k]*ddrprev[k];
      }
      if (ss1/sqrt(ss3)/sqrt(ss5)>ss2/sqrt(ss4)/sqrt(ss5))  
        dlambda=l1;
      else                          
        dlambda=l2;
      break;
    }
    case linearizedmeth:{
      nom=0.0;
      denom=0.0;
      for (k=0;k<n;k++){
        nom-=ddrprev[k]*(ddr[k]-ddrprev[k]);
        denom+=ddrprev[k]*v[k];
      }
      denom+=psi*psi*norfp*norfp*ddlambda;
    
      dlambda=nom/denom;
      break;
    }
    case fullmethod:
      //                  - a0 - ddrprev*K^{-1}*(fa-fi)
      // dlambda =  ----------------------------------------------
      //            2*(ddrprev*K^{-1}*fp + ddlambda*psi^2 * fp*fp)
      aux = displincr(lcid, istep, nlman, ddrprev, u, n);
      dlambda = -a0 - aux;
      aux = displincr(lcid, istep, nlman, ddrprev, v, n);
      dlambda /= 2.0*(aux + ddlambda*psi*psi*norfp*norfp);
      break;
    default:{
      print_err("unknown type of detemination of lambda parameter is required",__FILE__,__LINE__,__func__);
    }
  }
}


double dissipation_increment(double *fp, double *ddr, double *r, double lambda, long n)
{
  double ret = 0.0;
  ret += scprd(fp, ddr, n)*lambda;
  ret -= scprd(fp, r, n);

  return ret;
}



/**
  The function modifies the psi coefficient according to the
  ddlambda increment adaptively. If the increment decreases comparing to
  the initial/reference ddlambda increment (ddlamda0) than the value should be also decreased
  and thus the arclength is controlled rather by the displacements than by the load.
  The initial/reference value of the load increment ddlambda0 and proportional coefficient psi0
  are updated at the beginning of the main iteration loop and if the ddlambda changes
  the sign comparing to the ddlambda 0 (peak point of loading curve).
  
  @param i         - step id
  @param ddlambda  - actual increment of the load coefficeint lambda
  @param ddlambda0 - initial/reference value of the load coefficeint lambda increment (input/output parameter)
  @param psi0      - initial/reference value of the proportional coefficeint (input/output parameter)
  @param psi       - actual value of the proportional coefficeint (input/output parameter)

  @return The function returns the actual value of psi in the correspondign
          parameter. Also the parameters ddlambda0 and psi0 can be changed.

  Created by Tomas Koudelka 07.2011
*/
void adapt_psi_coeff(long i, double ddlambda, double &ddlambda0, double &psi0, double &psi)
{
  // at this point, ddlambda = dlambda
  // store initial values of load coefficient increment
  if (i==0)
    ddlambda0 = ddlambda;

  // load coefficient increment changed sign -> peak was attained ->
  // actualize values of stored initail load coefficient increment and initial proportional coefficient
  if ((sgn(ddlambda0) != sgn(ddlambda)) && (sgn(ddlambda) != 0.0))
  {
    psi0= psi;
    ddlambda0 = ddlambda;
  }
	
  // automatic correction of proportional coefficient
  // with respect to actual and initial load coefficient increment
  if (ddlambda0 != 0.0)
    psi = psi0*fabs(ddlambda/ddlambda0);

  return;
}
