#include "newtonraph.h"
#include "nssolver.h"
#include "backupsol.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "globmat.h"
#include "mechprint.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "mechprint.h"
#include "mathem.h"
#include "vector.h"
#include "matrix.h"
#include "node.h"
// for evaluation of hypoplasticity ???!!!
#include "intpoints.h"
//#include "hypoplunsatexptherm.h"
#include "hypoplunsatexptherm2.h"
#include "elemswitch.h"
//#include "hypoplunsatexptherm2.h"
#include "elemswitch.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>




/**
  Function solves system of nonlinear algebraic equations by
  Newton-Raphson method for the given load case.
  Solved system does not contain time variable.
  
  @param lcid  - load case id
  @param nlman - pointer to structure conatining setup of the solver
  @param ra    - %vector of attained displacements
  @param fa    - attained load %vector
  @param fc    - %vector of constant load
  @param fp    - %vector of proportional load
  @param li      - initial value of step id (default is 0)
  @param ilambda - initial value of load coefficient (default is 0)
  @param outres  - flag for performing of output of results (if yes -> print_step procedure is called for the each step)

  @return The function returns reached lambda parameter.

  Created by JK,  16.8.2001
  Rewritten by Tomas Koudelka, 08.2011
*/
double gnewton_raphson (long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp,
                        long li, double ilambda, answertype outres)
{
  long i,j,n,ni,ini,modif;
  double lambda;      // load coefficient
  double blambda;     // backup of load coefficient
  double dlambda;     // load coefficient increment
  double dlambdamin;  // minimum value of load coefficient increment
  double dlambdamax;  // maximum value of load coefficient increment
  double ierr;        // required error of residua
  double norf;        // norm of vector of unbalanced forces
  double norfa;       // norm of attained load vector (+ reaction vector eventually)
  //  double norfa;       // norm of attained load vector
  double *rb;         // backup of displacement vector
  double *dr;         // vector of displacement increments
  double *fi;         // vector of internal forces
  double *fb;         // vector of load increment
  resnormt normt;     // type of residual vector norm
  matrix lsm_a(3,3); // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  flagsw check_rv_fp = nlman->check_lambdar;  // flag for checking of attainment of required value for load coefficient lambdar 


  //  number of rows of the matrix
  n = Ndofm;
  // initial value of load coefficient
  lambda = ilambda;
  //  maximum number of increments
  ni = nlman->ninr;
  //  maximum number of iterations in one increment
  ini = nlman->niilnr;
  //  required error in the inner loop
  ierr = nlman->errnr;
  //  increment size
  dlambda=nlman->incrnr;
  //  minimum increment size
  dlambdamin=nlman->minincrnr;
  //  maximum increment size
  dlambdamax=nlman->maxincrnr;
  //type of residual vector norm
  normt = nlman->rnormtnr;

  norf = 0.0;

  //  initialization phase
  rb  = new double [n];
  fb  = new double [n];
  fi  = new double [n];
  dr  = new double [n];
  memset (rb,0,n*sizeof(double));
  memset (fb,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));
  memset (dr,0,n*sizeof(double));
  
  modif=0;
  
  // ***************************
  //  main iteration loop  ****
  // ***************************

  // assemble initial load vector, fa[j]=fc[j]+lambda*fp[j]
  // in case that no restorage from backup was performed (lambda == 0.0)
  if (lambda == 0.0)
    assemble_attained_load_vector(fa, fc, fp, n, lambda, nlman);

  for (i=li;i<ni;i++)
  {
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

    //  backup of left hand side vector
    //  rb[j]=ra[j];
    copyv(ra, rb, n);

    //  backup of reached lambda parameter
    blambda=lambda;
    Mp->lambda = lambda;
   
    fprintf (stdout,"\n increment %ld   lambda %e,    dlambda %e",i,lambda,dlambda);
    
    
    // check required value of load coefficient lambdar
    if ((check_rv_fp == on) && (dlambda+lambda > nlman->lambdar))
    {
      check_rv_fp = off;
      dlambda = nlman->lambdar - lambda;      
    }
    // update attained load vector
    // fa[j] += dlambda*fp[j] 
    addmultv(fa, fp, dlambda, n);

    // increment lambda for calculation of contributions from prescribed displacements
    // Mp->lambda = lambda+dlambda
    Mp->lambda += dlambda;
    Mp->dlambda = dlambda;

    // internal forces caused by the increased prescribed displacements
    internal_forces(lcid, fi);

    //  vector of load increment
    //  fb[j] = fa[j] - fi[j] = fc[j] + (lambda+dlambda)*fp[j] - fi[j]
    //  vector fi contains internal forces from the load vector of the previous time step and
    //  increment of internal forces from the new increment of prescribed displacements
    // 
    //  i.e. vector fb contains load increment caused force load + load caused by increment of 
    //  prescribed displacements
    subv(fa, fi, fb, n);

    
    //  assembling of tangent stiffness matrix
    assemble_stiffness_matrix(lcid,i,-1,li,no);
    
    //  solution of K(r).v=F
    Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

    //  update total displacement vector
    //  ra[j]+=dr[j];
    addv(ra, dr, n);
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    
    //  vector of unbalanced forces
    //  fb[j] = fa[j] - fi[j];
    subv(fa, fi, fb, n);

    // compute norm of the residual vector
    norf = compute_res_norm_nr(lcid, normt, fb, fa, n, norfa);
    
    fprintf (stdout,"\n    norm %e",norf);
    //    fprintf (Out,"\n norf=%e",norf);

    if (norf<ierr)
    {
      lambda+=dlambda;
      Mm->updateipval();
      if (outres == yes)
      {
        compute_req_val(lcid);
        print_step(lcid, i+1, lambda, fa, fb);
        print_flush();
      }
      modif++;
      if (modif>1)
      {
	dlambda*=2.0;
	if (dlambda>dlambdamax)
	  dlambda=dlambdamax;
	if (Mespr==1)
	  fprintf (stdout,"\n increment must be modified (dlambda increase) dlambda=%e",dlambda);
	modif=0;
      }
      // required value of lambda was attained for both load vectors -> end of newton-raphson procedure
      if ((nlman->check_lambdar == on) && (check_rv_fp == off))
        break;

      continue;
    }
    
    //  internal iteration loop
    fillm(0.0, lsm_a);
    fillv(0.0, lsm_r);
    for (j=0;j<ini;j++)
    {
      Mp->jstep = j;      

      // re-assemble stifness matrix if it is required
      assemble_stiffness_matrix(lcid,i,j,li,no);

      //  solution of K(r).v=F
      Mp->ssle->solve_system(Gtm,Smat,dr,fb,Out);
    
      // ra[k]+=dr[k];
      addv(ra, dr, n);
      
      //  computation of internal forces
      internal_forces(lcid,fi);
      
      //  vector of unbalanced forces
      //  fb[k]=fa[k]-fi[k]
      subv(fa, fi, fb, n);

      // compute norm of the residual vector
      norf = compute_res_norm_nr(lcid, normt, fb, fa, n, norfa);

      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf %.15e  dlambda %e",i,j,norf,dlambda);
      //fprintf (Out,"\n j=%ld     norf=%e",j,norf);//debug??!!
      
      if (norf<ierr) 
      // equilibrium was attained
        break;

      // divergence detection with help of the least square method
      if (check_divergency(nlman, lsm_a, lsm_r, lsm_l, j, norf))
        break;
    }
    
    modif=0;

    if (j==ini || norf>ierr)
    {
      // go back to the previous value of the attained load vector
      // fa[j] -= dlambda*fp[j]
      addmultv(fa, fp, -dlambda, n);
      //  ra[j]=rb[j];
      copyv(rb, ra, n);
      lambda=blambda;
      Mp->lambda = blambda;

      if (dlambda == dlambdamin)
      {
        if (Mespr==1)  fprintf(stdout,"\n increment cannot be decreased dlambda=%e    norf %e",dlambda,norf);
        break;
      }
      dlambda/=2.0;
      if (dlambda<dlambdamin)
      {
	dlambda=dlambdamin;
//	break;
      }
      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();
      if (Mespr==1)
      {
	fprintf(stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf:%e",dlambda,norf);
//	fprintf (Out,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
      }      
    }
    else
    {
      lambda+=dlambda;
      Mm->updateipval();
      if (outres == yes)
      {
        compute_req_val(lcid);
        print_step(lcid, i+1, lambda, fa, fb);
        print_flush();
      }
      // required value of lambda was attained for both load vectors -> end of arclength
      if ((nlman->check_lambdar == on) && (check_rv_fp == off))
        break;

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
    Mm->updateipval();
    compute_req_val (lcid);
    print_step_forced(lcid, i+1, lambda, fa, fb);
    print_flush();
  }

  if (Mp->hdbcont.save_stat())
  {
    //  backup of the attained state is required
    if (Mespr==1)
      fprintf (stdout,"\n Creating backup file\n");
    solver_save (ra, fp, i, lambda, dlambda, NULL, n);
  }  

  delete [] rb;
  delete [] dr;  
  delete [] fi;  
  delete [] fb;  

  return lambda;
}



/**
  Function solves system of nonlinear algebraic equations by
  Newton-Raphson method for the given load case.
  Solved system does not contain time variable.
  
  @param lcid  - load case id
  @param nlman - pointer to structure conatining setup of the solver
  @param ra    - %vector of attained displacements
  @param fa    - attained load %vector
  @param fc    - %vector of constant load
  @param fp    - %vector of proportional load
  @param flc -  constant component of load %vector due to forces
  @param flp -  proportional component of load %vector due to forces
  @param li      - initial value of load/time step id (default is 0)
  @param ilambda - initial value of load coefficient (default is 0)
  @param outres  - flag for performing of output of results (if yes -> print_step procedure is called for the each step)

  @return The function returns reached lambda parameter.

  Created by JK,  16.8.2001
  Rewritten by Tomas Koudelka, 11.2012
*/
double gnewton_raphson2 (long lcid, nonlinman *nlman, double *ra, double *fa, double *fc, double *fp, 
                         double *flc, double *flp, long li, double ilambda, answertype outres)
{
  long i,j,n,ni,ini,modif;
  double lambda;      // load coefficient
  double blambda;     // backup of load coefficient
  double dlambda;     // load coefficient increment
  double dlambdamin;  // minimum value of load coefficient increment
  double dlambdamax;  // maximum value of load coefficient increment
  double dtr;         // time step size coefficient required by material models
  double ierr;        // required error of residua
  double norf;        // norm of vector of unbalanced forces
  double norra;       // norm of vector of attained displacements
  double *rb;         // backup of displacement vector
  double *dr;         // vector of displacement increments
  double *fi;         // vector of internal forces
  double *fb;         // vector of load increment
  matrix lsm_a(3,3); // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  flagsw check_rv_fp = nlman->check_lambdar;  // flag for checking of attainment of required value for load coefficient lambdar 


  //  number of rows of the matrix
  n = Ndofm;
  // initial value of load coefficient
  lambda = ilambda;
  //  maximum number of increments
  ni = nlman->ninr;
  //  maximum number of iterations in one increment
  ini = nlman->niilnr;
  //  required error in inner loop
  ierr = nlman->errnr;
  //  increment size
  dlambda=nlman->incrnr;
  //  minimum increment size
  dlambdamin=nlman->minincrnr;
  //  maximum increment size
  dlambdamax=nlman->maxincrnr;

  //  initialization phase
  rb  = new double [n];
  fb  = new double [n];
  fi  = new double [n];
  dr  = new double [n];
  memset (rb,0,n*sizeof(double));
  memset (fb,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));
  memset (dr,0,n*sizeof(double));
  
  modif=0;

  // just for evaluation of hypoplasticity model ???!!!  
  Neval = 0.0;

  // ***************************
  //  main iteration loop  ****
  // ***************************

  // norm of attained displacement vector
  // norra>0.0 indicates start of the procedure from the nonzero stage and
  // therefore the constant load fc is not involved in the initial load increment vector fb
  norra = normv(ra, n);

  // assemble initial load vector, fa[j]=flc[j]+lambda*flp[j] due to forces
  // in case that no restorage from backup was performed (lambda == 0.0)
  if (lambda == 0.0)
    assemble_attained_load_vector(fa, flc, flp, n, lambda, nlman);

  for (i=li;i<ni;i++)
  {
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

    //  backup of left hand side vector
    //  rb[j]=ra[j];
    copyv(ra, rb, n);

    //  backup of reached lambda parameter
    blambda=lambda;
    Mp->lambda = lambda;
   
    fprintf (stdout,"\n increment %ld   lambda %e,    dlambda %e",i,lambda,dlambda);
    //fprintf (Out,"\n increment %ld   lambda %e,    dlambda %e",i,lambda,dlambda);debug??!!
    
    // check required value of load coefficient lambdar
    if ((check_rv_fp == on) && (dlambda+lambda > nlman->lambdar))
    {
      check_rv_fp = off;
      dlambda = nlman->lambdar - lambda;      
    }

    // update attained load vector due to forces
    // fa[j] += dlambda*flp[j] 
    addmultv(fa, flp, dlambda, n);

    // increment lambda for calculation of contributions from prescribed displacements
    // Mp->lambda = lambda+dlambda
    Mp->lambda += dlambda;
    Mp->dlambda = dlambda;

    // assemble load increment due to forces and prescribed displacements
    // fb[j] =  dlambda*fp[j]
    cmulv(dlambda, fp, fb, n);
    // for the initial increment (no attained displacements) add constant load vector
    // fb[j] += fc[j]
    if (norra == 0.0)
      addv(fb, fc, n);

    //  i.e. vector fb contains load increment caused force load + load caused by increment of 
    //  prescribed displacements
    
    // perform one step of Newton-Raphson iterative procedure
    norf = gnewton_raphson_one_step(lcid, nlman, fa, ra, fb, dr, fi, dtr, i, j, li, no);

    // update status of norm of attained dispalcement vector
    if (norra == 0.0)
      norra = normv(ra, n);
    // equilibrium was attanied in the performed step
    if (norf<ierr)
    {
      if ((nlman->check_lambdar == on) && (check_rv_fp == off))
        lambda = nlman->lambdar;
      else
        lambda+=dlambda;
      Mm->updateipval();
      if (outres == yes)
      {
        compute_req_val(lcid);
        print_step(lcid, i+1, lambda, fa, fb);
        print_flush();
      }
/*
      double ep = 0.0;
      for (long ii=0; ii<Mt->nn; ii++)
      {
        ivector cn;
        vector r;
        long ndofn = Mt->give_ndofn(ii);
        reallocv(ndofn, r);
        reallocv(ndofn, cn);
        noddispl (lcid, r.a, ii);
        Mt->give_node_code_numbers(ii, cn.a);
        for (long jj=0; jj<ndofn; jj++)
        {
          if (cn[jj] > 0)
          {
            ep += fi[cn[jj]]*r[jj];
            continue;
          }
          if ((cn[jj] < 0) && Mt->nodes[ii].react)
          {
            ep += Mt->nodes[ii].r[jj]*r[jj];
          }
        }
      }
     
      fprintf(stdout, "\nPotential energy Ep=%le\n", ep);
*/
      if (Mp->homog == 9){
        long j, ncomp = Mt->max_ncompstr;
        vector tifor(ASTCKVEC(ncomp)), ifor(ASTCKVEC(ncomp));
        
        for(j=0; j<Mm->tnip; j++)
          Mm->computenlstresses(j,Mm->ip[j]);
        for (j=0; j<Mt->ne; j++){
          nullv(ifor);
          elem_volintegration_quant(j, locstress, lcid, ifor);
          addv(tifor, ifor, tifor);
        }
        fprintf(stdout, "\nMacro-stress components:\n");
        printv(stdout, tifor);
      }

      // set step length modification flag
      if (j < 0.25*ini) // no inner loop was necessary or a few steps => set flag for step length increasing
        modif++;
      else           // inner loop was necessary with too much steps => no step length increasing
        modif=0;

      // modification of load step length 
      if (modif>1)   // at least two successive steps were performed without inner loop
      {
	dlambda*=2.0; // double step length
	if (dlambda>dlambdamax)
	  dlambda=dlambdamax;
	if (Mespr==1)
	  fprintf (stdout,"\n increment was modified (dlambda increase) dlambda=%e",dlambda);
	modif=0;
      }

      // required value of lambda was attained for both load vectors -> end of newton-raphson procedure
      if ((nlman->check_lambdar == on) && (lambda == nlman->lambdar))
      {
        // just for evaluation of hypoplasticity model ???!!!  
        if (Mm->hypoplustherm)
        {
          long i, ncompstr, nstatev = Mm->hypoplustherm[0].nstatev;
          for(i=0; i<Mm->tnip; i++)
          {
            ncompstr = Mm->ip[i].ncompstr;
            Neval += Mm->ip[i].eqother[2*ncompstr+3+nstatev+5];
          }
        }
        break;
      }
      continue;
    }
    
    
    // equilibrium was NOT attanied in the performed step
    if (norf>ierr)
    {
      // go back to the previous value of the attained load vector
      // fa[j] -= dlambda*flp[j]
      addmultv(fa, flp, -dlambda, n);
      //  ra[j]=rb[j];
      copyv(rb, ra, n);
      lambda=blambda;
      Mp->lambda = blambda;

      // copy state variables from the last attained equilibrium state (eqother array) to the actual state variables (other array)
      Mm->updateother();

      if (dlambda == dlambdamin) // load increment has been already set to minimum value and cannot be reduced
      {
        if (Mespr==1)
        {
          fprintf(stdout,"\n load increment cannot be decreased, dlambda=%e    norf %e",dlambda,norf);
          print_err("the step length was decreased to the minimum required value\n"
                    " but the equilibrium could not be attained.\n"
                    " FORCED output of the attained results was performed in this step.", 
                    __FILE__, __LINE__, __func__);
          // reaction should not be calculated because residual vector is sent
          // to the output instead of attained load vector and they decay the vector
          // plot of residual forces in GiD
          if (outres == yes)
          {
            long rc = Mp->reactcomp;
            Mp->reactcomp = 0;
            Mm->updateipval();
            compute_req_val (lcid);
            print_step_forced(lcid, i+1, lambda, fa, fb);
            Mp->reactcomp = rc;
            print_flush();
          }
        }
        break; // terminate main loop (for i)
      }
      if (dtr < 1.0)
        dlambda *= dtr;
      else
        dlambda/=2.0;  // decrease load coefficient increment because equilibrium was not attained 
      if (dlambda<dlambdamin)
      {
	dlambda=dlambdamin;
//	break;
      }
      if (Mespr==1)
      {
	fprintf(stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf:%e",dlambda,norf);
//	fprintf (Out,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
      }      
    }

    // just for evaluation of hypoplasticity model ???!!!  
    if (Mm->hypoplustherm)
    {
      long i, ncompstr, nstatev = Mm->hypoplustherm[0].nstatev;
      for(i=0; i<Mm->tnip; i++)
      {
        ncompstr = Mm->ip[i].ncompstr;
        Neval += Mm->ip[i].eqother[2*ncompstr+3+nstatev+5];
      }
    }
    
    // ------------------------------------
    //  finish of main iteration loop  ----
    // ------------------------------------
  }
  
  fprintf(stdout, "\n The total attained lambda =% le\n", lambda);

  if (Mp->hdbcont.save_stat())
  {
    //  backup of the attained state is required
    if (Mespr==1)
      fprintf (stdout,"\n Creating backup file\n");
    solver_save (ra, fa, i, lambda, dlambda, NULL, n);
  }  

  delete [] rb;
  delete [] dr;  
  delete [] fi;  
  delete [] fb;  

  return lambda;
}






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
  @param fi    - %vector of internal forces
  @param dtr   - the minimum ratio of required time step size to the actual one required by material models (output)
  @param istep - time/load step id
  @param j     - inner loop step id (output parameter, set to -1 if no inner loop was performed) 
  @param li    - initial value of time/load step id
  @param fusm  - flag for the force update of stiffnes %matrix
                 (fusm=yes - the stiffness %matrix will be updated regardless of other settings
                  fusm=no  - the stiffness %matrix will be updated as usual)

  @return The function returns reached norm of residual %vector.

  Created by Tomas Koudelka, 11.2012
*/
double gnewton_raphson_one_step(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, double *dr, double *fi, 
                                double &dtr, long istep, long &j, long li, answertype fusm)
{
  long n = Ndofm;    // number of rows of the matrix or vectors
  long ini;          // maximum number of inner iterations
  double ierr;       // required error of the residual vector
  double norf;       // norm of residual vector (vector of unbalanced forces)
  double norfa;       // norm of attained load vector (+ reaction vector eventually)
  resnormt normt;    // type of residual vector norm
  matrix lsm_a(3,3); // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  
  //  maximum number of iterations in one increment
  ini = nlman->niilnr;
  //  required error in the inner loop
  ierr = nlman->errnr;
  // type of residual vector norm
  normt = nlman->rnormtnr;

  //if(Mp->time >= 5.825400e+05)
  //fprintf (stdout,"\n zde cas!!"); //debug

  //  assembling of tangent stiffness matrix
  assemble_stiffness_matrix(lcid,istep,-1,li,fusm);

  /* if (Mp->time >= 5.375400e+05 && Mp->time<= 5.385400e+05){//debug??!!
     fprintf (Out,"\n\n Matice K ulohy pred, time %le\n",Mp->time);
     Smat->printmat (Out); 
     fflush(Out);//debug
     }
  */

  //  solution of K(r).\Delta r=\Delta F
  Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

  //  update total displacement vector
  //  ra[j]+=dr[j];
  addv(ra, dr, n);
    
  //  computation of internal forces
  internal_forces (lcid,fi);

  //double ep_fi = scprd(fi, ra, n);  // potential energy of internal forces
  //double ep_fe = scprd(fa, ra, n);  // potential energy of external load

  // check time step size requirements which originate in material models, 
  // after check, the dtr will contain the minimum ratio of required time step size to the actual one 
  // ( dtr \in (0.0; 1.0) ) collected from all integration points or dtr = 1.0 for no change 
  // in time step size
  dtr = Mm->dstep_red_mat();
  if (dtr < 1.0) // reduction of time step size was required by some material model
  {
    j=-1;
    return (2.0*ierr);
  }

  //  vector of unbalanced forces
  //  fb[j] = fa[j] - fi[j];
  subv(fa, fi, fb, n);

  // compute norm of the residual vector
  norf = compute_res_norm_nr(lcid, normt, fb, fa, n, norfa);
    
  //fprintf(Outdm->outf, "\n****** Load step %ld, norf=%le:\n\n", istep+1, norf);
  //Outdm->eo.print_out(Outdm->outf, lcid);

  fprintf (stdout,"\n  Initial residual norm norf=%.15le, norm of load vector norfa=%.15le",norf,norfa);
  //fprintf (Out,"\n norf=%e\n",norf);

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
    assemble_stiffness_matrix(lcid,istep,j,li,fusm);

    //fprintf (Out,"\n\n Matice K ulohy v iteraci c. %ld, time %le\n",j,Mp->time);
    //Smat->printmat (Out); 
    //fflush(Out);//debug

    // solution of K(r).v=F
    Mp->ssle->solve_system(Gtm,Smat,dr,fb,Out);

    // ra[k]+=dr[k];
    addv(ra, dr, n);
      
    // computation of internal forces
    internal_forces(lcid,fi);
    //ep_fi = scprd(fi, ra, n);

    // check time step size requirements which originate in material models, 
    // after check, the dtr will contain the minimum ratio of required time step size to the actual one 
    // ( dtr \in (0.0; 1.0) ) collected from all integration points or dtr = 1.0 for no change 
    // in time step size
    dtr = Mm->dstep_red_mat();
    if (dtr < 1.0) // reduction of time step size was required by some material model
    {
      j=-1;
      return (2.0*ierr);
    }
      
    // vector of unbalanced forces
    // fb[k]=fa[k]-fi[k]
    subv(fa, fi, fb, n);
      
    // compute norm of the residual vector
    norf = compute_res_norm_nr(lcid, normt, fb, fa, n, norfa);

    //fprintf(Outdm->outf, "\n****** Internal loop %ld, norf=%le:\n\n", j+1, norf);
    //Outdm->eo.print_out(Outdm->outf, lcid);
      
    fprintf (stdout,"\n increment %ld   inner loop j %ld     norm %.15le,   norfa %.15le", istep, j, norf, norfa);
    //fprintf (Out,"\n j=%ld     norf=%e",j,norf);
    /*
    vector pe(ASTCKVEC(1));
    double totpe = 0.0;
    long k;
    for(k=0; k<Mt->ne; k++){
      elem_volintegration_quant(k, penergydens, lcid, pe);
      totpe += pe(0);
      fprintf(Out, "\nistep=%ld, jstep=%ld, eid=%ld, pe=%le", Mp->istep, Mp->jstep, k+1, pe(0));
    }
    fprintf(Out, "\nistep=%ld, jstep=%ld, eid=%ld, totpe=%le\n", Mp->istep, Mp->jstep, k+1, totpe);
    compute_req_val (lcid);
    print_step_forced(lcid, istep*100+j, double(istep*100+j), fb);
    print_flush();
    */  
    if (norf<ierr) 
      // equilibrium was attained
      return norf;

    // divergence detection with help of the least square method
    if (check_divergency(nlman, lsm_a, lsm_r, lsm_l, j, norf))
      return norf;
  }

  // equilibrium was not attained in this step
  return norf;  
}


/**
  The function computes norm of the residual vector with respect to the given norm type.
  The norm is computed according to the argument normt which specifies whether the norm is either
  relative to the attained load vector fb or reactions or it may be absolute.

  @param[in] lcid - load case id
  @param[in] normt - type of the residual norm
  @param[in] fb - residual vector
  @param[in] fa - load vector
  @param[in] n  - dimension of vectors fb and fa
  @param[out] norfa - attained norm of the load vector

  @return The function returns norm of the residual vector fb with respect to settings given by normt.
*/
double compute_res_norm_nr(long lcid, resnormt normt, double *fb, double *fa, long n, double &norfa)
{
  double norf=0.0;
  norfa=0.0;
  const double zero = 1.0e-9;

  //  norm of vector of unbalanced forces (residual norm)
  norf=normv(fb,n);
  switch(normt){
    case rel_load_norm: // norm is relative to attained load vector
      // compute norm of attained load vector
      norfa=normv(fa,n);
      if (norfa > zero)
        norf /= norfa;
      break;
    case rel_react_norm: // norm is relative to attained reactions
      compute_reactions(lcid);
      norfa = Mt->compute_react_norm();
      if (norfa > zero)
        norf /= norfa;
      break;
    case rel_loadreact_norm: // norm is relative to attained reactions and load
      compute_reactions(lcid);
      norfa = Mt->compute_react_norm();
      norfa *= norfa;
      // compute norm of attained load vector
      norfa += scprd(fa, fa, n);
      norfa = sqrt(norfa);
      if (norfa > zero)
        norf /= norfa;
      break;
    case absol_norm: //absolute norm of the residual vector
      break;
    default:
      print_err("unknown type of residual norm %d is required", __FILE__, __LINE__, __func__, normt);
      abort();
  }

  return norf;
}
