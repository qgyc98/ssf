#include "slopesol.h"
#include "arclength.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "mechprint.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "intpoints.h"
#include "vector.h"
#include "matrix.h"
#include "mathem.h"
#include "gmatrix.h"
#include "tensor.h"
#include "vecttens.h"
#include "nssolver.h"

#include <math.h>
#include <string.h>


/**
  The function solves system of nonlinear algebraic equations
  by the arc-length method. This method was modified to reach the given value rlambda 
  of the lambda parameter. Rlambda value is reached with required error rlerr. 
  Only one right hand side vector is supported with respect to nonlinearity and absence of superposition. 
  Up to rlambda value the constant load vector is treated as proportional and then as constant. Proportional 
  vector grows until the nial number of steps is reached.

  @param lcid    - load case id
  @param rlambda - required value of lambda parameter
  @param rlerr   - required error between lambda and rlambda

  @return The function returns reached lambda parameter.

  Created by Tomas Koudelka,
*/
double arclengthrv1 (long lcid, double rlambda, double rlerr)
{
  long i,j,k,n,ni,ini,stop,modif,li,numr;
  double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,zero,ierr;
  double lambdao,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi;
  const char *mode;
  long back_dl = 1;
  long check_rv = 1;
  long nodiv_dl = 0;
  double bdl;
  
  
  
  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->nlman->nial;
  //  maximum number of iterations in one increment
  ini = Mp->nlman->niilal;
  //  computer zero
  zero = Mp->zero;
  //  required error in inner loop
  ierr = Mp->nlman->erral;
  //  length of arc
  dl = Mp->nlman->dlal;
  //  maximum length of arc
  dlmax = Mp->nlman->dlmaxal;
  //  minimum length of arc
  dlmin = Mp->nlman->dlminal;
  //  displacement-loading driving switch
  psi = Mp->nlman->psial;
  
  //  allocation of auxiliary arrays
  r = new double [n];
  ddr = new double [n];
  u = new double [n];
  v = new double [n];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  
  //  initialization phase
  ra = Lsrs->give_lhs (lcid);
  fp = Lsrs->give_rhs (lcid*2);
  fc = Lsrs->give_rhs (lcid*2+1);
  
  if (Mp->nlman->hdbackupal==2){
    arclopen (li,n,lambda,dl,ra,fp);    
  }
  else{
    lambda=0.0;
    lambdao=0.0;
    li=0;
  }
  
  if (li == 0)
    mode = "wt";
  else
    mode = "at";  
  
  
  for (j=0;j<n;j++)
  {
    if (lambda <= rlambda)
      fa[j]=lambda*(fc[j]+fp[j]);
    else
      fa[j]=rlambda*fc[j]+lambda*fp[j];
  }
  //  norm of proportionality vector
  norfp = ss(fp,fp,n);
  modif=0;


  // ***************************
  //  main iteration loop   ****
  // ***************************
  for (i=li;i<ni;i++) {

    //  backup of left hand side vector
    copyv (ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;
    if (back_dl)
      bdl = dl;

    fprintf (stdout,"\n arc-length: increment %ld   lambda %e  dl %e",i,lambda,dl);
/*    if (Mespr==1){
      fprintf (Out,"\n\n *******************************************************************");
      fprintf (Out,"\n arc-length: increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
    }*/
    
    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == li))
      stiffness_matrix (lcid);
    
    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    if (check_rv)
      addv(fc, fp, f, n);
    else
      copyv (fp, f, n);
    //  solution of K(r).v=F
    //Smat->solve_system (Gtm,v,f);
    Mp->ssle->solve_system (Gtm,Smat,v,f,Out);
    //  generalized norm of displacement increments
    norv = displincr (lcid, i, Mp->nlman, v, v, n);
    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    for (j=0;j<n;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      if (check_rv)
        fa[j]=(lambda+dlambda)*(fc[j]+fp[j]);
      else
        fa[j]=rlambda*fc[j]+(lambda+dlambda)*fp[j];
    }
    ddlambda=dlambda;
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    subv(fa, fi, f, n);    
    norf=normv(f,n);
    norfa = normv(fa, n);
    norf /= norfa;

    if (Mespr==1)  fprintf (stdout,"\n %e %e norf %e",lambda,dl,norf);
    
    if (norf<ierr){
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      modif++;

      if (modif>1) {
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  
          dl=dlmax;
	modif=0;
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl increase) dl=%e    norf=%e",dl,norf);
//	  fprintf (Out,"\n arc length must be modified (dl increase) dl=%e    norf=%e",dl,norf);
	}
      }
      lambda+=dlambda;
      if (check_rv && ((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
      {
        if (back_dl)
	{
          bdl = dl;
          back_dl = 0;
	}
        //  modification of the arc length
        dl/=2.0;
        if (dl<dlmin)
        {
          dl=dlmin;
          if (nodiv_dl)     
            break; 
          else
            nodiv_dl = 1;
        }
        if (Mespr==1)
        {
          fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
//          fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
        }
        //  restoring of left hand side vector
        for (j=0;j<n;j++)
          ra[j]=r[j];
        //  restoring of lambda parameter
        lambda=blambda;
      }
      if (fabs(lambda-rlambda) <= rlerr)
      {
        dl = bdl;
        check_rv = 0;
      }
      Mm->updateipval ();
        compute_req_val (lcid);
        print_step(lcid, i, lambda, fa);
        print_flush();
    }
    else{ //else
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){ //loop j
	//  back substitution
	//Smat->solve_system (Gtm,u,f);
	Mp->ssle->solve_system (Gtm,Smat,u,f,Out);
        copyv (ddr, f, n);
        addv(ddr, u, n);
	//  coefficient of quadratic equation
	quadeqcoeff (lcid, i, Mp->nlman, ddr, v, n, ddlambda, psi, norfp, dl, a0, a1, a2);
	/*
	  fprintf (Out,"\n\n\n Kontrola kvadraticke rovnice");
	  fprintf (Out,"\n norfp     %15.10le",norfp);
	  fprintf (Out,"\n norv      %15.10le",norv);
	  fprintf (Out,"\n (ddr,v)   %15.10le",ss(ddr,v,n));
	  fprintf (Out,"\n (ddr,ddr) %15.10le",ss(ddr,ddr,n));
	  fprintf (Out,"\n dl        %15.10le",dl);
	  fprintf (Out,"\n ddlambda  %15.10le",ddlambda);
	  fprintf (Out,"\n psi       %15.10le",psi);
	  fprintf (Out,"\n a2        %15.10le",a2);
	  fprintf (Out,"\n a1        %15.10le",a1);
	  fprintf (Out,"\n a0        %15.10le",a0);
	  fprintf (Out,"\n discrim   %15.10le",a1*a1-4.0*a2*a0);
	  fprintf (Out,"\n\n");
	*/
	
	//  solution of quadratic equation
        numr = solv_polynom_2(a2, a1, a0, l1, l2);
        switch (numr)
	{
	  case -1:
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
            break;
	  case 0:
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
            break;
    	  case 1:
            dlambda = l1;
            break;
	  default:
            break;
	}
	ss1=0.0;  ss2=0.0;
	ss3=0.0;  ss4=0.0;  ss5=0.0;
	for (k=0;k<n;k++){
	  ss1+=(ddr[k]+l1*v[k])*f[k];
	  ss2+=(ddr[k]+l2*v[k])*f[k];
	  ss3+=(ddr[k]+l1*v[k])*(ddr[k]+l1*v[k]);
	  ss4+=(ddr[k]+l2*v[k])*(ddr[k]+l2*v[k]);
	  ss5+=f[k]*f[k];
	}
	
	if (ss1/ss3/ss5>ss2/ss4/ss5)  
          dlambda=l1;
	else                          
          dlambda=l2;

	for (k=0;k<n;k++){
	  ddr[k]+=dlambda*v[k];
	  ra[k]+=u[k]+dlambda*v[k];
         if (check_rv)
            fa[k]+=dlambda*(fp[k]+fc[k]);
          else
            fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	
	//  computation of internal forces
	internal_forces (lcid,fi);
        subv(fa, fi, f, n);	
	norf=normv(f,n);
	norfa = normv(fa, n);
	norf /= norfa;

	if (Mespr==1){
	  fprintf (stdout,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
//	  fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	}
	if (norf<ierr){
	  lambda+=ddlambda;
          compute_req_val (lcid);
          print_step(lcid, i, lambda, fa);
          print_flush();
	  stop=1;
          if (check_rv && ((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
            stop=0;
          if (check_rv && (fabs(lambda-rlambda) <= rlerr))
	  {
            dl = bdl;
            check_rv = 0;
	  }
          if (stop > 0)
            Mm->updateipval();
          break; 
	}
      }

      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          if (nodiv_dl)     
            break; 
          else
            nodiv_dl = 1;
        }
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
//	  fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
	}
	//  restoring of left hand side vector
        copyv (r, ra, n);
	//  restoring of lambda parameter
	lambda=blambda;
      }
    }
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  if (Mp->nlman->hdbackupal==1)
    arclsave (i,n,lambda,dl,ra,fp);

  
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
  return (i);
}



/**
  Function solves system of nonlinear algebraic equations by Newton-Raphson method.
  Solved system does not contain time variable. Required value of the reached load coefficient lambda 
  can be specified.

  @param lcid - load case id
  @param rlambda - required value of reached load coefficient
  @param rlerr - tolerance for the reached rlambda

  @return The function does not return anything.

  Created by Tomas Koudelka, 16.8.2001
*/
void newton_raphsonrv1 (long lcid, double rlambda, double rlerr)
  //  function solves system of nonlinear algebraic
  //  equations by Newton-Raphson method
  //  solved system does not contain time variable
  //
  //  16.8.2001
{
  long i,j,k,n,ni,ini;
  double lambda,blambda,dlambda,dlambdamin,dlambdamax,zero,err,norf,norfa;
  double *r,*rb,*dr,*f,*fi,*fb,*fc,*fp;
  long back_incr = 1;
  long check_rv = 1;
  double binc;
  matrix lsm_a(3,3);
  vector lsm_r(3), lsm_l(3);
  Mp->nlman->divc_step = new double[10];
  Mp->nlman->divc_err  = new double[10];


  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->nlman->ninr;
  //  maximum number of iterations in one increment
  ini = Mp->nlman->niilnr;
  //  computer zero
  zero = Mp->zero;
  //  required error in inner loop
  err = Mp->nlman->errnr;
  //  increment size
  dlambda=Mp->nlman->incrnr;
  //  minimum increment size
  dlambdamin=Mp->nlman->minincrnr;
  //  maximum increment size
  dlambdamax=Mp->nlman->maxincrnr;

  rb = new double [n];
  f  = new double [n];
  fb = new double [n];
  fi = new double [n];
  dr = new double [n];
  memset (rb,0,n*sizeof(double));
  memset (f,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));
  memset (dr,0,n*sizeof(double));
  
  //  initialization phase
  r  = Lsrs->give_lhs (lcid);
  fp = Lsrs->give_rhs (lcid*2);
  fc = Lsrs->give_rhs (lcid*2+1);
  lambda=0.0;
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  print_init(-1, "wt");
  print_step(lcid, 0, 0.0, f);
  print_flush();
  for (i=0;i<ni;i++){
    
    //  indicator of strain computation
    //  it means, strains at integration points have not been computed
    Mp->strainstate=0;
    //  indicator of stress computation
    //  it means, stresses at integration points have not been computed
    Mp->stressstate=0;
    //  indicator of computation of other array
    //  it means, stresses at integration points have not been computed
    Mp->otherstate=0;


    //  backup of left hand side vector
    for (j=0;j<n;j++){
      rb[j]=r[j];
    }
    //  backup of reached lambda parameter
    blambda=lambda;
    if (back_incr)
      binc = dlambda;
   
    fprintf (stdout,"\n increment %ld   lambda %e,    dlambda %e",i,lambda,dlambda);
    
    //  vector of maximum load and vector of load increment
    for (j=0;j<n;j++){
/*      if (check_rv)
      {
        f[j]=lambda*(fp[j]+fc[j]);
        fb[j]=dlambda*(fp[j]+fc[j]);
      }
      else
      {
        f[j]=rlambda*fc[j]+lambda*fp[j];
        fb[j]=dlambda*fp[j];
      }*/
      if (check_rv)
      {
        f[j]=lambda*(fc[j]);
        fb[j]=dlambda*(fc[j]);
      }
      else
      {
        f[j]=rlambda*fc[j]+(lambda-rlambda)*fp[j];
        fb[j]=dlambda*fp[j];
      }
    }
    
    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == 0))
      stiffness_matrix (lcid);
    
    //  solution of K(r).v=F
    //Smat->solve_system (Gtm,dr,fb);
    Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

    for (j=0;j<n;j++){
      r[j]+=dr[j];
    }
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    
    //  vector of unbalanced forces
    for (j=0;j<n;j++){
/*      if (check_rv)
        fb[j]=f[j]+dlambda*(fp[j]+fc[j]);
      else
        fb[j]=f[j]+dlambda*fp[j];*/
      if (check_rv)
        fb[j]=f[j]+dlambda*(fc[j]);
      else
        fb[j]=f[j]+dlambda*fp[j];
    }
    norfa=normv(fb,n);
    for (j=0;j<n;j++){
      fb[j] -= fi[j];
    }
    //  norm of vector of unbalanced forces
    norf=normv(fb,n);
    norf /= norfa;
    
    if (norf<err){
      lambda+=dlambda;
      if (check_rv && ((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
      {
        if (back_incr)
	{
          binc = dlambda;
          back_incr = 0;
	}
        //  modification of the newton-rhapson
        dlambda/=2.0;
        if (dlambda<=dlambdamin)
        {
          if (dlambda == dlambdamin)
            break;
          dlambda=dlambdamin;
        }
        if (Mespr==1)
          fprintf (stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
        //  restoring of left hand side vector
        for (j=0;j<n;j++)
          r[j]=rb[j];
        //  restoring of lambda parameter
        lambda=blambda;
      }
      if (fabs(lambda-rlambda) <= rlerr)
      {
        dlambda = binc;
        check_rv = 0;
      }
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, f);
      print_flush();
      continue;
    }
    
    //  internal iteration loop
    fillm(0.0, lsm_a);
    fillv(0.0, lsm_r);
    for (j=0;j<ini;j++) {
      
      //  solution of K(r).v=F
      //Smat->solve_system (Gtm,dr,fb);
      Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

      for (k=0;k<n;k++){
	r[k]+=dr[k];
      }
      
      //  computation of internal forces
      internal_forces (lcid,fi);
      
      //  vector of unbalanced forces
      for (k=0;k<n;k++){
/*        if (check_rv)
          fb[k]=f[k]+dlambda*(fp[k]+fc[k])-fi[k];
        else
          fb[k]=f[k]+dlambda*fp[k]-fi[k];*/
        if (check_rv)
          fb[k]=f[k]+dlambda*(fc[k])-fi[k];
        else
          fb[k]=f[k]+dlambda*fp[k]-fi[k];
      }
      
      //  norm of vector of unbalanced forces
      norf=normv(fb,n);
      norf /= norfa;
      
      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf=%e  dlambda %e",i,j,norf,dlambda);
      
      if (norf<err)  break;

      // divergence detection with help of least square method
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
    }
    
    if (j==ini || norf>err){
      if (dlambda == dlambdamin)
      {
	fprintf (stdout,"\n increment of lambda cannot be decreased, iteration will be stopped\n");
	break;
      }
        
      dlambda/=2.0;
      if (dlambda<dlambdamin)  dlambda=dlambdamin;

      if (Mespr==1)
	fprintf (stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
      
      for (j=0;j<n;j++){
      	r[j]=rb[j];
      }
      lambda=blambda;
    }
    else{
      lambda+=dlambda;
      if (check_rv && ((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
      {
        if (back_incr)
	{
          binc = dlambda;
          back_incr = 0;
	}
        //  modification of the newton-rhapson
        dlambda/=2.0;
        if (dlambda<=dlambdamin)
        {
          if (dlambda == dlambdamin)
            break;
          dlambda=dlambdamin;
        }
        if (Mespr==1)
          fprintf (stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
        //  restoring of left hand side vector
        for (j=0;j<n;j++)
          r[j]=rb[j];
        //  restoring of lambda parameter
        lambda=blambda;
      }
      if (fabs(lambda-rlambda) <= rlerr)
      {
        if (back_incr == 0)
          dlambda = binc;
        check_rv = 0;
      }
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, f);
      print_flush();
    }
    
  }
  
  delete [] dr;  delete [] fi;  delete [] fb;  delete [] f;
  print_close();
}
