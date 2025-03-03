#include "nfssolver.h"
#include "nssolver.h"
#include "global.h"
#include "probdesc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "loadcase.h"
#include "gmatrix.h"
#include "mechprint.h"
#include "flsubdom.h"
#include "mathem.h"


/**
   function solves problem with floating subdomains
   incremental analysis is used, compliance can vary
   
   it may be used for pull-out tests
   
   JK, 23. 11. 2012
*/
void solve_incremental_floating_subdomains ()
{
  long i,n,nm,lcid,nincr;
  double proppar;
  double *lhs,*rhs,*d,*lambda;
  
  //  load case id, must be zero
  lcid=0;
  //  the number of increments
  nincr=Mp->nincr;
  //  the number of primal DOFs
  n=Ndofm;
  //  the number of dual DOFs - Lagrange multipliers
  nm=Mp->ssle->feti.ndofcp;
  //  proportional parameter
  proppar=1.0/nincr;
  
  //  array for total displacements
  d = new double [n];
  fillv (0.0,d,n);
  //  array for total Lagrange multipliers
  lambda = new double [nm];
  fillv (0.0,lambda,nm);
  //  array for increments of displacements, array is allocated in Lsrs
  lhs = Lsrs->give_lhs (lcid);
  //  array for increments of load, array is allocated in Lsrs
  rhs = Lsrs->give_rhs (lcid);

  //  assembling of the right hand side
  mefel_right_hand_side (lcid,rhs);
  //  computation of load increment
  cmulv (proppar,rhs,n);
  
  //  assembling of the stiffness matrix
  stiffness_matrix (lcid);
  
  Mp->ssle->feti.matrices_assembl (Smat,Out);
  Mp->ssle->feti.vectors_assembl (rhs,Out);
  
  print_init(-1, "wt");    
  
  //  loop over the number of increments
  for (i=0;i<nincr;i++){
    
    Mp->ssle->feti.define_b (Out);
    Mp->ssle->feti.define_h (Out);
    
    //  solution of equation system
    Mp->ssle->feti.mpcg (Out);
    Mp->ssle->feti.nodalunknowns (d,Out);
    
    //  total displacements are computed
    addv (lhs,d,n);
    
    print_init(-1, "wt");    
    //  computes and prints required quantities
    compute_req_val (lcid);
    print_step(lcid, i+1, 0.0, NULL);    
    
  }//  end of the loop over the number of increments
  
  for (i=0;i<Lsrs->nlc;i++){
    //  computes and prints required quantities
    //compute_ipstrains (i);
    //compute_ipstresses (i);
    compute_req_val (i);
    print_step(i, 0, 0.0, NULL);    
  }
  print_close();

}


/**
   function solves linear problem with floating subdomains
   it may be used for pull-out tests
   
   JK, 5.11.2005
*/
  /*

void solve_incremental_floating_subdomains ()
{
  long i,lcid,n,nincr;
  double *d,*lhs,*rhs;
  
  //  load case id, must be zero
  lcid=0;
  //  the number of unknowns
  n=Ndofm;
  //  the number of increments
  nincr=Mp->nincr;

  Fsd = new flsubdom;
  
  //  array of nodal displacements
  lhs=Lsrs->give_lhs (lcid);
  //  array of nodal forcess
  rhs=Lsrs->give_rhs (lcid);
  //  array of increments of nodal displacements
  d = new double [n];
  
  //  assembling of the right hand side
  mefel_right_hand_side (lcid,rhs);
  
  //  assembling of the stiffness matrix
  stiffness_matrix (lcid);
  
  Fsd->initialization (Ndofm,Mp->ense,rhs);
  
  print_init(-1, "wt");    
  for (i=0;i<nincr;i++){
    
    Fsd->solve_lin_alg_system (lhs,rhs);
    
    addv (lhs,d,n);
    
    Fsd->add_mult (1.0);
    Fsd->mult_correction ();
    
    for (i=0;i<Lsrs->nlc;i++){
      //compute_ipstrains (i);
      //compute_ipstresses (i);
      compute_req_val (i);
      print_step(i, 0, 0.0, NULL);    
    }
    
  }
  
  print_close();
  

}
  */  

/**
   function solves nonlinear problem of floating subdomain
   the nonlinearity is hidden in contact stresses and in
   material models, like plasticity, etc.
   
   function solves the problem by the arc-length method (continuation method)
   this function is modification of the function arclength
   from the file nssolver.cpp

   d stands for delta
   dd stands for capital delta
   
   fc - nonproportional %vector
   fp - proportional %vector
   n - number of unknowns
   
   JK, 24.6.2006
*/
/*
void solve_nonlinear_floating_subdomains ()
{
  long i,j,k,n,l,ni,ini,stop,modif,li,newmesh,numr,lcid,nlm;
  double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,zero,ierr;
  double lambdao,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi,*rhs,*blm;
  
  //  allocation of object floating subdomains
  Fsd = new flsubdom;
  
  //  only one load case can be solved
  //  it is nonlinear computation and the superposition is useless
  lcid=0;
  
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
  
  //  array for right hand side
  rhs=Lsrs->give_rhs (0);
  //  assembling of right hand side (load vector)
  mefel_right_hand_side (lcid,rhs);
  
  dlambda=0.0;
  
  //  allocation of auxiliary arrays
  r = new double [n];
  ddr = new double [n];
  u = new double [n];
  v = new double [n];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  
  //  initialization phase
  //  vector of nodal displacements
  ra = Lsrs->give_lhs (lcid);
  //  vector of proportional load
  fp = Lsrs->give_rhs (lcid*2);
  //  vector of nonproportional load
  fc = Lsrs->give_rhs (lcid*2+1);
  
  
  if (Mp->nlman->hdbackupal == 2){            // termit
    arclopen (li,n,lambda,dl,ra,fp);
  }
  else{
    lambda=0.0;
    lambdao=0.0;
    li=0;
  }
  
  
  //  norm of proportionality vector
  norfp = ss(fp,fp,n);
  modif=0;
  
  
  // assembling reached load vector
  for (j=0;j<n;j++){
    fa[j]=fc[j]+(lambda+dlambda)*fp[j];
  }
  if (li)
    print_init(-1, "at");
  else
    {
      print_init(-1, "wt");
      print_step(lcid, li, lambda, fa);
    }
  
  
  //  assembling of stiffness matrix
  stiffness_matrix (lcid);
  
  //  initialization of FETI method
  //  determination of number of Lagrange multipliers
  //  computation of rigid body motions
  Fsd->initialization (Ndofm,Mp->ense,f);

  //  number of Lagrange multipliers
  nlm = Fsd->nlm;
  
  
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  for (i=li;i<ni;i++){
    
    fprintf (Out,"\n\n arc-length  prirustek %ld",i);
    
    stop=1; // termit
    //  backup of left hand side vector
    copyv (ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;
    //  backup of reached Lagrange multipliers
    copyv (Fsd->tw,blm,nlm);

    //fprintf (stdout,"\n\n arc-length: increment %ld   lambda %e  dl %e",i,lambda,dl);
//    if (Mespr==1){
//  fprintf (Out,"\n\n *******************************************************************");
//  fprintf (Out,"\n arc-length: increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
//  }
  
    
    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == li)){
      //stiffness_matrix (lcid);
      //Fsd->initialization (Ndofm,Mp->ense,f);
    }
    
    stiffness_matrix (lcid);
    

    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    copyv(fp, f, n);
    
    Fsd->initialization (Ndofm,Mp->ense,f);
    
    Fsd->solve_lin_alg_system (v,f);
    
    //  solution of K(r).v=F
    //Smat->solve_system (v,f);
    
    //  generalized norm of displacement increments
    norv = displincr (v,n);
    
    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    Fsd->add_mult (dlambda);
    Fsd->mult_correction ();
    
    //  novinka
    if (Fsd->nch!=0){
      for (j=0;j<ini;j++){

	fprintf (stdout,"\n cyklime v casti jedna  %ld",j);

	//  restoring of Lagrange multipliers
	copyv (blm,Fsd->tw,nlm);
	
	Fsd->solve_lin_alg_system (v,f);
	
	//  generalized norm of displacement increments
	norv = displincr (v,n);
	
	//  compute new dlambda increment
	dlambda = dl/sqrt(norv+psi*psi*norfp);
	
	Fsd->add_mult (dlambda);
	Fsd->mult_correction ();
	
	if (Fsd->nch==0)
	  break;
      }
    }
    //  konec novinky
    
    
    for (j=0;j<n;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
    }
    ddlambda=dlambda;
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    subv(fa, fi, f, n);    
    norf=sizev(f,n);
    norfa = sizev(fa, n);
    //if (norfa<1.0e-8){}
    //else norf /= norfa;
    norf /= norfa;
    
    
    //if (Mespr==1)  fprintf (stdout,"\n ddlambda %e    dlambda %e   norf %e  ierr %e",ddlambda,dlambda,norf,ierr);
    if (Mespr==1)  fprintf (stdout,"\n increment %ld     norv %e      norf %e",i,norv,norf);
    
    if (norf<ierr){
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      modif++;
      
      if (modif>1){
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  
          dl=dlmax;
	modif=0;
	if (Mespr==1){
	  fprintf (stdout,"\n arc-length must be modified (dl increase) dl=%e",dl);
	  //	  fprintf (Out,"\n arc-length must be modified (dl increase) dl=%e",dl);
	}
      }
      lambda+=dlambda;
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, fa);
      print_flush();
    }
    else{
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){
	if (Mp->nlman->stmat==tangent_stiff){
	  stiffness_matrix (lcid);
	  Fsd->initialization (Ndofm,Mp->ense,f);
	}
	
	fprintf (Out,"\n arc-length  vnitrni smycka  %ld %ld",i,j);
	
	
        copyv(f, buf, n);
	
	Fsd->solve_lin_alg_system (u,f);
	
	//  back substitution
	//Smat->solve_system (u,f);
	
	//  backup of the vector ddr
        copyv(ddr, buddr, n);
	
	
        copyv(ddr, f, n);
        addv(ddr, u, n);
	//  coefficient of quadratic equation
	quadeqcoeff (ddr,v,n,ddlambda,psi,norfp,dl,a0,a1,a2);
	
	//  solution of quadratic equation
        numr = solv_polynom_2(a2, a1, a0, l1, l2);
        switch (numr)
	{
	  case -1 :
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
            break;
	  case 0 :
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
            break;
    	  case 1 :
            dlambda = l1;
            break;
	  default :
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
	
	if (fabs(l1)>fabs(l2))  dlambda=l2;
	else dlambda=l1;
	
	copyv (Fsd->tw,bulm,nlm);

	Fsd->add_mult (dlambda);
	Fsd->mult_correction ();
	
	//  novinka
	if (Fsd->nch!=0){
	  for (l=0;l<ini;l++){
	    
	    fprintf (stdout,"\n cyklime v casti dve  %ld",j);
	    
	    
	    //  restoring of Lagrange multipliers
	    copyv (bulm,Fsd->tw,nlm);
	    copyv (buddr,ddr,n);
	    copyv (buf,f,n);
	    
	    Fsd->solve_lin_alg_system (u,f);
	    
	    copyv(ddr, f, n);
	    addv(ddr, u, n);
	    //  coefficient of quadratic equation
	    quadeqcoeff (ddr,v,n,ddlambda,psi,norfp,dl,a0,a1,a2);
	    
	    //  solution of quadratic equation
	    numr = solv_polynom_2(a2, a1, a0, l1, l2);
	    switch (numr)
	      {
	      case -1 :
		fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
		break;
	      case 0 :
		fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
		break;
	      case 1 :
		dlambda = l1;
		break;
	      default :
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
	    
	    if (fabs(l1)>fabs(l2))  dlambda=l2;
	    else dlambda=l1;
	    
	    
	    Fsd->add_mult (dlambda);
	    Fsd->mult_correction ();
	    
	    if (Fsd->nch==0)
	      break;
	  }
	}
	//  konec novinky




	for (k=0;k<n;k++){
	  ddr[k]+=dlambda*v[k];
	  ra[k]+=u[k]+dlambda*v[k];
	  fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	
	//fprintf (stdout,"   ddlambda %e",ddlambda);
	
	//fprintf (Out,"\n ddlambda %e     dlambda %e",ddlambda,dlambda);
	//  computation of internal forces
	internal_forces (lcid,fi);
        subv(fa, fi, f, n);	
	norf=sizev(f,n);
	norfa = sizev(fa, n);
	norf /= norfa;
	//if (norfa<1.0e-8){}
	//else norf /= norfa;
	
	if (Mespr==1){
	  fprintf (stdout,"\n    norf=%e  ierr=%e",norf,ierr);
	  //fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	}

	if (norf<ierr){
	  lambda+=ddlambda;
	  Mm->updateipval ();
	  stop=1; 
	  compute_req_val (lcid);
          print_step(lcid, i+1, lambda, fa);
          print_flush();
	  break;
	}
      }
      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          break; 
        }
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e",dl);
	  //fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e",dl);
	}
	//  restoring of left hand side vector
        copyv(r, ra, n);
	//  restoring of lambda parameter
	lambda=blambda;
	//  restoring of Lagrange multipliers
	copyv (blm,Fsd->tw,nlm);
      }
    }
    
    fprintf (stdout,"\n increment %ld    total lambda %e",i,lambda);

    if (stop==0)
      continue;
    
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  
  fprintf (stdout,"\n\n multiplikatory   %lf",Fsd->tw[0]);
  
  if (Mp->nlman->hdbackupal==1)
    arclsave (i,n,lambda,dl,ra,fp);

  print_close();
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
}
*/


/**
   function solves nonlinear problem of floating subdomain
   the nonlinearity is hidden in contact stresses and in
   material models, like plasticity, etc.
   
   function solves the problem by the arc-length method (continuation method)
   this function is modification of the function arclength
   from the file nssolver.cpp

   d stands for delta
   dd stands for capital delta
   
   fc - nonproportional %vector
   fp - proportional %vector
   n - number of unknowns
   
   JK, 20.6.2006 (25.7.2001)
*/
/*
void solve_nonlinear_floating_subdomains_24_6 ()
{
  long i,j,k,n,l,ni,ini,stop,modif,li,numr,lcid,nlm;
  double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,zero,ierr;
  double lambdao,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi,*rhs,*blm;
  
  Fsd = new flsubdom;
  
  lcid=0;

  rhs=Lsrs->give_rhs (0);
  mefel_right_hand_side (lcid,rhs);
  
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
  
  dlambda=0.0;
  
  //  allocation of auxiliary arrays
  r = new double [n];
  ddr = new double [n];
  u = new double [n];
  v = new double [n];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  
  //  initialization phase
  //  vector of nodal displacements
  ra = Lsrs->give_lhs (lcid);
  //  vector of proportional load
  fp = Lsrs->give_rhs (lcid*2);
  //  vector of nonproportional load
  fc = Lsrs->give_rhs (lcid*2+1);
  
  
  if (Mp->nlman->hdbackupal == 2){            // termit
    arclopen (li,n,lambda,dl,ra,fp);
  }
  else{
    lambda=0.0;
    lambdao=0.0;
    li=0;
  }
  
  
  //  norm of proportionality vector
  norfp = ss(fp,fp,n);
  modif=0;
  
  
  // assembling reached load vector
  for (j=0;j<n;j++){
    fa[j]=fc[j]+(lambda+dlambda)*fp[j];
  }
  if (li)
    print_init(-1, "at");
  else
    {
      print_init(-1, "wt");
      print_step(lcid, li, lambda, fa);
    }
  
  
  
  stiffness_matrix (lcid);
  Fsd->initialization (Ndofm,Mp->ense,f);
  //  number of Lagrange multipliers
  nlm = Fsd->nlm;
  
  blm = new double [nlm];

  double *buddr;
  buddr = new double [n];
  double *buf;
  buf = new double [n];
  double *bulm;
  bulm = new double [nlm];

  // ***************************
  //  main iteration loop  ****
  // ***************************
  for (i=li;i<ni;i++){
    
    fprintf (Out,"\n\n arc-length  prirustek %ld",i);
    
    stop=1; // termit
    //  backup of left hand side vector
    copyv (ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;
    //  backup of reached Lagrange multipliers
    copyv (Fsd->tw,blm,nlm);

    //fprintf (stdout,"\n\n arc-length: increment %ld   lambda %e  dl %e",i,lambda,dl);
        if (Mespr==1){
	  fprintf (Out,"\n\n *******************************************************************");
	  fprintf (Out,"\n arc-length: increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
	  }
    
    
    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == li)){
      //stiffness_matrix (lcid);
      //Fsd->initialization (Ndofm,Mp->ense,f);
    }
    
    stiffness_matrix (lcid);
    

    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    copyv(fp, f, n);
    
    Fsd->initialization (Ndofm,Mp->ense,f);
    
    Fsd->solve_lin_alg_system (v,f);
    
    //  solution of K(r).v=F
    //Smat->solve_system (v,f);
    
    //  generalized norm of displacement increments
    norv = displincr (v,n);
    
    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    Fsd->add_mult (dlambda);
    Fsd->mult_correction ();
    
    //  novinka
    if (Fsd->nch!=0){
      for (j=0;j<ini;j++){

	fprintf (stdout,"\n cyklime v casti jedna  %ld",j);

	//  restoring of Lagrange multipliers
	copyv (blm,Fsd->tw,nlm);
	
	Fsd->solve_lin_alg_system (v,f);
	
	//  generalized norm of displacement increments
	norv = displincr (v,n);
	
	//  compute new dlambda increment
	dlambda = dl/sqrt(norv+psi*psi*norfp);
	
	Fsd->add_mult (dlambda);
	Fsd->mult_correction ();
	
	if (Fsd->nch==0)
	  break;
      }
    }
    //  konec novinky
    
    
    for (j=0;j<n;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
    }
    ddlambda=dlambda;
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    subv(fa, fi, f, n);    
    norf=sizev(f,n);
    norfa = sizev(fa, n);
    //if (norfa<1.0e-8){}
    //else norf /= norfa;
    norf /= norfa;
    
    
    //if (Mespr==1)  fprintf (stdout,"\n ddlambda %e    dlambda %e   norf %e  ierr %e",ddlambda,dlambda,norf,ierr);
    if (Mespr==1)  fprintf (stdout,"\n increment %ld     norv %e      norf %e",i,norv,norf);
    
    if (norf<ierr){
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      modif++;
      
      if (modif>1){
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  
          dl=dlmax;
	modif=0;
	if (Mespr==1){
	  fprintf (stdout,"\n arc-length must be modified (dl increase) dl=%e",dl);
	  //	  fprintf (Out,"\n arc-length must be modified (dl increase) dl=%e",dl);
	}
      }
      lambda+=dlambda;
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, fa);
      print_flush();
    }
    else{
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){
	if (Mp->nlman->stmat==tangent_stiff){
	  stiffness_matrix (lcid);
	  Fsd->initialization (Ndofm,Mp->ense,f);
	}
	
	fprintf (Out,"\n arc-length  vnitrni smycka  %ld %ld",i,j);
	
	
        copyv(f, buf, n);
	
	Fsd->solve_lin_alg_system (u,f);
	
	//  back substitution
	//Smat->solve_system (u,f);
	
	//  backup of the vector ddr
        copyv(ddr, buddr, n);
	
	
        copyv(ddr, f, n);
        addv(ddr, u, n);
	//  coefficient of quadratic equation
	quadeqcoeff (ddr,v,n,ddlambda,psi,norfp,dl,a0,a1,a2);
	
	//  solution of quadratic equation
        numr = solv_polynom_2(a2, a1, a0, l1, l2);
        switch (numr)
	{
	  case -1 :
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
            break;
	  case 0 :
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
            break;
    	  case 1 :
            dlambda = l1;
            break;
	  default :
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
	
	if (fabs(l1)>fabs(l2))  dlambda=l2;
	else dlambda=l1;
	
	copyv (Fsd->tw,bulm,nlm);

	Fsd->add_mult (dlambda);
	Fsd->mult_correction ();
	
	//  novinka
	if (Fsd->nch!=0){
	  for (l=0;l<ini;l++){
	    
	    fprintf (stdout,"\n cyklime v casti dve  %ld",j);
	    
	    
	    //  restoring of Lagrange multipliers
	    copyv (bulm,Fsd->tw,nlm);
	    copyv (buddr,ddr,n);
	    copyv (buf,f,n);
	    
	    Fsd->solve_lin_alg_system (u,f);
	    
	    copyv(ddr, f, n);
	    addv(ddr, u, n);
	    //  coefficient of quadratic equation
	    quadeqcoeff (ddr,v,n,ddlambda,psi,norfp,dl,a0,a1,a2);
	    
	    //  solution of quadratic equation
	    numr = solv_polynom_2(a2, a1, a0, l1, l2);
	    switch (numr)
	      {
	      case -1 :
		fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
		break;
	      case 0 :
		fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
		break;
	      case 1 :
		dlambda = l1;
		break;
	      default :
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
	    
	    if (fabs(l1)>fabs(l2))  dlambda=l2;
	    else dlambda=l1;
	    
	    
	    Fsd->add_mult (dlambda);
	    Fsd->mult_correction ();
	    
	    if (Fsd->nch==0)
	      break;
	  }
	}
	//  konec novinky




	for (k=0;k<n;k++){
	  ddr[k]+=dlambda*v[k];
	  ra[k]+=u[k]+dlambda*v[k];
	  fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	
	//fprintf (stdout,"   ddlambda %e",ddlambda);
	
	//fprintf (Out,"\n ddlambda %e     dlambda %e",ddlambda,dlambda);
	//  computation of internal forces
	internal_forces (lcid,fi);
        subv(fa, fi, f, n);	
	norf=sizev(f,n);
	norfa = sizev(fa, n);
	norf /= norfa;
	//if (norfa<1.0e-8){}
	//else norf /= norfa;
	
	if (Mespr==1){
	  fprintf (stdout,"\n    norf=%e  ierr=%e",norf,ierr);
	  //fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	}

	if (norf<ierr){
	  lambda+=ddlambda;
	  Mm->updateipval ();
	  stop=1; 
	  compute_req_val (lcid);
          print_step(lcid, i+1, lambda, fa);
          print_flush();
	  break;
	}
      }
      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          break; 
        }
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e",dl);
	  //fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e",dl);
	}
	//  restoring of left hand side vector
        copyv(r, ra, n);
	//  restoring of lambda parameter
	lambda=blambda;
	//  restoring of Lagrange multipliers
	copyv (blm,Fsd->tw,nlm);
      }
    }
    
    fprintf (stdout,"\n increment %ld    total lambda %e",i,lambda);

    if (stop==0)
      continue;
    
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  
  fprintf (stdout,"\n\n multiplikatory   %lf",Fsd->tw[0]);
  
  if (Mp->nlman->hdbackupal==1)
    arclsave (i,n,lambda,dl,ra,fp);

  print_close();
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
}
*/


/**
   JK, 26.6.2006
*/
/*
void solve_nonlinear_floating_subdomains_29_6 ()
{
  long i,j,k,n,lcid,ni,ini,stop,modif,li,numr,nlm;
  double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,zero,ierr;
  double lambdao;//,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi,*d,*rhs;
  
  
  //  load case id must be equal to zero
  //  only one load case can be solved
  //  one load case may contain several proportional vectors
  lcid = 0;
  
  //  assembling of stiffness matrix
  stiffness_matrix (lcid);
  
  //  number of unknown displacements
  n = Ndofm;
  
  //  auxiliary assigment of pointer
  rhs=Lsrs->give_rhs (lcid);
  
  //  assembling of right hand side vector (vector of nodal forces)
  mefel_right_hand_side (lcid,rhs);
  
  //  vector of nodal displacements
  d   = Lsrs->give_lhs (lcid);
  //  vector of proportional load
  fp = Lsrs->give_rhs (lcid*2);
  //  vector of constant load
  fc = Lsrs->give_rhs (lcid*2+1);
  
  
  
  //  allocation of object floating subdomains
  Fsd = new flsubdom;
  
  //  computation of rigid body motions
  //  list of dependent equations
  //  determination of number of Lagrange multipliers
  Fsd->initialization (Ndofm,Mp->ense,fp);
  
  //  number of Lagrange multipliers
  nlm = Fsd->nlm;
  
  
  
  // *******************************
  //  data about arc-length solver
  // *******************************

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
  
  
  // *********************************
  //  allocation of auxiliary arrays
  // *********************************
  r = new double [nlm];
  ra = new double [nlm];
  ddr = new double [nlm];
  u = new double [nlm];
  v = new double [nlm];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  
  
  lambda=0.0;
  lambdao=0.0;
  dlambda=0.0;
  li=0;
  modif=0;  
  
  
  
  
  
  //  norm of proportionality vector
  norfp = ss(fp,fp,n);
  
  
  
  // assembling reached load vector
  for (j=0;j<n;j++){
    fa[j]=fc[j]+(lambda+dlambda)*fp[j];
  }
  if (li)
    print_init(-1, "at");
  else
    {
      print_init(-1, "wt");
      print_step(lcid, li, lambda, fa);
    }
  
  
  
  // ***************************
  //  main iteration loop   ****
  // ***************************
  for (i=li;i<ni;i++){
    
    fprintf (Out,"\n\n arc-length  prirustek %ld",i);
    

    stop=1; // termit
    //  backup of left hand side vector
    copyv (ra, r, nlm);
    //  backup of reached lambda parameter
    blambda=lambda;

    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    copyv(fp, f, n);

    //  solution of K(r).v=F
    Fsd->pmcg (f,Out);
    copyv (Fsd->w,v,nlm);

    //  generalized norm of displacement increments
    norv = displincr (v,nlm);

    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    for (j=0;j<nlm;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
    }
    for (j=0;j<n;j++){
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
    }
    ddlambda=dlambda;
    
    Fsd->add_mult (dlambda);
    Fsd->mult_correction ();

    Fsd->nonlinlagrmultdispl (d,fa);

    
    //  computation of internal forces
    internal_forces (lcid,fi);
    subv(fa, fi, f, n);    
    norf = sizev(f,n);
    norfa = sizev(fa, n);
    norf /= norfa;
    
    
    if (Mespr==1)  fprintf (stdout,"\n increment %ld     norv %e      norf %e",i,norv,norf);
    
    if (norf<ierr){
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      modif++;

      if (modif>1){
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  
          dl=dlmax;
	modif=0;
	if (Mespr==1){
	  fprintf (stdout,"\n arc-length must be modified (dl increase) dl=%e",dl);
//	  fprintf (Out,"\n arc-length must be modified (dl increase) dl=%e",dl);
	}
      }
      lambda+=dlambda;
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, fa);
      print_flush();
    }
    else{
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){
	
	fprintf (Out,"\n arc-length  vnitrni smycka  %ld %ld",i,j);
	
	//  solution of K(r).v=F
	Fsd->pmcg (f,Out);
	copyv (Fsd->w,u,nlm);
	
	
        //copyv(ddr, f, nlm);
        addv(ddr, u, nlm);
	//  coefficient of quadratic equation
	//quadeqcoeff (ddr,v,n,ddlambda,psi,norfp,dl,a0,a1,a2);
	quadeqcoeff (ddr,v,nlm,ddlambda,psi,norfp,dl,a0,a1,a2);
	
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

	
	//  solution of quadratic equation
        numr = solv_polynom_2(a2, a1, a0, l1, l2);
        switch (numr)
	{
	  case -1 :
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
            break;
	  case 0 :
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
            break;
    	  case 1 :
            dlambda = l1;
            break;
	  default :
            break;
	}

	if (fabs(l1)>fabs(l2))  dlambda=l2;
	else dlambda=l1;

	for (k=0;k<nlm;k++){
	  ddr[k]+=dlambda*v[k];
	  ra[k]+=u[k]+dlambda*v[k];
	}
	for (k=0;k<n;k++){
	  fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	
	Fsd->add_mult (dlambda);
	Fsd->mult_correction ();
	
	Fsd->nonlinlagrmultdispl (d,fa);
	
	//  computation of internal forces
	internal_forces (lcid,fi);
        subv(fa, fi, f, n);	
	norf=sizev(f,n);
	norfa = sizev(fa, n);
	norf /= norfa;
	
	if (Mespr==1){
	  fprintf (stdout,"\n    norf=%e  ierr=%e",norf,ierr);
	  //fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	}

	if (norf<ierr){
	  lambda+=ddlambda;
	  Mm->updateipval ();
	  stop=1; 
	  compute_req_val (lcid);
          print_step(lcid, i+1, lambda, fa);
          print_flush();
	  break;
	}
      }
      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          break; 
        }
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e",dl);
	  //fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e",dl);
	}
	//  restoring of left hand side vector
        copyv(r, ra, nlm);
	//  restoring of lambda parameter
	lambda=blambda;
      }
    }
    
    fprintf (stdout,"\n increment %ld    total lambda %e",i,lambda);
    
    if (stop==0)
      continue;
    
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  
  print_close();
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
}
*/

/*
void solve_nonlinear_floating_subdomains ()
{
  long i,j,k,n,lcid,nlm,ni,ini,stop,modif,li,newmesh, numr;
  double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,zero,ierr;
  double lambdao,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi,*rhs;
  
  
  //  load case id must be equal to zero
  //  only one load case can be solved
  //  one load case may contain several proportional vectors
  lcid = 0;
  
  //  assembling of stiffness matrix
  stiffness_matrix (lcid);
  
  //  number of unknown displacements
  n = Ndofm;
  
  //  auxiliary assigment of pointer
  rhs=Lsrs->give_rhs (lcid);
  
  //  assembling of right hand side vector (vector of nodal forces)
  mefel_right_hand_side (lcid,rhs);
  
  //  vector of nodal displacements
  ra = Lsrs->give_lhs (lcid);
  //  vector of proportional load
  fp = Lsrs->give_rhs (lcid*2);
  //  vector of constant load
  fc = Lsrs->give_rhs (lcid*2+1);
  
  
  
  //  allocation of object floating subdomains
  Fsd = new flsubdom;
  
  //  computation of rigid body motions
  //  list of dependent equations
  //  determination of number of Lagrange multipliers
  Fsd->initialization (Ndofm,Mp->ense,fp);
  
  //  number of Lagrange multipliers
  nlm = Fsd->nlm;
  
  
  
  // *******************************
  //  data about arc-length solver
  // *******************************

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
  
  
  dlambda=0.0;
  lambda=0.0;
  lambdao=0.0;
  li=0;

  
  //  norm of proportionality vector
  norfp = ss(fp,fp,n);
  modif=0;


  // assembling reached load vector
  for (j=0;j<n;j++){
    fa[j]=fc[j]+(lambda+dlambda)*fp[j];
  }
  if (li)
    print_init(-1, "at");
  else
  {
    print_init(-1, "wt");
    print_step(lcid, li, lambda, fa);
  }
  

  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  for (i=li;i<ni;i++){

    fprintf (Out,"\n\n arc-length  prirustek %ld",i);


    stop=1; // termit
    //  backup of left hand side vector
    copyv (ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;

    
    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    copyv(fp, f, n);
    
    
    stiffness_matrix (lcid);
    Fsd->factorize ();


    //  solution of K(r).v=F
    Fsd->pmcg (f,Out);
    
    Fsd->lagrmultdispl (v,f);
    
    
    //  generalized norm of displacement increments
    norv = displincr (v,n);

    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    //if (norv+psi*psi*norfp < 1.0e-8)  dlambda=0.0;
    //else dlambda = dl/sqrt(norv+psi*psi*norfp);
    
    for (j=0;j<n;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
    }
    ddlambda=dlambda;
    
    Fsd->add_mult (dlambda);
    Fsd->mult_correction ();
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    subv(fa, fi, f, n);    
    norf=sizev(f,n);
    norfa = sizev(fa, n);
    //if (norfa<1.0e-8){}
    //else norf /= norfa;
    norf /= norfa;

    
    //if (Mespr==1)  fprintf (stdout,"\n ddlambda %e    dlambda %e   norf %e  ierr %e",ddlambda,dlambda,norf,ierr);
    if (Mespr==1)  fprintf (stdout,"\n increment %ld     norv %e      norf %e",i,norv,norf);
    
    if (norf<ierr){
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      modif++;

      if (modif>1){
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  
          dl=dlmax;
	modif=0;
	if (Mespr==1){
	  fprintf (stdout,"\n arc-length must be modified (dl increase) dl=%e",dl);
//	  fprintf (Out,"\n arc-length must be modified (dl increase) dl=%e",dl);
	}
      }
      lambda+=dlambda;
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, fa);
      print_flush();
    }
    else{
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){


	fprintf (Out,"\n arc-length  vnitrni smycka  %ld %ld",i,j);
	fprintf (stdout,"\n arc-length  vnitrni smycka  %ld %ld",i,j);
	
	
	stiffness_matrix (lcid);
	Fsd->factorize ();
	
	
	//  back substitution
	Fsd->pmcg (f,Out);
	Fsd->lagrmultdispl (u,f);
	
        copyv(ddr, f, n);
        addv(ddr, u, n);
	//  coefficient of quadratic equation
	quadeqcoeff (ddr,v,n,ddlambda,psi,norfp,dl,a0,a1,a2);
	
	//  solution of quadratic equation
        numr = solv_polynom_2(a2, a1, a0, l1, l2);
        switch (numr)
	{
	  case -1 :
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
            break;
	  case 0 :
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
            break;
    	  case 1 :
            dlambda = l1;
            break;
	  default :
            break;
	}

	if (fabs(l1)>fabs(l2))  dlambda=l2;
	else dlambda=l1;

	//fprintf (stdout,"\n increment  %ld     inner loop  %ld   x1=%e x2=%e   dlambda %e",i,j,l1,l2,dlambda);
	//fprintf (stdout,"\n increment  %ld     inner loop  %ld   x1=%e x2=%e   ddlambda %e",i,j,l1,l2,ddlambda);

	for (k=0;k<n;k++){
	  ddr[k]+=dlambda*v[k];
	  ra[k]+=u[k]+dlambda*v[k];
	  fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	
	Fsd->add_mult (dlambda);
	Fsd->mult_correction ();
	
	
	//  computation of internal forces
	internal_forces (lcid,fi);
        subv(fa, fi, f, n);	
	norf=sizev(f,n);
	norfa = sizev(fa, n);
	norf /= norfa;
	
	if (Mespr==1){
	  fprintf (stdout,"\n    norf=%e  ierr=%e",norf,ierr);
	  //fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	}

	if (norf<ierr){
	  lambda+=ddlambda;
	  Mm->updateipval ();
	  stop=1; 
	  compute_req_val (lcid);
          print_step(lcid, i+1, lambda, fa);
          print_flush();
	  break;
	}
      }
      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          break; 
        }
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e",dl);
	  //fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e",dl);
	}
	//  restoring of left hand side vector
        copyv(r, ra, n);
	//  restoring of lambda parameter
	lambda=blambda;
      }
    }
    
    fprintf (stdout,"\n increment %ld    total lambda %e",i,lambda);

    if (stop==0)
      continue;
    
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  
  if (Mp->nlman->hdbackupal==1)
    arclsave (i,n,lambda,dl,ra,fp);

  print_close();
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
}
*/


/**
   function solves system of nonlinear algebraic equation
   by arc-length method (continuation method)
   only one right hand side vector is supported with respect
   to nonlinearity and absence of superposition

   d stands for delta
   dd stands for capital delta
   
   fc - nonproportional %vector
   fp - proportional %vector
   n - order of the system
   
   @param lcid - load case id
   
   25.7.2001
*/
/*
void solve_nonlinear_floating_subdomains ()
{
  long i,j,k,n,nlm,ni,ini,li,numr,stop,modif,lcid;
  double dl,dlmax,dlmin,psi,ierr,zero,lambda,dlambda,ddlambda,blambda,norf,norfa;
  double a0,a1,a2,l1,l2,norv;
  double *r,*ddr,*fc,*fp,*f,*dv,*du,*dnu,*rhs;
  double *br,*bmu,*av,*av1,*av2,*fa,*fi,*aa1,*aa2,*a,*b,*c;

  
  //long i,j,k,n,ni,ini,stop,modif,li,newmesh, numr;
  //double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,zero,ierr;
  //double lambdao,ss1,ss2,ss3,ss4,ss5;
  //double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi;

  
  
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
  
  //  load case id must be equal to zero
  //  only one load case can be solved
  //  one load case may contain several proportional vectors
  lcid = 0;
  
  //  assembling of stiffness matrix
  stiffness_matrix (lcid);


  //  auxiliary assigment of pointer
  rhs=Lsrs->give_rhs (lcid);
  
  //  assembling of right hand side vector (vector of nodal forces)
  mefel_right_hand_side (lcid,rhs);
  
  //  vector of nodal displacements
  r  = Lsrs->give_lhs (lcid);
  //  vector of proportional load
  fp = Lsrs->give_rhs (lcid*2);
  //  vector of constant load
  fc = Lsrs->give_rhs (lcid*2+1);
  
  
  
  //  allocation of object floating subdomains
  Fsd = new flsubdom;
  
  //  computation of rigid body motions
  //  list of dependent equations
  //  determination of number of Lagrange multipliers
  Fsd->initialization (Ndofm,Mp->ense,fp);
  
  //  number of Lagrange multipliers
  nlm = Fsd->nlm;
  
  
  //  allocation of auxiliary arrays
  //  array for backup of nodal displacements
  br = new double [n];
  //  array for backup of Lagrange multipliers
  bmu = new double [nlm];

  du = new double [n];
  dv = new double [n];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  av = new double [n];
  av1 = new double [n];
  av2 = new double [n];
  aa1 = new double [n];
  aa2 = new double [n];

  dnu = new double [nlm];
  
  a = new double [n];
  b = new double [n];
  c = new double [n];
  

  //  proportional factor
  lambda = 0.0;
  dlambda = 0.0;
  ddlambda = 0.0;
  
  //  attained level in arclength method
  //  it is 0 by default
  //  otherwise it is equal to performed number of steps in previous computations
  li=0;
  
  modif=0;
  
  
  
  //  initiation of load vector
  for (j=0;j<n;j++){
    fa[j]=fc[j]+lambda*fp[j];
  }
  
  
  if (li)
    print_init(-1, "at");
  else
  {
    print_init(-1, "wt");
    print_step(lcid, li, lambda, fa);
  }
  
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  for (i=li;i<ni;i++){

    fprintf (Out,"\n\n arc-length increment %ld",i);
    fprintf (stdout,"\n\n arc-length increment   %ld",i);


    //stop=1; // termit



    //  backup of nodal displacements
    copyv (r, br, n);
    //  backup of Lagrange multipliers
    copyv (Fsd->muo,bmu,nlm);
    //  backup of reached lambda parameter
    blambda=lambda;



    
    //  assembling of tangent stiffness matrix
    if ((Mp->nlman->stmat==tangent_stiff) || (i == li)){
      stiffness_matrix (lcid);
      Fsd->factorize ();
    }
    
    //  solution of K(r).dv = fp-B.dnu
    Fsd->pmcg (fp,Out);

    for (k=0;k<nlm;k++){
      fprintf (Out,"\n krok %ld                 lambda %ld   %lf",i,k,Fsd->w[k]);
    }
    

    Fsd->lagrmultdispl (dv,fp);
    copyv (Fsd->w,dnu,nlm);
    
    
    //  generalized norm of displacement increments
    norv = displincr (dv,n);
    
    Gtm->flsub.coarse_local (dnu,av);
    
    norf=0.0;
    for (k=0;k<n;k++){
      norf+=(fp[k]-av[k])*(fp[k]-av[k]);
    }
    
    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norf);
    
    
    for (k=0;k<n;k++){
      ddr[k]=dlambda*dv[k];
      r[k]+=ddr[k];
      fa[k]=fc[k]+(lambda+dlambda)*fp[k];
    }
    ddlambda=dlambda;
    for (k=0;k<nlm;k++){
      Fsd->ddmu[k]=dlambda*dnu[k];
      Fsd->mu[k]+=Fsd->ddmu[k];
    }
    Fsd->mult_correction ();
    
    
    //  computation of internal forces
    internal_forces (lcid,fi);
    subv(fa, fi, f, n);    
    norf=sizev(f,n);
    norfa = sizev(fa, n);
    norf /= norfa;

    
    //if (Mespr==1)  fprintf (stdout,"\n ddlambda %e    dlambda %e   norf %e  ierr %e",ddlambda,dlambda,norf,ierr);
    if (Mespr==1)  fprintf (stdout,"\n increment %ld     norv %e      norf %e",i,norv,norf);
    
    if (norf<ierr){
      // ******************************************
      //  no inner iteration loop is required  ***
      // ******************************************
      modif++;

      if (modif>1){
	//  arc length modification
	dl*=2.0;
	if (dl>dlmax)  
          dl=dlmax;
	modif=0;
	if (Mespr==1){
	  fprintf (stdout,"\n arc-length must be modified (dl increase) dl=%e",dl);
//	  fprintf (Out,"\n arc-length must be modified (dl increase) dl=%e",dl);
	}
      }
      lambda+=dlambda;
      
      for (k=0;k<nlm;k++){
	Fsd->muo[k]+=Fsd->ddmu[k];
      }
      
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, fa);
      print_flush();
    }
    else{
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){
	if (Mp->nlman->stmat==tangent_stiff)
	  stiffness_matrix (lcid);
	
	Fsd->factorize ();
	
	//  back substitution
	Fsd->pmcg (f,Out);
	Fsd->lagrmultdispl (du,f);

	for (k=0;k<nlm;k++){
	  fprintf (Out,"\n krok %ld iterace %ld     lambda %ld   %lf",i,j,k,Fsd->w[k]);
	}
	
	fprintf (Out,"\n arc-length  vnitrni smycka  %ld %ld",i,j);
	
	
	Gtm->flsub.coarse_local (Fsd->ddmu,av1);
	Gtm->flsub.coarse_local (Fsd->w,av2);
	
	for (k=0;k<n;k++){
	  a[k]=ddr[k]+du[k];
	}
	
	for (k=0;k<n;k++){
	  b[k]=ddlambda*fp[k]-av1[k]-av2[k];
	}
	
	
	Gtm->flsub.coarse_local (dnu,av);
	for (k=0;k<n;k++){
	  c[k]=fp[k]-av[k];
	}
	
	a2=0.0;  a1=0.0;  a0=0.0;
	for (k=0;k<n;k++){
	  a2+=c[k]*c[k]*psi*psi+dv[k]*dv[k];
	  a1+=2.0*b[k]*c[k]*psi*psi+2.0*dv[k]*a[k];
	  a0+=a[k]*a[k]+b[k]*b[k]*psi*psi;
	}
	a0-=dl*dl;
	
	//  solution of quadratic equation
        numr = solv_polynom_2(a2, a1, a0, l1, l2);
        switch (numr)
	{
	  case -1 :
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
            break;
	  case 0 :
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
            break;
    	  case 1 :
            dlambda = l1;
            break;
	  default :
            break;
	}

	if (fabs(l1)>fabs(l2))  dlambda=l2;
	else dlambda=l1;





	for (k=0;k<n;k++){
	  ddr[k]+=du[k]+dlambda*dv[k];
	  r[k]+=du[k]+dlambda*dv[k];
	  fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	for (k=0;k<nlm;k++){
	  Fsd->ddmu[k]+=Fsd->w[k]+dlambda*dnu[k];
	  Fsd->mu[k]+=Fsd->w[k]+dlambda*dnu[k];
	}
	Fsd->mult_correction ();
	
	//fprintf (stdout,"   ddlambda %e",ddlambda);
	
	//fprintf (Out,"\n ddlambda %e     dlambda %e",ddlambda,dlambda);
	//  computation of internal forces
	internal_forces (lcid,fi);
        subv(fa, fi, f, n);	
	norf=sizev(f,n);
	norfa = sizev(fa, n);
	norf /= norfa;
	
	if (Mespr==1){
	  fprintf (stdout,"\n    norf=%e  ierr=%e",norf,ierr);
	  //fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
	}

	if (norf<ierr){
	  lambda+=ddlambda;
	  
	  for (k=0;k<nlm;k++){
	    Fsd->muo[k]+=Fsd->ddmu[k];
	  }
	  
	  Mm->updateipval ();
	  stop=1; 
	  compute_req_val (lcid);
          print_step(lcid, i+1, lambda, fa);
          print_flush();
	  break;
	}
      }
      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          break; 
        }
	if (Mespr==1){
	  fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e",dl);
	  //fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e",dl);
	}
	//  restoring of nodal displacements
        copyv(br, r, n);
	//  restoring of Lagrange multipliers
        copyv(bmu,Fsd->muo,nlm);
	//  restoring of lambda parameter
	lambda=blambda;
      }
    }
    
    fprintf (stdout,"\n increment %ld    total lambda %e",i,lambda);

    if (stop==0)
      continue;
    
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  
  
  
  //if (Mp->nlman->hdbackupal==1)
  //arclsave (i,n,lambda,dl,ra,fp);

  print_close();
  

  delete [] br;
  delete [] bmu;
  delete [] dv;
  delete [] f;
  delete [] fa;  
  delete [] fi;
  delete [] av;
  delete [] av1;
  delete [] av2;
  delete [] aa1;
  delete [] aa2;
  delete [] dnu;

}
*/




/**
   function solves system of nonlinear algebraic equations
   load is applied by prescribed displacements
   
   JK, 19.1.2007
*/
/*
void solve_nonlinear_floating_subdomains_19_1 ()
{
  long i,j,ni,ini,lcid;
  double lambda,dlambda,mindlambda,err,norfi;
  double *f,*lhs,*rhs,*arhs,*fi;
  
  //  only one load case is assumed
  lcid=0;
  //  number of iteration steps/increments
  ni=Mp->nlman->ninr;
  //  maximum number of iterations in one increment
  ini = Mp->nlman->niilnr;
  //  total magnitude of the right hand side
  lambda=0.0;
  //  magnitude of the step increment
  dlambda=Mp->nlman->incrnr;
  //  minimum magnitude of the step increment
  mindlambda=Mp->nlman->minincrnr;
  //  required error
  err=Mp->nlman->errnr;
  
  //
  // fi = new double [?]
  
  Fsd = new flsubdom;
  
  //  nodal displacements
  lhs=Lsrs->give_lhs (lcid);
  //  proportional vector of the prescribed forces
  rhs=Lsrs->give_rhs (lcid);
  //  auxiliary array
  f = new double [Ndofm];
  //  auxiliary vector
  arhs = new double [Ndofm];
  
  
  //  proportional vector of the prescribed forces
  mefel_right_hand_side (lcid,rhs);
  
  //  stiffness matrix
  stiffness_matrix (lcid);
  
  Fsd->initialization (Ndofm,Mp->ense,rhs);
  
  //  loop over iteration steps/increments
  for (i=0;i<ni;i++){
    
    //  increment of the right hand side
    for (j=0;j<Ndofm;j++){
      f[j]=dlambda*rhs[j];
    }
    
    //  solution of the system of algebraic equations
    //  increment of the Lagrange multipliers are obtained
    Fsd->pmcg (f,Out);
    
    //  update of Lagrange multipliers
    Fsd->add_mult (1.0);
    
    //  correction of the Lagrange multipliers
    Fsd->mult_correction ();
    
    
    for (j=0;j<Ndofm;j++){
      arhs[j]=lambda*rhs[j];
    }
    
    //  computation of displacements from Lagrange multipliers
    Fsd->lagrmultdispl (lhs,arhs);
    
    
    //  actual internal forces
    internal_forces (lcid,fi);
    
    //  norm of the vector of internal forces
    norfi = ss (fi,fi,Ndofm);
    
    if (norfi<err){
      //  no inner loop is required
    }
    else{
      //  inner loop is required
      for (j=0;j<ini;j++){
	
	Fsd->solve_lin_alg_system (lhs,fi);
	
	
	internal_forces (lcid,fi);
	
      }
      
    }
    
    
    print_init(-1, "wt");    
    for (i=0;i<Lsrs->nlc;i++){
      //compute_ipstrains (i);
      //compute_ipstresses (i);
      //compute_req_val (i);
      print_step(i, 0, 0.0, NULL);    
    }

    
  }
  
  print_close();
  

}
*/
