#include "epsolver.h"
#include "arclength.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "mechprint.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "nssolver.h"
#include "springel.h"
#include "intpoints.h"
#include "vector.h"
#include "matrix.h"
#include "mathem.h"
#include "gmatrix.h"
#include "tensor.h"
#include "vecttens.h"
#include <math.h>
#include <string.h>



/**
  The function controls run of earth pressure problem computation.

  @retval The function does not return anything.

  Created by Tomas Koudelka
*/
void solve_epressure ()
{
  long i;
  double *rhs, *lhs;


  for (i=0;i<Lsrs->nlc;i++){

    rhs=Lsrs->give_rhs (i);
    lhs=Lsrs->give_lhs (i);
    //  load vectors
    nullv (rhs+i*2*Ndofm, 2*Ndofm);
    if (Mp->tepsol < epplast_sol)
    {
      Mb->lc[i].assemble (i,rhs+i*Ndofm,NULL,1.0);
      Mb->dlc[i].assemble (rhs+i*Ndofm, lhs+i*Ndofm, NULL);
    }
    else
    {
      Mb->lc[i*2+0].assemble (i*2+0,rhs+(i*2+0)*Ndofm,NULL,1.0);
      Mb->lc[i*2+1].assemble (i*2+1,rhs+(i*2+1)*Ndofm,NULL,1.0);
    }
    //  solver
    epressure_solver (i);
  }
}



/**
  The function calls appropriate solver for the earth pressure problem.

  @param lcid - load case id.

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void epressure_solver (long lcid)
{
  switch (Mp->tepsol)
  {
    case gep_sol:
      general_epressure (lcid);
      break;
    case gepvarsup_sol:
      general_epressure_varsup (lcid);
      break;
    case epplast_sol:
      femplast_epressure (lcid);
      break;
    default:
      print_err("unknown solver of nonlinear equation system is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function solves erath pressure problem using GLPT. Iteration of load is used.

  @param lcid - load case id.

  @return The function does not return anything.
  
  Created by Tomas Koudelka
*/
void general_epressure (long lcid)
{
  long i, j, n, ni;
  double err, norf, dl;
  double *r,*f,*fi,*fb, *fc, lambda, rlam;


  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->niep;
  // required norm of unbalanced forces
  err = Mp->errep;
  // required dynamic load step
  dl = Mp->stepep;

  fb = new double [n];
  fi = new double [n];
  memset (fb,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));

  //  initialization phase
  r  = Lsrs->give_lhs (lcid);
  f  = Lsrs->give_rhs (2*lcid);
  fc = Lsrs->give_rhs (2*lcid+1);

  //***************************
  //  main iteration loop  ****
  //***************************
  lambda = 0;
  for (i=0;i<ni;i++){

    print_init(i, "wt");
    //  vectors of loading
    nullv (f+lcid*n, n);
    Mb->lc [lcid].assemble (0,f+2*lcid*n,NULL,1.0);
    Mb->dlc[lcid].assemble (f+2*lcid*n,r+lcid*n, NULL);
    //  backup of load vector
    for (j=0;j<n;j++){
      f[j] = (f[j] - fb[j]); // difference between new and backupped load vector 
      if (lambda < 1.0) // in case that the load is not fully applied
        f[j] *= dl/(1.0-lambda);     // use only required increment
      fc[j] = fb[j]; // store previously reached load
      fb[j] += f[j]; // add new increment to previously reached load
    }
    lambda += dl;
    if (lambda > 1.0)
      lambda = 1.0;
    //    memset (r,0,n*sizeof(double));
    rlam = arclengthrv(lcid, 1.0, Mp->rlerrep);    
    if (fabs(rlam-1.0) > Mp->rlerrep)
    {
      print_err("arclength could not reach required lambda", __FILE__, __LINE__, __func__);
      abort();
    }
    //  computation of internal forces, updating strains in the integration points
//    internal_forces (lcid,fi);
    print_step(lcid, i+1, lambda, fb);

    nullv (f+lcid*n, n);
    Mb->lc [lcid].assemble (0,f+2*lcid*n,NULL,1.0);
    Mb->dlc[lcid].assemble (f+2*lcid*n,r+lcid*n, NULL);
    for (j = 0; j < n; j++)
      fi[j] = f[j] - fb[j];
    norf = ss(fi, fi, n);
    norf /= ss(f, f, n);
    fprintf (stdout, "\n\n MAIN ITERATION LOOP: NORF=%e, lambda = %e\n", norf, lambda);
    fprintf (stdout, "********************************************\n\n");
    if ((norf < err) && (lambda == 1.0))
      break;
    print_close();
  }
  print_close();
  delete [] fi;  delete [] fb;
}


/**
  The function solves erath pressure problem using GLPT. Temporary supports are used for
  better algorithm convergence.

  @param lcid - load case id.

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void general_epressure_varsup (long lcid)
{
  long i, j, n, ni;
  double err, norf;
  double *r,*rb,*f,*fi,*fb, lambda;
  double **br;


  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->niep;
  // required norm of unbalanced forces
  err = Mp->errep;

  rb = new double [n];
  fb = new double [n];
  fi = new double [n];
  br = new double*[Gtm->nn];
  memset (rb,0,n*sizeof(double));
  memset (fb,0,n*sizeof(double));
  memset (fi,0,n*sizeof(double));
  memset (br,0,Gtm->nn*sizeof(*br));
  for (i = 0; i < Gtm->nn; i++)
    br[i] = new double[Gtm->give_ndofn(i)];

  //  initialization phase
  r  = Lsrs->give_lhs (lcid);
  f  = Lsrs->give_rhs (lcid);

  //***************************
  //  main iteration loop  ****
  //***************************
  for (i=0;i<ni;i++) 
  {
    print_init(i, "wt");
    //  vectors of loading
    nullv (f+lcid*n, n);
    Mb->lc [lcid].assemble (0,f+lcid*n,NULL,1.0);
    Mb->dlc[lcid].assemble (f+lcid*n, NULL,r+lcid*n);
    //  backup of load vector
    for (j=0;j<n;j++){
      fb[j] = f[j];
    }

    memset (r,0,n*sizeof(double));
    lambda = arclengthrv(lcid, 1.0, Mp->rlerrep);

    //  computation of internal forces, updating strains in the integration points
    internal_forces (lcid,fi);
    for (j = 0; j < n; j++)
      f[j] = fb[j] - fi[j];
    print_step(lcid, i+1, lambda, fb);
    norf = ss(f, f, n);
    fprintf (stdout, "\n MAIN ITERATION LOOP: NORF=%e", norf);
    if ((norf < err) && (i >= Mp->nselnodep))
      break;
    if (i < Mp->nselnodep)
    {
      // restore original support state
      Gtm->restore_codnum();
      fprintf(stdout, "\n removing support in the node number %ld", Mp->selnodep[i]+1);
      for (j = 0; j <= i; j++)
      {
        // remove support in the x direction from the selected nodes depending current iteration step
        Gtm->gnodes[Mp->selnodep[j]].cn[0] = 1;
      }
      // generation new code numbers
      n = Ndofm = Gtm->gencodnum();
      // backup of displacements in the nodes
      for (j = 0; j < Mt->nn; j++)
        noddispl (lcid, br[j], j);
      //  initialization phase for next step
      //  left and right hand sides should be deleted
      delete Lsrs;
      delete Smat;  // size of global matrix will grow -> we can delete it
      Smat = NULL;
      delete [] fi;
      delete [] fb;
      delete [] rb;
      Lsrs = new lhsrhs;
      Lsrs->alloc ();
      rb = new double [n];
      fb = new double [n];
      fi = new double [n];
      memset (rb,0,n*sizeof(double));
      memset (fb,0,n*sizeof(double));
      memset (fi,0,n*sizeof(double));

      r  = Lsrs->give_lhs (lcid);
      f  = Lsrs->give_rhs (lcid);
      restore_displ(lcid, br);
    }
    print_close();
  }
  print_close();
  delete [] fi;  delete [] fb; delete [] rb;
  for (i = 0; i < Gtm->nn; i++)
    delete [] br[i];
  delete [] br;
}

/**
  The function solves erath pressure problem using FEM plasticity models. 
  Iteration of rigid spring support is used. Spring supports are incrementally 
  removed during each itreation step. Removed supports are expected to be ordered
  as last elements, each step is removed the last element and number of elements is 
  decremented. So the user should order removed spring elements so the first removed
  spring element is the last element (element with the greatest number) and so on.
  Newton-Raphson method is used for given plasticity problem solution. This solver supposes
  the static load equal to dead weight of excaved soil is applied to the bottom of 
  excaved part of the soil body and this static load is decremented in the each step when 
  the spring support is removed.

  @param lcid - load case id.

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void femplast_epressure (long lcid)
{
  long i,j,idxp, idxs,ni, *nod, ndofn, *cnum;
  double *rhs, *rb, *react, *fpb, *fsb, *forig, lambda;
  vector ifor;
  ivector cn;

  vector sig;
  vector sigt(6);
  //matrix sigt(3,3);
  double *p = new double[Mm->tnip];
  double *porig = new double[Mm->tnip];
  char fn[1001];
  char fn2[1001];
  FILE *femcad;

  ni = Mp->niep+1;
  //***************************
  //  main iteration loop  ****
  //***************************
  rb    = new double[Ndofm];
  fpb   = new double[Ndofm];
  fsb   = new double[Ndofm];
  forig = new double[Ndofm];
  nod   = new long  [Mp->niep];
  react = new double[Mp->niep];
  cnum  = new long  [Mp->niep]; 
  for (i=0;i<ni;i++)
  {
    print_init(i+1, "wt");
    if (i < 1)
      fprintf (stdout, "\n MAIN ITERATION LOOP: step %ld\n", i+1);
    else
      fprintf (stdout, "\n MAIN ITERATION LOOP: step %ld, elem. %ld was removed\n", i+1, Mt->ne+1);
    fprintf(stdout, " ---------------------------------------\n\n");
    rhs=Lsrs->give_rhs (lcid);
    idxs = (lcid*2+1)*Ndofm;
    idxp = (lcid*2+0)*Ndofm;
    if (i == 0)
    {
      copyv (rhs+idxp, fpb, Ndofm);
      copyv (rhs+idxs, fsb, Ndofm);
      addv(rhs+idxp, fsb, Ndofm);
      nullv (rhs+idxs, Ndofm);
      nullv (forig, Ndofm);
    }
    if (i > 0)
    {      
      nullv (rhs+idxp, Ndofm);
      for (j = 0; j < Ndofm; j++)
      {
        if (i == 1)
          rhs[idxp+j] = (Mp->loadcoef[i-1]-1.0)*fsb[j];
        else
          rhs[idxp+j] = (Mp->loadcoef[i-1]-Mp->loadcoef[i-2])*fsb[j];
      }
/*      if (cnum[i-1] > 0)
	rhs[idxp+cnum[i-1]-1] += react[i-1];*/
    }
    lambda = newton_raphsonep(lcid, forig);
//    nonlinear_solver (lcid);
    for (j=0; j < Ndofm; j++)   
      forig[j] += lambda * (rhs+idxp)[j];
    
    Mm->computenlstresses(0,Mm->ip[0]);
    for (j=0; j<Mm->tnip; j++)
    {
      reallocv(RSTCKVEC(Mm->ip[j].ncompstr, sig));
      Mm->givestress (0,j,sig);
      //vector_tensor (sig, sigt, Mm->ip[j].ssst, stress);
      give_full_vector(sigt,sig,Mm->ip[j].ssst);
      p[j] = first_invar (sigt)/3;
    }
    if (i == 0)
      copyv(p, porig,Ndofm);
    else
    {
      for (j=0; j<Mm->tnip; j++)
      {
        p[j] -= porig[j];      
      }
      sprintf(fn,"%s%s.jk0", Mp->path, Mp->filename);
      sprintf(fn2,"%s%s.jk0.%ld", Mp->path, Mp->filename,i);
      femcad = fopen(fn, "at");
      char dlcid[50];
      sprintf(dlcid, "%ld", i+1); 
      write_elemscalar(femcad, p, "p", dlcid);
      fclose(femcad);
      rename (fn, fn2);
    }
	     

    if (i == Mp->niep)
      break;
    nod[i] = Gtm->gelements[Gtm->ne-1].nodes[0];
    ndofn = Mt->give_ndofn(nod[i]);
    reallocv(RSTCKVEC(ndofn, ifor));
    reallocv(RSTCKIVEC(ndofn, cn));
    Spring->res_internal_forces (lcid, Mt->ne-1, ifor);
    react[i] = ifor[0];
    Mt->give_node_code_numbers (nod[i],cn.a);
    cnum[i] = cn[0];
    // removing currently last spring support element.
    Mt->ne--;
    Gtm->ne--;
    print_close();
  }
  print_close();
  delete [] rb;
  delete [] fpb;
  delete [] fsb;
  delete [] nod;
  delete [] react;
  delete [] cnum;
  delete [] forig;
  delete [] p;
  delete [] porig;
}



/**
  The function extracts nodal displacements bckr to the displacement %vector of Lsrs.

  @param lcid  - number of load case
  @param bckr  - array with displacement in the each node

  @return The function does not return anything.

  Created by Tomas Koudelka  3.12.2002
*/
void restore_displ(long lcid, double **bckr)
{
  long i, j, ii;

  for (i = 0; i < Gtm->nn; i++)
  {
    for (j = 0; j < Gtm->gnodes[i].ndofn; j++)
    {
      ii=Gtm->gnodes[i].cn[j];
      if (ii >  0)   Lsrs->lhs[lcid*Ndofm+ii-1] = bckr[i][j];
    }
  }
}

/**
  The function solves system of nonlinear algebraic equations
  by the arc-length method. This method was modified to stop iteration
  on the given value rlambda of the lambda parameter. Rlambda value is reached with
  required error rlerr. Only one right hand side vector is supported with respect
  to nonlinearity and absence of superposition.

  @param lcid    - load case id
  @param rlambda - required value of lambda parameter
  @param rlerr   - required error between lambda and rlambda

  @return The function returns reached lambda parameter.

  Cretated by Tomas Koudelka
*/
double arclengthrv (long lcid, double rlambda, double rlerr)
{
  long i,j,k,n,ni,ini,stop,modif,li, numr;
  double a0,a1,a2,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norv,norfa,ierr;
  double ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi;
  long cn, ecn[6];
  char file[255];
//  FILE *gr;
  FILE *grout;
  
  
  
  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->nlman->nial;
  //  maximum number of iterations in one increment
  ini = Mp->nlman->niilal;
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
    li=0;
  }
  
  sprintf(file, "%s%s.epg", Mp->path, Mp->filename);
  grout = fopen (file, "wt");  
  
  for (j=0;j<n;j++)
    fa[j]=fc[j]+(lambda)*fp[j];

  //  norm of proportionality vector
  norfp = loadincr(lcid, Mp->nlman, lambda, fp, n);
  modif=0;

  // pouze docasne pro jeden typ ulohy
  long kk;
  matrix sm(6,6);

  // ***************************
  //  main iteration loop   ****
  // ***************************
  for (i=li;i<ni;i++) {

    fprintf(grout, "Lambda %e\n", lambda);
    fprintf(grout, "Posuny ux --------\n");
    for (kk = 0; kk < Mt->nn; kk++)
    {
      cn = Mt->give_dof(kk,0)-1;
      fprintf(grout, "   %2ld %e\n", kk+1, ra[cn]);
    }
    fprintf(grout, "Zatizeni fx --------\n");
    for (kk = 0; kk < Mt->nn; kk++)
    {
      cn = Mt->give_dof(kk,0)-1;
      double fpv = fp[cn];
      fprintf(grout, "   %2ld %le\n", kk+1, fpv);
    }
    fprintf(grout, "Tuhosti a sily--------\n");
    for (kk = 0; kk < Mt->ne; kk++)
    {
      if (Mt->give_elem_type(kk) == spring_1)
      {
        Spring->res_stiffness_matrix (kk,sm);
        Mt->give_code_numbers (kk, ecn);
        if (ecn[0] > 0)
          fprintf(grout, "   %2ld %e %e\n", kk+1, sm[0][0], sm[0][0]*ra[ecn[0]-1]);
      }
    }
    fprintf(grout, "\n");
    fflush(grout);
    //  backup of left hand side vector
    copyv (ra, r, n);
    //  backup of reached lambda parameter
    blambda=lambda;

    fprintf (stdout,"\n arc-length: increment %ld   lambda %e  dl %e",i,lambda,dl);
/*    if (Mespr==1){
      fprintf (Out,"\n\n *******************************************************************");
      fprintf (Out,"\n arc-length: increment %ld,   lambda=%e   dl=%e",i,lambda,dl);
    }*/
    
    //  assembling of tangent stiffness matrix
    stiffness_matrix (lcid);
    
    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    copyv (fp, f, n);
    //  solution of K(r).v=F
    //Smat->solve_system (Gtm,v,f);
    Mp->ssle->solve_system (Gtm,Smat,v,f,Out);
    //  generalized norm of displacement increments
    norv = displincr (lcid, i, Mp->nlman, v, v, n);
    //  compute new dlambda increment
    dlambda = dl/sqrt(norv+psi*psi*norfp*norfp);
    
    for (j=0;j<n;j++){
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
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

      lambda+=dlambda;
      if (((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
      {
        //  modification of the arc length
        dl/=2.0;
        if (dl<dlmin)
        {
          dl=dlmin;
          break;
        }
        modif = 0;
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
      if (fabs(lambda-rlambda) <= rlerr)
        break;
      Mm->updateipval ();
//      aux_mech_nonlin_print (gr,ra,lambda);
//      print_output(lcid);
      print_step(lcid, i, lambda, fa);
      print_flush();
      stop = 0;
    }
    else{
      // ****************************
      //  inner iteration loop  ****
      // ****************************
      stop=0;
      for (j=0;j<ini;j++){
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
            print_err("infinite number of solution of constrained condition in function arclength", __FILE__, __LINE__, __func__);
            break;
	  case 0:
            print_err("nonsolvable constrained condition in function arclength", __FILE__, __LINE__, __func__);
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
	  fa[k]+=dlambda*fp[k];
	}
	ddlambda+=dlambda;
	
//	fprintf (Out,"\n ddlambda %e     dlambda %e",ddlambda,dlambda);
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
//	  aux_mech_nonlin_print (gr,ra,lambda);
          print_step(lcid, i, lambda, fa);
          print_flush();
	  stop=1;
          if (((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
            stop=0;
          if (fabs(lambda-rlambda) <= rlerr)
	    stop=2;
          if (stop > 0)
            Mm->updateipval();
          break; 
	}
      }

      // limit value of lambda was reached
      if (stop == 2)
        break;

      modif=0;
      if (stop==0){
	//  modification of the arc length
	dl/=2.0;
	if (dl<dlmin){
          dl=dlmin;  
          break; 
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
    if (stop==0)
      continue;

    // termit
    if (i%10 && Mm->plast)
      continue;
    // timret
//    print_output(lcid);
  }
  
  // ------------------------------------
  //  finish of main iteration loop  ----
  // ------------------------------------
  
  if (Mp->nlman->hdbackupal==1)
    arclsave (i,n,lambda,dl,ra,fp);

  
  //  if (i==ni){   sprintf (file,"%s%s.dx.%ld.finito",Mp->path,Mp->filename,i);   print_default_2_dx (Gtm,Mp,Mt,0,file);  }
  // timret
  
  delete [] r;		    
  delete [] fi;  
  delete [] fa;  
  delete [] f;
  delete [] v;  
  delete [] u;  
  delete [] ddr;
  
  //fclose (gr); 
  return (lambda);
}














/* -------------------------------------------------------------------------------------- */







/**
  This function solves system of nonlinear algebraic equations
  by the arc-length method. This method was modified to stop iteration
  on the given value rlambda of the lambda parameter. Rlambda value is reached with
  required error rlerr. Only one right hand side vector is supported with respect
  to nonlinearity and absence of superposition.


  @param lcid    - load case id
  @param rlambda - required value of lambda parameter
  @param rlerr   - required error between lambda and rlambda

  @return The function returns reached lambda parameter.

double arclengthrv (long lcid, double rlambda, double rlerr)
  //
  //  d odpovida delte
  //  dd odpovida DELTE (trojuhelniku)
  //
  //  fc - konstantni vektor
  //  fp - vektor proporcionality
  //
  //  n - pocet radku matice
{
  long i,j,k,n,ni,ini,stop,modif,li, kk;
  double a0,a1,a2,d,l1,l2,dl,dlmax,dlmin,psi,lambda,blambda,dlambda,ddlambda,norf,norfp,norfa,norv,zero,ierr;
  double lambdao,ss1,ss2,ss3,ss4,ss5;
  double *r,*ra,*ddr,*u,*v,*f,*fa,*fp,*fc,*fi;
  FILE *grout;

  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->nial;
  //  maximum number of iterations in one increment
  ini = Mp->niilal;
  //  computer zero
  zero = Mp->zero;
  //  required error in inner loop
  ierr = Mp->erral;
  //  length of arc
  dl = Mp->dlal;
  //  maximum length of arc
  dlmax = Mp->dlmaxal;
  //  minimum length of arc
  dlmin = Mp->dlminal;
  //  displacement-loading driving switch
  psi = Mp->psial;

  //  allocation of auxiliary arrays
  r = new double [n];
  ddr = new double [n];
  u = new double [n];
  v = new double [n];
  f = new double [n];
  fa = new double [n];
  fi = new double [n];
  fc = new double [n];
  memset(fc, 0, sizeof(*fc)*n);

  //  allocation of backup arrays for state variables in integration points
  //  Mm->alloc_backup ();

  //  initialization phase
  ra = Lsrs->give_lhs (lcid);
  fp = Lsrs->give_rhs (lcid);
//  fc = Lsrs->give_rhs (lcid*2+1);

  //  array containing selected DOFs
  seldofinit ();

  lambda=0.0;
  lambdao=0.0;
  li=0;

  //  norm of proportionality vector
  norfp = ss (fp,fp,n);
  matrix sm(6,6);
  modif=0;
  grout = fopen("gr.dat", "wt");
  // ***************************
  //   main iteration loop  ****
  // ***************************
  for (i=li;i<ni;i++){

    fprintf(grout, "%e\n", lambda);
    for (kk = 0; kk < Mt->nn; kk++)
      fprintf(grout, "   %e\n", ra[3*kk]);
    fprintf(grout, "Tuhosti --------\n");
    for (kk = 15; kk < Mt->ne; kk++)
    {
      Spring->res_stiffness_matrix (kk,sm);
      fprintf(grout, "   %e\n", sm[0][0]);
    }
    fprintf(grout, "\n");
    fflush(grout);
    //  backup of left hand side vector
    for (j=0;j<n;j++)
      r[j]=ra[j];
    //  backup of reached lambda parameter
    blambda=lambda;
    //  backup of state variables
//    Mm->backup_state_var ();
//    aux_nonlin_print (gr,ra,lambda);

    fprintf (stdout,"\n arc-length: increment %ld   lambda %e  dl %e",i,lambda,dl);
//    if (Mespr==1)
//    {
//      if (Mp->adp==1)
//        Mp->ado->print (Out,i,lambda,0);
//    }

    //  assembling of tangent stiffness matrix
    stiffness_matrix (lcid);

    //  backup of the fp, in ldl_sky the right hand side will be destroyed
    for (j=0;j<n;j++)
      f[j]=fp[j];
    //  solution of K(r).v=F
    Smat->solve_system (v,f);

    //  generalized norm of displacement increments
    norv = displincr (v,n);

    dlambda = dl/sqrt(norv+psi*psi*norfp);

    for (j=0;j<n;j++)
    {
      ddr[j]=dlambda*v[j];
      ra[j]+=ddr[j];
      fa[j]=fc[j]+(lambda+dlambda)*fp[j];
    }

    ddlambda=dlambda;

    //  computation of internal forces
    internal_forces (lcid,fi);

    for (k=0;k<n;k++)
      f[k]=fa[k]-fi[k];

    norf=ss(f,f,n);
    norf=sqrt(norf);
    norfa = ss(fa,fa,n);
    norfa = sqrt(norfa);
    norf /= norfa;

    if (Mespr==1)  fprintf (stdout,"\n %e %e norf %e",lambda,dl,norf);

    if (norf<ierr)
    {
      // ******************************************
      //   no inner iteration loop is required  ***
      // ******************************************
      modif++;
      if (modif>1)
      {
        //  arc length modification
        dl*=2.0;
        if (dl>dlmax)  dl=dlmax;
        modif=0;

        if (Mespr==1)
        {
          fprintf (stdout,"\n arc length must be modified (dl increase) dl=%e    norf=%e",dl,norf);
          fprintf (Out,"\n arc length must be modified (dl increase) dl=%e    norf=%e",dl,norf);
        }
      }

      lambda+=dlambda;
      if (((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
      {
        //  modification of the arc length
        dl/=2.0;
        if (dl<dlmin)
        {
          dl=dlmin;
          break;
        }

        if (Mespr==1)
        {
          fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
          fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
        }

        //  restoring of left hand side vector

        for (j=0;j<n;j++)
          ra[j]=r[j];

        //  restoring of lambda parameter
        lambda=blambda;
      }
      if (fabs(lambda-rlambda) <= rlerr)
        break;
      Mm->updateipval();
      continue;
    }

    // ****************************
    //   inner iteration loop  ****
    // ****************************
    stop=0;
    for (j=0;j<ini;j++)
    {
//Smat->solve_system (u,f);
    Mp->ssle.solve_system (Gtm,Smat,u,f,Out);

      for (k=0;k<n;k++)
      {
        f[k]=ddr[k];
        ddr[k]+=u[k];
      }

      //  coefficient of quadratic equation
      a2=norv+psi*psi*norfp;
      a1=2.0*ss(ddr,v,n)+2.0*ddlambda*psi*psi*norfp;
      a0=ss(ddr,ddr,n)+ddlambda*ddlambda*psi*psi*norfp-dl*dl;

      //  solution of quadratic equation
      if (fabs(a2)<zero)
      {
        if (fabs(a1)<zero)
        {
          if (fabs(a0)>zero)
            fprintf (stderr,"\n\n nonsolvable constrained condition in function arclength");
          else
            fprintf (stderr,"\n\n infinite number of solution of constrained condition in function arclength");
        }
        else
         dlambda=(0.0-a0)/a1;
      }
      else
      {
        d=a1*a1-4.0*a2*a0;
        if (d<0.0)
        {
          fprintf (stderr,"\n negative discriminant in function arclength");
          break;
        }
        l1=(0.0-a1+sqrt(d))/2.0/a2;
        l2=(0.0-a1-sqrt(d))/2.0/a2;
      }

      //  zacatek novinky

      ss1=0.0;  ss2=0.0;
      ss3=0.0;  ss4=0.0;  ss5=0.0;
      for (k=0;k<n;k++)
      {
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
      //  konec novinky
//      fprintf(Out, "Dlambda %g, i = %ld, j = %ld\n", dlambda, i, j);

      for (k=0;k<n;k++)
      {
        ddr[k]+=dlambda*v[k];
        ra[k]+=u[k]+dlambda*v[k];
        fa[k]+=dlambda*fp[k];
      }
      ddlambda+=dlambda;

      //  computation of internal forces
      internal_forces (lcid,fi);

      for (k=0;k<n;k++)
        f[k]=fa[k]-fi[k];

      norf=ss(f,f,n);
      norf=sqrt(norf);
      norfa = ss(fa,fa,n);
      norfa=sqrt(norfa);
      norf /= norfa;

      if (Mespr==1)
      {
        fprintf (stdout,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
        fprintf (Out,"\n increment  %ld     inner loop  %ld    norf=%e",i,j,norf);
      }

      if (norf<ierr)
      {
        lambda+=ddlambda;

        stop=1;
        if (((lambda - rlambda) > 0.0) && (fabs(lambda-rlambda) > rlerr))
          stop=0;
        if (fabs(lambda-rlambda) <= rlerr)
	  stop=2;
        if (stop > 0)
          Mm->updateipval();
        break;
      }
    }

    // limit value of lambda was reached
    if (stop == 2)
      break;

    modif=0;

    if (stop==0)
    {
      //  modification of the arc length
      dl/=2.0;
      if (dl<dlmin)
      {
        dl=dlmin;
        break;
      }

      if (Mespr==1)
      {
        fprintf (stdout,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
        fprintf (Out,"\n arc length must be modified (dl decrease) dl=%e    norf=%e",dl,norf);
      }

      //  restoring of left hand side vector

      for (j=0;j<n;j++)
        ra[j]=r[j];

      //  restoring of lambda parameter
      lambda=blambda;
      //  restoring of state variables
//      Mm->restore_state_var ();
    }
  }
  fclose(grout);
  delete [] r;
  delete [] fi;
  delete [] fa;
  delete [] fc;
  delete [] f;
  delete [] v;
  delete [] u;
  delete [] ddr;
  return lambda;
}
*/



/**
  Function solves system of nonlinear algebraic equations by Newton-Raphson method.
  Solved system does not contain time variable, it takes into account
  previous reached displacement/stress state and starts with nonzero 
  displacement vector. pfp is proportional load vector of the previous 
  reached initial displacements.

  @param lcid - loadcase id
  @param pfp - proportional load %vector of the prvious reached initial displacements.

  @retrun The function returns reached load coefficient lambda.

  Created by Tomas Koudelka, 16.8.2001
*/
double newton_raphsonep (long lcid, double *pfp)
{
  long i,j,k,n,ni,ini;
  double lambda,blambda,dlambda,dlambdamin,err,norf,norfa;
  double *r,*rb,*dr,*f,*fi,*fb,*fc,*fp;

  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->nlman->ninr;
  //  maximum number of iterations in one increment
  ini = Mp->nlman->niilnr;
  //  required error in inner loop
  err = Mp->nlman->errnr;
  //  increment size
  dlambda=Mp->nlman->incrnr;
  //  minimum increment size
  dlambdamin=Mp->nlman->minincrnr;

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
  
/*  char fn[1001];
  FILE *femcad;
  sprintf(fn,"%s%s.jk0", Mp->path, Mp->filename);
  femcad = fopen(fn, "wt");
  export_femcad(femcad, 0);
  write_displ(femcad, 0, 0);
  Mm->compute_nodeothers(0);
  for (j=0; j<Mm->ip[0].ncompother; j++)
    write_nodscalar(femcad, other, 0, j, 0);
  Mm->compute_nodestresses(0);
  for (j=0; j<Mm->ip[0].ncompstr; j++)
    write_nodscalar(femcad, stress, 0, j, 0);*/
  // ***************************
  //  main iteration loop  ****
  // ***************************
  for (i=0;i<ni;i++){
    
    //  backup of left hand side vector
    for (j=0;j<n;j++){
      rb[j]=r[j];
    }
    //  backup of reached lambda parameter
    blambda=lambda;
    //  backup of state variables
//    Mm->backup_state_var ();
   
    fprintf (stdout,"\n lambda %e,    dlambda %e",lambda,dlambda);
//    fprintf (Out,"\n\n i=%ld     lambda %f",i,lambda);
    //print_displacements (Out,0);

    
    //  vector of maximum load and vector of load increment
    for (j=0;j<n;j++){
      f[j]=fc[j]+lambda*fp[j];
      fb[j]=dlambda*fp[j];
    }
    
    //  assembling of tangent stiffness matrix
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
      fb[j]=pfp[j]+f[j]+dlambda*fp[j];
    }
    norfa=ss(fb,fb,n);
    for (j=0;j<n;j++){
      fb[j] -= fi[j];
    }
    //  norm of vector of unbalanced forces
    norf=ss(fb,fb,n);
    norf /= norfa;
    
//    fprintf (Out,"\n norf=%e",norf);

    if (norf<err){
      lambda+=dlambda;
      Mm->updateipval ();
/*      write_displ(femcad, 0, i+1);
      Mm->compute_nodeothers(0);
      for (j=0; j<Mm->ip[0].ncompother; j++)
        write_nodscalar(femcad, other, 0, j, i+1);
      Mm->compute_nodestresses(0);
      for (j=0; j<Mm->ip[0].ncompstr; j++)
        write_nodscalar(femcad, stress, 0, j, i+1);*/
      print_step(lcid, i, lambda, f);
      print_flush();
      continue;
    }
    
    //  internal iteration loop
    for (j=0;j<ini;j++){
      
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
	fb[k]=pfp[k]+f[k]+dlambda*fp[k]-fi[k];
      }
      
      //  norm of vector of unbalanced forces
      norf=ss(fb,fb,n);
      norf /= norfa;
      
      fprintf (stdout,"\n j=%ld     norf=%e",j,norf);
//      fprintf (Out,"\n j=%ld     norf=%e",j,norf);
      
      if (norf<err)  break;
    }
    
    if (j==ini || norf>err){
      dlambda/=2.0;
      if (dlambda<dlambdamin)  dlambda=dlambdamin;

      if (Mespr==1){
	fprintf (stdout,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
//	fprintf (Out,"\n increment must be modified (dlambda decrease) dlambda=%e    norf=%e",dlambda,norf);
      }
      
      for (j=0;j<n;j++){
      	r[j]=rb[j];
      }
      lambda=blambda;
      //      Mm->restore_state_var ();
    }
    else{
      lambda+=dlambda;
      Mm->updateipval ();
      print_step(lcid, i, lambda, f);
      print_flush();
/*      write_displ(femcad, 0, i+1);
      Mm->compute_nodeothers(0);
      for (j=0; j<Mm->ip[0].ncompother; j++)
        write_nodscalar(femcad, other, 0, j, i+1);
      Mm->compute_nodestresses(0);
      for (j=0; j<Mm->ip[0].ncompstr; j++)
        write_nodscalar(femcad, stress, 0, j, i+1);*/
    }
    
  }
  
  delete [] dr;  delete [] fi;  delete [] fb;  delete [] f;  delete [] rb;
  
//  fclose (femcad);
  return lambda;
}
