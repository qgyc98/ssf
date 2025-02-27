#include "nssolver.h"
#include "arclength.h"
#include "newtonraph.h"
#include "slopesol.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "backupsol.h"
#include "mechprint.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "intpoints.h"
#include "mathem.h"
#include "vector.h"
#include "meshtransfer.h"
#include "elemswitch.h"

#include "j2flow.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>



/**
  Function solves non-linear statics problems.
   
  Created by JK, 2001
  Rewritten by TKo, 09.2011
*/
void solve_nonlinear_statics ()
{
  long i, li;
  double *rhs, *ra, *fa, *fp, *fc, *fl, *flc, *flp;
  double ilambda, dl;
  
  // vector of attained load
  fa = new double[Ndofm];
  // load vector due to forces for both proportional and constant loads
  fl = new double[2*Ndofm];

  switch (Mp->nlman->tnlinsol)
  {
    case arcl:
    case arclrv1:
    case arclrv:
    case arclg:
      seldofinit();
      break;
    default:
      break;
  }
  
  Mp->lambda = 0.0;  
  
  for (i=0;i<Lsrs->nlc;i++){

    // It is necessary to set Mp->lambda to one in order to
    // calculate of contribution from the unit prescribed displacements 
    Mp->lambda = 1.0; 
    // pointer to the right hand side vectors of the actual load case
    rhs=Lsrs->give_rhs(i*2);
    // assembling of load vectors
    mefel_right_hand_side (i, rhs, fl);

    //  attained displacements
    ra = Lsrs->give_lhs(i);
    // clear the vector of attained displacements
    memset(ra, 0, sizeof(*ra)*Ndofm);
    // clear the vector of attained load
    memset(fa, 0, sizeof(*fa)*Ndofm);
    //  vector of proportional load including contributions from prescribed displacements
    fp = Lsrs->give_rhs(i*2);
    //  vector of constant load including contributions from prescribed displacements
    fc = Lsrs->give_rhs(i*2+1);
    //  vector of proportional load due to forces
    flp = fl;
    //  vector of constant load due to forces
    flc = fl+Ndofm;


    if ((Mp->nlman->hdbackupal == 2) || (Mp->hdbcont.restore_stat()))
    { 
      //  initiation of mechanical material models
      Mm->initmaterialmodels(i, true);
      if (Mespr==1)
        fprintf (stdout,"\n Reading of backup file\n");
      // restorage from the backup is performed
      solver_restore(ra, fa, li, ilambda, dl, NULL, Ndofm);
      Mm->updateother();
      Mp->lambda = ilambda;
      print_init(-1, "at");
    }
    else 
    {
      Mp->lambda = ilambda = 0.0;
      li=0;
      //  initiation of mechanical material models
      Mm->initmaterialmodels(i, false);
      print_init(-1, "wt");
      print_step(i, li, Mp->lambda, fa);
      print_flush();
    }
    
    nonlinear_solver(i, ra, fa, fc, fp, flc, flp, ilambda, li);
    print_close();
    
  }
  delete [] fa;
  delete [] fl;
}



/**
  Function calls selected solver of systems of nonlinear equation 
  for the given load case.

  @param lcid - load case id
  @param ra  - %vector of attained displacements
  @param fa  - %vector of attained load
  @param fc  - constant component of load %vector (due to prescribed displacements + forces)
  @param fp  - proportional component of load %vector (due to prescribed displacements + forces)
  @param flc -  constant component of load %vector due to forces
  @param flp -  proportional component of load %vector due to forces
  @param ilambda - initial value of load coefficient lambda (used for backup restorage)
  @param li - initial step id (used for backup restorage)


  @return The function does not return anything.

  Created by Tomas Koudelka and JK.
  Rewritten by Tomas Koudelka, 09.2011
*/
void nonlinear_solver (long lcid, double *ra, double *fa, double *fc, double *fp, double *flc, double *flp, 
                       double ilambda, long li)
{
  switch (Mp->nlman->tnlinsol){
  case arcl:
  case arclg:{
    garclength(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  case newton:
  case newtong:{
    // newton_raphson (lcid);
    gnewton_raphson2(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  case arclrv1:{
    // arclengthrv1 (lcid,1.1,1.0e-3);
    Mp->nlman->check_lambdar = on;
    Mp->nlman->lambdar       = 1.1;
    Mp->nlman->errl         = 1.0e-3;
    garclength(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  case newtonrv1:{
    // newton_raphsonrv1 (lcid,Mp->nlman->rvlam,1.0e-3);
    Mp->nlman->check_lambdar = on;
    Mp->nlman->errl          = 1.0e-3;
    gnewton_raphson2(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  case newtonrestart:{
    // special case for Hinkley
    newton_raphson_restart (lcid);
    break;
  }
  case displctrl:{
    // displ_control (lcid);
    gnewton_raphson2(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  case displctrlrv:{
    // displ_controlrv (lcid, 0.0, 1.0, Mp->nlman->rvlam, Mp->nlman->errl, 1);
    Mp->nlman->check_lambdar = on;
    gnewton_raphson2(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  case arclrv:{
    // arclengthrv1 (lcid, Mp->nlman->rvlam, Mp->nlman->errl);
    Mp->nlman->check_lambdar = on;
    garclength(lcid, Mp->nlman, ra, fa, fc, fp, flc, flp, li, ilambda, yes);
    break;
  }
  default:
    print_err("unknown solver of nonlinear equations system is required", __FILE__, __LINE__, __func__);
  }
}



/**
  The function assembles the stiffness %matrix
  it depends on the variable Mp->nlman->stmat
  for possible values (see alias.h).
  
  @param lcid - load case id
  @param i    - index in outer loop
  @param j    - index in inner loop
  @param li   - lower index of the outer loop = number of the first outer increment
                li=0 if the computation starts from the initial state
	        li>0 if the computation starts from backup
  @param fusm - flag for forced update of the stiffness %matrix 
                fu=yes => assemble stiffness %matrix regardless of other settings
                fu=no  => assemble stiffness %matrix as usual
  @return The function returns 1 if the stiffness_matrix function has been called and zero 
          if the content of the stiffness %matrix has not been changed.
          Additionally, it calls functions which changes the content of 
          the global stiffness %matrix.

  JK, TKo 30.8.2010	    
*/
long assemble_stiffness_matrix (long lcid,long i,long j,long li, answertype fusm)
{
  long ret = 0;

  if ((Smat == NULL) || (fusm == yes)) // matrix has not yet been assembled (changed number of DOFs) or forced update is required
  {
    stiffness_matrix (lcid);
    return 1;
  }
  
  switch (Mp->nlman->stmat){
  case initial_stiff:{
    if (i==li && j<0){
      //  this condition describes the first step in the computation
      stiffness_matrix (lcid);
      ret = 1;      
    }
    break;
  }
  case tangent_stiff:{
    long tmp = Mespr;

    if ((j >= 0) || (i>0))
      Mespr=0;

    stiffness_matrix (lcid);

    if ((j >= 0) || (i>0))
      Mespr=tmp;

    ret = 1;      
    break;
  }
  case secant_stiff:{
    stiffness_matrix (lcid);
    ret = 1;      
    break;
  }
  case incr_tangent_stiff:{
    if (j<0){
      //  the stiffness matrix is updated only after outer step
      //  this fact is indicated by the variable j which is set to negative number
      stiffness_matrix (lcid);
      ret = 1;      
    }
    break;
  }
  case ijth_tangent_stiff:{
    if (j<0){
      //  this is modification in outer loop
      if ((i%Mp->nlman->ithstiffchange == 0) || (i==li))
      {
	stiffness_matrix (lcid);
        ret = 1;
      }
    }
    else{
      //  this is modification in inner loop
      if ((j+1)%Mp->nlman->jthstiffchange == 0)
      {
	stiffness_matrix (lcid);
        ret = 1;      
      }
    }
    break;
  }
  default:{
    print_err("unknown type of stiffness matrix is required",__FILE__,__LINE__,__func__);
  }
  }
  return ret;
}



/**
  The function assembles the attained load %vector fa. It depends on the
  two load vectors fc, fp, attained load coefficient lambda and the setup specified 
  in the nlman.

  @param fa - array of components of attained load %vector
  @param fc - array of the constant load %vector which is added without lambda multiplier.
  @param fp - array of proportional load %vector which is added multiplied by lambda multiplier usually.
              In certain cases, the %vector fp has to be applied until the required value lambdar is 
              attained and than the increasing of load level stops.
  @param n      - number of components of fa, fc and fp
  @param lambda - attained value of the load coefficient
  @param nlman  - structure with setup of handling of lambda and lambda tresholds

  @return The function returns updated attained load %vectro in the parameter fa.

  Created by TKo, 08.2011
*/
void assemble_attained_load_vector(double *fa, double *fc, double *fp, long n, double lambda, nonlinman *nlman)
{
  double aux = lambda;

  if ((nlman->check_lambdar == on) && (lambda > nlman->lambdar))
    aux = nlman->lambdar;
  addmultv(fc, fp, aux, fa, n);
}



/**
  The function performs the divergency check of the inner loop if it is required in the
  nlman setup. The divergency is detected with help of the least square method where
  the quadratic approximaton of the error/step_id curve is performed. If the minimum of the
  approximing quadratic function is exceeded the divergency of the iteration process is assumed 
  and signalized by the nonzero return value. The minimum number of steps of the inner loop is 
  performed depending on the nlman->div_min_steps setup but at least three steps has to be 
  performed for the correct approximation.

  @param nlman - pointer to the structure with the setup of the nonlinear solver
  @param lsm_a - matrix used for the least square method. It will be used for the storage of contributions
                 from the particular inner loop steps.
  @param lsm_r - right hand side for the least square method. It will be used for the storage of contributions
                 from the particular inner loop steps.
  @param lsm_l - left hand side vector (unknown parameters of quadratic function)
  @param j     - inner loop step id
  @param norf  - error norm (norm of unbalanced forces)
                  

  @retval 0 - no divergency detected
  @retval 1 - divergency was detected

  Created by TKo, 08.2011
*/
long check_divergency(nonlinman *nlman, matrix &lsm_a, vector &lsm_r, vector &lsm_l, long j, double norf)
{
  double ret;
  double *divc_step = nlman->divc_step;
  double *divc_err  = nlman->divc_err;
  long min_steps = nlman->div_min_steps;
  long i;

  if (nlman->check_div == on)
  {
    if (j >= min_steps-1)
    { 
      i = j%min_steps;
      if (j == min_steps-1)       // accumulation of values + calculation of parameters   
        ret = lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero, 1);
      else                        // accumulation/removal of values + calculation of parameters   
        ret = lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, divc_step[i], divc_err[i], Mp->zero, 2); 
      divc_step[i] = double(j+1);
      divc_err[i]  = norf;
      i = (j-1)%min_steps;
      //fprintf(stdout, " % le*x*x%+le*x%+le", lsm_l[2], lsm_l[1], lsm_l[0]);
      if ((ret > 0.0) && (divc_err[i] < norf)) 
      { // peak of quadratic approximation was passed and error norm is increasing
        fprintf (stdout,"\n\n Divergence of inner iteation loop is detected. Inner loop is stopped.\n\n");
        return(1);
      }
    }
    else
    {
      // accumulation of values only
      lsm_quad(lsm_a, lsm_r, lsm_l, double(j+1), norf, 0.0, 0.0, Mp->zero, 0);
      i = j%min_steps;
      divc_step[i] = double(j+1);
      divc_err[i]  = norf;
    }
  }
  return(0);
}



/**  
  Function solves system of nonlinear algebraic equations by 
  Newton-Raphson method. Solved system does not contain time variable.
  The method retarts analysis from the saved data of the time dependent 
  analysis. Internal data rearangement is performed with respect of 
  solved problem in Hinkley nuclear power plant.
  
  @param lcid      - required load case id

  @return The function returns reached load coefficient lambda.

  Created by JK, Tomas Koudelka
*/
void newton_raphson_restart (long lcid)
  //  function solves system of nonlinear algebraic
  //  equations by Newton-Raphson method
  //  solved system does not contain time variable
  //
  //  16.8.2001
{
  long i,j,k,n,ni,ini,nii;
  double lambda,blambda,dlambda,dlambdamin,err,norf,norfa;
  double *r,*rb,*dr,*f,*fi,*fb,*fc,*fp;


  //  number of rows of the matrix
  n = Ndofm;
  //  maximum number of increments
  ni = Mp->nlman->ninr;
  // number of initial increments
  nii = Mp->nlman->nienr;
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




  long ncompstr;
  double tmp;
  FILE *aux;

  // if restart from backup file is required
  if (Mp->nlman->hdbr)
  {
    aux = fopen (Mp->nlman->backupfname, "rt");
    if (aux == NULL)
    {
      fprintf(stderr, "\n\nError - unable to open backup file %s\n", Mp->nlman->backupfname);
      fprintf(stderr, "file %s, line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Program will be terminated\n");
      abort();
    }
    // searching required backup
    fprintf(stdout, "\n");
    solver_restore(r,fc,i,tmp,tmp,NULL,n);
    /*
    for (k=0; k<Mp->nlman->hdbid; k++)
    {
      fprintf(stdout, "\r reading backup record %ld", k+1);
      fflush(stdout);
      solver_restore(aux, r, fc, i, tmp, j, 39, 8);
    }
    fprintf(stdout, "\n backup time: %e\n", tmp);
    if (j != Ndofm)
    {
      fprintf(stderr, "\n\nError - incorrect number of DOFs is stored in backup file\n");
      fprintf(stderr, "file %s, line %d\n", __FILE__, __LINE__);
    }
    */
/*
  Begining of Hinkley modification
*/
    for (i=0; i<Mm->tnip; i++)
    {
      if (Mm->ip[i].ncompeqother >= 8)
      {
        ncompstr = Mm->ip[i].ncompstr;
        for (j=0; j<ncompstr; j++)
          Mm->ip[i].stress[j] = Mm->ip[i].eqother[ncompstr+j];
      }
    }
  }

  Mp->strcomp = 0;
  internal_forces(lcid, fc);
  for (i=0; i < Ndofm; i++)
    fc[i] *= 1.0e+6;
  Mp->strcomp = 1;
/*
 End of Hinkley modification
*/


  // compute and print values initial values
  compute_req_val (lcid);
  print_init(-1, "wt");
  print_step(lcid, 0, 0, fc);
  print_flush();


  //  assembling of tangent stiffness matrix
  stiffness_matrix (lcid);
  // norm of load vector
  norfa=normv(fc,n);
  //  internal iteration loop
  for (j=0;j<nii;j++)
  {
    //  computation of internal forces
    internal_forces (lcid,fi);
    //  vector of unbalanced forces
    for (k=0;k<n;k++)
      fb[k]=fc[k]-fi[k];
    //  norm of vector of unbalanced forces
    norf=normv(fb,n);
    norf /= norfa;
    
    if (norf<err)
    {
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, j+1, -0.0000001, fc);
      print_flush();
      break;
    }
    //Smat->solve_system (Gtm,dr,fb);
    Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

/*
    double tmp1, tmp2, tmp3, tmp4;
    long ncomp;
    tmp1=tmp2=tmp3=tmp4=0.0;
    fprintf(Out, "\n\n*Krok %ld\n", j);
    fprintf(Out,"\n\nn\tfc\tfi\tfc-fi\tr\tr+dr\tddr\n");
    for(k=0; k<n; k++)
    {
      tmp1+=fc[k]*fc[k];
      tmp2+=fi[k]*fi[k];
      tmp3+=r[k]*r[k];
      r[k]+=dr[k];
      tmp4+=r[k]*r[k];
      if (k%100==0)
        fprintf(Out, "%ld\t%e\t%e\t%e\t%e\t%e\t%e\n",k, tmp1, tmp2, tmp1-tmp2, tmp3, tmp4, tmp3-tmp4);
    }
    fprintf(Out, "%ld\t%e\t%e\t%e\t%e\t%e\t%e\n",k, tmp1, tmp2, tmp1-tmp2, tmp3, tmp4, tmp3-tmp4);
    fprintf(Out,"\n\nIpp\t omega\n");
   
    for (k=0; k < Mm->tnip; k++)
    {
      if (Mm->ip[k].ncompother > 8)
      {
        ncomp = 2*Mm->ip[k].ncompstr;
        fprintf(Out, "%ld\t%e\t%e\n", k, Mm->ip[k].other[ncomp],Mm->ip[k].other[ncomp+1]);
      }
    }
    fflush(Out);*/
    for(k=0; k<n; k++)
      r[k]+=dr[k];
/*    swap_other();
    compute_req_val (lcid);
    print_step(lcid, j+1, j+1, fc);
    print_flush();
    swap_other();*/
    fprintf (stdout,"\n Equilibrium interation loop j %ld     norf=%e", j, norf);
  }
  if ((nii > 0) && (norf > err))
  {
    fprintf(stdout, "\n Initial equilibrium could not be reached\n");
    fprintf(stdout, " Program will be terminated\n");
    print_close();
    delete [] rb;
    delete [] f;
    delete [] fb;
    delete [] fi;
    delete [] dr;
    return;
  }
  else
    fprintf(stdout, "\n\n Starting incremental loading . . .\n");


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
   
    fprintf (stdout,"\n increment %ld   lambda %e,    dlambda %e",i,lambda,dlambda);
    
    //  vector of maximum load and vector of load increment
    for (j=0;j<n;j++){
      f[j]=fc[j]+lambda*fp[j];
      fb[j]=dlambda*fp[j];
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
      fb[j]=f[j]+dlambda*fp[j];
    }
    norfa=normv(fb,n);
    for (j=0;j<n;j++){
      fb[j] -= fi[j];
    }
    //  norm of vector of unbalanced forces
    norf=normv(fb,n);
    norf /= norfa;
    
    
    
//    fprintf (Out,"\n norf=%e",norf);

    if (norf<err){
      lambda+=dlambda;
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, f);
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
	fb[k]=f[k]+dlambda*fp[k]-fi[k];
      }
      
      //  norm of vector of unbalanced forces
      norf=normv(fb,n);
      norf /= norfa;
      
      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf=%e  dlambda %e",i,j,norf,dlambda);
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
    }
    else{
      lambda+=dlambda;
      Mm->updateipval ();
      compute_req_val (lcid);
      print_step(lcid, i+1, lambda, f);
      print_flush();
    }
    
  }
  
  delete [] rb; delete [] dr;  delete [] fi;  delete [] fb;  delete [] f;
  print_close();
}



/**
  The function saves actual step of arclength method to the text file.

  @param fi - id of saved step
  @param n - number of r components
  @param blambda - attained load coefficient lambda
  @param dl - legth of arc
  @param r - %vector of attained displacements
  @param fp - %vector of proportional load

  @return The function does not return anything.

  Created by JK,
*/
void arclsave (long fi,long n,double blambda,double dl,double *r,double *fp)
{
  long i,j;
  char *file;
  FILE *aux;
  file = new char[strlen(Mp->path)+11+1];
  sprintf (file,"%sarcl.backup",Mp->path);
  aux = fopen (file,"w");

  //aux = fopen ("arcl.backup","w");
  
  //  attained number of steps
  fprintf (aux,"%ld\n",fi);
  //  number of rows
  fprintf (aux,"%ld\n",n);
  //  attained coefficient of proportionality
  fprintf (aux,"%e\n",blambda);
  //  arc-length
  fprintf (aux,"%e\n",dl);
  //  number of selected dofs
  fprintf (aux,"%ld\n",Mp->nlman->nsdofal);
  //  selected dofs
  for (i=0;i<Mp->nlman->nsdofal;i++){
    fprintf (aux,"%ld\n",Mp->nlman->seldofal[i]);
  }
  //  left hand side vector
  for (i=0;i<n;i++){
    fprintf (aux,"%e\n",r[i]);
  }
  //  proportionality vector
  for (i=0;i<n;i++){
    fprintf (aux,"%e\n",fp[i]);
  }
  //  number of integration points
  fprintf (aux,"%ld\n",Mm->tnip);
  //  integration points
  for (i=0;i<Mm->tnip;i++){
    fprintf (aux,"%ld %ld\n",Mm->ip[i].ncompstr,Mm->ip[i].ncompother);
    for (j=0;j<Mm->ip[i].ncompstr;j++){
      fprintf (aux,"%e\n",Mm->ip[i].strain[j]);
    }
    for (j=0;j<Mm->ip[i].ncompstr;j++){
      fprintf (aux,"%e\n",Mm->ip[i].stress[j]);
    }
    for (j=0;j<Mm->ip[i].ncompother;j++){
      fprintf (aux,"%e\n",Mm->ip[i].eqother[j]);
    }
  }
  delete [] file;
  fclose (aux);
}



/**
  The function restores saved step of arclength method from the text file.

  @param fi - id of saved step
  @param n - number of r components
  @param blambda - attained load coefficient lambda
  @param dl - legth of arc (output)
  @param r - %vector of attained displacements (output)
  @param fp - %vector of proportional load (output)

  @return The function returns restored values in the passed parameters.

  Created by JK,
*/
void arclopen (long &fi,long &n,double &blambda,double &dl,double *r,double *fp)
{
  long i,j;
  char *file;
  FILE *aux;

  file = new char[strlen(Mp->path)+11+1];
  sprintf (file,"%sarcl.backup",Mp->path);
  aux = fopen (file,"r");
  
  //  attained number of steps
  fscanf (aux,"%ld",&fi);
  //  number of rows of the matrix = number of unknowns
  fscanf (aux,"%ld",&n);
  //  attained coefficient of proportionality
  fscanf (aux,"%le",&blambda);
  //  arc-length
  fscanf (aux,"%le",&dl);
  //  number of selected dofs
  fscanf (aux,"%ld",&Mp->nlman->nsdofal);
  //  selected dofs
  for (i=0;i<Mp->nlman->nsdofal;i++){
    fscanf (aux,"%ld",&Mp->nlman->seldofal[i]);
  }
  //  left hand side vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",&r[i]);
  }
  //  proportionality vector
  for (i=0;i<n;i++){
    fscanf (aux,"%le",&fp[i]);
  }
  //  number of integration points
  fscanf (aux,"%ld",&Mm->tnip);
  //  integration points
  for (i=0;i<Mm->tnip;i++){
    fscanf (aux,"%ld %ld",&Mm->ip[i].ncompstr,&Mm->ip[i].ncompother);
    for (j=0;j<Mm->ip[i].ncompstr;j++){
      fscanf (aux,"%le",&Mm->ip[i].strain[j]);
    }
    for (j=0;j<Mm->ip[i].ncompstr;j++){
      fscanf (aux,"%le",&Mm->ip[i].stress[j]);
    }
    for (j=0;j<Mm->ip[i].ncompother;j++){
      fscanf (aux,"%le",&Mm->ip[i].eqother[j]);
    }
  }

  delete [] file;
  fclose (aux);
}
