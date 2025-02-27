#include <errno.h>

#include "slesolv.h"
#include "mathem.h"
#include "iotools.h"
#include "precond.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "seqselnodes.h"
#include "saddpoint.h"
#include "aggregator.h"
#include "fllopsif.h"

slesolv::slesolv ()
{
  //  type of solver of system of linear equations
  tlinsol=gauss_elim;
  tsol=gauss_elim;
  //  computer zero
  zero=1.0e-15;
  
  //  maximum number of iterations
  ni=0;
  //  required norm of residual vector
  res=0.0;
  //  performed number of iterations
  ani=0;
  //  attained norm of residual vector
  ares=0.0;
  
  //  initial values in iterative methods
  //  iv - the input array is multiplied by zero
  iv = 0;
  
  //  number of rows which are not condensed
  nbdof=0;

  //prec=NULL;
  
  bsize=0;
  sp = NULL;
  fllopptr = NULL;  
}

slesolv::~slesolv ()
{
  delete sp;
  delete fllopptr;
}

/**
   function reads information about solver of linear equations
   
   @param top - pointer to the general topology
   @param in - input stream
   @param mespr - type of message printing
   
   6.4.2003, JK
*/
void slesolv::read (gtopology *top,XFILE *in,long mespr)
{
  //  type of solver of system of linear equations
  xfscanf (in,"%k%m","typelinsol",&linsolvertype_kwdset,(int*)&tlinsol);
  
  switch (tlinsol){
  case gauss_elim:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the gaussian elimination method");
    break;
  }
  case ldl:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the LDL method");
    break;
  }
  case lu:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the LU method");
    break;
  }
  case ll:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the LL method (Cholesky)");
    break;
  }
  case ill:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the incomplete LL method (Cholesky)");
    break;
  }
  case conden:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by condensation method");
    xfscanf (in,"%k%m","typesecondsol",&linsolvertype_kwdset,(int*)&tsol);
    xfscanf (in,"%k%ld","nbdof",&nbdof);
    break;
  }
  case cg:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the CG method");
    //  maximum number of iterations in conjugate gradient method
    //  required norm of vector of residuals
    xfscanf (in,"%k%ld","number_of_iterations",&ni);
    xfscanf (in,"%k%lf","norm_of_residual",&res);
    xfscanf (in,"%k%ld","initial_values",&iv);
    break;
  }
  case bicg:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equationes will be solved by the bi-conjugate gradient method");
    //  maximum number of iterations in bi-conjugate gradient method
    //  required norm of vector of residuals
    xfscanf (in,"%k%ld","number_of_iterations",&ni);
    xfscanf (in,"%k%lf","norm_of_residual",&res);
    xfscanf (in,"%k%ld","initial_values",&iv);
    break;
  }
    
  case lapack_sol:{
    break;
  }
    
  case spdirldl:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by sparse direct method (LDL)");
    break;
  }
  case spdirlu:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by sparse direct method (LU)");
    break;

  }
  case spdirll:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by sparse direct method (LL)");
    break;
  }

  case pardisolver:{
    //if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by PARDISO");
    print_err("PARDISO solver is not supported at this time",__FILE__,__LINE__,__func__);
    break;
  }
    
  case sschur:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the Schur complement method");
    
    schur.read (top,mespr,in);
    
    top->rst=1;
    top->cngen=0;
    
    break;
  }
  case sfeti:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the FETI method");
    
    feti.read (top,in);
    
    break;
  }
  case saddle_point:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the saddle point method");
    
    //sp = new saddpoint ();
    
    //xfscanf (in,"%ld",&sp->ns);
    
    break;
  }

  case jacobi:{
    xfscanf (in,"%k%ld","number_of_iterations",&ni);
    xfscanf (in,"%k%lf","norm_of_residual",&res);
    break;
  }
  case gauss_seidel:{
    xfscanf (in,"%k%ld","number_of_iterations",&ni);
    xfscanf (in,"%k%lf","norm_of_residual",&res);
    break;
  }
  case permon_solver:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the FLLOP library");
    
    fllopptr = new fllopsif ();
    break;
  }
  case cholmod_solver:{
    if (mespr==1)  fprintf (stdout,"\n system of linear equations will be solved by the CHOLMOD library");    
    break;
  }

  default:{
    print_err("unknown type of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
  
  fprintf (stdout,"\n computer zero in slesolv is defined as %le",zero);  
  
  switch (tlinsol){
  case gauss_elim:
  case ldl:
  case lu:
  case ll:
  case ill:
  case conden:
  case lapack_sol:
  case spdirldl:
  case spdirlu:
  case spdirll:
  case pardisolver:
  case permon_solver:
  case cholmod_solver:{
    //  direct methods do not need preconditioners
    break;
  }
  case jacobi:
  case gauss_seidel:
  case cg:
  case bicg:{
    prec.read (top,in,mespr);
    break;
  }
  case saddle_point:
  case sschur:
  case sfeti:{
    break;
  }
  default:{
    print_err("unknown type of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }

  
  if (prec.pt==boss){
    if (prec.agg->impl==2){
      //  sequential topology will be read in function transtop::read
      top->rst=1;
    }
  }
  
  
}

/**
   function prints
   
   @param out - output stream
   
   6.4.2003, JK
*/
void slesolv::print (FILE *out)
{
  fprintf (out,"%d\n",(int)tlinsol);

  switch (tlinsol){
  case gauss_elim:
  case ldl:
  case lu:
  case ll:
  case ill:
  case lapack_sol:
  case spdirldl:
  case spdirlu:
  case spdirll:
  case pardisolver:
  case permon_solver:
  case cholmod_solver:{
    break;
  }
  case conden:{
    fprintf (out,"%ld\n",nbdof);
    break;
  }
  case sfeti:{
    feti.print (out);
    break;
  }
    
  case cg:{
    fprintf (out,"%ld %e %ld\n", ni,res,iv);
    break;
  }
  case bicg:{
    fprintf (out,"%ld %e %ld\n",ni,res,iv);
    break;
  }
  case saddle_point:{
    break;
  }

  case jacobi:{
    fprintf (out,"%ld %e\n", ni,res);
    break;
  }
  case gauss_seidel:{
    fprintf (out,"%ld %e\n", ni,res);
    break;
  }
  default:{
    print_err("unknown type of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
  switch (tlinsol){
  case gauss_elim:
  case ldl:
  case lu:
  case ll:
  case ill:
  case lapack_sol:
  case spdirldl:
  case spdirlu:
  case spdirll:
  case pardisolver:
  case permon_solver:
  case cholmod_solver:{
    //  direct methods do not need preconditioners
    break;
  }
  case jacobi:
  case gauss_seidel:
  case cg:
  case bicg:{
    prec.print (out);
    break;
  }
  case sfeti:
  case saddle_point:{
    break;
  }
  default:{
    print_err("unknown type of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
  
  fprintf (out,"\n");
}

/**
   function initiates values from part of parallel code
   
   @param out - output stream
   
   10.4.2003, JK
*/
void slesolv::initiate (slesolv *ssle)
{
  tlinsol = ssle->tlinsol;
  tsol = ssle->tsol;

  bsize = ssle->bsize;
  zero = ssle->zero;

  ni = ssle->ni;
  res = ssle->res;
  ani = ssle->ani;
  ares = ssle->ares;
  iv = ssle->iv;
}

/**
   function prepares all necessary data for some solvers,
   namely for the Schur complement method and for the FETI method
   
   @param top - pointer to the general topology
   @param out - output file
   
   JK, 8.8.2010
*/
void slesolv::prepare_data (long *ndof,gtopology *top,FILE *out)
{
  long i,j;
  
  if (tlinsol == sschur || tlinsol == sfeti){
    if (top->stop->md==glob_glued){
      top->auto_subdom_detection (out);
    }
  }
  
  // **********************************************************
  //  single processor variant of the Schur complement method
  // **********************************************************
  if (tlinsol == sschur){
    top->stop->coarse_local_map (top,out);

    long **schurnodes;
    schurnodes = new long* [top->stop->ns];
    for (i=0;i<top->stop->ns;i++){
      schurnodes[i] = new long [top->stop->nnsd[i]];
    }
    
    //  selection of boundary/interface nodes
    for (i=0;i<top->stop->ns;i++){
      for (j=0;j<top->stop->nnsd[i];j++){
	if (top->stop->nodmultip[i][j]>1)
	  schurnodes[i][j]=1;
	else
	  schurnodes[i][j]=-1;
      }
    }
    
    seqselnodes *selnodschur;
    selnodschur = new seqselnodes (top->stop->ns,top->stop->nnsd,schurnodes,
				   top->stop->nodmultip,top->stop->ggnbn,top->stop->icnbnmas,
				   top->stop->tnbn,top->stop->icmultip,top->stop->lnbncn,top->stop->ggnbncn,top->stop->sid,
				   top->stop->md,out,1);
    
    //  the following two functions have to be called in this order
    //  otherwise, array cn in the class gnode is rewritten
    //  and wrong results are obtained
    selnodschur->prepare_schur (top,out);
    ndof[0]=top->stop->schur_ordering (top,out);
    
    schur.initiate (selnodschur);
    schur.assemble_subdom_unknowns (top,out);
    
    delete selnodschur;
    
    for (i=0;i<top->stop->ns;i++){
      delete [] schurnodes[i];
    }
    delete [] schurnodes;

  }
  
  // **********************************************
  //  single processor variant of the FETI method
  // **********************************************
  if (tlinsol == sfeti){
    top->stop->coarse_local_map (top,out);
    
    long **fetinodes;
    fetinodes = new long* [top->stop->ns];
    for (i=0;i<top->stop->ns;i++){
      fetinodes[i] = new long [top->stop->nnsd[i]];
    }
    
    for (i=0;i<top->stop->ns;i++){
      for (j=0;j<top->stop->nnsd[i];j++){
	if (top->stop->nodmultip[i][j]>1)
	  fetinodes[i][j]=1;
	else
	  fetinodes[i][j]=-1;
      }
    }
    
    seqselnodes *selnodfeti;
    
    selnodfeti = new seqselnodes (top->stop->ns,top->stop->nnsd,fetinodes,
				  top->stop->nodmultip,top->stop->ggnbn,top->stop->icnbnmas,
				  top->stop->tnbn,top->stop->icmultip,top->stop->lnbncn,top->stop->ggnbncn,top->stop->sid,
				  top->stop->md,out,1);
    
    ndof[0]=selnodfeti->prepare_feti (feti.fetiimpl,top,out);
    
    feti.initiate (selnodfeti,top,out);

    feti.assemble_subdom_unknowns (top,out);
    
    //delete [] selnodfeti;

    for (i=0;i<top->stop->ns;i++){
      delete [] fetinodes[i];
    }
    delete [] fetinodes;
  }
  
}

/**
   function solves system of linear algebraic equations
   
   @param gt - general topology
   @param gm - %matrix of the system
   @param lhs - %vector of unknowns
   @param rhs - %vector of right hand side
   
   JK, 17.12.2006
*/
void slesolv::solve_system (gtopology *gt,gmatrix *gm,double *lhs,double *rhs,FILE *out)
{
  time_t bt,et;
  long hod,min;
  double sec = clock();
 
  //  initiation of preconditioners
  //  if necessary
  prec.initiation (gt,gm,out);
  
  switch (tlinsol){
  case gauss_elim:
  case ldl:
  case lu:
  case ll:
  case ill:
  case jacobi:
  case gauss_seidel:
  case cg:
  case bicg:
  case lapack_sol:
  case spdirldl:
  case spdirlu:
  case spdirll:
  case pardisolver:
  case cholmod_solver:{
    bt = time (NULL);

    gm->solve_system (gt,prec,lhs,rhs,out);

    et = time (NULL);
    sec = (clock() - sec) / (double)CLOCKS_PER_SEC;
    hod = (long)sec/3600;  sec -= hod*3600;
    min = (long)sec/60;    sec -= min*60;
    //fprintf (stdout,"\n ---------------------------------------");
    //fprintf (stdout,"\n Consumed time by solver %ld:%ld:%5f    %ld",hod,min,sec,et-bt);
    //fprintf (stdout,"\n ---------------------------------------\n");

    break;
  }
  case conden:{
    gt->gnodes[5].ai=0;
    gt->nidof=gm->n-nbdof;
    gm->nbdof=nbdof;
    gm->solve_system (gt,prec,lhs,rhs,out);
    break;
  }
  case sschur:{

    schur.solve_system (gt,gm,lhs,rhs,out);

    break;
  }
  case sfeti:{
    
    /*
    //ptop->initiation (top,ltg);
    //ptop->numbers_of_all_nodes_on_subdomains (domproc,out);


    gt->stop->compute_multiplicity (out);
    gt->stop->find_boundary_nodes (out);
    gt->stop->rewrite_ltg ();
    //gt->stop->codnum_renumb (gt);
    
    //feti.subdomain_ordering (gt,out);
    
    long **fetinodes;
    fetinodes = new long* [feti.ns];
    for (i=0;i<feti.ns;i++){
      fetinodes[i] = new long [gt->stop->nnsd[i]];
    }
    
    for (i=0;i<feti.ns;i++){
      for (j=0;j<gt->stop->nnsd[i];j++){
	if (gt->stop->ltg[i][j]>-1)
	  fetinodes[i][j]=gt->stop->ltg[i][j];
	else
	  fetinodes[i][j]=-1;
      }
    }
    
    selnodfeti = new seqselnodes (feti.ns,gt->stop->nnsd,fetinodes,gt->stop->md,out,1);

    selnodfeti->nodes_on_master (gt->nn,out);
    selnodfeti->node_multiplicity (out);

    selnodfeti->number_all_dofs (gt,out);
    selnodfeti->ndofn_on_master (gt,out);
    selnodfeti->dof_indicators (gt,out);
    selnodfeti->group_local_nodes (out);

    selnodfeti->dof_feti (out);

    selnodfeti->number_contrib (out);
    selnodfeti->contrib_dofs (gt,out);
    
    feti.initiate (selnodfeti,gt,out);

    feti.assemble_subdom_unknowns (gt,out);

    //ndof=f1->ndof;
    */
    feti.solve_system (gt,gm,lhs,rhs,out);
    
    break;
  }
    
  case saddle_point:{
    
    //sp->solve_system (gt,gm,lhs,rhs,out);
    break;
  }
  case permon_solver:{
    fllopptr->solver_fllop(gm->cr->n, gm->cr->ci, gm->cr->adr, gm->cr->a, rhs, lhs, 1, 5, NULL, NULL, NULL, NULL);
    break;
  }
  default:{
    print_err("unknown type of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
  check_math_err();
}
