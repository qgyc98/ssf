#include "mpi.h"
#include "psolver.h"
#include <string.h>
#include <math.h>
#include <limits.h>

#include "selnodes.h"



psolver::psolver (int np,int mr, char *nameproc, int nl)
{
  //  number of processors
  nproc=np;
  //  my rank = processor id
  myrank=mr;
  // name of processor
  strcpy (procName,nameproc);
  nameLength=nl;
  
  //  number of executed domain
  ndom=0;
  
  
  //  type of domain decomposition method
  tdd = no_par_solver;
  //  type of solver of reduced system of equations
  trssol = (redsystsolver) 0;
  //  type of storage of the matrix of the reduced system
  rsmstor = (storagetype) 0;
  //  type of mesh description
  md = (meshdescription) 0;

  
  zero=1.0e-15;
  
  condfixing = -1;
  nmember = -1;
  //  local to global array
  ltg = NULL;
  //  correspondence between subdomains and processors
  domproc = NULL;
  
  //  list of numbers of nodes on subdomain
  nndom = NULL;
  //  list of numbers of elements on subdomain
  nedom = NULL;
  //  list of numbers of stored matrix entries on subdomain
  negmdom = NULL;

  
  //  object containing the Schur complement method
  schcom = NULL;
  //  object containing the one-level FETI method
  f1 = NULL;
  //  object containing the FETI-DP method
  dpf = NULL;
  //  object containing the layered plate problem
  lpp = NULL;
  //  object containing the parallel conjugate gradient method
  //parcg = NULL;
  //  object containing the boundary version of the parallel conjugate gradient method
  //boundparcg = NULL;
  
  selnodschur = NULL;
  selnodfeti = NULL;

  //  parallel topology
  ptop=NULL;

  //fixing node selection class pointer
  fixnodes = NULL;
  methodcondcor = nomethod;
  typecondcur = notype;
  nmembercur = -1;  
  nuserdefnod = 0;
  userdefnod = NULL;
  typecondsurf = notype;
  nmembersurf = -1;
  nring = -1;
  ring = NULL;
  

  ptopjk = NULL;
  gquantfluxres = NULL;
  
  ssle = new slesolv();
}



psolver::~psolver()
{
  delete [] ltg;
  delete [] domproc;
  delete [] nndom;
  delete [] nedom;
  delete [] negmdom;

  delete schcom;
  delete f1;
  delete dpf;
  delete lpp;
  //delete parcg;
  //delete boundparcg;

  delete selnodschur;
  delete selnodfeti;

  delete ptop;
  delete ssle;
  delete ptopjk;
  if (gquantfluxres){
    for(long i=0; i<ptopjk->tnnp; i++)
      delete [] gquantfluxres[i];
    delete [] gquantfluxres;
  }
}



/**
   function initiates parallel solver
   
   @param top - pointer to the general topology
   
   JK
*/
void psolver::initiate (gtopology *top,int /*argc*/,const char**/*argv*/)
{
  //  allocation of parallel topology
  ptop = new partop (nproc,myrank,ndom,md,procName,nameLength);
  ptopjk = new partopjk (nproc,myrank,md,procName,nameLength);
  
  switch (tdd){
    // *********************
    //  NO PARALLEL SOLVER
    // *********************
  case no_par_solver:{
    break;
  }
 
  case schurcompldd:{
    schcom = new schurcompl (nproc,myrank,ndom);
        
    schcom->trssol    = trssol;
    schcom->rsmstor   = rsmstor;
    schcom->ssle->initiate (ssle);
    schcom->dmstor    = dmstor;
    schcom->nicgsch   = nicg;
    schcom->errcgsch  = errcg;
    schcom->zero      = zero;

    break;
  }
  case fetidd:{
    f1 = new feti1 (nproc,myrank,ndom);
    
    f1->enrbm = enrbm;
    f1->lithr = lithr;
    f1->nicg  = nicg;
    f1->errcg = errcg;
    f1->zero  = zero;
    f1->fetiprecond = fetiprecond;
    break;
  }
  case dpfetidd:{
    dpf = new dpfeti (nproc,myrank,ndom);

    dpf->rsmstor = rsmstor;
    dpf->ssle->initiate (ssle);
    dpf->nicgdpfeti  = nicg;
    dpf->errcgdpfeti = errcg;
    dpf->zero  = zero;
    
    if(fixnodes != NULL){
      delete fixnodes;
    }
    fixnodes = new fixnodesel (nproc,myrank,ndom,md,domproc,procName,nameLength,mespr);
    fixnodes->condfixing = condfixing;
    fixnodes->methodcondcor = methodcondcor;
    fixnodes->typecondcur = typecondcur;
    fixnodes->nmembercur = nmembercur;  
    fixnodes->nuserdefnod  = nuserdefnod;
    fixnodes->userdefnod = userdefnod;
    fixnodes->typecondsurf = typecondsurf;
    fixnodes->nmembersurf = nmembersurf;
    fixnodes->nring = nring;
    fixnodes->ring = ring;
    break;
  }
  case parconjuggrad:{
    //parcg = new parcongrad (nproc,myrank,ndom,mespr);
    
    //parcg->ni = nicg;
    //parcg->err = errcg;
    break;
  }
  case boundparconjuggrad:{
    //boundparcg = new boundparcongrad (nproc,myrank,ndom,mespr);
    //boundparcg->ni = nicg;
    //boundparcg->err = errcg;
    //boundparcg->prec = prec;
    //if(argc > 2){
    //boundparcg->Argc=argc-3;
    // boundparcg->Argv=new char*[boundparcg->Argc];
    //long i;
    //for(i = 0; i < boundparcg->Argc; i++){
    //boundparcg->Argv[i] = new char[1000];
    //strcpy(boundparcg->Argv[i],argv[i+3]);
    //}
    //}
    //else{
    //boundparcg->Argc = 0;
    //boundparcg->Argv = NULL;
    //}
    break;
  }
  case layered_plate:{
    lpp = new lplate (nproc,myrank,ndom,top);
    
    lpp->nicg = nicg;
    lpp->errcg = errcg;
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }
  
  
}

/**
   function reads basic facts about parallel solver of algebraic equations
   
   @param in - input stream
   @param nd - number of subdomain
   @param tsm - storage type of subdomain %matrix
   @param mespr - message printing (yes if 1)

   9.10.2003, JK
*/
void psolver::read (XFILE *in,gtopology *gt,long nd,storagetype tsm,long mes)
{
  //  message printing
  mespr=mes;
  
  //  number of subdomain
  ndom=nd;
  
  //  type of storage of matrix of subdomains
  dmstor = tsm;
  
  //  type of mesh description
  //
  //  all_nodes=1
  //  bound_nodes=2
  //  neg_bound_nodes=3
  //  glob_glued=4
  //  metis=11
  xfscanf (in,"%k%m","mesh_description",&meshdescription_kwdset,(int*)&md);

  //  type of domain decomposition solver
  //
  //  no_par_solver=0
  //  schurcompldd=1
  //  fetidd=5
  //  dpfetidd=10
  //  parconjuggrad=40
  //  boundparconjuggrad=45
  //  layered_plate=50
  //
  xfscanf (in,"%k%m","type_of_domain_decomposition_solver",&domdectype_kwdset,(int*)&tdd);
  
  
  switch (tdd){
    // *********************
    //  NO PARALLEL SOLVER
    // *********************
  case no_par_solver:{
    //  this case is used for parallel homogenization
    //  the master processor computes the macroscale problem
    //  while the slaves computes the microscale problems
    //  therefore no parallel solver is needed
    break;
  }

    // ************************************************************************************
    //  SCHUR COMPLEMENT METHOD - PRIMAL DOMAIN DECOMPOSITION METHOD - CONDENSATION METHOD
    // ************************************************************************************
  case schurcompldd:{
    if (mespr==1)  fprintf (stdout,"\n primal domain decomposition (Schur complement method) is used");
    
    //  type of reduced solver
    //
    //  master_sol=1
    //  paral_sol=2
    xfscanf (in,"%k%m","type_of_reduced_solver",&redsystsolver_kwdset,(int*)&trssol);
    switch (trssol){
    case master_sol:{
      if (mespr==1)  fprintf (stdout,"\n reduced system of equations is solved sequentially on master processor");
      xfscanf (in,"%k%m","storage_of_reduced_system_matrix",&storagetype_kwdset,(int*)&rsmstor);
      
      ssle->read (gt,in,mespr);
      
      break;	
    }
    case paral_sol:{
      if (mespr==1)  fprintf (stdout,"\n reduced system of equations is solved in parallel");
      xfscanf (in,"%k%ld","number_of_iterations",&nicg);
      xfscanf (in,"%k%le","error_of_computation",&errcg);
      break;
    }
    default:{
      par_print_err(myrank+1, procName, "unknown reduced system solver is required",__FILE__, __LINE__, __func__);
    }
    }
    break;
  }
    // **************************************************************************************
    //  FINITE ELEMENT TEARING AND INTERCONNECTING METHOD - DUAL DOMAIN DECOMPOSITION METHOD
    // **************************************************************************************
  case fetidd:{
    if (mespr==1)  fprintf (stdout,"\n dual domain decomposition (Finite Element Tearing and Interconnecting method) is used");
    xfscanf (in,"%k%ld","est_number_of_RBM",&enrbm);
    xfscanf (in,"%k%lf","treshold_lin_dep",&lithr);
    //fprintf (stdout,"\n proc %d  lithr %e",myrank,lithr);
    xfscanf (in,"%k%ld","number_of_iterations",&nicg);
    xfscanf (in,"%k%le","error_of_computation",&errcg);

    //  type of preconditioner
    //  nofetiprecond=0
    //  lumped=1
    //  dirichlet=2
    xfscanf (in,"%k%m","type_of_preconditioning",&fetiprecond_kwdset,&fetiprecond);
    switch (fetiprecond){
    case nofetiprecond:{
      fprintf (stdout,"\n No preconditioning is used for FETI method");
      break;}
    case lumped:{
      fprintf (stdout,"\n Lumped preconditioning is used for FETI method");
      break;
    }
    case dirichlet:{
      fprintf (stdout,"\n Dirichlet preconditioning is used for FETI method");
      break;
    }
    }
    break;
  }
    // ***************************************************************
    //  DUAL-PRIMAL FINITE ELEMENT TEARING AND INTERCONNECTING METHOD
    // ***************************************************************
  case dpfetidd:{
    long i;
    if (mespr==1)  fprintf (stdout,"\n dual-primal finite element tearing and interconnecting method is used");
    
    //  type of reduced solver
    xfscanf (in,"%k%m","type_of_reduced_solver",&redsystsolver_kwdset,(int*)&trssol);
    switch (trssol){
    case master_sol:{
      if (mespr==1)  fprintf (stdout,"\n reduced system of equations is solved sequentially on master processor");
      xfscanf (in,"%k%m","storage_of_reduced_system_matrix",&storagetype_kwdset,(int*)&rsmstor);
      ssle->read (gt,in,mespr);
      xfscanf (in,"%k%ld","number_of_iterations",&nicg);
      xfscanf (in,"%k%le","error_of_computation",&errcg);
      break;
    }
    default:{
      par_print_err(myrank+1, procName, "unknown reduced system solver is required",__FILE__, __LINE__, __func__);
      break;
    }
    }
    xfscanf (in,"%k%m","condensation_of_fixing_nodes",&condfixing_kwdset,&condfixing);
    //fprintf (stdout,"condfixing je %ld\n",condfixing);    
    switch(condfixing){
    case nocondconer:{
      if (mespr==1)  fprintf (stdout,"\n Fixing nodes will not be condensed");
      methodcondcor = nomethod;
      break;
    }
    case automatic:{
      if (mespr==1)  fprintf (stdout,"\n Fixing nodes will be condensed automatically");
      break;
    }
    case userdef:{
      if (mespr==1)  fprintf (stdout,"\n Fixing nodes will be condensed by user definition");
      xfscanf (in,"%k%m","place_of_condensation",&methodcond_kwdset,&methodcondcor);
      if (methodcondcor == cursurf ||  methodcondcor == curvecond){
	if(methodcondcor == curvecond){
	  if (mespr==1)  fprintf (stdout,"\n Additional fixing nodes will be added on boundary curve only");
	
	}
	else{
	  if (mespr==1)  fprintf (stdout,"\n Additional fixing nodes will be added on boundary curve and boundary surfaces");
	}
	xfscanf (in,"%k%m","method_of_condensation_curve",&typecondfixing_kwdset,&typecondcur);
	switch(typecondcur){
	case centroid_fix:{
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into centroid between two fixing nodes");
	  break;
	}
	case nth_memb:{
	  xfscanf (in,"%ld",&nmembercur);
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into %ld nodes between two fixing nodes",nmembercur);  
	  break;
	}
	case all_memb:{
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into all nodes between two fixing nodes");
	  break;
	}
	case rand_memb:{
	  xfscanf (in,"%ld",&nmembercur);
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into %ld nodes between two fixing nodes on random positions",nmembercur);
	  break;
	}
	case n_part_curve:{
	  xfscanf (in,"%ld",&nmembercur);
	  if (mespr==1)  fprintf (stdout,"\n Boundary curve between two fixing nodes will be cut into %ld parts",nmembercur);
	  break;
	}
	case userposdef:{
	  xfscanf (in,"%ld",&nuserdefnod);
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into %ld nodes which user specified in input\n",nuserdefnod);
	  userdefnod = new long[nuserdefnod];
	  for(i = 0; i < nuserdefnod; i++){
	    xfscanf (in,"%ld",&userdefnod[i]);
	  }
	  break;
	}
	}
      }
      if (methodcondcor == cursurf ||  methodcondcor == surfacecond){
	xfscanf (in,"%k%m","method_of_condensation_surface",&typecondfixing_kwdset,&typecondsurf);
	switch(typecondsurf){
	case centroid_fix:{
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into centroid of boundary surface");
	  break;
	}
	case rand_memb:{
	  xfscanf (in,"%ld",&nmembersurf);
	  if (mespr==1)  fprintf (stdout,"\n %ld nodes will be added on boundary surface",nmembersurf);
	  break;
	}
	case all_memb:{
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into all nodes between two fixing nodes");
	  break;
	}
	case n_mark:{
	  xfscanf (in,"%ld",&nmembersurf);
	  break;
	}
	case chose_ring:{
	  xfscanf (in,"%ld",&nring);
	  ring = new long[nring];
	  for(i = 0; i < nring; i++){
	    xfscanf (in,"%ld",&ring[i]);
	  }
	  break;
	}
	case userposdef:{
	  xfscanf (in,"%ld",&nuserdefnod);
	  if (mespr==1)  fprintf (stdout,"\n New fixing nodes will be added into %ld nodes which user specified in input\n",nuserdefnod);
	  userdefnod = new long[nuserdefnod];
	  for(i = 0; i < nuserdefnod; i++){
	    xfscanf (in,"%ld",&userdefnod[i]);
	  }
	  break;
	}
	}
      }
      break;
    }
    }
    break;
  }
    
    // *************************************
    //  PARALLEL CONJUGATE GRADIENT METHOD
    // *************************************
  case parconjuggrad:{
    if (mespr==1)  fprintf (stdout,"\n parallel conjugate gradient method is used");
    
    xfscanf (in,"%k%ld","number_of_iterations",&nicg);
    xfscanf (in,"%k%le","error_of_computation",&errcg);
    
    break;
  }
    // *******************************************************
    // BOUNDARY VERSION OF  PARALLEL CONJUGATE GRADIENT METHOD
    // *******************************************************
  case boundparconjuggrad:{
    if (mespr==1)  fprintf (stdout,"\n boundary version of parallel conjugate gradient method is used");
    
    xfscanf (in,"%k%ld","number_of_iterations",&nicg);
    xfscanf (in,"%k%le","error_of_computation",&errcg);

    //  type of preconditioner
    //  noprec=0
    //  petscilu=1
    //  pardiagprec=2
    //
    xfscanf (in,"%k%m","type_of_preconditioning",&parcgprec_kwdset,&prec);
    switch(prec){
    case noprec:{
      if (mespr==1) fprintf (stdout,"\n boundary version of parallel conjugate gradient method will not be preconditioned");
      break;
    }
    case petscilu:{
      if (mespr==1) fprintf (stdout,"\n boundary version of parallel conjugate gradient method will be preconditioned by PETSC ILU");
      break;
    }
    case pardiagprec:{
      if (mespr==1) fprintf (stdout,"\n boundary version of parallel conjugate gradient method will be preconditioned by Jacobi preconditioning");
    }
    }
    break;
  }
   


  case layered_plate:{
    
    if (mespr==1)  fprintf (stdout,"\n layered linear plate problem is solved with help of orthonormalization of constraints");
    
    //  data about solver of system of linear equations
    ssle->read (gt,in,mespr);
    
    //  data about conjugate gradient method
    xfscanf (in,"%k%ld","number_of_iterations",&nicg);
    xfscanf (in,"%k%le","error_of_computation",&errcg);
    
    break;
  }
    
    // ****************************************************************
    // ****************************************************************
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }
  
}

/**
   function prints facts about parallel solver
   
   @param out - output file
   
   9.10.2003, JK
*/
void psolver::print (FILE *out)
{
  //  type of mesh description
  fprintf (out,"%d\n",md);
  //  type of domain decomposition solver
  fprintf (out,"%d\n", tdd);
  
  switch (tdd){
  case no_par_solver:{
    //  this case is used for parallel homogenization
    //  the master processor computes the macroscale problem
    //  while the slaves computes the microscale problems
    //  therefore no parallel solver is needed
    break;
  }
    // ************************************************************************************
    //  SCHUR COMPLEMENT METHOD - PRIMAL DOMAIN DECOMPOSITION METHOD - CONDENSATION METHOD
    // ************************************************************************************
  case schurcompldd:{
    
    //  type of reduced solver
    fprintf (out,"%d ",(int)trssol);
    
    switch (trssol){
    case master_sol:{
      fprintf (out,"%d\n",(int)rsmstor);
      
      ssle->print (out);
      fprintf (out,"\n");
      
      break;	
    }
    case paral_sol:{
      fprintf (out,"%ld %e\n",nicg,errcg);
      break;
    }
    default:{
      par_print_err(myrank+1, procName, "unknown reduced system solver is required",__FILE__, __LINE__, __func__);
    }
    }
    
    break;
  }
    // **************************************************************************************
    //  FINITE ELEMENT TEARING AND INTERCONNECTING METHOD - DUAL DOMAIN DECOMPOSITION METHOD
    // **************************************************************************************
  case fetidd:{
    fprintf (out,"%ld %e ",enrbm,lithr);
    fprintf (out,"%ld %e\n",nicg,errcg);
    fprintf (out,"%ld\n",fetiprecond);
    switch (fetiprecond){
    case nofetiprecond:{
      fprintf (stdout,"\n No preconditioning is used for FETI method");
      break;}
    case lumped:{
      fprintf (stdout,"\n Lumped preconditioning is used for FETI method");
      break;
    }
    case dirichlet:{
      fprintf (stdout,"\n Dirichlet preconditioning is used for FETI method");
      break;
    }
    }
    break;
  }
    // ***************************************************************
    //  DUAL-PRIMAL FINITE ELEMENT TEARING AND INTERCONNECTING METHOD
    // ***************************************************************
  case dpfetidd:{
    //  type of reduced solver
    fprintf (out,"%d\n",(int)trssol);
    switch (trssol){
    case master_sol:{
      fprintf (out,"%d\n",(int)rsmstor);
      ssle->print (out);
      fprintf (out,"%ld %e\n",nicg,errcg);
      break;    
    }
    default:{
      par_print_err(myrank+1, procName, "unknown reduced system solver is required",__FILE__, __LINE__, __func__);
    }
    }
    long i;
    fprintf (out,"%ld\n",condfixing);
    if(condfixing == userdef){
      fprintf (out,"%ld\n",methodcondcor);
      if (methodcondcor == cursurf ||  methodcondcor == curvecond){
	fprintf (out,"%ld\n",typecondcur);
	if(typecondcur == nth_memb || typecondcur == rand_memb || typecondcur == n_part_curve){
	  fprintf (out,"%ld\n",nmembercur);
	}
	if(typecondcur == userposdef){
	  fprintf (out,"%ld\n",nuserdefnod);
	  for(i = 0; i < nuserdefnod; i++){
	    fprintf (out,"%ld   ",userdefnod[i]);
	  }
	  fprintf (out,"\n");
	}
	break;
      }
      if (methodcondcor == cursurf ||  methodcondcor == surfacecond){
	fprintf (out,"%ld",typecondsurf);
	switch(typecondsurf){
	case rand_memb:{
	  fprintf (out,"   %ld\n",nmembersurf);
	  break;
	}
	case n_mark:{
	  fprintf (out,"   %ld",nmembersurf);
	  break;
	}
	case chose_ring:{
	  fprintf (out,"   %ld    ",nring);
	  for(i = 0; i < nring; i++){
	    fprintf (out,"   %ld",ring[i]);
	  }
	  fprintf (out,"\n");
	  break;
	}
	case userposdef:{
	  fprintf (out,"   %ld\n",nuserdefnod);
	  for(i = 0; i < nuserdefnod; i++){
	    fprintf (out,"   %ld   ",userdefnod[i]);
	  }
	  fprintf (out,"\n");
	  break;
	}
	}
      }
    }
    
    
    break;
  }
    // *************************************
    //  PARALLEL CONJUGATE GRADIENT METHOD
    // *************************************
  case parconjuggrad:{
    
    fprintf (out,"%ld %le",nicg,errcg);
    
    break;
  }
    // *******************************************************
    // BOUNDARY VERSION OF  PARALLEL CONJUGATE GRADIENT METHOD
    // *******************************************************
  case boundparconjuggrad:{
    
    fprintf (out,"%ld %le",nicg,errcg);
    fprintf (out,"  %ld",prec);
    
    break;
  }
    
  case layered_plate:{
    //  data about solution strategy
    ssle->print (out);
    fprintf (out,"\n");
    
    //  data about conjugate gradient method
    fprintf (out,"%ld %e\n",nicg,errcg);
    
    break;
  }
    
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }
  
  fprintf (out,"\n");
  
}


/**
   function reads local to global correspondence of nodes
   
   @param in - input file
   @param top - pointer to the sequential general toplogy
   @param out - output file (used for auxiliary uotputs)
   
   JK, 20.1.2002
*/
void psolver::read_ltg (XFILE *in,gtopology *top,FILE */*out*/)
{
  long i;
  
  if (tdd != no_par_solver){
    
    ltg = new long [top->nn];
    
    for (i=0;i<top->nn;i++){
      xfscanf (in,"%ld",ltg+i);
      
      if (ltg[i]<0)
	ltg[i]=0-ltg[i];
      
      ltg[i]--;
    }
  }
  
}

/**
   function establishes correspondence among processors and subdomains
   
   domproc[i]=j - the i-th processor contains the j-th subdomain

   JK
*/
void psolver::procdomcorr ()
{
  long i,j;
  MPI_Status stat;
  
  if (myrank==0){
    domproc = new long [nproc];
    domproc[0]=ndom;
    for (i=1;i<nproc;i++){
      MPI_Recv (&j,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      domproc[stat.MPI_TAG]=j;
    }
  }
  else{
    MPI_Send (&ndom,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
}

/*
void psolver::nodesplit (gtopology *top,FILE *out)
{
  ptop->initiation (top,ltg);
  ptop->numbers_of_all_nodes_on_subdomains (domproc,out);
  ptop->compute_multiplicity (ltg,domproc,out);
  ptop->find_boundary_nodes (ltg,domproc,out);
  ptop->boundary_nodes_on_master (top,domproc,out);
}
*/

/**
   function generates ordering of unknowns with respect
   to the selected type of domain decomposition method
   
   @param top - pointer to the sequential general topology
   @param out - output file (used for auxiliary output)
   @param proc_name - processor name
   
   JK, 11.7.2008 - revised
*/
long psolver::ordering (gtopology *top,FILE *out,char *proc_name)
{
  long i,ndof;
  long *schurnodes,*fetinodes;
  
  top->ordering = new long [top->nn];
  
  //  switch over renumbering type
  switch (top->nodren){
  case no_renumbering:{
    for (i=0;i<top->nn;i++){
      top->ordering[i]=i;
    }
    break;
  }
  case cuthill_mckee:
  case rev_cuthill_mckee:{
    top->cuthill_mckee_renumb (out);
    break;
  }
  case sloan:{
    top->sloan_renumb(out);
    break;
  }
  default:{
    print_err("unknown type of renumbering is required",__FILE__,__LINE__,__func__);
  }
  }


  switch (tdd){
    // *********************
    //  NO PARALLEL SOLVER
    // *********************
  case no_par_solver:{
    ndof = top->codenum_generation (out);
    break;
  }
    
    // **************************
    //  SCHUR COMPLEMENT METHOD
    // **************************
  case schurcompldd:{
    /*
    ptop->prepare_data (domproc,ltg,top,out,proc_name);
    
    schurnodes = new long [top->nn];
    
    for (i=0;i<top->nn;i++){
      if (ptop->nodmultip[i]>1)
	schurnodes[i]=1;
      else
	schurnodes[i]=-1;
    }

    selnodschur = new selnodes (nproc,myrank,md,domproc,ptop->nnsd,schurnodes,
				ptop->tnbn,ptop->icmultip,ptop->nodmultip,ptop->icnbn,ptop->gnbndom,
				ptop->gnbncn,ptop->sid,
				out,proc_name,mespr);

    selnodschur->node_coarse_numbers (out);
    selnodschur->number_all_dofs (top,domproc,out);

    selnodschur->ndofn_on_master (top,domproc,out);
    selnodschur->dof_indicators (top,domproc,out);

    if (md==all_nodes){
      selnodschur->schur_ordering (top,ptop->dofind,out);
    }
    if (md==bound_nodes){
      selnodschur->schur_ordering (top,out);
    }

    ndof = ptop->schur_ordering (top,domproc,out,proc_name);

    
    schcom->initiate (ptop,selnodschur);
    
    //  this must be used for sparse direct solver
    for (i=0;i<top->nn;i++){
      top->gnodes[i].ai=ltg[i];
    }

    delete [] schurnodes;
    */
    
    //  25.8.2011
    ndof = ptopjk->vse (ltg,top,domproc,out,proc_name);
    schcom->initiate (ptopjk,out);
    
    break;
  }
    
    // **************
    //  FETI METHOD
    // **************
  case fetidd:{
    
    ptop->initiation (top,ltg);
    ptop->numbers_of_all_nodes_on_subdomains (domproc,out);
    ptop->compute_multiplicity (ltg,domproc,out);
    ptop->find_boundary_nodes (ltg,domproc,out);
    ptop->rewrite_ltg (ltg);

    
    //  nova verze
    //ptop->initiation (top,ltg);
    //ptop->numbers_of_all_nodes_on_subdomains (domproc,out);
    //ptop->assemble_multip (ltg,domproc,out,proc_name);
    //ptop->node_coarse_numbers (domproc,out,proc_name);
    //  konec nove verze

    f1->subdomain_ordering (top,ltg,out);
    fetinodes = new long [top->nn];

    fprintf (out,"\n\n kontrola fetinodes");
    for (i=0;i<top->nn;i++){
      if (ltg[i]>-1)
	fetinodes[i]=ltg[i];
      else
	fetinodes[i]=-1;
      
      fprintf (out,"\n fetinodes %ld   %ld    ltg %ld",i,fetinodes[i],ltg[i]);
    }
    

    selnodfeti = new selnodes (nproc,myrank,ndom,top->nn,fetinodes,md,out,mespr);
    selnodfeti->number_of_selected_nodes (domproc,out);
    selnodfeti->nodes_on_master (domproc,out);
    selnodfeti->node_multiplicity (out);
    selnodfeti->number_all_dofs (top,domproc,out);
    selnodfeti->ndofn_on_master (top,domproc,out);
    selnodfeti->dof_indicators (top,domproc,out);
    selnodfeti->group_local_nodes (out);
    selnodfeti->dof_feti (out);
    selnodfeti->number_contrib (domproc,out);
    selnodfeti->contrib_dofs (top,domproc,out);
    
    // dofmultip for scaling
    selnodfeti->dof_multiplicity (out);
      
    f1->initiate (selnodfeti,out);
    
    ndof=f1->ndof;

    delete [] fetinodes;
    break;
  }
    //
    //  FETI-DP METHOD
    //
  case dpfetidd:{
    
    ptop->prepare_data (domproc,ltg,top,out,proc_name);
    bool selection;
    // checking of array ltg - fixing nodes
    selection=fixnodes->check_ltg(ltg,ptop->nn,out);
    switch(selection){
    case true:{
      // top - pointer to class gtopology
      // ptop - pointer to class partop
      fixnodes->initiate(top,ptop,out);
      // array ltg - local to global corespondence
      fixnodes->fixing_detection(ltg);
      // rewrite mesh descripton
      md = bound_nodes;
      delete fixnodes;
      break;
    }
    case false:{
      delete fixnodes;
      break;
    }
    }
   
    
    //  this must be used for sparse direct solver
    for (i=0;i<top->nn;i++){
      if (ltg[i]>-2)
    	top->gnodes[i].ai=-1;
      else
    	top->gnodes[i].ai = 0-ltg[i]-2;
    }
    
    fprintf (out,"\n\n kontrola pole ltg");
    for (i=0;i<top->nn;i++){
      fprintf (out,"\n ltg %ld   %ld",i+1,ltg[i]);
    }
    fprintf (out,"\n");
    
    dpf->dpfetiordering (top,ltg);
    ndof=dpf->ndof;
    
    fetinodes = new long [top->nn];
    
    //fprintf (out,"\n\n kontrola fetinodes");
    for (i=0;i<top->nn;i++){
      if (ltg[i]>-1)
    	fetinodes[i]=ltg[i];
      else
    	fetinodes[i]=-1;
      
      //fprintf (out,"\n fetinodes %ld   %ld    ltg %ld",i,fetinodes[i],ltg[i]);
    }
    
    selnodfeti = new selnodes (nproc,myrank,ndom,top->nn,fetinodes,md,out,mespr);
    selnodfeti->number_of_selected_nodes (domproc,out);
    selnodfeti->nodes_on_master (domproc,out);
    selnodfeti->node_multiplicity (out);
    selnodfeti->number_all_dofs (top,domproc,out);
    selnodfeti->ndofn_on_master (top,domproc,out);
    selnodfeti->dof_indicators (top,domproc,out);
    selnodfeti->group_local_nodes (out);
    selnodfeti->dof_feti (out);
    selnodfeti->number_contrib (domproc,out);
    selnodfeti->contrib_dofs (top,domproc,out);
    
    // dofmultip for scaling
    selnodfeti->dof_multiplicity (out);
    
    fprintf (out,"\n\n\n KONEC FETI, ZACATEK SCHURA \n\n\n");
    fprintf (out,"\n\n kontrola schurnodes");
    
    schurnodes = new long [top->nn];
    for (i=0;i<top->nn;i++){
      if (ltg[i]<-1)
    	schurnodes[i]=0-ltg[i]-2;
      else
    	schurnodes[i]=-1;
      //fprintf (out,"\n schurnodes %ld   %ld    ltg %ld",i,schurnodes[i],ltg[i]);
    }
    selnodschur = new selnodes (nproc,myrank,ndom,top->nn,schurnodes,md,out,mespr);
    selnodschur->number_of_selected_nodes (domproc,out);
    selnodschur->nodes_on_master (domproc,out);
    selnodschur->node_multiplicity (out);
    selnodschur->number_all_dofs (top,domproc,out);

    selnodschur->ndofn_on_master (top,domproc,out);
    selnodschur->dof_indicators (top,domproc,out);


//     if (md==all_nodes){
//       selnodschur->schur_ordering (top,ptop->dofind,out);
//     }
    //if (md==bound_nodes){
      //selnodschur->schur_ordering (top,out);
    //}
    
//     dpf->initiate (selnodschur,selnodfeti);
//     delete [] schurnodes;
    delete [] fetinodes;
    break;
  }
  case parconjuggrad:{
    
    //ptop->prepare_data (domproc,ltg,top,out,proc_name);
    //ptop->assemble_gcnd(top,domproc,out);
    
    //parcg->initiate(ptop,top,out);
    //ndof=parcg->ndof;
    break;
  }
    

  case boundparconjuggrad:{
    ptop->prepare_data (domproc,ltg,top,out,proc_name);


    schurnodes = new long [top->nn];


    for (i=0;i<top->nn;i++){
      if (ptop->nodmultip[i]>1)
	schurnodes[i]=1;
      else
	schurnodes[i]=-1;
    }

    selnodschur = new selnodes (nproc,myrank,md,domproc,ptop->nnsd,schurnodes,
				ptop->tnbn,ptop->icmultip,ptop->nodmultip,ptop->icnbn,ptop->gnbndom,
				ptop->gnbncn,ptop->sid,
				out,proc_name,mespr);

    selnodschur->node_coarse_numbers (out);
    selnodschur->number_all_dofs (top,domproc,out);

    selnodschur->ndofn_on_master (top,domproc,out);
    selnodschur->dof_indicators (top,domproc,out);

    
    if (md==all_nodes){
      selnodschur->schur_ordering (top,ptop->dofind,out);
    }
    if (md==bound_nodes){
      selnodschur->schur_ordering (top,out);
    }

    ndof = ptop->schur_ordering (top,domproc,out,proc_name);

    //boundparcg->initiate (selnodschur,ptop,top,out);

    delete [] schurnodes;
    break;
  }
  case layered_plate:{
    lpp->ndof = top->gencodnum ();
    lpp->globcnnum_lpp (top,domproc,out);
    ndof=lpp->ndof;
    break;
  }
  default:{
    par_print_err(myrank,"unknown type of domain decomposition solver is required",__FILE__,__LINE__,__func__);
  }
  }
  
  /*
  long i;
  for (i=0;i<top->nn;i++){
    
    switch (tdd){
    case schurcompldd:{
      top->gnodes[i].ai = ltg[i]+1;
      break;
    }
    case fetidd:{
      top->gnodes[i].ai = ltg[i]+1;
      break;
    }
    case dpfetidd:{
      if (ltg[i]+1>=0)   top->gnodes[i].ai = 0;
      if (ltg[i]+1<0)    top->gnodes[i].ai = 0-ltg[i]-1;
      break;
    }
    case parconjuggrad:{
      top->gnodes[i].ai = ltg[i]+1;
      break;
      }
    case layered_plate:{
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of domain decomposition is required in function read (file %s, line %d).\n",__FILE__,__LINE__);
      abort ();
    }
    }
    
  }
  */
  
  
  //ptop->dof_multiplicity (top,out);

  return ndof;
}

/**
   function assembles constraint matrix
   function is used in layered static analysis
   
   @param th - array containing thicknesses in nodes
   @param out - output stream
   
   JK, 4.3.2003
*/
void psolver::constr_mat (double *th,FILE *out)
{
  switch (tdd){
  case layered_plate:{
    lpp->constraintmat (th,domproc,out);
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }

}

/**
   function solves system of linear algebraic equations in parallel
   
   @param top - pointer to the general topology describing appropriate subdomain
   @param gm - pointer to the %matrix of subdomain
   @param lhs - array containing solution
   @param rhs - array containing right hand side
   @param out - output file (for auxiliary output)
   @param mespr - message printing indicator
   
   JK, revised 24.7.2008
*/
void psolver::par_linear_solver (gtopology *top,gmatrix *gm,
				 double *lhs,double *rhs,FILE *out,long /*mespr*/)
{
  switch (tdd){
  case schurcompldd:{
    schcom->solve_system (top,gm,domproc,lhs,rhs,out);
    break;
  }
  case fetidd:{
    f1->solve_system (top,gm,domproc,lhs,rhs,out);
    break;
  }
  case dpfetidd:{
    /*long i,j;
    fprintf (out,"matice podoblasti\n");
    fprintf (out,"%ld \n",gm->sky->negm);
    for (i=0;i<gm->sky->n;i++){
      for (j=gm->sky->adr[i];j<gm->sky->adr[i+1];j++){
	fprintf (out,"%ld ",i-j+gm->sky->adr[i]+1);
	fprintf (out,"%ld ",i+1);
	fprintf (out,"%le \n",gm->sky->a[j]);
      }
    }
    fprintf (out,"konec matice podoblasti\n");*/
    //dpf->solve_system (top,gm,domproc,lhs,rhs,out,mespr);
    break;
  }
  case parconjuggrad:{
    //parcg->solve_system (ptop,top,gm,domproc,lhs,rhs,out,0);
    break;
  }
  case boundparconjuggrad:{
    //boundparcg->solve_system (ptop,top,gm,domproc,lhs,rhs,out,0);
    break;
  }
  case layered_plate:{
    lpp->solve_system (top,gm,domproc,lhs,rhs,out);
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }

}

/**
   function gathers contributions of boundary parts of local vectors
   (defined on subdomains) into one global vector (defined on boundaries)
   the global vector must be deleted outside of this subroutine
   
   @param lv - local vector
   @param gv - global vector, it is allocated only on the master processor
   
   JK, 6.11.2004
*/
void psolver::gather_bound_vect (double *lv,double *gv)
{
  switch (tdd){
  case schurcompldd:{
    schcom->gather_bound_vect (lv,gv,domproc);
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }

}

/**
   function scatters contributions from global vector (defined on boundaries)
   into local vectors (defined on subdomains)
   function adds contributions to the local vectors
   
   @param lv - local vector
   @param gv - global vector, it is allocated only on the master processor
   
   JK, 6.11.2004
*/
void psolver::scatter_bound_vect (double *lv,double *gv)
{
  switch (tdd){
  case schurcompldd:{
    schcom->scatter_bound_vect (lv,gv,domproc);
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }

}

/**
   function computes norm of vector of unbalanced values
   function gathers vector of unbalanced values on boundaries
   function scatters correct values into subdomains
   function returns norm of vector of unbalanced values
   
   @param lv - local vector
   
   JK, 6.11.2004
*/
double psolver::unbalanced_values (double *lv,FILE *out)
{
  double norm;
  
  switch (tdd){
  case schurcompldd:{
    norm = schcom->unbalanced_values (lv,domproc,out);
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }
  
  return norm;
}

/**
   function computes scalar product (dot product) in parallel
   
   @param lv1,lv2 - arrays containing vectors of appropriate subdomain
   @param out - output file (for auxiliary output)
   
   JK, 24.7.2008
*/
double psolver::pss (double *lv1,double *lv2,FILE */*out*/)
{
  long i,n;
  double dp=0.0;
  double *gv1=NULL,*gv2=NULL;
  MPI_Status stat;
  
  switch (tdd){
  case schurcompldd:{
    if (myrank==0){
      //  number of components of the global vector
      n=schcom->ndofcp;
      //  allocation of the global vectors
      gv1 = new double [n];
      gv2 = new double [n];
      nullv (gv1,n);
      nullv (gv2,n);
      //  assembling of the global vector
      //  dp contains contributions from the internal DOFs
      dp=schcom->pss_gather_bound_vect (lv1,gv1,lv2,gv2,domproc);
      //  computation of contributions from the interface DOFs
      for (i=0;i<n;i++){
	dp+=gv1[i]*gv2[i];
      }
      delete [] gv1;
      delete [] gv2;
    }
    else{
      //  assembling of the global vector
      //  dp is useless on slaves, it is equal to zero
      dp=schcom->pss_gather_bound_vect (lv1,gv1,lv2,gv2,domproc);
    }
    //  scattering of scalar product of vector to the slaves
    if (myrank==0){
      for (i=1;i<nproc;i++){
        MPI_Send(&dp,1,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (&dp,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    break;
  }
  default:{
    par_print_err(myrank,"unknown type of domain decomposition is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return dp;
}



/**
  The function computes norm of quantity flux resultants of the whole problem. 
  It may used for the scaling of residual %vector norm in Newton-Raphson procedure.

  @param lv[in] - %vector of local values on thegiven subdomain, the first componet must be
                  the local norm contribution due to internal nodes while quantity flux resultant components 
                  at boundary/interface nodes follows
  @param out[in] - poimter to the opened text file for eventual log output

  @return The function returns the resulting norm of %vector of quantity flux resultant components. 

  Created by Tomas Koudelka, 08.2021
*/
double psolver::compute_quantfluxres_norm(vector &lv, FILE */*out*/)
{
  double norm = 0.0;

  long i, j, k, l, m,nid;
  MPI_Status stat;
  vector gv;
  // total number of all nodes in the problem
  long tnnp = ptopjk->tnnp;
  //  number of components of the local contribution vector from subdomains
  long n = ptopjk->maxnbnwcd*ptopjk->maxndofbn+1;

  switch (tdd){
    case schurcompldd:{
      if (myrank==0){
        if (gquantfluxres == NULL){
          gquantfluxres = new double*[tnnp];
          memset(gquantfluxres, 0, sizeof(*gquantfluxres)*tnnp);
        }
        else{
          for (i=0; i<tnnp; i++){
            if (gquantfluxres[i])
              memset(gquantfluxres[i], 0, sizeof(*gquantfluxres[i])*ptopjk->ndofnm[i]);
          }
        }
        // allocation of global vector
        reallocv(n, gv);
        // assemble contributions from master
        l = 0;
        norm += lv(l);
        l++;
        for (j=0; j<ptopjk->nbnwc; j++){
          nid = ptopjk->mltg[0][ptopjk->bnwcd[0][j]];
          if (gquantfluxres[nid] == NULL){
            gquantfluxres[nid] = new double[ptopjk->ndofnm[nid]];
            memset(gquantfluxres[nid], 0, sizeof(*gquantfluxres[nid])*ptopjk->ndofnm[nid]);
          }
          for(k=0; k<ptopjk->ndofnm[nid]; k++){
            gquantfluxres[nid][k] += lv[l];
            l++;
          }
        }
        // assemble contributions from slaves
        for(i=1; i<nproc; i++){
          nullv(gv);        
          MPI_Recv (gv.a, n, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
          l = 0;
          norm += gv(l); // contribution from internal DOFs
          l++;
          m = domproc[stat.MPI_TAG];
          // collect contribution from boundary/interface nodes
          for (j=0; j<ptopjk->nbnwcd[m]; j++){
            nid = ptopjk->mltg[i][ptopjk->bnwcd[m][j]];
            if (gquantfluxres[nid] == NULL){
              gquantfluxres[nid] = new double[ptopjk->ndofnm[nid]];
              memset(gquantfluxres[nid], 0, sizeof(*gquantfluxres[nid])*ptopjk->ndofnm[nid]);
            }
            for(k=0; k<ptopjk->ndofnm[nid]; k++){
              gquantfluxres[nid][k] += gv[l];
              l++;
            }
          }
        }
        // compute norm contributions from interface nodes
        for (i=0; i<tnnp; i++){
          if (gquantfluxres[i]){
            for(k=0; k<ptopjk->ndofnm[i]; k++){
              norm += gquantfluxres[i][k]*gquantfluxres[i][k];
            }
          }
        }
      }
      else{
        // send partial norm from internal nodes and interface node values to master      
        MPI_Send(lv.a, lv.n, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      
      //  scattering of norm of quantity flux resultant vector to the slaves
      if (myrank==0){
        for (i=1;i<nproc;i++){
          MPI_Send(&norm, 1, MPI_DOUBLE, i, myrank, MPI_COMM_WORLD);
        }
      }
      else{
        MPI_Recv (&norm, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
      }
      MPI_Barrier (MPI_COMM_WORLD);
      break;
    }
    default:{
      par_print_err(myrank,"unknown type of domain decomposition is required",__FILE__,__LINE__,__func__);
      abort();
    }
  }
  
  return norm;
}



/**
   function computes norm of selected components of a %vector
   function computes norm of the appropriate components of the right hand side
   function compares the norm of unbalanced components with required residual
   
   @param cid - component id (id of required component)
   @param err - required residual
   @param thresh - threshold (computer zero)
   @param top - pointer to the subdomain topology
   @param lv1,rhs - local vectors
   @param gv1,grhs - global vectors, they are allocated only on the master processor
                     the arrays must be deleted outside of the subroutine, dimension must be ndofcp
   @param n[in] - dimension of vectors lv1 and rhs
   @param out[in] - pointer to the opened text file for logging messages
   @param norfv[out] - resulting norm of lv %vector of whole domain
   
   JK, 14.7.2011
*/
long psolver::selected_norm_calculation (long cid,double err,double thresh,long sc,gtopology *top,double *lv1,double *gv1,double *rhs,double *grhs,long n,FILE *out, double &norfv)
{
  long stop;
  
  switch (tdd){
  case schurcompldd:{
    stop = schcom->selected_norm_calculation (cid,err,thresh,sc,ptopjk,top,lv1,gv1,rhs,grhs,domproc,n,out, norfv);
    break;
  }
  default:{
    par_print_err(myrank+1, procName, "unknown domain decomposition solver is required",__FILE__, __LINE__, __func__);
  }
  }
  
  return stop;
}


















/**
   function collects informations about computation
   it collects number of nodes, elements and stored %matrix entries on subdomains
   it searches for minimum and maximum of the previously mentioned numbers
   
   JK, 16.8.2007
*/
void psolver::computation_stat (gtopology *top,gmatrix *gm)
{
  long i,j,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize = 3;
  buff = new long [3];

  //  number of nodes on subdomain
  buff[0]=top->nn;
  //  number of finite elements on subdomain
  buff[1]=top->ne;
  //  number of matrix entries
  buff[2]=gm->give_negm ();
  
  if (myrank==0){
    nndom = new long [nproc];
    nedom = new long [nproc];
    negmdom = new long [nproc];
    
    j=domproc[0];
    nndom[j]=buff[0];
    nedom[j]=buff[1];
    negmdom[j]=buff[2];
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      nndom[j]=buff[0];
      nedom[j]=buff[1];
      negmdom[j]=buff[2];
    }
    
    minnn=LONG_MAX;
    maxnn=0;
    for (i=0;i<nproc;i++){
      if (minnn>nndom[i])
	minnn=nndom[i];
      if (maxnn<nndom[i])
	maxnn=nndom[i];
    }
    
    minne=LONG_MAX;
    maxne=0;
    for (i=0;i<nproc;i++){
      if (minne>nedom[i])
	minne=nedom[i];
      if (maxne<nedom[i])
	maxne=nedom[i];
    }
    
    minnegm=LONG_MAX;
    maxnegm=0;
    for (i=0;i<nproc;i++){
      if (minnegm>negmdom[i])
	minnegm=negmdom[i];
      if (maxnegm<negmdom[i])
	maxnegm=negmdom[i];
    }
    
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
}


void psolver::computation_stat_print (FILE *out)
{
  long i;
  
  if (myrank==0){
    fprintf (out,"\n\n\n informations about computation \n");
    
    fprintf (out,"\n i   NN   NE    NEGM");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %3ld %9ld  %9ld  %12ld",i+1,nndom[i],nedom[i],negmdom[i]);
    }

    fprintf (out,"\n");

    fprintf (out,"\n minimum number of nodes           %ld",minnn);
    fprintf (out,"\n maximum number of nodes           %ld",maxnn);

    fprintf (out,"\n minimum number of elements        %ld",minne);
    fprintf (out,"\n maximum number of elements        %ld",maxne);

    fprintf (out,"\n minimum number of matrix entries  %ld",minnegm);
    fprintf (out,"\n maximum number of matrix entries  %ld",maxnegm);
  }
}



/**
   Function compares two double values on master computer and scatters result to 
   slaves. It is used in nonlinear solvers for constraining common decision whether the
   error is in tolerance or not.
   
   @param val - compared value
   @param lim - limit to which the val is compared
   
   @retval  1 - in case that val > lim
   @retval  0 - in case that val == lim
   @retval -1 - in case that val < lim
   TKo, 31.7.2008
*/
long psolver::compare_on_master (double val, double lim)
{
  long ret = -2;
  long i;
  MPI_Status stat;

  if (myrank==0)
  {
    if (val > lim)
      ret = 1;
    else
    {
      if (val < lim)
        ret = -1;
      else
        ret = 0;
    }
    for (i=1;i<nproc;i++)
      MPI_Send(&ret,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
  }
  else
    MPI_Recv(&ret,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

  MPI_Barrier (MPI_COMM_WORLD);
  
  return ret;
}
