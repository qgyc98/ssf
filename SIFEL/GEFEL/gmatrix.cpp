#include "gmatrix.h"
#include "precond.h"

gmatrix::gmatrix ()
{
  //  number of rows/columns = number of unknowns
  n=0;
  //  computer zero
  zero=1.0e-30;
  fprintf (stdout,"\n computer zero in gmatrix is defined as %e",zero);

  //  maximum number of iterations in an iterative method
  ni=0;
  //  required norm of the residual vector in an iterative method
  res=0.0;
  //  actual number of iterations performed in an iterative method
  ani=0;
  //  attained norm of the residual vector in an iterative method
  ares=0.0;

  //  parameter for SSOR preconditioning
  omega=0.0;
  //  parameter gamma for incomplete factorization
  indegamma=0.0;
  //  threshold for selection of matrix entries
  limit=0.0;
  
  
  //  dense matrix storage
  dm=NULL;
  //  diagonal matrix storage (e.g. for explicit dynamics)
  diagm=NULL;
  //  skyline storage
  sky=NULL;
  //  double skyline storage
  dsky=NULL;
  //  compressed row storage
  cr=NULL;
  //  symmetric compressed row storage
  scr=NULL;
  //  compressed column storage
  cc=NULL;
  //  symmetric compressed column storage
  scc=NULL;
  //  element by element storage
  em=NULL;

  //  sparse solver
  sdirect = NULL;
  //  auxiliary array for sparse solver
  auxsd=NULL;
  //  block size for sparse solver
  bsize=0;
  //  indicator of data transfer from cr to sdirect
  pmat=0;
  

  //  defined for LAPACK
  //Matrix = NULL;
  
  //iss = NULL;

  //  PETSC
  //petscmat = NULL;
}

gmatrix::~gmatrix ()
{
  delete dm;
  delete diagm;
  delete sky;
  delete dsky;
  delete cr;
  delete scr;
  delete cc;
  delete scc;
  delete em;

  delete sdirect;
  delete auxsd;
  
  //  defined for LAPACK
  //delete Matrix;

  //delete iss;
}

/**
   function allocates appropriate storage scheme
   
   JK
*/
void gmatrix::alloc ()
{
  switch (ts){
  case dense_matrix:{
    if (dm==NULL)
      dm = new densemat ();
    break;
  }
  case diag_mat:{
    if (diagm==NULL)
      diagm = new diagmat ();
    break;
  }
  case skyline_matrix:{
    if (sky==NULL)
      sky = new skyline ();
    break;
  }
  case double_skyline:{
    if (dsky==NULL)
      dsky = new dskyline ();
    break;
  }
  case compressed_rows:{
    if (cr==NULL)
      cr = new comprow ();
    break;
  }
  case symm_comp_rows:{
    if (scr==NULL)
      scr = new symcomprow ();
    break;
  }
  case compressed_columns:{
    if (cc==NULL)
      cc = new compcol ();
    break;
  }
  case symm_comp_columns:{
    if (scc==NULL)
      scc = new symcompcol ();
    break;
  }
  case element_matrices:{
    if (em==NULL)
      em = new elemmat ();
    break;
  }
    
    
  case lapack_stor:{
    //if (cr==NULL)  cr = new comprow ();
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case petsc:{
    //if (cr==NULL)  cr = new comprow ();
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case spdirect_stor_cr:{
    if (cr==NULL)
      cr = new comprow ();
    if (sdirect==NULL)
      sdirect = new DSSolver();
    //if (iss==NULL) iss = new spasol ();
    break;
  }
    
  case spdirect_stor_scr:{
    if (scr==NULL)
      scr = new symcomprow ();
    if (sdirect==NULL)
      sdirect = new DSSolver();
    break;
  }
    
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
  //  arrays are newly allocated, data were not sent from cr/scr to sparse solver
  pmat=0;
}

/**
   function deallocates storage scheme
*/
void gmatrix::dealloc ()
{
  switch (ts){
  case dense_matrix:{
    delete dm;
    break;
  }
  case skyline_matrix:{
    delete sky;
    break;
  }
  case double_skyline:{
    delete dsky;
    break;
  }
  case compressed_rows:{
    delete cr;
    break;
  }
  case symm_comp_rows:{
    delete scr;
    break;
  }
  case compressed_columns:{
    delete cc;
    break;
  }
  case symm_comp_columns:{
    delete scc;
    break;
  }
  case element_matrices:{
    delete em;
    break;
  }
    
    
  case lapack_stor:{
    //delete cr;
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case petsc:{
    //delete cr;
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case spdirect_stor_cr:{
    delete cr;
    delete sdirect;
    
    //if (iss==NULL) iss = new pardiso ();
    break;
  }
    
  case spdirect_stor_scr:{
    delete scr;
    delete sdirect;
    break;
  }
    
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   initialization of storage schemes
   
   @param top - pointer to the general topology
   @param ndof - number of degrees of freedom = number of unknowns
   @param ats - type of storage of the %matrix
   @param mespr - indicator of message printing
   
   JK
*/
void gmatrix::initiate (gtopology *top,long ndof,storagetype ats,long mespr,FILE *out)
{
  //  definition of the number of rows/columns
  n=ndof;
  //  type of matrix storage
  ts = ats;
  //  allocation of storage schemes
  alloc ();
  
  switch (ts){
  case dense_matrix:{
    dm->initiate (n,mespr);
    break;
  }
  case skyline_matrix:{
    sky->initiate (top,n,mespr);
    break;
  }
  case double_skyline:{
    dsky->initiate (top,n,mespr);
    break;
  }
  case compressed_rows:{
    cr->initiate (top,n,mespr);
    break;
  }
  case symm_comp_rows:{
    scr->initiate (top,n,mespr);
    break;
  }
  case compressed_columns:{
    cc->initiate (top,n,mespr);
    break;
  }
  case symm_comp_columns:{
    scc->initiate (top,n,mespr);
    break;
  }
  case element_matrices:{
    em->initiate (top,n,mespr);
    break;
  }
    
  case lapack_stor:{
    //cr->initiate (top,ndof,mespr);
    //int Error;
    //Matrix = spCreate(ndof, &Error);
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case petsc:{
    //cr->initiate (top,ndof,mespr);
    //PetscInitialize(&argc,&args,(char *)0,help);
    //ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,ndof,ndof,PETSC_DECIDE,PETSC_NULL,&petscmat); CHKERRA(ierr);
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }

  case spdirect_stor_cr:{
    cr->initiate (top,n,mespr);
    break;
  }
  case spdirect_stor_scr:{
    scr->initiate (top,n,mespr);
    break;
  }
  case diag_mat:{
    diagm->initiate (n,mespr);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
  if (ts == spdirect_stor_cr || ts == spdirect_stor_scr){
    if (sdirect->IsInitialized () == 0){
      // redirect the output to default console MathTracer
      sdirect->SetMT();
      sdirect->SetMTMtx();
      //  object sdirect was allocated but it was not initialized
      //  the initialization is performed now
      
      // Makes necessary intialization
      // 0 - produces sparse symmetrical LDL^T factorization (default)
      // 1 - produces sparse symmetrical LLT Cholesky factorization
      // 2 - produces sparse LU factorization on symmetrical pattern -> undirected graph
      // 3 - produces sparse incomplete LDL^T factorization 
      switch (tlinsol){
      case spdirldl:{
	sdirect->Initialize(1,eDSSFactorizationLDLT);
	break;
      }
      case spdirlu:{
	sdirect->Initialize(1,eDSSFactorizationLU);
	break;
      }
      case spdirll:{
	sdirect->Initialize(1,eDSSFactorizationLLT);
	break;
      }
      case cg:{
	if (tprec == sparseindec)
	  sdirect->Initialize(1,eDSSFactorizationLDLTIncomplete);
	break;
      }
      case conden:{
	switch (tsol){
	case spdirldl:{
	  sdirect->Initialize(1,eDSSFactorizationLDLT);
	  break;
	}
	case spdirlu:{
	  sdirect->Initialize(1,eDSSFactorizationLU);
	  break;
	}
	case spdirll:{
	  sdirect->Initialize(1,eDSSFactorizationLLT);
	  break;
	}
	case cg:{
	  if (tprec == sparseindec)
	    sdirect->Initialize(1,eDSSFactorizationLDLTIncomplete);
	  break;
	}
	default:{
	  print_err("unknown type of sparse solver is required",__FILE__,__LINE__,__func__);
	}
	}
	
	break;
      }
	
      default:{
	print_err("unknown type of sparse solver is required",__FILE__,__LINE__,__func__);
      }
      }
      
      
      //  Block size for sparse direct solver
      if (top->gnodes == NULL){
	bsize = (unsigned char) 6;
      }
      else{
	bsize = (unsigned char) top->gnodes[0].ndofn;
      }
      if (bsize<1){
	print_err("wrong block size",__FILE__,__LINE__,__func__);
	abort ();
      }
      
      //  n - number of unknowns
      //  cr->a - array containing the matrix
      //  cr->ci - array of column indices
      //  cr->adr - array of addresses of the first row entries
      //  0,0 - initial indices in arrays cr->a and cr->ci
      //  1 - SIFEL indicator (true)
      //  0 - nonsymmetric matrix
      if (ts == spdirect_stor_cr){
	SparseMatrixF sm((unsigned long)n,cr->a,(unsigned long*)cr->ci,(unsigned long*)cr->adr,0,0,1,0);
	//  function copies pointers from the object sm to sdirect
	sdirect->SetMatrixPattern (&sm,bsize);
      }
      else{
	SparseMatrixF sm((unsigned long)n,scr->a,(unsigned long*)scr->ci,(unsigned long*)scr->adr,0,0,1,1);
	//  function copies pointers from the object sm to sdirect
	sdirect->SetMatrixPattern (&sm,bsize);
      }
      
      if (tlinsol==conden){
	//  if the condensation is required, additional informations have to be sent
	auxdatsparsesolver (top,out);
	
	if (top->nn != 0){
	  sdirect->LoadMCN (top->nn,bsize,auxsd,1);
	}
	
      }
      
      
    }
    else{
      //  function cleans all arrays in the object sdirect
      sdirect->LoadZeros ();
    }
  }
  
}


/**
   function defines necessary parameters for the class %gmatrix from the class %slesolv
   
   @param ssle - object containing informations about solution of system of linear equations
   
   JK
*/
void gmatrix::setval (slesolv *ssle)
{
  //  type of solver of linear system
  tlinsol=ssle->tlinsol;
  //  if condensation is used, additional type of solver is needed
  tsol=ssle->tsol;
  //  type of preconditioner
  tprec=ssle->prec.pt;
  
  //  maximum number of iterations in conjugate gradient method
  ni=ssle->ni;
  //  required residual
  res=ssle->res;
  //  initial vector
  iv=ssle->iv;
  
  //omega=ssle.ssoromega;
  //indegamma=ssle.indegamma;
}


/**
   function localizes local %matrix lm to the global one
   %matrix is stored in the object of the class %matrix

   @param lm - local chracteristic %matrix of the element
   @param cn - array containing code numbers
   @param eid - element id
   @param ndofe - the number of DOfs 
   16.7.2002
*/
void gmatrix::localize (matrix &lm,ivector &cn,long eid)
{
  if (lm.m != cn.n){

  }

  switch (ts){
  case dense_matrix:{
    dm->localize (lm,cn.a);
    break;
  }
  case skyline_matrix:{
    sky->localize (lm,cn.a);
    break;
  }
  case double_skyline:{
    dsky->localize (lm,cn.a);
    break;
  }
  case compressed_rows:{
    cr->localize (lm,cn.a);
    break;
  }
  case symm_comp_rows:{
    scr->localize (lm,cn.a);
    break;
  }
  case compressed_columns:{
    cc->localize (lm,cn.a);
    break;
  }
  case symm_comp_columns:{
    scc->localize (lm,cn.a);
    break;
  }
  case element_matrices:{
    em->localize (lm,cn.a,eid);
    break;
  }

  case lapack_stor:{
    //cr->localize (lm,cn);
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
  case petsc:{
    //cr->localize (lm,cn);
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }

  case spdirect_stor_cr:{
    cr->localize (lm,cn.a);
    break;
  }
  case spdirect_stor_scr:{
    scr->localize (lm,cn.a);
    break;
  }
  case diag_mat:{
    diagm->localize (lm,cn.a);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function localizes local %matrix lm to the global one
   %matrix is stored in array of type double
   
   @param lm - local chracteristic %matrix of the element
   @param cn - array containing code numbers
   @param nc - number of rows and columns
   @param eid - element id
   
   16.7.2002
*/
void gmatrix::localized (double *lm,long *cn,long nc,long eid)
{
  switch (ts){
  case dense_matrix:{
    dm->localized (lm,cn,nc);
    break;
  }
  case skyline_matrix:{
    sky->localized (lm,cn,nc);
    break;
  }
  case double_skyline:{
    dsky->localized (lm,cn,nc);
    break;
  }
  case compressed_rows:{
    cr->localized (lm,cn,nc);
    break;
  }
  case symm_comp_rows:{
    scr->localized (lm,cn,nc);
    break;
  }
  case element_matrices:{
    em->localized (lm,cn,nc,eid);
    break;
  }
  case lapack_stor:{
    //cr->localized (lm,cn,nc);
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
  case petsc:{
    //cr->localized (lm,cn,nc);
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
  case spdirect_stor_cr:{
    cr->localized (lm,cn,nc);
    break;
  }
  case spdirect_stor_scr:{
    scr->localized (lm,cn,nc);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function localizes local %matrix lm to the global one
   %matrix km can be nonsquare, it means with different number
   of rows and columns
   
   @param lm - local chracteristic %matrix of the element
   @param rcn - array containing row code numbers
   @param ccn - array containing column code numbers
   
   14.8.2002
*/
void gmatrix::glocalize (matrix &lm,ivector &rcn,ivector &ccn)
{
  switch (ts){
  case dense_matrix:{
    dm->glocalize (lm,rcn.a,ccn.a);
    break;
  }
  case skyline_matrix:{
    sky->glocalize (lm,rcn.a,ccn.a);
    break;
  }
  case double_skyline:{
    dsky->glocalize (lm,rcn.a,ccn.a);
    break;
  }
  case compressed_rows:{
    //cr->localize (lm,cn);
    break;
  }
  case symm_comp_rows:{
    //scr->localize (lm,cn);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function localizes contributions from Lagrange multipliers to the general %matrix
   
   @param nm - number of Lagrange multipliers
   @param ncn1 - nodal code numbers of the first node (on one side of interface)
   @param ncn2 - nodal code numbers of the second node (on the other side of interface)
   @param mcn - code numbers of Lagrange multipliers defined between the previous nodes

   JK, 8.8.2008
*/
void gmatrix::mult_localize (long nm,long *ncn1,long *ncn2,long *mcn)
{
  switch (ts){
  case dense_matrix:{
    dm->mult_localize (nm,ncn1,ncn2,mcn);
    break;
  }
  case skyline_matrix:{
    sky->mult_localize (nm,ncn1,ncn2,mcn);
    break;
  }
  case double_skyline:{
    dsky->mult_localize (nm,ncn1,ncn2,mcn);
    break;
  }
  default:{
    print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function initializes general %matrix
   
   @param limit - threshold for selection of %matrix entries
   @param mespr - message printing indicator
   
   JK
*/
void gmatrix::prepmat (double /*limit*/,long mespr)
{
  switch (ts){
  case dense_matrix:{
    dm->decompid=0;
    break;
  }
  case skyline_matrix:{
    sky->setnotfact ();
    break;
  }
  case double_skyline:{
    dsky->setnotfact ();
    break;
  }
  case compressed_rows:{
    cr->decompid=0;
    long rejected=0;
    //rejected = cr->minimize (limit);
    if (mespr==1)  fprintf (stdout,"\n number of rejected entries  %ld",rejected);
    if (mespr==1)  fprintf (stdout,"\n number of matrix entries  %ld",cr->negm);
    break;
  }
  case symm_comp_rows:{
    scr->decompid=0;
    long rejected=0;
    //rejected = scr->minimize (limit);
    if (mespr==1)  fprintf (stdout,"\n number of rejected entries  %ld",rejected);
    if (mespr==1)  fprintf (stdout,"\n number of matrix entries  %ld",scr->negm);
    break;
  }
  case compressed_columns:{
    cc->decompid=0;
    if (mespr==1)  fprintf (stdout,"\n number of matrix entries  %ld\n",cc->negm);
    break;
  }

  case symm_comp_columns:{
    scc->decompid=0;
    if (mespr==1)  fprintf (stdout,"\n number of matrix entries  %ld\n",scc->negm);
    break;
  }

  case lapack_stor:{
    //long i,j;
    //ElementPtr pElement; 
    //for (i=0;i<cr->n;i++){
    //for (j=cr->adr[i];j<cr->adr[i+1];j++){
    //pElement = spGetElement(Matrix, i+1, cr->ci[j]+1);
    //pElement->Real = cr->a[j];
    //}
    //}
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case petsc:{
    /*
    long i,j,k,l;
    int ierr,*col,*row;
    Scalar *val;

    row = (int *) PetscMalloc(1*sizeof(int)); CHKPTRA(col);
    for (i=0;i<n;i++){
      l=cr.adr[i+1]-cr.adr[i];

      val = (Scalar *) PetscMalloc(l*sizeof(Scalar)); CHKPTRA(val);
      col = (int *) PetscMalloc(l*sizeof(int)); CHKPTRA(col);
      
      k=0;
      for (j=cr.adr[i];j<cr.adr[i+1];j++){
	val[k]=cr.a[j];
	col[k]=cr.ci[j];
	k++;
      }
      row[0]=i;
      
      ierr = MatSetValues(petscmat,1,&row,l,&col,&val,INSERT_VALUES); CHKERRA(ierr);    
      
      PetscFree(col);  PetscFree(val);
    }
    PetscFree(row);
    
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
    */
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
    
  case spdirect_stor_cr:{
    
    break;
  }
  case spdirect_stor_scr:{
    
    break;
  }
    
    
  case element_matrices:{
    em->decompid=0;
    break;
  }
  case diag_mat:{
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}



/**
   function initializes some sparse solvers
   it copies data stored in the compressed row or
   symmetric compressed row format to a sparse solver
   
   @param top - pointer to general topology
   @param out - output file for auxiliary output
   
*/
void gmatrix::prepmat2 (gtopology */*top*/,FILE */*out*/)
{
  switch (ts){
  case dense_matrix:{
    break;
  }
  case skyline_matrix:{
    break;
  }
  case double_skyline:{
    break;
  }
  case compressed_rows:{
    break;
  }
  case symm_comp_rows:{
    break;
  }
  case compressed_columns:{
    break;
  }
  case symm_comp_columns:{
    break;
  }

  case lapack_stor:{
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
  case petsc:{
    print_err("PETSC solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
    
  case spdirect_stor_cr:{
    if (sdirect->IsAllocated()==0){
      //  arrays in sparse direct solver are not allocated
      
      //  memory for matrix entries is allocated at this place
      //  it computes the symbolic factorization
      sdirect->PreFactorize ();
    }

    //  n - number of unknowns
    //  cr->a - array containing the matrix
    //  cr->ci - array of column indices
    //  cr->adr - array of addresses of the first row entries
    //  0,0 - initial indices in arrays cr->a and cr->ci
    //  1 - SIFEL indicator (true)
    //  0 - nonsymmetric matrix
    SparseMatrixF sm((unsigned long)n,cr->a,(unsigned long*)cr->ci,(unsigned long*)cr->adr,0,0,1,0);
    sdirect->LoadNumbers (&sm);
    
    //delete cr;
    //cr = NULL;
    break;
  }
    
  case spdirect_stor_scr:{
    if (sdirect->IsAllocated()==0){
      //  arrays in sparse direct solver are not allocated

      //  memory for matrix entries is allocated at this place
      //  it computes the symbolic factorization
      sdirect->PreFactorize ();
    }
    
    //  n - number of unknowns
    //  cr->a - array containing the matrix
    //  cr->ci - array of column indices
    //  cr->adr - array of addresses of the first row entries
    //  0,0 - initial indices in arrays cr->a and cr->ci
    //  1 - SIFEL indicator (true)
    //  0 - nonsymmetric matrix
    SparseMatrixF sm((unsigned long)n,scr->a,(unsigned long*)scr->ci,(unsigned long*)scr->adr,0,0,1,1);
    sdirect->LoadNumbers (&sm);

    //delete scr;
    //scr = NULL;
    break;
  }
    
    
  case element_matrices:{
    break;
  }
  case diag_mat:{
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
  //  data from cr/scr were sent to sparse solver
  pmat=1;
}


/**
   function assembles data for sparse direct solver
   especially data for condensation are assembled
   
   @param top - pointer to the general topology
   @param out - output file
   
   JK
*/
void gmatrix::auxdatsparsesolver (gtopology *top,FILE *out)
{
  long i,j,k,m,nn,ndofn,rndofn;
  
  
  fprintf (out,"\n\n\n\n kontrola Richarda \n\n");
  for (i=0;i<top->nn;i++){
    fprintf (out,"\n %ld  %ld",i,top->gnodes[i].ai);
  }
  
  /*
  fprintf (stdout,"\n\n\n\n kontrola Richarda \n\n");
  for (i=0;i<top->nn;i++){
    fprintf (stdout,"\n %ld  %ld",i,top->gnodes[i].ai);
  }
  fflush(stdout);
  */

  if (top->gnodes != NULL){
    
    nn = top->nn;
    ndofn = top->gnodes[0].ndofn;
    
    //fprintf (stderr,"\n\n pocet nodes na domene %d",nn);
    //fprintf (stderr,"\n pocet elements na domene %d\n\n",top->ne);
    
    auxsd = new long [nn*ndofn];
    
    k=0;
    for (i=0;i<nn;i++){
      rndofn=top->gnodes[i].ndofn;
      if (rndofn!=ndofn){
	print_err("wrong number of degrees of freedom of node is required",__FILE__,__LINE__,__func__);
	abort ();
      }
      
      if (top->gnodes[i].ai>-1){
	//  boundary nodes
	for (j=0;j<ndofn;j++){
	  m=top->give_dof (i,j);
	  if (m>-1){
	    //  unknowns
	    auxsd[k] = -1 - top->give_dof (i,j);
	    k++;
	  }
	  if (m<0){
	    //  prescribed values
	    auxsd[k] = -1;
	    k++;
	  }
	}
      }
      else{
	//  inner nodes
	for (j=0;j<ndofn;j++){
	  m=top->give_dof (i,j);
	  if (m>-1){
	    //  unknowns
	    auxsd[k] = m-1;
	    k++;
	  }
	  if (m<0){
	    //  prescribed values
	    auxsd[k] = -1;
	    k++;
	  }
	}
      }
    }
    
  }
  
  fprintf (out,"\n\n\n\n kontrola Richarda \n\n");
  k=0;
  for (i=0;i<top->nn;i++){
    fprintf (out,"\n uzel %5ld",i);
    for (j=0;j<ndofn;j++){
      fprintf (out,"  %ld",auxsd[k]);
      k++;
    }
  }
  fprintf (out,"\n\n\n"); 
}



/**
   function solves system of linear algebraic equations
   
   @param top - pointer to the general topology
   @param prec - preconditioner
   @param lhs - array containing solution of the system
   @param rhs - array containing the right hand side
   @param out - output file for auxiliary output
   
   JK
*/
void gmatrix::solve_system (gtopology *top,precond &prec,double *lhs,double *rhs,FILE *out)
{
  //time_t bt = time (NULL), et;
  
  if (pmat==0)
    prepmat2 (top,out);

  switch (tlinsol){
    
  case gauss_elim:{
    switch (ts){
    case dense_matrix:{
      dm->gemp (lhs,rhs,1,zero,2);
      break;
    }
    case diag_mat:{
      diagm->gemp (lhs,rhs,zero);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case diagsolv:{
    switch (ts){
    case diag_mat:{
      diagm->gemp (lhs,rhs,zero);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }

  case ll:{
    switch (ts){
    case dense_matrix:{
      dm->ll (lhs,rhs,zero,1);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case ill:{
    switch (ts){
    case dense_matrix:{
      dm->ill (lhs,rhs,zero,zero,1);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }

  case ldl:{
    switch (ts){
    case dense_matrix:{
      if (dm->decomp()==0){
        dm->ll (NULL,NULL,zero,2);
      }
      dm->ll (lhs,rhs,zero,3);
      break;
    }
    case skyline_matrix:{
      if (sky->decomp()==0){
        sky->ldl_sky (NULL,NULL,zero,2);
      }
      sky->ldl_sky (lhs,rhs,zero,3);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case lu:{
    switch (ts){
    case dense_matrix:{
      if (dm->decomp()==0){
        dm->lu (lhs,rhs,zero,2);
      }
      dm->lu (lhs,rhs,zero,3);
      break;
    }
    case double_skyline:{
      if (dsky->decomp()==0){
        dsky->lu_dsky (NULL,NULL,zero,2);
      }
      dsky->lu_dsky (lhs,rhs,zero,3);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case conden:{
    double *condvect,*av;
    densemat adm;
    
    condvect=new double[nbdof];
    av=new double[nbdof];
    
    adm.alloc(nbdof);

    condense (top,adm.a,condvect,lhs,rhs,nbdof,1,out);
    
    adm.gemp (av,condvect,1,1.0e-10,1);

    copyv (av,condvect,nbdof);

    condense (top,adm.a,condvect,lhs,rhs,nbdof,2,out);

    delete [] condvect;
    delete [] av;
    adm.dealloc ();

    break;
  }

  case cg:{
    switch (ts){
    case dense_matrix:{
      switch (tprec){
      case noprecond:{
	dm->cg (lhs,rhs,ni,res,ani,ares,zero,iv);
	break;
      }
      case incomdec:{
        fprintf (stdout,"\n incomplete decomposition of matrix");
	densemat dmid;
	dmid.copy (dm);
	
        if (dmid.decomp()==0){
          dmid.ill (lhs,rhs,zero,zero,2);
          //dmid.changedecomp ();
        }
        
        //dm->cg_prec (dmid,lhs,rhs,nicg,errcg,anicg,aerrcg,zero,0,tprec,0.0);
        dm->cg_prec_new (prec,lhs,rhs,ni,res,ani,ares,zero,iv);
        break;
      }
      default:{
	print_err("unknown preconditioner is required",__FILE__,__LINE__,__func__);
      }
      }

      break;
    }
    case compressed_rows:{
      cr->cg_prec (prec,lhs,rhs,ni,res,ani,ares,zero,iv);
      //cr->cg (lhs,rhs,nicg,errcg,anicg,aerrcg,zero,0);
      //cr->cg_cr_rev (lhs,rhs,nicg,errcg,anicg,aerrcg,zero,0);
      //cr->cg_new (lhs,rhs,nicg,errcg,anicg,aerrcg,zero,0);
      break;
    }
    case symm_comp_rows:{
      switch (tprec){
      case noprecond:{
        scr->cg (lhs,rhs,ni,res,ani,ares,zero,iv);
        break;
      }
      case diagprec:{
        scr->cg_prec (lhs,rhs,ni,res,ani,ares,zero,iv,tprec,0.0,sdirect);
        break;
      }
      case ssorprec:{
        scr->cg_prec (lhs,rhs,ni,res,ani,ares,zero,iv,tprec,omega,sdirect);
        break;
      }
      case incomdec:{
        fprintf (stdout,"\n incomplete decomposition of matrix");
        if (scr->decomp()==0){
          scr->incomplete_ldl (indegamma);
          scr->changedecomp ();
        }
        
        scr->cg_prec (lhs,rhs,ni,res,ani,ares,zero,iv,tprec,0.0,sdirect);
        break;
      }
      case sparseindec:{
        fprintf (stdout,"\n incomplete sparse decomposition of matrix");
	
        if (sdirect->IsFactorized()==0){
          sdirect->ReFactorize ();
        }
	
        scr->cg_prec (lhs,rhs,ni,res,ani,ares,zero,iv,tprec,0.0,sdirect);
	break;
      }

      default:{
        fprintf (stderr,"\n\n unknown preconditioner type is required in function cgsol (file %s, line %d).\n",__FILE__,__LINE__);
      }
      }
      
      break;
    }
    case element_matrices:{
      em->cg (lhs,rhs,ni,res,ani,ares,zero,iv);
      break;
    }

    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }

  case bicg:{
    switch (ts){
    case double_skyline:{
      dsky->bicg (lhs,rhs,ni,res,ani,ares,zero,iv);
      break;
    }
    case compressed_rows:{
      cr->bicg (lhs,rhs,ni,res,ani,ares,zero,iv);
      //cr->bicg_new (lhs,rhs,nicg,errcg,anicg,aerrcg,zero,0);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case lapack_sol:{
    //int Error;
    //Error = spFactor(Matrix);
    //spSolve( Matrix, rhs, lhs);
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case spdirldl:
  case spdirlu:
  case spdirll:{
    if (sdirect->IsFactorized()==0){
      sdirect->ReFactorize();
    }
    sdirect->Solve(lhs,rhs);
    mute();
    break;
  }
    
    
  case pardisolver:{
    //  pozor na cislovani od 1
    //    printf("JSEM TU 1");
    //    fflush(stdout);
    //    iss->symbfact (cr->a,cr->ci,cr->adr,cr->n);
    //    printf("JSEM TU 2");
    //    fflush(stdout);
    
    //    iss->numfact (cr->a,cr->ci,cr->adr,cr->n);
    //    printf("JSEM TU 3");
    //    fflush(stdout);
    
    //    iss->backsubst (cr->a,cr->ci,cr->adr,cr->n,lhs,rhs);
    //    printf("JSEM TU 4");
    //    fflush(stdout);
    
    print_err("PARDISO solver is not supported at this time",__FILE__,__LINE__,__func__);
    
    break;
  }
    
  case jacobi:{
    switch (ts){
    case dense_matrix:{
      dm->jacobi (lhs,rhs,ni,res,ani,ares,out);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }    
    break;
  }
    
    
  case gauss_seidel:{
    switch (ts){
    case dense_matrix:{
      dm->gauss_seidel (lhs,rhs,ni,res,ani,ares,out);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case cholmod_solver:{
    switch (ts){
      case symm_comp_columns:
        scc->solve(lhs, rhs);
        break;
      case compressed_columns:
        cc->solve(lhs, rhs);
        break;
      default:{
        print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
      }
    }
    break;
  }
  default:{
    print_err("unknown type of linear solveris required",__FILE__,__LINE__,__func__);
  }
  }

  //et = time (NULL);
  //fprintf (stdout,"\n time of solution of linear equations    %ld",et-bt);
  fflush(stdout);
}

/**
   function decomposes (factorizes) %matrix to the LU or LDL or LL form
   
   JK
*/
void gmatrix::decompose_matrix ()
{
  switch (tlinsol){
    
  case ldl:{
    switch (ts){
    case skyline_matrix:{
      if (sky->decomp()==0){
        sky->ldl_sky (NULL,NULL,zero,2);
      }
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case lu:{
    switch (ts){
    case dense_matrix:{
      if (dm->decomp()==0){
        dm->lu (NULL,NULL,zero,2);
      }
      break;
    }
    case double_skyline:{
      if (dsky->decomp()==0){
        dsky->lu_dsky (NULL,NULL,zero,2);
      }
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case lapack_sol:{
    //int Error;
    //Error = spFactor(Matrix);
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case spdirldl:
  case spdirlu:
  case spdirll:{
    if (sdirect->IsFactorized()==0){
      sdirect->ReFactorize();
    }
    break;
  }
    
    
  default:{
    print_err("unknown type of linear solver is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes back substitution
   
   @param lhs - array containing %vector of solution
   @param rhs - array containing %vector of the right hand side
   
   JK
*/
void gmatrix::back_substitution (double *lhs,double *rhs)
{
  switch (tlinsol){
    
  case ldl:{
    switch (ts){
    case skyline_matrix:{
      if (sky->decomp()==0){
        sky->ldl_sky (NULL,NULL,zero,2);
      }
      sky->ldl_sky (lhs,rhs,zero,3);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case lu:{
    switch (ts){
    case dense_matrix:{
      if (dm->decomp()==0){
        dm->lu (lhs,rhs,zero,2);
      }
      dm->lu (lhs,rhs,zero,3);
      break;
    }
    case double_skyline:{
      if (dsky->decomp()==0){
        dsky->lu_dsky (NULL,NULL,zero,2);
      }
      dsky->lu_dsky (lhs,rhs,zero,3);
      break;
    }
    default:{
      print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
    }
    }
    break;
  }
    
  case lapack_sol:{
    //spSolve( Matrix, rhs, lhs);
    print_err("LAPACK solver is not supported at this time",__FILE__,__LINE__,__func__);
    abort ();
    break;
  }
    
  case spdirldl:
  case spdirlu:
  case spdirll:{
    if (sdirect->IsFactorized()==0){
      sdirect->ReFactorize();
    }
    sdirect->Solve(lhs,rhs);
    break;
  }
  case pardisolver:{
    //  pozor na cislovani od 1
    //printf("chybi u PARDISO");
    //fflush(stdout);
    print_err("PARDISO solver is not supported at this time",__FILE__,__LINE__,__func__);
    
    break;
  }    
  default:{
    print_err("unknown type of linear solver is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes incomplete factorization of the %matrix
   
   @param incompltresh - treshold for incomplete factorization
   
   JK, 15.3.2007
*/
void gmatrix::incomplete_fact (double incompltresh)
{
  double *x,*y;
  x=NULL;
  y=NULL;
  
  switch (ts){
  case dense_matrix:{
    dm->ill (x,y,zero,incompltresh,2);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function computes incomplete factorization of the %matrix
   
   @param incompltresh - treshold for incomplete factorization
   
   JK, 15.3.2007
*/
void gmatrix::back_incomplete_fact (double *x,double *y,double incompltresh)
{
  
  switch (ts){
  case dense_matrix:{
    dm->ill (x,y,zero,incompltresh,3);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function condenses %matrix
   
   @param top - pointer to the general topology
   @param condmat - reduced %matrix (stored as dense %matrix)
   @param condvect - reduced %vector
   @param lhs - the left hand side (non-reduced)
   @param rhs - the right hand side (non-reduced)
   @param nrdof - number of condensed unknowns
   @param tc - computation type
   @param tc=1 - condensation (forward reduction) only
   @param tc=2 - back-substitution only
   
   JK
*/
void gmatrix::condense (gtopology *top,double *condmat,double *condvect,double *lhs,double *rhs,long nrdof,long tc,FILE *out)
{
  long i,j;
  
  if (pmat==0)
    prepmat2 (top,out);
  
  switch (ts){
  case dense_matrix:{
    dm->gempkon (condmat,condvect,lhs,rhs,nrdof,zero,tc);
    break;
  }
  case skyline_matrix:{
    sky->ldlkon_sky (condmat,condvect,lhs,rhs,nrdof,tc);
    break;
  }
  case double_skyline:{
    dsky->lukon_dsky (condmat,condvect,lhs,rhs,zero,nrdof,tc);
    break;
  }
  case spdirect_stor_cr:
  case spdirect_stor_scr:{
    

    if (tc==2){
      j=0;
      for (i=top->nidof;i<n;i++){
	lhs[i]=condvect[j];
	j++;
      }
    }

    sdirect->condense (condmat,lhs,rhs,tc);
    
    if (tc==1){
      j=0;
      for (i=top->nidof;i<n;i++){
	condvect[j]=rhs[i];
	j++;
      }
    }

    break;
  }
  default:{
    print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes kernel of %matrix
   
   @param rbm - array containing the rigid body modes
   @param nse - the number of linearly dependent rows
   @param se - array of linearly dependent rows
   @param ense - estimated number of linearly dependent rows
   @param limit - the threshold for linear dependency determination
   @param tc - computation type
   
   JK
*/
void gmatrix::kernel (double *rbm,long &nse,long *se,long ense,double limit,long tc)
{
  switch (ts){
  case dense_matrix:{
    dm->ker (rbm,nse,se,ense,limit);
    break;
  }
  case skyline_matrix:{
    sky->ker (rbm,nse,se,ense,limit,tc);
    break;
  }
  default:{
    print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes 
   
   
   JK
*/
void gmatrix::ldl_feti (double *lhs,double *rhs,long nse,long *se,double zero)
{
  switch (ts){
  case skyline_matrix:{
    sky->ldl_feti_sky (lhs,rhs,nse,se,zero);
    break;
  }
  default:{
    print_err("wrong storage type of matrix is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function computes %matrix-%vector multiplication
   
   @param a,b - array containing %vector components
   
   JK,14.6.2002
*/
void gmatrix::gmxv (double *a,double *b)
{
  switch (ts){
  case dense_matrix:{
    dm -> mxv_dm(a,b);
    break;
  }
  case skyline_matrix:{
    sky -> mxv_sky(a,b);
    break;
  }
  case double_skyline:{
    dsky -> mxv_dsky(a,b);
    break;
  }
  case compressed_rows:{
    cr -> mxv_cr(a,b);
    break;
  }
  case symm_comp_rows:{
    scr -> mxv_scr(a,b);
    break;
  }
  case compressed_columns:{
    cc -> mxv_cc(a,b);
    break;
  }
  case symm_comp_columns:{
    scc -> mxv_scc(a,b);
    break;
  }
  case element_matrices:{
    em -> mxv_em(a,b);
    break;
  }
  case spdirect_stor_scr:{
    scr -> mxv_scr(a,b);
    break;
  }
  case spdirect_stor_cr:{
    cr -> mxv_cr(a,b);
    break;
  }
  case diag_mat:{
    diagm->mxv_diag (a,b);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }

}

/**
   function computes %matrix-%vector multiplication
   %matrix is decomposed
   
   @param a,b - array containing %vector components
   
   JK, 14.6.2002
*/
void gmatrix::decompgmxv (double *a,double *b)
{
  switch (ts){
  case skyline_matrix:{
    sky -> ldlmxv_sky (a,b);
    break;
  }
  case spdirect_stor_cr:{
    sdirect -> MulMatrixByVector (a,b);
    break;
  }
  case spdirect_stor_scr:{
    sdirect -> MulMatrixByVector (a,b);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function adds to general %matrix another general %matrix
   premultiplied by scalar a
*/
void gmatrix::addgm (double a,gmatrix &gm)
{
  switch (ts){
  case dense_matrix:{
    dm -> addmat_dm (a,*gm.dm);
    break;
  }
  case skyline_matrix:{
    sky -> addmat_sky (a,*gm.sky);
    break;
  }
  case double_skyline:{
    dsky -> addmat_dsky (a,*gm.dsky);
    break;
  }
  case compressed_rows:{
    cr -> addmat_cr (a,*gm.cr);
    break;
  }
  case symm_comp_rows:{
    scr -> addmat_scr (a,*gm.scr);
    break;
  }
  case compressed_columns:{
    cc->addmat_cc(a,*gm.cc);
    break;
  }
  case symm_comp_columns:{
    scc ->addmat_scc (a,*gm.scc);
    break;
  }
  case element_matrices:{
    em -> addmat_em (a,*gm.em);
    break;
  }
  case spdirect_stor_scr:{
    scr -> addmat_scr (a,*gm.scr);
    //sdirect -> AddNumbers (a,gm.sdirect->GetSparseMatrix());
    break;
  }
  case spdirect_stor_cr:{
    cr -> addmat_cr (a,*gm.cr);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function multiplies general %matrix by a variable a
   
   @param a - multiplication constant
   
   JK, 17.3.2007
*/
void gmatrix::scalgm (double a)
{
  switch (ts){
  case dense_matrix:{
    dm -> scalmat_dm (a);
    break;
  }
  case skyline_matrix:{
    sky -> scalmat_sky (a);
    break;
  }
  case double_skyline:{
    dsky -> scalmat_dsky (a);
    break;
  }
  case compressed_rows:{
    cr -> scalmat_cr (a);
    break;
  }
  case symm_comp_rows:{
    scr -> scalmat_scr (a);
    break;
  }
  case compressed_columns:{
    cc -> scalmat_cc (a);
    break;
  }
  case symm_comp_columns:{
    scc -> scalmat_scc (a);
    break;
  }
  case element_matrices:{
    em   -> scalmat_em (a);
    break;
  }
  case spdirect_stor_scr:{
    scr -> scalmat_scr (a);
    //sdirect -> ScaleMatrix (a);
    break;
  }
  case spdirect_stor_cr:{
    cr -> scalmat_cr (a);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function copies a general %matrix  to another general %matrix gm
   the argument is filled
   
   @param gm - general %matrix where the copy will be stored
   
   JK, 17.3.2007
*/
void gmatrix::copygm (gmatrix &gm)
{
  //  type of matrix storage
  gm.ts=ts;
  
  switch (ts){
  case dense_matrix:{
    if (gm.dm==NULL){
      gm.dm = new densemat ();
    }
    dm -> copy_dm (*gm.dm);
    break;
  }
  case skyline_matrix:{
    if (gm.sky==NULL){
      gm.sky = new skyline ();
    }
    sky -> copy_sky (*gm.sky);
    break;
  }
  case double_skyline:{
    if (gm.dsky==NULL){
      gm.dsky = new dskyline ();
    }
    dsky -> copy_dsky (*gm.dsky);
    break;
  }
  case compressed_rows:{
    if (gm.cr==NULL){
      gm.cr = new comprow ();
    }
    cr -> copy_cr (*gm.cr);
    break;
  }
  case symm_comp_rows:{
    if (gm.scr==NULL){
      gm.scr = new symcomprow ();
    }
    scr -> copy_scr (*gm.scr);
    break;
  }
    /* asi toto: ale pada to
       case spdirect_stor_scr:{
       if (gm.scr==NULL){
       gm.scr = new symcomprow ();
       }
       scr -> copy_scr (*gm.scr);
       break;
       }
       case spdirect_stor_cr:{
       if (gm.cr==NULL){
       gm.cr = new comprow ();
       }
       cr -> copy_cr (*gm.cr);
       break;
       }
    */
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

void gmatrix::a12block (gtopology *top,double *block,long nrdof,FILE *out)
{
  if (pmat==0)
    prepmat2 (top,out);

  switch (ts){
  case skyline_matrix:{
    sky->ldl_a12block (block,nrdof);
    break;
  }
  case spdirect_stor_scr:{
    sdirect->GetA12block(block);
    break;
  }
  default:{
    print_err("wrong type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function prints general %matrix
   
   @param out - output stream
   
   JK, 17.3.2007
*/
void gmatrix::printmat (FILE *out)
{
  switch (ts){
  case dense_matrix:{
    dm->printmat (out);
    break;
  }
  case skyline_matrix:{
    sky->printmat (out);
    break;
  }
  case double_skyline:{
    dsky->printmat (out);
    break;
  }
  case compressed_rows:{
    cr->printmat (out);
    break;
  }
  case symm_comp_rows:{
    scr->printmat (out);
    break;
  }
  case compressed_columns:{
    cc->printmat (out);
    break;
  }
  case symm_comp_columns:{
    scc->printmat (out);
    break;
  }
  case element_matrices:{
    em->printmat (out);
    break;
  }
  case spdirect_stor_scr:{
    scr->printmat (out);
    break;
  }
  case spdirect_stor_cr:{
    cr->printmat (out);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function prints diagonal entries of general %matrix
   
   @param out - output stream
   
   JK, 17.3.2007
*/
void gmatrix::printdiag (FILE *out)
{
  switch (ts){
  case dense_matrix:{
    dm->printdiag (out);
    break;
  }
  case skyline_matrix:{
    sky->printdiag (out);
    break;
  }
  case double_skyline:{
    dsky->printdiag (out);
    break;
  }
  case compressed_rows:{
    cr->printdiag (out);
    break;
  }
  case symm_comp_rows:{
    scr->printdiag (out);
    break;
  }
  case spdirect_stor_cr:{
    cr->printdiag (out);
    break;
  }
  case spdirect_stor_scr:{
    scr->printdiag (out);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
}

void gmatrix::changedecomp ()
{
  switch (ts){
  case dense_matrix:{
    dm->changedecomp ();
    break;
  }
  case skyline_matrix:{
    sky->changedecomp ();
    break;
  }
  case double_skyline:{
    dsky->changedecomp ();
    break;
  }
  case compressed_rows:{
    //cr->changedecomp ();
    break;
  }
  case symm_comp_rows:{
    //scr->changedecomp ();
    break;
  }
  case compressed_columns:{
    cc->changedecomp ();
    break;
  }
  case symm_comp_columns:{
    scc->changedecomp ();
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }

}


/*
void gmatrix::setfact ()
{
  switch (ts){
  case dense_matrix:{
    dm->setfact ();
    break;
  }
  case skyline_matrix:{
    sky->setfact ();
    break;
  }
  case double_skyline:{
    dsky->setfact ();
    break;
  }
  case spdirldl:
  case spdirlu:
  case spdirll:{
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }

}

void gmatrix::setnotfact ()
{
  switch (ts){
  case dense_matrix:{
    dm->setnotfact ();
    break;
  }
  case skyline_matrix:{
    sky->setnotfact ();
    break;
  }
  case double_skyline:{
    dsky->setnotfact ();
    break;
  }
  case spdirldl:
  case spdirlu:
  case spdirll:{
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }

}

*/

long gmatrix::decomp ()
{
  long d=0;

  switch (ts){
  case dense_matrix:{
    d=dm->decomp ();
    break;
  }
  case skyline_matrix:{
    d=sky->decomp ();
    break;
  }
  case double_skyline:{
    d=dsky->decomp ();
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  return d;
}

/**
   function returns required matrix entry
   
   @param ri,ci - row and column indices
   
   JK, 24.7.2005
*/
double gmatrix::give_entry (long ri,long ci)
{
  double e;
  
  switch (ts){
  case dense_matrix:{
    e=dm->give_entry (ri,ci);
    break;
  }
  case skyline_matrix:{
    e=sky->give_entry (ri,ci);
    break;
  }
  case double_skyline:{
    //dsky->changedecomp ();
    e = 0.0;
    print_err(" give_entry function is not supported on the double_skyline storage.", __FILE__, __LINE__, __func__);
    break;
  }
  case compressed_rows:{
    e=cr->give_entry (ri,ci);
    break;
  }
  case symm_comp_rows:{
    e=scr->give_entry (ri,ci);
    break;
  }
  case diag_mat:{
    e = diagm->give_entry (ri);
    break;
  }
  default:{
    e = 0.0;
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return e;
}

/**
   function adds matrix entry to required position
   
   @param e - matrix entry
   @param ri,ci - row and column indices
   
   JK, 25.7.2005
*/
void gmatrix::add_entry (double e,long ri,long ci)
{
  switch (ts){
  case dense_matrix:{
    dm->add_entry (e,ri,ci);
    break;
  }
  case skyline_matrix:{
    sky->add_entry (e,ri,ci);
    break;
  }
  case double_skyline:{
    //dsky->changedecomp ();
    break;
  }
  case compressed_rows:{
    cr->add_entry (e,ri,ci);
    break;
  }
  case symm_comp_rows:{
    scr->add_entry (e,ri,ci);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }

}

/**
   function returns number of stored %matrix entries
   
   JK, 16.8.2007
*/
long gmatrix::give_negm ()
{
  long negm = 0;
  
  switch (ts){
  case dense_matrix:{
    negm=dm->give_negm ();
    break;
  }
  case skyline_matrix:{
    negm=sky->give_negm ();
    break;
  }
  case double_skyline:{
    negm=dsky->give_negm ();
    break;
  }
  case compressed_rows:{
    negm=cr->give_negm ();
    break;
  }
  case symm_comp_rows:{
    negm=scr->give_negm ();
    break;
  }
  case spdirect_stor_cr:{
    negm=cr->give_negm ();
    break;
  }
  case spdirect_stor_scr:{
    negm=scr->give_negm ();
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return negm;
}

/**
   function scales %matrix by its diagonal elements
   
   JK, 23.5.2008
*/
void gmatrix::diag_scale (double *d)
{
  switch (ts){
  case dense_matrix:{
    dm->diag_scale (d);
    break;
  }
  case skyline_matrix:{
    sky->diag_scale (d);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function checks diagonal entries
   
   the function is used in some nonlinear nostationary problems
   where high jumps in coefficients occur
   some of element matrices are zero matrices and this
   function puts nonzero values on the diagonal
   
   JK, 14.7.2008
*/
void gmatrix::diag_check (double thr,double *rhs)
{
  switch (ts){
  case skyline_matrix:{
    sky->diag_check (thr);
    break;
  }
  case double_skyline:{
    dsky->diag_check (thr,rhs);
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
}


/**
   function performes the power method
   the method computes the largest eigenvalue of %matrix
   
   @param v - array containing the eigenvector
   @param ni - maximum number of iterations
   @param err - required error
   
   JK, 28.5.2008
*/
double gmatrix::power_method (double *v,long ni,double err)
{
  long i;
  double ev,evo,nom,denom;
  double *p;
  
  p = new double [n];
  
  for (i=0;i<n;i++){
    v[i]=1.0;
  }
  
  evo=0.0;
  
  for (i=0;i<ni;i++){
    gmxv (v,p);
    
    //  nominator
    nom=ss(v,p,n);
    //  denominator
    denom=ss(v,v,n);
    
    ev=nom/denom;
    
    fprintf (stdout,"\n power method, step %6ld,  quotient %e",i,ev);

    copyv(p,v,n);
    
    cmulv(1.0/sqrt(denom),v,n);
    
    if (fabs(ev-evo)/ev<err && i>0)
      break;
    evo=ev;
  }

  delete [] p;
  
  return ev;
}

/**
   function performes the inverse iteration method
   the method computes the smallest eigenvalue of %matrix
   
   @param v - array containing the eigenvector
   @param ni - maximum number of iterations
   @param err - required error
   
   JK, 28.5.2008
*/
double gmatrix::inverse_iteration (double *v,long ni,double err)
{
  long i;
  double ev,evo,nom,denom;
  double *p;
  
  p = new double [n];
  
  for (i=0;i<n;i++){
    v[i]=1.0;
  }
  
  evo=0.0;

  dm->lu (p,v,zero,2);
 
  for (i=0;i<ni;i++){
    dm->lu (p,v,zero,3);
    
    //  nominator
    nom=ss(v,p,n);
    //  denominator
    denom=ss(v,v,n);
    
    ev=nom/denom;
    
    fprintf (stdout,"\n inverse iteration, step %6ld,  quotient %e",i,ev);

    copyv(p,v,n);
    
    cmulv(1.0/sqrt(denom),v,n);
    
    if (fabs(ev-evo)/ev<err && i>0){
      break;
    }
    
    evo=ev;
  }

  delete [] p;
  
  return ev;

}


/**
   function estimates spectral radius

   the estimates are based on the Gershgorin circle theorem
   for details see G.H. Golub, C.F. Van Loan: Matrix computations.
   The Johns Hopkins University Press, 3rd edition, page 320, 1996
   
   JK, 27.8.2008
*/
double gmatrix::estim_spect_radius ()
{
  double sr=0.0;
  
  switch (ts){
  case dense_matrix:{
    sr=dm->estim_spect_radius ();
    break;
  }
  case compressed_rows:{
    sr=cr->estim_spect_radius ();
    break;
  }
  default:{
    print_err("unknown type of matrix storage is required",__FILE__,__LINE__,__func__);
  }
  }
  
  return sr;
}


/**
  The function turns off output messages from equation system solvers.

  @return The function does not return anything but it may change the status of solver flags for output.

  Created by Tomas Koudelka, 10.5.2016
*/
void gmatrix::mute()
{
  if (ts == spdirect_stor_cr || ts == spdirect_stor_scr){
    // redirect solver output messeges to the NULL device
    sdirect->SetMTN();
    sdirect->SetMTMtx();
  }
}
