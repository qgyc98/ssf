#include "seqfeti.h"
#include "gmatrix.h"
#include <limits.h>

seqfeti::seqfeti ()
{
  //  type of FETI implementation
  //  this means no type is selected
  fetiimpl=no_impl;
  //  indicator of the presence of discontinuities
  //  discont=0 - the classical FETI method
  //  discont=1 - the FETI method with interface discontinuities
  discont=1;
  //  type of storage of subdomain matrices
  smst=skyline_matrix;
  //smst=dense_matrix;
  
  //  indicator of assembled matrices
  matassem=0;
  //  indicator of assembled right hand side vectors
  vecassem=0;
  
  
  //  the number of subdomains
  ns=0;
  //  estimated number of rigid body modes (estimated dimension of the kernel)
  ense=0;
  

  //  the maximum number of degrees of freedom on one subdomain
  ndofmax=0;
  //  the number of DOFs (unknowns) in coarse problem
  ndofcp=0;
  //  threshold for kernel detection
  thresh=0.0;
  //  the maximum number of iterations in conjugate gradient method
  nicg=0;
  //  the number of performed iterations in conjugate gradient method
  anicg=0;
  //  required error
  errcg=0.0;
  //  attained error
  aerrcg=0.0;
  //  computer zero
  zero=1.0e-100;
  
  //  type of preconditioner
  prec=noprecond;
  //  size of matrix G
  gsize=0;
  //  the maximum number of degrees of freedom on one subdomain
  ndofmax=0;
  //  the number of DOFs (unknowns) in coarse problem
  ndofcp=0;
  
  
  //  the numbers of DOFs on subdomains
  ndofmas=NULL;
  
  //  node-subdomain correspondence
  nsid=NULL;
  
  //  list of DOFs on subdomains
  //  it is used in connection with preconditioning
  cndom=NULL;
  
  //  array of numbers of unknowns (DOFs) contributing to the coarse problem
  ncdofd=NULL;
  
  //  array containing code numbers contributing to the coarse problem
  //  extracted values from subdomains to the coarse problem
  edofs=NULL;
  
  //  array containing numbers of RBM on subdomains
  nrbmdom=NULL;

  //  rigid body modes / kernel
  rbmdom=NULL;
  
  //  array containing addresses of first RBM in coarse matrix
  rbmadr=NULL;

  //  list of linearly dependent equations
  se=NULL;

  //  the %matrix G
  g=NULL;
  
  //  the (G^T G)^{-1}
  invgg = NULL;
  
  //  the arryas with the right hand sides
  ff = NULL;
  
  //  e vector
  e = NULL;
  
  //  array for the nodal unknowns
  d = NULL;
  
  //  array of Lagrange multipliers
  lambda = NULL;
  
  //  array for the %vector b
  //  it is a constant %vector containing prescribed discontinuities
  b = NULL;
  
  //  array for the compliances in the %matrix H
  h = NULL;
  
  //  the number of contributions in the arrays booldatar, booldatac and booldata
  ncbool=NULL;

  //  array containing row indices for construction of the Boolean matrices
  booldatar=NULL;

  //  array containing column indices for construction of the Boolean matrices
  booldatac=NULL;

  //  array containing %matrix entries of the Boolean matrices
  booldata=NULL;

  //  array of coarse code numbers
  ccn=NULL;


  //  subdomain matrices stored in the %skyline storage
  smsky=NULL;
  //  subdomain matrices stored in the %dense format
  smdm=NULL;

}

seqfeti::~seqfeti ()
{
  long i;
  
  // ************************
  //  twodimensional arrays
  // ************************
  for (i=0;i<ns;i++){
    //  array of coarse code numbers
    if (ccn!=NULL)
      delete [] ccn[i];
    //  array containing %matrix entries of the Boolean matrices
    if (booldata!=NULL)
      delete [] booldata[i];
    //  array containing column indices for construction of the Boolean matrices
    if (booldatac!=NULL)
    delete [] booldatac[i];
    //  array containing row indices for construction of the Boolean matrices
    if (booldatar!=NULL)
    delete [] booldatar[i];
    //  list of linearly dependent equations
    if (se!=NULL)
    delete [] se[i];
    //  rigid body modes / kernel
    if (rbmdom!=NULL)
    delete [] rbmdom[i];
    //  array containing code numbers contributing to the coarse problem
    if (edofs!=NULL)
    delete [] edofs[i];
    //  list of DOFs on subdomains
    if (cndom!=NULL)
    delete [] cndom[i];
    //  the arrays with the right hand sides
    if (ff!=NULL)
    delete [] ff[i];
    //  array for the nodal unknowns
    if (d!=NULL)
    delete [] d[i];
  }

  //  array of coarse code numbers
  delete [] ccn;
  //  array containing %matrix entries of the Boolean matrices
  delete [] booldata;
  //  array containing column indices for construction of the Boolean matrices
  delete [] booldatac;
  //  array containing row indices for construction of the Boolean matrices
  delete [] booldatar;
  //  list of linearly dependent equations
  delete [] se;
  //  rigid body modes / kernel
  delete [] rbmdom;
  //  array containing code numbers contributing to the coarse problem
  delete [] edofs;
  //  list of DOFs on subdomains
  delete [] cndom;
  //  the arrays with the right hand sides
  delete [] ff;
  //  array for the nodal unknowns
  delete [] d;
 
  
  // ************************
  //  onedimensional arrays
  // ************************
  //  number of contributions in the arrays booldatar, booldatac and booldata
  delete [] ncbool;
  //  array containing addresses of first RBM in coarse matrix
  delete [] rbmadr;
  //  array containing numbers of RBM on subdomains
  delete [] nrbmdom;
  //  array of numbers of unknowns (DOFs) contributing to the coarse problem
  delete [] ncdofd;
  //  node-subdomain correspondence
  delete [] nsid;
  //  numbers of DOFs on subdomains
  delete [] ndofmas;
  
  //  the %matrix G
  delete [] g;
  //  the (G^T G)^{-1}
  delete [] invgg;
  //  e vector
  delete [] e;
  //  array of Lagrange multipliers
  delete [] lambda;
  //  array for the %vector b
  delete [] b;
  //  array for the compliances in the %matrix H
  delete [] h;
  
  //  subdomain matrices stored in the %skyline storage
  delete [] smsky;
  //  subdomain matrices stored in the %dense format
  delete [] smdm;

}

/**
   function reads basic data
   
   @param top - pointer to general topology
   @param in - input file
   
   JK, 20.11.2007
*/
void seqfeti::read (gtopology *top,XFILE *in)
{
  //  type of implementation
  //  fetiimpl=1 - boolean matrices
  //  fetiimpl=2 - nonredundant definition of multipliers
  //  fetiimpl=3 - redundant definition of multipliers
  xfscanf (in,"%m",&fetiimplem_kwdset,&fetiimpl);
  
  if (fetiimpl!=boolean_matrices && fetiimpl!=nonredundant && fetiimpl!=redundant){
    print_err("unknown type of sequential FETI implementation",__FILE__,__LINE__,__func__);
  }
  
  //  the number of subdomains
  //  the number of estimated rigid body modes on one subdomain
  //  threshold for zero detection on the diagonal (threshold for detection of dependent row)
  //  the maximum number of iterations in the conjugate gradient method
  //  required norm of the residual
  xfscanf (in,"%ld %ld %le %ld %le",&ns,&ense,&thresh,&nicg,&errcg);
  //  type of preconditioner
  xfscanf (in,"%d",(int*)&prec);
  
  
  // **********************************
  //  special requirements of solvers
  // **********************************
  if (fetiimpl==nonredundant || fetiimpl==redundant){
    //  sequential topology will be read in function mechtop or transtop :: read
    top->rst=1;
  }
  if (fetiimpl==boolean_matrices){
    //  sequential topology will be read in function mechtop or transtop :: read
    //  Boolean matrices will be read too
    top->rst=2;
  }
  if (prec==feti_dirichlet){
    //  code numbers will be generated with respect to the Schur complement method
    //  internal DOFs are ordered first, interface DOFs are ordered afterwards
    top->cngen=2;
  }
  
  //  the code numbers are not generated
  top->cngen=0;
  
}

/**
   function prints basic data
   
   @param out - output file
   
   JK, 24.11.2008
*/
void seqfeti::print (FILE *out)
{
  //  type if implementation
  fprintf (out,"%d",fetiimpl);
  
  //  the number of subdomains
  //  the number of estimated rigid body modes on one subdomain
  //  threshold for zero detection on the diagonal (threshold for detection of dependent row)
  //  the maximum number of iterations in the conjugate gradient method
  //  required norm of residual
  fprintf (out,"%ld %ld %le %ld %le\n",ns,ense,thresh,nicg,errcg);
  //  type of preconditioner
  fprintf (out,"%d\n",prec);
}

/**
   function reads data for construction of Boolean matrices
   
   @param in - input file

   JK, 19.2.2008
*/
void seqfeti::read_booldata (XFILE *in)
{
  long i,j;
  
  //  allocation of array
  ncbool = new long [ns];
  booldatar = new long* [ns];
  booldatac = new long* [ns];
  booldata = new double* [ns];
  
  for (i=0;i<ns;i++){
    //  loop over the number of subdomains
    
    //  the number of contributions on the i-th subdomain
    xfscanf (in,"%ld",ncbool+i);
    
    //  allocation of array
    booldatar[i] = new long [ncbool[i]];
    booldatac[i] = new long [ncbool[i]];
    booldata[i] = new double [ncbool[i]];
    
    for (j=0;j<ncbool[i];j++){
      xfscanf (in,"%ld %ld %lf",&booldatar[i][j],&booldatac[i][j],&booldata[i][j]);
    }
  }

}

/**
   function assembles the following arrays
   ncdofd - array of numbers of unknowns (DOFs) contributing to the coarse problem
   wscalmat - 
   edofs - array containing code numbers contributing to the coarse problem
   ccn - array of coarse code numbers
   
   @param selnodfeti - pointer to the object containing selected nodes
   @param top - pointer to object with topology
   @param out - output file
   
   JK, 17.9.2007
*/
void seqfeti::initiate (seqselnodes *selnodfeti,gtopology *top,FILE *out)
{
  long i,j,k,ndofn,dof,nid;
  long *red;
  red = new long [ns+1];
  
  //  the number of DOFs (unknowns) in coarse problem = the total number of DOFs in selected nodes
  ndofcp = selnodfeti->tndofsn;
  
  //  array of numbers of unknowns (DOFs) contributing to the coarse problem
  ncdofd = new long [ns];
  for (i=0;i<ns;i++){
    ncdofd[i] = selnodfeti->snndofmas[i];
  }
  
  //  array containing local code numbers contributing to the coarse problem
  edofs = new long* [ns];
  for (i=0;i<ns;i++){
    edofs[i] = new long [ncdofd[i]];
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<ncdofd[i];j++){
      edofs[i][j] = selnodfeti->lndofmas[i][j];
    }
  }
  //  at this moment, edofs contains DOFs in global glued ordering
  //  nodes are taken into account in increasing order and code numbers
  //  are generated, the ordering has to be recalculated from the global
  //  glued ordering to local ordering, this is done with the help of
  //  the arry red
  
  red[0]=0;
  nid=0;
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    red[i+1]=0;
    //  loop over the number of nodes in the i-th subdomain
    for (j=0;j<top->stop->nnsd[i];j++){
      ndofn = top->give_ndofn (nid);
      for (k=0;k<ndofn;k++){
	dof = top->give_dof (nid,k);
	if (dof>red[i+1]){
	  red[i+1]=dof;
	}
      }
      nid++;
    }
  }
  
  //  array red contains the minimum numbers used in subdomains
  //  these numbers have to be subtracted from all code numbers
  
  for (i=0;i<ns;i++){
    for (j=0;j<ncdofd[i];j++){
      if (edofs[i][j]>0)
	edofs[i][j] -= red[i];
      else
	edofs[i][j] += red[i];
    }
  }
  
  
  fprintf (out,"\n\n\n check of the array edofs \n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n domain %ld",i);
    for (j=0;j<ncdofd[i];j++){
      fprintf (out,"\n edofs  %6ld  %6ld",j,edofs[i][j]);
    }
  }
  
  
  


  // scaling matrix
  wscalmat = new double [ndofcp]; 
  for (i = 0; i < ndofcp; i++){
    //wscalmat[i] = double(selnodfeti->dofmultip[i]);
  }
  




  ccn = new long* [ns];
  for (i=0;i<ns;i++){
    ccn[i] = new long [ncdofd[i]];
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<ncdofd[i];j++){
      ccn[i][j] = selnodfeti->cndofmas[i][j];
    }
  }
  
  delete [] red;

  /*
  // lokalni wscalmat
  lwscalmat = new double* [ns]; 
  nlwscalmat = new long [ns];
  for (i = 0; i < ns; i++){
    nlwscalmat[i] = selnodfeti->ndofmas[i];
    lwscalmat[i] = new double [ nlwscalmat[i]]; 
    for (j = 0; j <  nlwscalmat[i]; j++){
      lwscalmat[i][j] = double(selnodfeti->ldofmultip[i][j]);
    }
  }
  */
}


/**
   function determines ndofmax
   the variable is used for simple allocation in other functions
   
   JK, 16.9.2007, checked 20.9.2009
*/
void seqfeti::det_ndofmax ()
{
  long i,j;
  
  if (fetiimpl==nonredundant || fetiimpl==redundant){
    ndofmax=0;
    for (i=0;i<ns;i++){
      if (ndofmax<ndofmas[i])
	ndofmax=ndofmas[i];
    }
  }
  
  if (fetiimpl==boolean_matrices){
    ndofmax=0;
    ndofcp=0;
    for (i=0;i<ns;i++){
      for (j=0;j<ncbool[i];j++){
	if (ndofmax<booldatac[i][j])
	  ndofmax=booldatac[i][j];
	if (ndofcp<booldatar[i][j])
	  ndofcp=booldatar[i][j];
      }
    }
  }
}

/**
   function assembles list of unknowns belonging to the subdomains
   
   function assembles:
   sid - node-subdomain correspondence
   ndofmas - numbers of DOFs on subdomains
   cndom - list of DOFs on subdomains
   
   @param top - pointer to topology
   @param out - output file
   
   JK, 21.8.2007
*/
void seqfeti::assemble_subdom_unknowns (gtopology *top,FILE *out)
{
  long i,j,k,l,m,ndofn;

  //  node-subdomain correspondence
  nsid = new long [top->nn];
  l=0;
  for (i=0;i<ns;i++){
    for (j=0;j<top->stop->nnsd[i];j++){
      nsid[l]=i;
      l++;
    }
  }
  
  
  //  number of DOFs on subdomains
  ndofmas = new long [ns];
  for (i=0;i<ns;i++){
    ndofmas[i]=0;
  }
  
  for (i=0;i<top->nn;i++){
    //  number of DOFs in node
    ndofn=top->give_ndofn (i);
    //  subdomain id
    l=nsid[i];
    for (j=0;j<ndofn;j++){
      //  DOF id
      k=top->give_dof (i,j);
      if (k>0)
	ndofmas[l]++;
    }
  }
  
  //  list of DOFs on subdomains
  cndom = new long* [ns];
  for (i=0;i<ns;i++){
    cndom[i] = new long [ndofmas[i]];
    ndofmas[i]=0;
  }
  
  for (i=0;i<top->nn;i++){
    //  number of DOFs in node
    ndofn=top->give_ndofn (i);
    //  subdomain id
    l=nsid[i];
    for (j=0;j<ndofn;j++){
      //  DOF id
      k=top->give_dof (i,j);
      if (k>0){
	cndom[l][ndofmas[l]]=k;
	ndofmas[l]++;
      }
    }
  }
  
  

  
  
  if (prec == feti_dirichlet){
    for (i=0;i<ns;i++){
      for (j=0;j<ndofmas[i];j++){
	l=cndom[i][j];
	m=j;
	for (k=j;k<ndofmas[i];k++){
	  if (l>cndom[i][k]){
	    l=cndom[i][k];
	    m=k;
	  }
	}
	cndom[i][m]=cndom[i][j];
	cndom[i][j]=l;
      }
    }
  }

  
  fprintf (out,"\n\n kontrola pole ndofmas");
  for (i=0;i<ns;i++){
    fprintf (out,"\n ndofmas %6ld   %6ld",i,ndofmas[i]);
  }
  fprintf (out,"\n\n kontrola pole cndom");
  for (i=0;i<ns;i++){
    fprintf (out,"\n domain %ld",i);
    for (j=0;j<ndofmas[i];j++){
      fprintf (out,"\n dom %6ld  dof %6ld   %6ld",i,j,cndom[i][j]);
    }
  }  
  
}


/**
   function assembles subdomain matrices from the %matrix of the whole system
   
   @param gm - pointer to the %matrix of the system
   @param out - output file
   
   JK, 21.8.2007
*/
void seqfeti::subdomain_matrices (gmatrix *gm,FILE *out)
{
  long i,j,k;
  
  if (smst==dense_matrix){
    smdm = new densemat [ns];
    for (i=0;i<ns;i++){
      //smdm[i].assemble_dense_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,ndofmas[i],cndom[i]);
    }
  }
  if (smst==skyline_matrix){
    smsky = new skyline [ns];
    for (i=0;i<ns;i++){
      smsky[i].assemble_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,ndofmas[i],cndom[i]);
    }
  }
  
    
  if (fetiimpl==nonredundant || fetiimpl==redundant){
    long *aux;
    ndofprec = new long [ns];
    cnprec = new long* [ns];
    cpreccn = new long* [ns];
    
    for (i=0;i<ns;i++){
      aux = new long [ndofmas[i]];
      for (j=0;j<ndofmas[i];j++){
	aux[j]=0;
      }
      
      for (j=0;j<ncdofd[i];j++){
	k=edofs[i][j];
	if (k<0)
	  k=0-k-1;
	else
	  k--;
	if (k>-1)
	  aux[k]++;
      }
      
      ndofprec[i]=0;
      for (j=0;j<ndofmas[i];j++){
	if (aux[j]>0)
	  ndofprec[i]++;
      }
      
      fprintf (out,"\n\n aux  domena %ld\n",i);
      for (j=0;j<ndofmas[i];j++){
	fprintf (out," %ld",aux[j]);
      }
      
      
      cnprec[i] = new long [ndofprec[i]];
      cpreccn[i] = new long [ndofprec[i]];
      k=0;
      for (j=0;j<ndofmas[i];j++){
	if (aux[j]>0){
	  cnprec[i][k]=cndom[i][j];
	  cpreccn[i][k]=j;
	  k++;
	}
      }
      
      delete [] aux;
    }
  
    fprintf (out,"\n\n\n\n kontrola ve funkci subdomain_matrices (file %s, line %d)\n\n",__FILE__,__LINE__);
    for (i=0;i<ns;i++){
      fprintf (out,"\n domena %5ld     ndofprec %ld\n",i,ndofprec[i]);
      for (j=0;j<ndofprec[i];j++){
	fprintf (out,"  %ld",cnprec[i][j]);
      }
      fprintf (out,"\n");
      for (j=0;j<ndofprec[i];j++){
	fprintf (out,"  %ld",cpreccn[i][j]);
      }
    }
    fprintf (out,"\n\n\n\n\n");
  }
  

  if (fetiimpl==nonredundant || fetiimpl==redundant){
    if (prec==feti_lumped){
      smscr = new symcomprow [ns];
      for (i=0;i<ns;i++){
	//gm->scr->select_submatrix (&smscr[i],ndofmas[i],cndom[i]);
	gm->scr->select_submatrix (&smscr[i],ndofprec[i],cnprec[i]);
      }
    }
    
    if (prec==feti_dirichlet){
      psmdm = new densemat [ns];
      
      double *x,*y,*condvect;
      for (i=0;i<ns;i++){
	x = new double [ndofmas[i]];
	y = new double [ndofmas[i]];
	psmdm[i].a = new double [ndofprec[i]*ndofprec[i]];
	psmdm[i].negm=ndofprec[i]*ndofprec[i];
	psmdm[i].n=ndofprec[i];
	condvect=new double [ndofprec[i]];
	
	smsky[i].ldlkon_sky (psmdm[i].a,condvect,x,y,ndofprec[i],1);
	smsky[i].~skyline ();
	
	delete [] x;
	delete [] y;
	delete [] condvect;
      }
      
      smsky = new skyline [ns];
      for (i=0;i<ns;i++){
	smsky[i].assemble_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,ndofmas[i],cndom[i]);
      }
      
    }
  }
  
  if (fetiimpl==boolean_matrices){
    if (prec==feti_lumped){
      smskyprec = new skyline [ns];
      for (i=0;i<ns;i++){
	smskyprec[i].assemble_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,ndofmas[i],cndom[i]);
      }
    }
    if (prec==feti_dirichlet){
      
      ndofprec = new long[ns];
      for (i=0;i<ns;i++){
	ndofprec[i] = 1;
	for(j = 1; j < ncbool[i]; j++){
	  if(booldatac[i][j] != booldatac[i][j-1]){
	    ndofprec[i]++;
	  }
	}
      }
      
      psmdm = new densemat [ns];
      
      double *x,*y,*condvect;
      for (i=0;i<ns;i++){
	x = new double [ndofmas[i]];
	y = new double [ndofmas[i]];
	psmdm[i].a = new double [ndofprec[i]*ndofprec[i]];
	psmdm[i].negm=ndofprec[i]*ndofprec[i];
	psmdm[i].n=ndofprec[i];
	condvect=new double [ndofprec[i]];
	
	smsky[i].ldlkon_sky (psmdm[i].a,condvect,x,y,ndofprec[i],1);
	smsky[i].~skyline ();
	
	delete [] x;
	delete [] y;
	delete [] condvect;
      }
      
      smsky = new skyline [ns];
      for (i=0;i<ns;i++){
	smsky[i].assemble_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,ndofmas[i],cndom[i]);
      }
    }
  }  
  
}

/**
   function computes kernel (rigid body modes) of the %matrix
   
   JK, 28.8.2007
*/
void seqfeti::kernel (FILE */*out*/)
{
  long i,j;
  
  //  array for numbers of rigid body modes
  if (nrbmdom!=NULL){
    delete [] nrbmdom;
  }
  nrbmdom = new long [ns];
  
  //  array for linearly dependent equation numbers
  if (se!=NULL){
    for (i=0;i<ns;i++){
      delete [] se[i];
    }
    delete [] se;
  }
  se = new long* [ns];
  for (i=0;i<ns;i++){
    se[i] = new long [ense];
    for (j=0;j<ense;j++){
      se[i][j]=-1;
    }
  }
  
  //  array for rigid body modes
  if (rbmdom!=NULL){
    for (i=0;i<ns;i++){
      delete [] rbmdom[i];
    }
    delete [] rbmdom;
  }
  rbmdom = new double* [ns];
  for (i=0;i<ns;i++){
    rbmdom[i] = new double [ense*ndofmas[i]];
    for (j=0;j<ense*ndofmas[i];j++){
      rbmdom[i][j]=0.0;
    }
  }
  
  
  if (smst==dense_matrix){
    //  loop over subdomains
    for (i=0;i<ns;i++){
      smdm[i].ker (rbmdom[i],nrbmdom[i],se[i],ense,thresh);
    }
  }
  if (smst==skyline_matrix){
    //  loop over subdomains
    for (i=0;i<ns;i++){
      smsky[i].ker (rbmdom[i],nrbmdom[i],se[i],ense,thresh,3);
      
      fprintf (stdout,"\n the number of rigid body motions %ld",nrbmdom[i]);
      
      //smsky[i].printmat (out);
    }
  }

}
 













/**
   function assembles the Boolean %matrix of one subdomain
   
   @param sdid - subdomain id
   @param a - array storing the Boolean %matrix

   JK, 17.2.2009
*/
void seqfeti::boolean_matrix (long sdid,double *a)
{
  long i,ri,ci;
  
  for (i=0;i<ncbool[sdid];i++){
    ri=booldatar[sdid][i]-1;
    ci=booldatac[sdid][i]-1;
    a[ri*ndofmas[sdid]+ci]=booldata[sdid][i];
  }
}



/**
   function extracts contributions from local %vector to coarse %vector
   
   @param nd - domain id (number of required subdomain)
   @param lv - local %vector
   @param cv - coarse %vector
   
   JK, 14.9.2007
*/
void seqfeti::local_coarse (long nd,double *lv,double *cv)
{
  if (fetiimpl==nonredundant || fetiimpl==redundant){
    
    long i,j,k;
    
    for (i=0;i<ncdofd[nd];i++){
      j=edofs[nd][i];
      k=ccn[nd][i]-1;
      if (j<0){
	j=0-j-1;
	cv[k]-=lv[j];
      }
      else{
	j=j-1;
	cv[k]+=lv[j];
      }
    }
  }
  
  if (fetiimpl==boolean_matrices){
    
    double *a,*av;
    
    //  allocation of Boolean matrix
    a = new double [ndofcp*ndofmas[nd]];
    nullv (a,ndofcp*ndofmas[nd]);
    
    //  allocation of auxiliary vector
    av = new double [ndofcp];
    nullv (av,ndofcp);
    
    //  assembling of Boolean matrix
    boolean_matrix (nd,a);
    //  Bu
    mxv (a,lv,av,ndofcp,ndofmas[nd]);
    //
    addv (cv,av,ndofcp);
    
    delete [] a;
    delete [] av;
  }
  
}


/**
   function extracts contributions from coarse %vector to local %vector
   
   @param nd - domain id (number of required subdomain)
   @param lv - local %vector
   @param cv - coarse %vector

   JK, 14.9.2007
*/
void seqfeti::coarse_local (long nd,double *lv,double *cv)
{
  if (fetiimpl==nonredundant || fetiimpl==redundant){
    long i,j,k;
    
    for (i=0;i<ncdofd[nd];i++){
      k=ccn[nd][i]-1;
      j=edofs[nd][i];
      if (j<0){
	j=0-j-1;
	lv[j]-=cv[k];
      }
      else{
	j=j-1;
	lv[j]+=cv[k];
      }
    }
  }
  
  if (fetiimpl==boolean_matrices){
    
    double *a,*av;
    
    //  allocation of Boolean matrix
    a = new double [ndofcp*ndofmas[nd]];
    nullv (a,ndofcp*ndofmas[nd]);
    
    //  allocation of auxiliary vector
    av = new double [ndofmas[nd]];
    nullv (av,ndofmas[nd]);
    
    //  assembling of Boolean matrix
    boolean_matrix (nd,a);
    //  Bu
    mtxv (a,cv,av,ndofcp,ndofmas[nd]);
    //
    addv (lv,av,ndofmas[nd]);
    
    delete [] a;
    delete [] av;
  }
  
}







/**
   function computes size of %matrix G
   
   @param out - output stream
   
   JK, 16.9.2007
*/
void seqfeti::g_matrixsize (FILE */*out*/)
{
  long i;
  
  //  addresses computed from nrbmdom array
  rbmadr = new long [ns+1];
  for (i=0;i<ns+1;i++){
    rbmadr[i]=0;
  }
  
  //  computation of addresses
  rbmadr[0]=0;
  for (i=0;i<ns;i++){
    rbmadr[i+1]=rbmadr[i]+nrbmdom[i];
  }
  
  //  size of the matrix G
  gsize=rbmadr[ns];
  
}



/**
   function assembles %matrix G
   columns are RBM after localization
   
   JK, 16.9.2007
*/
void seqfeti::g_matrix (FILE */*out*/)
{
  long i,k,l;
  double *av;
  
  //  allocation of the array for G matrix
  if (g==NULL)
    g = new double [ndofcp*gsize];
  
  
  av = new double [ndofcp];
  
  for (i=0;i<ns;i++){
    for (k=0;k<nrbmdom[i];k++){
      nullv (av,ndofcp);
      local_coarse (i,rbmdom[i]+k*ndofmas[i],av);
      for (l=0;l<ndofcp;l++){
	g[l*gsize+rbmadr[i]+k]=0.0-av[l];
      }
    }
  }
  
  delete [] av;
}


/**
   function assembles %vector e
   
   @param f - right hand side
   
   JK, 16.9.2007
*/
void seqfeti::evector (FILE *out)
{
  long i,j,k;
  
  if (e==NULL)
    e = new double [gsize];
  
  k=0;
  for (i=0;i<ns;i++){
    for (j=0;j<nrbmdom[i];j++){
      e[k]=0.0-ss(&rbmdom[i][j*ndofmas[i]],ff[i],ndofmas[i]);
      k++;
    }
  }
  
  fprintf (out,"\n\n\n kontrola vektoru e \n");
  for (j=0;j<gsize;j++){
    fprintf (out,"\n %4ld   %f",j,e[j]);
  }
  
}

/**
   function computes the inverse %matrix to the %matrix G^T G
   
   JK, 29. 7. 2012
*/
void seqfeti::inverse_matrix_GG (FILE *out)
{
  long i,j;
  double *ee,*gg;
  
  if (invgg == NULL){
    invgg = new double [gsize*gsize];
  }
  gg = new double [gsize*gsize];
  ee = new double [gsize*gsize];
  
  for (i=0;i<gsize;i++){
    for (j=0;j<gsize;j++){
      ee[i*gsize+j]=0.0;
    }
    ee[i*gsize+i]=1.0;
  }
  
  mtxm (g,g,gg,ndofcp,gsize,gsize);
  
  fprintf (out,"\n\n\n kontrola matice G^T G \n");
  for (i=0;i<gsize;i++){
    fprintf (out,"\n %4ld ",i);
    for (j=0;j<gsize;j++){
      fprintf (out," %lf",gg[i*gsize+j]);
    }
  }
  
  gemp (gg,invgg,ee,gsize,gsize,zero,1);
  
  fprintf (out,"\n\n\n kontrola matice (G^T G)^{-1} \n");
  for (i=0;i<gsize;i++){
    fprintf (out,"\n %4ld ",i);
    for (j=0;j<gsize;j++){
      fprintf (out," %lf",invgg[i*gsize+j]);
    }
  }
  
  delete [] gg;
  delete [] ee;
}

/**
   function computes projection in FETI method
   v_new = v_old - G . (G^T . G)^{-1} . G^T . v_old
   function overwrites the %vector v by the projected %vector
   
   @param v - %vector
   
   JK, 16.9.2007
*/
void seqfeti::feti_projection (double *v)
{
  long i;
  double *p,*q,*w;
  
  p = new double [gsize];
  q = new double [gsize];
  w = new double [ndofcp];

  //  p = G^T.v
  mtxv (g,v,p,ndofcp,gsize);
  
  //  q = (G^T.G)^{-1}.p
  mxv (invgg,p,q,gsize,gsize);
  
  //  w = G.q
  mxv (g,q,w,ndofcp,gsize);
  
  //  v_new = v_old - w
  for (i=0;i<ndofcp;i++){
    v[i]-=w[i];
  }
  
  delete [] w;  delete [] q;  delete [] p;
}

/**
   function defines the compliances stored in the array h
   
   @param out - output stream
   
   JK, 3. 8. 2012
*/
void seqfeti::define_h (FILE */*out*/)
{
  long i;
  
  if (h==NULL){
    h = new double [ndofcp];
  }
  
  for (i=0;i<ndofcp;i+=2){
    h[i]=1.0e-5;
    h[i+1]=0.0;
  }
}

/**
   function defines the prescribed discontinuities stored in the array b
   
   @param out - output stream
   
   JK, 3. 8. 2012
*/
void seqfeti::define_b (FILE */*out*/)
{
  long i;
  
  if (b==NULL){
    b = new double [ndofcp];
  }
  
  for (i=0;i<ndofcp;i++){
    b[i]=0.0;
  }
}


void seqfeti::scaling(double *invect,double *outvect,long n,FILE */*out*/)
{
  double *aux;
  long i;
  aux = new double[n];
  
  for(i = 0; i < n; i++){
    aux[i] = invect[i]/wscalmat[i];
  }
  
  /*
  fprintf (out,"\n\n\n kontrola skalovani \n");
  for(i = 0; i < n; i++){
    fprintf (out,"%ld  pred %le po %le wscalmat %lf\n",i+1,invect[i],aux[i],wscalmat[i]); 
  }
  */
  
  for(i = 0; i < n; i++){
    outvect[i] = aux[i];
  }
  
  delete []aux;
}

void seqfeti::lscaling(double *invect,double *outvect,long ndom,FILE *out)
{
  double *aux;
  long i;
  
  aux = new double[nlwscalmat[ndom]];
  
  for(i = 0; i < nlwscalmat[ndom]; i++){
    aux[i] = invect[i]/lwscalmat[ndom][i];
  }
  
  
  fprintf (out,"\n\n\n kontrola skalovani \n");
  for(i = 0; i <  nlwscalmat[ndom]; i++){
    fprintf (out,"%ld  pred %le po %le wscalmat %lf ccn %ld\n",i+1,invect[i],aux[i],lwscalmat[ndom][i],ccn[ndom][i]); 
  }
  
  
  for(i = 0; i <  nlwscalmat[ndom]; i++){
    outvect[i] = aux[i];
  }
  
  delete []aux;
}



/**

*/
void seqfeti::lumpedprec (long nd,double *dd,double *pp)
{
  long i;
  double *d,*p;
  
  if (fetiimpl == nonredundant || fetiimpl==redundant){
    d = new double [ndofprec[nd]];
    p = new double [ndofprec[nd]];
    
    for (i=0;i<ndofprec[nd];i++){
      d[i]=0.0;
      p[i]=0.0;
    }
    
    for (i=0;i<ndofprec[nd];i++){
      d[i]=dd[cpreccn[nd][i]];
    }
    
    smscr[nd].mxv_scr (d,p);
    
    for (i=0;i<ndofprec[nd];i++){
      pp[cpreccn[nd][i]]=p[i];
    }
    delete [] p;
    delete [] d;
  }
  if(fetiimpl == boolean_matrices){
    
    /*for (i=0;i<ndofmas[nd];i++){
      pp[i]=dd[i];
    }
    */
    smskyprec[nd].mxv_sky (dd,pp);
    
    
  }
}

/**
   JK, 7.5.2008
*/
void seqfeti::dirichletprec (long nd,double *dd,double *pp)
{
  long i,j;
  double *d,*p;
  
  if (fetiimpl == nonredundant || fetiimpl==redundant){
    d = new double [ndofprec[nd]];
    p = new double [ndofprec[nd]];
    
    for (i=0;i<ndofprec[nd];i++){
      d[i]=0.0;
      p[i]=0.0;
    }
    
    for (i=0;i<ndofprec[nd];i++){
    d[i]=dd[cpreccn[nd][i]];
    }
    
    psmdm[nd].mxv_dm (d,p);
    
    
    for (i=0;i<ndofprec[nd];i++){
      pp[cpreccn[nd][i]]=p[i];
    }
    
    delete [] p;
    delete [] d;
  }
  if (fetiimpl == boolean_matrices){
    d = new double [ndofprec[nd]];
    p = new double [ndofprec[nd]];
    j = 0;
    for (i=0;i<ndofmas[nd];i++){
      if(dd[i] != 0){
	d[j]=dd[i];
	j++;
      }
    }
    psmdm[nd].mxv_dm (d,p);
    
    j = 0;
    for (i=0;i<ndofmas[nd];i++){
      if(dd[i] != 0){
	pp[i]=p[j];
	j++;
      }
      else{
	pp[i] = 0;
      }
    }
  }
}



/**
   function performs the modified conjugate gradient method
   
   @param out - output stream
   
   for details se the book J. Kruis: Domain Decomposition Methods
   for Distributed Computing. Saxe-Coburg Publications, 2006
   
   JK, 16.9.2007, revised 29. 7. 2012
*/
void seqfeti::mpcg (FILE *out)
{
  long i,j;
  double nom,denom,alpha,beta;
  double *s,*sd,*r,*p,*pp,*invgge;
  
  if (lambda==NULL)
    lambda = new double [ndofcp];
  
  if (ndofmax==0){
    print_err("the attribute ndofmax is not determined",__FILE__,__LINE__,__func__);
  }
  
  //  vectors defined on each subdomain
  sd = new double [ndofmax];
  pp = new double [ndofmax];
  
  //  direction vector allocation
  s = new double [ndofcp];
  //  residual vector allocation
  r = new double [ndofcp];
  //  auxiliary vector allocation
  p = new double [ndofcp];
  
  
  //  initiation of conjugate gradients
  invgge = new double [gsize];
  //  (G^T.G)^{-1} e
  mxv (invgg,e,invgge,gsize,gsize);
  //  lambda = G [(G^T.G)^{-1} e]
  mxv (g,invgge,lambda,ndofcp,gsize);
  delete [] invgge;
  
  nullv (r,ndofcp);
  for (j=0;j<ns;j++){
    nullv (sd,ndofmas[j]);
    coarse_local (j,sd,lambda);
    //  \lambda_0 - f
    //  vector dd is rewritten, ff cannot be rewritten
    //  therefore, lambda_0 - f is used instead of f - lambda_0
    subv (sd,ff[j],ndofmas[j]);
    
    //  K^+ (\lambda_0 - f) = -r
    if (smst==dense_matrix)
      break;
    if (smst==skyline_matrix)
      smsky[j].ldl_feti_sky (pp,sd,nrbmdom[j],se[j],zero);
    
    //  odecet skoku
    //subv (pp,jum[j],ndofmas[j]);
    
    local_coarse (j,pp,r);
  }
  //  the vector r contains gradient
  cmulv(-1.0,r,ndofcp);
  //  the vector r contains residual now
  
  if (discont==1){
    //  there are discontinuities
    
    //  loop over the number of DOFs in the coarse problem
    for (i=0;i<ndofcp;i++){
      r[i]-=b[i]+h[i]*lambda[i];
    }//  end of the loop over the number of DOFs in the coarse problem
  }//  end of the statement if (discont==1)
  
  //  P r
  feti_projection (r);
  
  //  initialization of direction vector
  for (i=0;i<ndofcp;i++){
    s[i]=r[i];
  }
  
  //  nominator evaluation
  nom = ss (r,r,ndofcp);
  
  
  // *************************
  //  main iteration loop   **
  // *************************
  for (i=0;i<nicg;i++){
    
    // ******************************************
    //   auxiliary vector computation F.s = p  **
    // ******************************************
    nullv (p,ndofcp);
    for (j=0;j<ns;j++){
      nullv (sd,ndofmax);
      coarse_local (j,sd,s);
      
      //  computation of K^+ d on master
      if (smst==dense_matrix)
	break;
      if (smst==skyline_matrix)
	smsky[j].ldl_feti_sky (pp,sd,nrbmdom[j],se[j],zero);
      
      local_coarse (j,pp,p);
    }

    if (discont==1){
      //  there are discontinuities
      
      //  loop over the number of DOFs in the coarse problem
      for (j=0;j<ndofcp;j++){
	//p[j]+=h[j]*lambda[j];
	p[j]+=h[j]*s[j];
      }//  end of the loop over the number of DOFs in the coarse problem
    }//  end of the statement if (discont==1)
    
    //  denominator of alpha
    denom = ss (s,p,ndofcp);
    
    fprintf (out,"\n MPCG in FETI, iterace %ld, kontrola citatele a jmenovatele alpha %e / %e",i,nom,denom);
    
    if (fabs(denom)<zero){
      print_err("zero denominator of alpha in the modified conjugate gradient method",__FILE__,__LINE__,__func__);
      break;
    }
    
    // *****************************
    //   determination of  alpha  **
    // *****************************
    alpha = nom/denom;
    
    // **************************************************************
    //  new residual and new approximation of Lagrange multipliers
    // **************************************************************
    for (j=0;j<ndofcp;j++){
      r[j]-=alpha*p[j];
      lambda[j]+=alpha*s[j];
    }
    
    feti_projection (r);
    
    
    // *******************************************************
    // *******************************************************
    //  vlozka s predpodminenim
    // *******************************************************
    // *******************************************************
    /*
    nullv (p,ndofcp);
    
    switch (prec){
    case noprecond:{
      for (j=0;j<ns;j++){
	nullv (dd,ndofmax);
	coarse_local (j,dd,g);
	nullv (pp,ndofmax);
	for (k=0;k<ndofmas[j];k++){
	  pp[k]=dd[k];
	}
	local_coarse (j,pp,p);
      }
      break;
    }
    case feti_lumped:{
      //scaling(g,g,ndofcp,out);
      for (j=0;j<ns;j++){
	nullv (dd,ndofmax);
	coarse_local (j,dd,g);
	//lscaling(dd,dd,j,out);
	nullv (pp,ndofmax);
	
	fprintf (out,"\n\nkontrola v lumpedprec\n");
	for (k=0;k<ndofmas[j];k++){
	  fprintf (out,"%lf\n",dd[k]);
	}
	
	lumpedprec (j,dd,pp);
		
	//lscaling(pp,pp,j,out);
	local_coarse (j,pp,p);
      }
      //scaling(p,p,ndofcp,out);
      break;
    }
    case feti_dirichlet:{
      //scaling(g,g,ndofcp,out);
      for (j=0;j<ns;j++){
	nullv (dd,ndofmax);
	coarse_local (j,dd,g);
	//lscaling(dd,dd,j,out);
	nullv (pp,ndofmax);
	
	fprintf (out,"\n\nkontrola v dirichlet\n");
	for (k=0;k<ndofmas[j];k++){
	  fprintf (out,"%lf\n",dd[k]);
	}
	dirichletprec (j,dd,pp);
	//for (k=0;k<ndofmas[j];k++){
	//pp[k]=dd[k];
	//}
	//lscaling(pp,pp,j,out);
	local_coarse (j,pp,p);
      }
      //scaling(p,p,ndofcp,out);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of preconditioner is required in function mpcg (file %s, line %d)\n",__FILE__,__LINE__);
    }
    }
    
    feti_projection (p,h,h1);
    */
    // *******************************************************
    // *******************************************************
    //   konec vlozky s predpodminenim
    // *******************************************************
    // *******************************************************
    
    
    denom = nom;
    if (fabs(denom)<zero){
      print_err("zero denominator of beta in the modified conjugate gradient method",__FILE__,__LINE__,__func__);
      break;
    }
    
    nom=ss(r,r,ndofcp);
    
    fprintf (stdout,"\n iteration  %4ld        norm r   %e",i,sqrt(nom));
    fprintf (out,"\n iteration  %4ld        norm r   %e",i,sqrt(nom));
    
    fprintf (out,"\n kontrola citatele a jmenovatele pred betou  %e / %e",nom,denom);
    
    if (sqrt(nom)<errcg){
      break;
    }
    
    
    //  computation of beta coefficient
    beta = nom/denom;
    
    //  new direction vector
    for (j=0;j<ndofcp;j++){
      s[j]=beta*s[j]+r[j];
    }
    
  }
  
  //  the number of iterations performed
  anicg=i;
  //  the attained norm of the residual
  aerrcg=nom;
  
  delete [] sd;
  delete [] pp;
  delete [] s;
  delete [] r;
  delete [] p;

}





/**
   function performs modified preconditioned reorthonormalized conjugate gradient method
   
   @param w - %vector of Lagrange multipliers
   @param ff - %vectors of right hand side
   @param h - %matrix H
   @param h1 - %matrix (H^T.H)^{-1}
   @param ndofcp,hsize - number of rows and columns of %matrix H
   
   ndofcp - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 5.5.2008
*/
 /*
void seqfeti::mprcg (gtopology *top,gmatrix *gm,
		     double *lambda,double **ff,double *e,double *g,double *g1,FILE *out)
{
  long i,j,k;
  double nom,denom,alpha,s;
  double *ge,*r,*lambdad,*rr,*w,*p,*y,*op,*fp,*pp,*fpd,*ofp,*coef,*wd,*yd,*beta;

    r = new double [ndofcp];
  lambdad = new double [ndofmax];
  rr = new double [ndofmax];
  w = new double [ndofcp];
  p = new double [ndofcp];
  y = new double [ndofcp];
  op = new double [ndofcp*101];
  fp = new double [ndofcp];
  pp = new double [ndofmax];
  fpd = new double [ndofmax];
  ofp = new double [ndofcp*101];
  coef = new double [101];
  wd = new double [ndofmax];
  yd = new double [ndofmax];
  beta = new double [101];


  //  vectors defined on each subdomain

  


  //  initiation of conjugate gradients
  
  ge = new double [hsize];
  //  (G^T.G)^{-1} e
  mxv (g1,e,ge,hsize,hsize);
  //  G [(G^T.G)^{-1} e]
  mxv (g,ge,lambda,ndofcp,hsize);
  delete [] ge;
  
  nullv (w,ndofcp);
  for (j=0;j<ns;j++){
    nullv (lambdad,ndofmas[j]);
    coarse_local (j,lambdad,lambda);
    // f - \lambda_0
    subv (lambdad,ff[j],ndofmas[j]);
    
    //  K^+ (f - \lambda_0) = r
    if (smst==dense_matrix)
      break;
    if (smst==skyline_matrix)
      smsky[j].ldl_feti_sky (wd,lambdad,nrbmdom[j],se[j],zero);
    
    
    //  zacatek - streda 28.11.2007
    //  odecet skoku
    //subv (pp,jum[j],ndofmas[j]);
    //  konec - streda 28.11.2007
    
    
    //nullvr (buff,maxlggl+1);
    local_coarse (j,wd,w);
  }
  
  for (j=0;j<ndofcp;j++){
    w[j]*=-1.0;
  }

  //  P r
  feti_projection (w,g,g1);
  
  //  initialization of direction vector
  for (i=0;i<ndofcp;i++){
    p[i]=w[i];
    y[i]=w[i];
  }

  for (j=0;j<ndofcp;j++){
    op[0*ndofcp+j]=p[j];
  }
  
  
  // *************************
  //  main iteration loop   **
  // *************************
  for (i=0;i<nicg;i++){
    
    // ******************************************
    //   auxiliary vector computation K.d = p  **
    // ******************************************
    nullv (fp,ndofcp);
    for (j=0;j<ns;j++){
      nullv (pp,ndofmax);
      coarse_local (j,pp,p);
      
      //  computation of K^+ d on master
      if (smst==dense_matrix)
	break;
      if (smst==skyline_matrix)
	smsky[j].ldl_feti_sky (fpd,pp,nrbmdom[j],se[j],zero);
      
      //nullvr (buff,maxlggl+1);
      local_coarse (j,fpd,fp);
    }
    
    for (j=0;j<ndofcp;j++){
      ofp[i*ndofcp+j]=fp[j];
    }
    
    nom = ss (p,w,ndofcp);

    //  denominator of alpha
    denom = ss (p,fp,ndofcp);
    
    coef[i]=denom;

    //fprintf (stdout,"\n denominator u alpha   %le",denom);
    
    fprintf (out,"\n\n kontrola citatele a jmenovatele pred alpha %e / %e",nom,denom);
    //fprintf (stderr,"\n\n kontrola citatele a jmenovatele pred alpha %lf / %lf",nom,denom);
    
    if (fabs(denom)<zero){
      fprintf (stderr,"\n\n zero denominator in modified conjugate gradient method (file %s, line %d).\n",__FILE__,__LINE__);
      break;
    }
    
    // *******************
    //   vypocet alpha  **
    // *******************
    alpha = nom/denom;
    
    // **************************************************************
    //  vypocet noveho gradientu g a nove aproximace neznamych x   **
    // **************************************************************
    nom=0.0;
    for (j=0;j<ndofcp;j++){
      lambda[j]+=alpha*p[j];
      w[j]-=alpha*fp[j];
      nom+=w[j]*w[j];
    }
    
    fprintf (stdout,"\n iteration  %4ld        norm g   %e",i,sqrt(nom));
    fprintf (out,"\n iteration  %4ld        norm g   %e",i,sqrt(nom));

    feti_projection (w,g,g1);
    
    // *******************************************************
    // *******************************************************
    //  vlozka s predpodminenim
    // *******************************************************
    // *******************************************************
    
    nullv (y,ndofcp);
    
    switch (prec){
    case noprecond:{
      for (j=0;j<ndofcp;j++){
	y[j]=w[j];
      }
      break;
    }
    case feti_lumped:{
      //scaling(w,w,ndofcp,out);
      for (j=0;j<ns;j++){
	nullv (wd,ndofmax);
	coarse_local (j,wd,w);
	//lscaling(dd,dd,j,out);
	nullv (pp,ndofmax);
	
	lumpedprec (j,wd,yd);
	
	//lscaling(pp,pp,j,out);
	local_coarse (j,yd,y);
      }
      //scaling(y,y,ndofcp,out);
      break;
    }
    default:{
      fprintf (stderr,"\n\n unknown type of preconditioner is required in function mpcg (file %s, line %d)\n",__FILE__,__LINE__);
    }
    }
    
    feti_projection (y,g,g1);
    
    // *******************************************************
    // *******************************************************
    //   konec vlozky s predpodminenim
    // *******************************************************
    // *******************************************************
    
    for (k=0;k<=i;k++){
      s=0.0;
      for (j=0;j<ndofcp;j++){
	s+=y[j]*ofp[k*ndofcp+j];
      }
      beta[k]=s/coef[k];
    }
    for (j=0;j<ndofcp;j++){
      s=0.0;
      for (k=0;k<=i;k++){
	s+=beta[k]*op[k*ndofcp+j];
      }
      p[j]=y[j]-s;
      op[(i+1)*ndofcp+j]=p[j];
    }
    
    
    if (sqrt(nom)<errcg){
      break;
    }
    
    
  }
  
  anicg=i;  aerrcg=nom;
  
  fprintf (out,"\n\n\n\n kontrola Lagrangeovych multiplikatoru \n");
  for (i=0;i<ndofcp;i++){
    //fprintf (out,"\n lagr. mult %4ld    %e",i,w[i]);
  }
  
}
 */

/**
   function computes the nodal unknowns on subdomains
   
   @param lhs - array containing the left hand side
   @param out - output file
   
   JK, 29. 7. 2012
*/
 void seqfeti::nodalunknowns (double *lhs, FILE *out)
{
  long i,j;
  
  if (d == NULL){
    d = new double* [ns];
    for (i=0;i<ns;i++){
      d[i] = new double [ndofmas[i]];
      for (j=0;j<ndofmas[i];j++){
	d[i][j]=0.0;
      }
    }
  }
  else{
    for (i=0;i<ns;i++){
      for (j=0;j<ndofmas[i];j++){
	d[i][j]=0.0;
      }
    }
  }
  
  lagrmultnodunknowns (out);
  
  
  for (i=0;i<ns;i++){
    for (j=0;j<ndofmas[i];j++){
      lhs[cndom[i][j]-1]=d[i][j];
    }
  }
  
}



/**
   function computes nodal unknowns from Lagrange multipliers
   function is used in the FETI method
   
   @param out - output stream
   
   JK, 16.9.2007
*/
void seqfeti::lagrmultnodunknowns (FILE *out)
{
  long i,j,k,l;
  double *r,*pp,*alpha,*u,*av;
  
  av = new double [6];
  pp = new double [ndofmax];
  r = new double [ndofcp];
  alpha = new double [gsize];
  u = new double [gsize];
  
  //  zeroing of the array r
  nullv (r,ndofcp);
  //  loop over the number of subdomains
  for (j=0;j<ns;j++){
    nullv (pp,ndofmax);
    //  multipliers connected to the j-th subdomain are selected
    coarse_local (j,pp,lambda);
    
    //  pp = f - lambda
    subv (pp,ff[j],ndofmas[j]);
    
    //  K^+ (f-\lambda) = d 
    if (smst==dense_matrix)
      break;
    if (smst==skyline_matrix)
      smsky[j].ldl_feti_sky (d[j],pp,nrbmdom[j],se[j],zero);
    
    //  vector d contains -1 multiple of K^+(f-B^T\lambda) at this moment
    local_coarse (j,d[j],r);
  }
  
  cmulv (-1.0,r,ndofcp);
  
  fprintf (out,"\n");
  for (i=0;i<ndofcp;i++){
    fprintf (out,"\n g  %20.15le",r[i]);
  }
  
  if (discont==1){
    //  there are discontinuities
    
    //  loop over the number of DOFs in the coarse problem
    for (i=0;i<ndofcp;i++){
      r[i]-=b[i]+h[i]*lambda[i];
    }//  end of the loop over the number of DOFs in the coarse problem
  }
  
  //  G^T (g-F \lambda -H\lambda-const)
  mtxv (g,r,u,ndofcp,gsize);
  //  (G^T G)^{-1} G^T (g-F \lambda)
  mxv (invgg,u,alpha,gsize,gsize);
  //cmulv (-1.0,alpha,gsize);
  
  for (i=0;i<ns;i++){
    //  vector d contains -1 multiple of K^+(f-B^T\lambda) at this moment
    cmulv (-1.0,d[i],ndofmas[i]);
    
    nullv (av,6);
    l=rbmadr[i];
    for (j=0;j<nrbmdom[i];j++){
      av[j] = alpha[l];  l++;
    }
    
    for (j=0;j<ndofmas[i];j++){
      for (k=0;k<nrbmdom[i];k++){
	d[i][j]+=rbmdom[i][k*ndofmas[i]+j]*av[k];
      }
    }
  }
  
  delete [] av;
  delete [] pp;
  delete [] r;
  delete [] alpha;
  delete [] u;
}




/**
   function assembles the right hand sides for subdomains
   
   @param rhs - array containing the whole %vector of the right hand side
   
   29. 7. 2012, JK
*/
void seqfeti::assemble_ff (double *rhs)
{
  long i,j;
  
  if (ff==NULL){
    ff = new double* [ns];
    for (i=0;i<ns;i++){
      ff[i] = new double [ndofmas[i]];
      for (j=0;j<ndofmas[i];j++){
	ff[i][j]=0.0;
      }
    }
  }else{
    for (i=0;i<ns;i++){
      for (j=0;j<ndofmas[i];j++){
	ff[i][j]=0.0;
      }
    }
  }
  
  //  allocation of vectors of right hand sides on subdomains
  for (i=0;i<ns;i++){
    globloc (rhs,ff[i],cndom[i],ndofmas[i]);
  }
  
}


/**
   function assembles the matrices G, (G^T.G)^{-1}
   
   @param gm - %matrix of the system
   @param out - output stream
   
   JK, 2. 8. 2012
*/
void seqfeti::matrices_assembl (gmatrix *gm,FILE *out)
{
  long i,j;
  time_t t1,t2,t3,t4,t5;

  //  determination of the ndofmax
  det_ndofmax ();

  t1 = time (NULL);

  //  assembling of subdomain matrices
  subdomain_matrices (gm,out);

  t2 = time (NULL);

  //  computation of kernels / rigid body modes
  kernel (out);

  //  investigation of the G matrix sizes
  g_matrixsize (out);
  
  if (gsize==0){
    print_err("gsize is 0",__FILE__,__LINE__,__func__);
    abort ();
  }
  
  t3 = time (NULL);

  //  assembling of G matrix
  g_matrix (out);

  t4 = time (NULL);

  fprintf (out,"\n\n\n kontrola matice G \n");
  for (i=0;i<ndofcp;i++){
    fprintf (out,"\n %4ld",i);
    for (j=0;j<gsize;j++){
      fprintf (out,"  %f",g[i*gsize+j]);
    }
  }

  //  computation of the matrix (G^T G)^{-1}
  inverse_matrix_GG (out);

  t5 = time (NULL);

  fprintf (stdout,"\n time of subdomain matrices assembling     %ld", long(t2-t1));
  fprintf (stdout,"\n time of kernel computation                %ld", long(t3-t2));
  fprintf (stdout,"\n time of G matrix assembling               %ld", long(t4-t3));
  fprintf (stdout,"\n time of inverse matrix G^T G computation  %ld", long(t5-t4));
  
  //  matrices of the FETI method are assembled
  matassem=1;
}

/**
   function assembles the matrices G, (G^T.G)^{-1}
   
   @param gm - %matrix of the system
   @param out - output stream
   
   JK, 2. 8. 2012
*/
void seqfeti::vectors_assembl (double *rhs,FILE *out)
{
  time_t t1,t2,t3;
  
  t1 = time (NULL);

  //  assembling of the vectors of the right hand side for each subdomain
  assemble_ff (rhs);

  t2 = time (NULL);

  //  assembling of the e vector
  evector (out);
  
  t3 = time (NULL);
  
  fprintf (stdout,"\n time of f assembling   %ld",long(t2-t1));
  fprintf (stdout,"\n time of e assembling   %ld",long(t3-t2));
}

/**
   function solves system of equations by the FETI method
   
   @param top - pointer to object with topology
   @param gm - %matrix of the system
   @param lhs - %vector of solution (left hand side)
   @param rhs - %vector of right hand side
   @param out - output stream
   
   JK, 16.9.2007
*/
void seqfeti::solve_system (gtopology */*top*/,gmatrix *gm,
			    double *lhs,double *rhs,FILE *out)
{
  long i,j;
  
  if (matassem==0)
    matrices_assembl (gm,out);
  
  vectors_assembl (rhs,out);

  fprintf (stdout,"\n\n gsize   %ld \n",gsize);
  fprintf (out,"\n ndofcp   %ld",ndofcp);
  fprintf (stdout,"\n ndofcp   %ld",ndofcp);
  
  define_b (out);
  define_h (out);

  //  modified conjugate gradient method
  mpcg (out);
  
  fprintf (out,"\n\n\n kontrola multiplikatoru \n");
  for (j=0;j<ndofcp;j++){
    fprintf (out,"\n %4ld   %20.15le",j,lambda[j]);
  }
  
  
  nodalunknowns (lhs,out);

  fprintf (out,"\n\n kontrola");
  for (i=0;i<ns;i++){
    fprintf (out,"\n");
    for (j=0;j<ndofmas[i];j++){
      fprintf (out,"\n displ %4ld %4ld     %20.15le",i,j,d[i][j]);
    }
  }
  
}
