#include "mpi.h"
#include "schurcompl.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "precond.h"
#include "iotools.h"

/**
   constructor
   
   @param np - number of processors
   @param mr - my rank (process id)
   @param nd - number of subdomains
   
   JK, revised 24.7.2008
*/
schurcompl::schurcompl(int np,int mr,long nd)
{
  //  number of processors
  nproc=np;
  //  my rank (process id)
  myrank=mr;
  //  number of subdomains
  ndom=nd;
  
  //  number of DOFs on subdomain
  ndof=0;
  //  number of internal DOFs on subdomain
  nidof=0;
  //  number of boundary/interface DOFs on subdomain
  nbdof=0;
  //  maximum number of boundary/interface DOFs on subdomain
  maxnbdof=0;

  // **************************************************************************
  //  the following variables are initialized in subroutine psolver::initiate
  // **************************************************************************
  //  type of solver of reduced/coarse problem
  trssol = (redsystsolver) 0;
  //  type of storage of matrix of reduced/coarse problem
  rsmstor = (storagetype) 0;
  //  maximum number of iterations in the conjugate gradient method used for solution of reduced problem
  nicgsch=0;
  //  required/prescribed norm of residuum
  errcgsch=0.0;
  // ***************************************************************
  //  end of variables initialized in subroutine psolver::initiate
  // ***************************************************************
  
  //  actual number of iterations (number of iterations performed during the solution)
  anicgsch=0;
  //  attained norm of residuum
  aerrcgsch=0.0;
  //  computer zero
  zero=1.0e-15;
  
  //  array containing numbers of boundary/interface DOFs on subdomains
  nbdofmas=NULL;
  
  
  gtop=NULL;
  arr=NULL;
  
  //  solver of system of linear eqautions
  ssle = new slesolv ();
  // deallocation of nbdofmas array in the destructor is not required by default
  destr_nbdofmas = 0;
}

schurcompl::~schurcompl()
{
  if (myrank==0){
    if (destr_nbdofmas)
      delete [] nbdofmas;
    
    delete gtop;
    delete arr;
  }
  delete ssle;
}



/**
  Function search prescribed common code numbers on boundary nodes in the given subdomain.
  @param top   - pointer to the sequential general topology
  @param ltg   - local to global array
  @param idgnn - array with node numbers of boundary nodes
  @param nbn   - number of boundary nodes
  @param gnnc  - array with boundary node numbers at positions of nodes with prescribed common code numbers 
  @param id    - searched node id
  @param ccn   - searched common code number id

*/
/*
void schurcompl::search_comcn(gtopology *top,long *ltg, long *idgnn, long nbn, long **gnnc, long id, long ccn)
{
  long i, j, k, ii, ndofn, flag;
  for (i=0; i<nbn; i++)
  {
    flag = 0;
    ii = idgnn[i];
    ndofn=top->give_ndofn(ii);
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k == ccn)
      {
        i = nbn;
        flag = 1;
        break;
      }
    }
  }
  if (flag)
  {
    for (i=id;i<top->nn;i++){
      for (j=0;j<ndofn;j++){
        k=top->give_dof (i,j);
        if (k == ccn)
        {
          gnnc[i][k]=ltg[ii];
        }
      }
    }
  }
}

*/

/**
   function generates ordering of unknowns on subdomains
   
   @param top - pointer to the sequential general topology
   @param ltg - local to global array
   
   JK, 7.8.2007
*/
/*
void schurcompl::subdomain_ordering (gtopology *top,long *ltg,FILE *out)
{
  long i,j,k,ii,ndofn,*aux,*idgnn,**gnnc, nbn;
  
  //  searching of maximum code number
  ndof=0;
  for (i=0;i<top->nn;i++){
    ndofn=top->give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>ndof)  ndof=k;
    }
  }
  ndof--;
  if (ndof<0)  ndof=0;
  aux = new long [ndof];
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }

  gnnc = new long*[top->nn];
  memset(gnnc, 0, sizeof(*gnnc)*top->nn);
  idgnn = new long[top->nn];
  memset(idgnn, 0, sizeof(*idgnn)*top->nn);
  j = 0;
  // array with indeces of boundary nodes
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1)
    {
      idgnn[j]=ltg[i];
      j++;
    }
  }
  nbn = j;
  for (i=0;i<top->nn;i++){
    ndofn=top->give_ndofn (i);
    gnnc[i] = new long[ndofn];
    for (j=0;j<ndofn;j++){
      gnnc[i][j]=-1;
    }
  }
  //  contributions from nodes
  ndof=1;
  for (i=0;i<top->nn;i++){
    ii=top->ordering[i];
    if (ltg[ii]==-1){
      //  inner nodes
      ndofn=top->give_ndofn (ii);
      for (j=0;j<ndofn;j++){
	k=top->give_dof(ii,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k==1){
	  top->save_dof (ii,j,ndof);  ndof++;
	}
        // common code numbers
        if (k>1){
          // node has not assigned code number 
	  if ((aux[k-2]==-1) && (gnnc[ii][j]==-1))
          {
            // search for the same common code number on the boundary nodes
            search_comcn(top,ltg,idgnn,nbn,gnnc,ii,k);
          }
          if ((aux[k-2]==-1) && (gnnc[ii][j]==-1))
          {
            // node has not assigned common code number and it is not coupled with the domain boundary
            top->save_dof (ii,j,ndof);
            aux[k-2]=ndof;
            ndof++;
          }
          if ((aux[k-2]>-1) && (gnnc[ii][j]==-1))
          {
            // node has already assigned common code number and  it is not coupled with the domain boundary
            top->save_dof (ii,j,aux[k-2]);
          }
        }
      }
    }
  }
  nidof=ndof-1;
  top->nidof=nidof;

  for (i=0;i<top->nn;i++){
    ii=top->ordering[i];
    if (ltg[ii]>-1){
      //  boundary nodes
      ndofn=top->give_ndofn (ii);
      for (j=0;j<ndofn;j++){
	k=top->give_dof (ii,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k==1){
	  top->save_dof (ii,j,ndof);  ndof++;
	}
        // common code numbers
        if (k>1)
        {
          if (aux[k-2]==-1){
            // node has not assigned common code number
            top->save_dof (ii,j,ndof);
            aux[k-2]=ndof;
            ndof++;
          }
          else{
            // node has already assigned common code number
            top->save_dof (ii,j,aux[k-2]);
          }
        }
      }
    }
  }
  // numbering of nodes which are not on boundary and they have not assigned
  // common code number with the boundary nodes
  for (i=0;i<top->nn;i++){
    ii=top->ordering[i];
    ndofn=top->give_ndofn (ii);
    for (j=0;j<ndofn;j++){
      k=top->give_dof(ii,j);
      if ((k>1) && (gnnc[ii][j]!=-1) && (ltg[ii]==-1)){
        // node has already assigned common code number from boundary node and it does not lie on boundary
        top->save_dof (ii,j,aux[k-2]);
      }
    }
  }
  ndof--;
  
  nbdof=ndof-nidof;
  top->nbdof=nbdof;

  delete [] aux;
  delete [] idgnn;
  for (i=0;i<top->nn;i++)
    delete [] gnnc[i];
  delete [] gnnc;
}

*/

/**
   function orders unknowns on subdomains with respect of Schur complement method
   inner variables are ordered at the beginning
   boundary variables are ordered at he end
   
   @param top - pointer to domain topology
   @param ptop - pointer to topology for paralle computation
   
   JK, 8.7.2005
*/
/*
void schurcompl::subdomain_ordering (gtopology *top,partop *ptop)
{
  long i,j,k,nn,nbn,nin,ndofn,nid;
  
  //  number of all nodes on subdomain
  nn = ptop->nn;
  //  number of internal nodes on subdomain
  nin = ptop->nin;
  //  number of boundary nodes on subdomain
  nbn = ptop->nbn;
  

  // ********************************
  //  ordering of internal unknowns
  // ********************************
  ndof=1;
  for (i=0;i<nin;i++){
    //  node id (in local ordering)
    nid=ptop->lnin[i];
    //  number of DOFs at required node
    ndofn = top->give_ndofn (nid);
    
    for (j=0;j<ndofn;j++){
      k=top->give_dof (nid,j);
      if (k>0){
	top->save_dof (nid,j,ndof);
	ndof++;
      }
    }
  }
  //  number of internal DOFs (unknowns)
  nidof = ndof-1;

  
  // ********************************
  //  ordering of boundary unknowns
  // ********************************
  for (i=0;i<nbn;i++){
    //  node id (in local ordering)
    nid=ptop->lnbn[i];
    //  number of DOFs at required node
    ndofn = top->give_ndofn (nid);
    
    for (j=0;j<ndofn;j++){
      k=top->give_dof (nid,j);
      if (k>0){
	top->save_dof (nid,j,ndof);
	ndof++;
      }
    }
  }
  ndof--;
  

  //  contributions from elements
  for (i=0;i<top->ne;i++){
    if (top->give_cne (i)==1){
      for (j=0;j<top->give_ndofe(i);j++){
	if (top->give_dof (i,j)>ndof)  ndof=top->give_dof (i,j);
      }
    }
  }
  
  //  number of boundary DOFs (unknowns)
  nbdof = ndof-nidof;
  
  //  state of code numbers is changed
  //  code numbers are generated
  top->cnstate=1;
  
}
*/

/**
   JK, 7.8.2007
*/
void schurcompl::initiate (partop *ptop,selnodes *selnodschur)
{
  //  number of degrees of freedom on subdomain
  ndof = ptop->ndof;
  //  number of internal degrees of freedom (unknowns) on subdomain
  nidof = ptop->nidof;
  //  number of boundary degrees of freedom (unknowns) on subdomain
  nbdof = ptop->nbdof;
  
  //  maximum number of reduced DOFs on subdomain
  maxnbdof = selnodschur->maxsnndof;
  
  if (myrank==0){
    //  number of global DOFs = total number of boundary DOFs
    //  number of DOFs of coarse problem
    ndofcp = selnodschur->tndofsn;
    
    //  array containing numbers of boundary DOFs on subdomains
    nbdofmas = selnodschur->snndofmas;
    
    gtop = new gtopology;
    gtop->initiate (selnodschur->cnmas,selnodschur->snndofmas,nproc);
    gtop->cnstate=1;
  }
}

void schurcompl::initiate (partopjk *ptop,FILE *out)
{
  long i;
  
  //  number of degrees of freedom on subdomain
  ndof = ptop->ndof;
  //  number of internal degrees of freedom (unknowns) on subdomain
  nidof = ptop->nidof;
  //  number of boundary degrees of freedom (unknowns) on subdomain
  nbdof = ptop->nbdof;
  
  //  maximum number of reduced DOFs on subdomain
  maxnbdof = ptop->maxnbdofd;
  
  fprintf (out,"\n\n jouda \n\n");
  fprintf (out,"\n ndof   %ld",ndof);
  fprintf (out,"\n nidof  %ld",nidof);
  fprintf (out,"\n nbdof  %ld",nbdof);
  fprintf (out,"\n maxnbdof  %ld",maxnbdof);
  
  if (myrank==0){
    //  number of global DOFs = total number of boundary DOFs
    //  number of DOFs of coarse problem
    ndofcp = ptop->ndofc;

    fprintf (out,"\n ndofcp  %ld",ndofcp);
    
    //  array containing numbers of boundary DOFs on subdomains
    nbdofmas = ptop->nbdofd;
    
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nbdofmas %ld  %ld",i,nbdofmas[i]);
    }
    
    gtop = new gtopology;
    gtop->initiate (ptop->bdofd,ptop->nbdofd,nproc);
    gtop->cnstate=1;
  }

}

/**
   function orders unknowns in coarse problem
   
   function initializates object gtop of the class gtopology
   it represents subdomains as superelements
   
   @param ptop - pointer to parallel topology
   
   JK, 8.7.2005
*/
/*
void schurcompl::coarse_problem_ordering (partop *ptop,FILE *out)
{
  long i,j,k,m,tnbn;
  long *nbnd,*ndofncn;
  long **cnbn,**nbdofnd,**cncn,**bcn;
  long ***lbcn;
  MPI_Status stat;
  
  if (myrank==0){


    //  array of numbers of boundary nodes on subdomains
    nbnd = ptop->nbnd;
    //  array containing coarse node numbers of boundary nodes
    cnbn = ptop->cnbn;
    //  array containing numbers of DOFs on boundary nodes
    nbdofnd = ptop->nbdofnd;
    //  total number of boundary nodes
    tnbn = ptop->tnbn;
    //  array containing local boundary code numbers
    lbcn = ptop->lbcn;
    
    if (nbnd==NULL){
      fprintf (stderr,"\n\n array containing local boundary code numbers (lbcn) is not allocated");
      fprintf (stderr,"\n in function coarse_problem_ordering (file %s, line %d)",__FILE__,__LINE__);
      fprintf (stderr,"\n function code_numbers_on_master (class partop) should be called first,");
    }
    
    //  number of DOFs at coarse nodes
    ndofncn = new long [tnbn];
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nbnd[i];j++){
	k=cnbn[i][j];
	ndofncn[k]=nbdofnd[i][j];
      }
    }
    

    //  assignment of constraint indicators
    cncn = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      cncn[i] = new long [ndofncn[i]];
    }


    //  number of boundary DOFs on subdomains
    //  it must not be deleted, it is used in further subroutines
    nbdofmas = new long [nproc];
    destr_nbdofmas = 1;  // flag for deallocation of nbdofmas in the destructor

    maxnbdof=0;
    for (i=0;i<nproc;i++){
      nbdofmas[i]=0;
      for (j=0;j<nbnd[i];j++){
	m=cnbn[i][j];
	for (k=0;k<nbdofnd[i][j];k++){
	  cncn[m][k]=lbcn[i][j][k];
	  if (lbcn[i][j][k]>0)
	    nbdofmas[i]++;
	}
      }
      if (maxnbdof<nbdofmas[i])
	maxnbdof=nbdofmas[i];
    }
    

    //  generation of code numbers (unknown ordering)
    ndofcp=1;
    for (i=0;i<tnbn;i++){
      for (j=0;j<ndofncn[i];j++){
	if (cncn[i][j]>0){
	  cncn[i][j]=ndofcp;
	  ndofcp++;
	}
      }
    }
    ndofcp--;
    
    
    bcn = new long* [nproc];
    for (i=0;i<nproc;i++){
      bcn[i] = new long [nbdofmas[i]];
      nbdofmas[i]=0;
    }
    
    //  assembling of lists of boundary code numbers on subdomains
    for (i=0;i<nproc;i++){
      for (j=0;j<nbnd[i];j++){
	m=cnbn[i][j];
	for (k=0;k<ndofncn[m];k++){
	  if (cncn[m][k]>0){
	    bcn[i][nbdofmas[i]]=cncn[m][k];
	    nbdofmas[i]++;
	  }
	}
      }
    }
    

    fprintf (out,"\n\n\n kontrola bcn \n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %ld",i+1);
      for (j=0;j<nbdofmas[i];j++){
	fprintf (out,"\n %ld  %ld",j+1,bcn[i][j]);
      }
    }


    gtop = new gtopology;
    gtop->initiate (bcn,nbdofmas,nproc);
    gtop->cnstate=1;
    
    
    for (i=0;i<nproc;i++){
      delete [] bcn[i];
    }
    delete [] bcn;
    
    for (i=0;i<tnbn;i++){
      delete [] cncn[i];
    }
    delete [] cncn;
    
    delete [] ndofncn;

  }
  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnbdof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&maxnbdof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
}
*/

/**
   function solves reduced system of equations from primal
   domain decomposition by conjugate gradient method
   %matrix of reduced system is not assembled
   iteration process is directed by the master processor
   
   @param domproc - domain-processor correspondence
   @param condmat - array containing reduced matrices (Schur complements)
   @param condvect - array containing reduced vectors
   @param out - output file (for auxiliary output)
   
   JK, 2.12.2001
*/
void schurcompl::solve_red_sys_iter (long *domproc,double *condmat,double *condvect,FILE */*out*/)
{
  int ind;
  long i,j,k,sizebuff,gndofe,*cn;
  double alpha,beta,nom,denom,norrhs;
  double *lhs,*rhs,*buff,*p,*d,*r;
  MPI_Status stat;
  rhs=NULL;
  lhs=NULL;
  sizebuff = maxnbdof+1;
  buff = new double [sizebuff];
  memset (buff,0,sizebuff*sizeof(double));
  
  //  right hand side assembling
  if (myrank==0){
    lhs = new double [ndofcp];
    memset (lhs,0,ndofcp*sizeof(double));
    rhs = new double [ndofcp];
    memset (rhs,0,ndofcp*sizeof(double));
    
    gndofe = gtop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    gtop->give_code_numbers (domproc[0],cn);
    locglob (rhs,condvect,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];  ind=0;
      
      gndofe = gtop->give_ndofe (j);
      cn = new long [gndofe];
      gtop->give_code_numbers (j,cn);
      locglob (rhs,buff,cn,gndofe);
      delete [] cn;

    }
  }
  else{
    copyv (condvect,buff,ndof-nidof);
    MPI_Send(buff,sizebuff,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  
  //  iteration loop block
  if (myrank==0){
    r = new double [ndofcp];
    d = new double [ndofcp];
    p = new double [ndofcp];
    
    //  iteration initialization
    norrhs=0.0;  nom=0.0;
    for (i=0;i<ndofcp;i++){
      lhs[i]=0.0;
      r[i]=rhs[i];
      d[i]=rhs[i];
      norrhs+=rhs[i]*rhs[i];
      nom+=r[i]*r[i];
    }
    
    //  iteration loop
    for (i=0;i<nicgsch;i++){
      
      //  direction vector scattering
      for (j=1;j<nproc;j++){
	memset (buff,0,sizebuff*sizeof(double));
	gndofe = gtop->give_ndofe (domproc[j]);
	cn = new long [gndofe];
	gtop->give_code_numbers (domproc[j],cn);
	globloc (d,buff,cn,gndofe);
	delete [] cn;
	
	MPI_Send(buff,sizebuff,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }
      
      memset (buff,0,sizebuff*sizeof(double));
      gndofe = gtop->give_ndofe (domproc[0]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[0],cn);
      globloc (d,buff,cn,gndofe);
      delete [] cn;
      
      //  matrix-vector multiplication
      mxv (condmat,buff,condvect,nbdof,nbdof);
      
      //  matrix-vector multiplication collection
      memset (p,0,ndofcp*sizeof(double));

      gndofe = gtop->give_ndofe (domproc[0]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[0],cn);
      locglob (p,condvect,cn,gndofe);
      delete [] cn;

      for (j=1;j<nproc;j++){
	MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	k=domproc[stat.MPI_TAG];
	
	gndofe = gtop->give_ndofe (k);
	cn = new long [gndofe];
	gtop->give_code_numbers (k,cn);
	locglob (p,buff,cn,gndofe);
	delete [] cn;
	
      }
      
      //  denominator
      denom=ss(d,p,ndofcp);
      
      if (fabs(denom)<zero){
	if (i!=nicgsch-1){
	  buff[maxnbdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send(buff,sizebuff,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	fprintf (stderr,"\n\n denominator in alpha expression in parallel conjugate gradient method is equal to zero.\n");
	break;
      }
      
      //  coefficient alpha
      alpha=nom/denom;
      
      //  new global vectors
      for (j=0;j<ndofcp;j++){
	r[j]-=alpha*p[j];
	lhs[j]+=alpha*d[j];
      }
      
      
      denom=nom;
      nom=ss(r,r,ndofcp);
      
      //if (i%100==0)  fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,nom/norrhs);
      fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,nom/norrhs);
      if (nom/norrhs<errcgsch){
        fprintf (stdout,"\n last iteration %ld    norres/norrhs %e",i,nom/norrhs);
	if (i!=nicgsch-1){
	  buff[maxnbdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send(buff,sizebuff,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      //  beta coefficient
      beta=nom/denom;
      
      //  new direction vector
      for (j=0;j<ndofcp;j++){
	d[j]=r[j]+beta*d[j];
      }
      
    }
    anicgsch=i;
    aerrcgsch=nom/norrhs;
    
    /*
    fprintf (out,"\n\n\n REDUCED SYSTEM SOLUTION (on master)\n");
    for (i=0;i<ndofcp;i++){
      fprintf (out,"\n lhs %ld     %e",i,lhs[i]);
    }
    */
    
  }
  else{
    //  iteration loop
    for (i=0;i<nicgsch;i++){
      //  direction vector receipt
      MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (buff[maxnbdof]>0.5)  break;
      
      //  matrix-vector multiplication
      mxv (condmat,buff,condvect,nbdof,nbdof);
      
      //  matrix-vector multiplication collection
      copyv (condvect,buff,nbdof);
      MPI_Send(buff,sizebuff,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
    }
    anicgsch=i;
  }
  
  //  global solution distribution
  if (myrank==0){
    for (i=1;i<nproc;i++){

      memset (buff,0,sizebuff*sizeof(double));
      gndofe = gtop->give_ndofe (domproc[i]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[i],cn);
      globloc (lhs,buff,cn,gndofe);
      delete [] cn;

      MPI_Send(buff,sizebuff,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }

    nullv (condvect,nbdofmas[domproc[0]]);

    gndofe = gtop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    gtop->give_code_numbers (domproc[0],cn);
    globloc (lhs,buff,cn,gndofe);
    delete [] cn;

  }
  else{
    MPI_Recv(buff,sizebuff,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  copyv (buff,condvect,nbdof);
  
  delete [] buff;
  if (myrank == 0){
    delete [] r;
    delete [] d;
    delete [] p;
    delete [] rhs;
    delete [] lhs;
  }
}

/**
   function solves reduced system of linear algebraic equations on the master
   processor by a direct method
   
   @param domproc - domain-processor correspondence
   @param condmat - array containing reduced matrices (Schur complements)
   @param condvect - array containing reduced vectors
   @param decomp - indicator for valid condensed matrix of reduced system.
                  (0=condensation have to be performed/1=only backward substitution is necessary) 
   @param out - output file (for auxiliary output)
   
   JK
*/
void schurcompl::solve_red_sys_fin (long *domproc,double *condmat,double *condvect,long decomp,FILE *out)
{
  long i,j,k,l,ii,gndofe,*cn;
  double *lhs,*rhs;
  double *buff;
  MPI_Status stat;
  precond prec;
  
  if (decomp==0)
    buff = new double [maxnbdof*maxnbdof+maxnbdof];
  else
    buff = new double [maxnbdof];

  ii=0;
  if (decomp==0)
  {  
    for (i=0;i<nbdof;i++){
      for (j=0;j<nbdof;j++){
        buff[ii]=condmat[i*nbdof+j];  ii++;
      }
    }
  }
  for (i=0;i<nbdof;i++){
    buff[ii]=condvect[i];  ii++;
  }
  
  if (myrank==0){
    lhs = new double [ndofcp];
    memset (lhs,0,ndofcp*sizeof(double));
    rhs = new double [ndofcp];
    memset (rhs,0,ndofcp*sizeof(double));
    if (decomp==0){
      if (arr){
	//        scanf("%*s");
        delete arr;
      }
      arr = new gmatrix;    
      arr->ts = rsmstor;
      arr->setval (ssle);
      arr->alloc (); // uz je v arr->initiate
      arr->initiate (gtop,ndofcp,rsmstor,1,out);    
      
      gndofe = gtop->give_ndofe (domproc[0]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[0],cn);
      arr->localized (condmat,cn,gndofe,domproc[0]);
      locglob (rhs,condvect,cn,gndofe);
      delete [] cn;
      
      
      for (i=1;i<nproc;i++){
	MPI_Recv(buff,maxnbdof*maxnbdof+maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	l=domproc[stat.MPI_TAG];
	
	ii=0;
	for (j=0;j<nbdofmas[l];j++){
	  for (k=0;k<nbdofmas[l];k++){
	    condmat[j*nbdofmas[l]+k]=buff[ii];  ii++;
	  }
	}
	for (j=0;j<nbdofmas[l];j++){
	  condvect[j]=buff[ii];  ii++;
	}
	
	gndofe = gtop->give_ndofe (l);
	cn = new long [gndofe];
	gtop->give_code_numbers (l,cn);
	arr->localized (condmat,cn,gndofe,l);
	locglob (rhs,condvect,cn,gndofe);
	delete [] cn;
	
      }
      
      arr->prepmat (0.0,1);
    }
    else
    {
      gndofe = gtop->give_ndofe (domproc[0]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[0],cn);
      locglob (rhs,condvect,cn,gndofe);
      delete [] cn;
      for (i=1;i<nproc;i++){
        MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        l=domproc[stat.MPI_TAG];
        
        ii=0;
        for (j=0;j<nbdofmas[l];j++){
          condvect[j]=buff[ii];  ii++;
        }
      
        gndofe = gtop->give_ndofe (l);
        cn = new long [gndofe];
        gtop->give_code_numbers (l,cn);
        locglob (rhs,condvect,cn,gndofe);
        delete [] cn;
      }
    }
    
    
    //fprintf (out,"\n\n kontrola matice redukovaneho problemu \n\n");
    //arr->printmat (out);
    //fprintf (out,"\n\n konec kontroly matice redukovaneho problemu \n\n");
    
    //  solution of the reduced problem
    ssle->solve_system (gtop,arr,lhs,rhs,out);
    
    
    //fprintf (out,"\n\n\n REDUCED SYSTEM SOLUTION (on master)\n");
    //for (i=0;i<ndofcp;i++){
    //fprintf (out,"\n lhs %ld     %e",i,lhs[i]);
    //}
    

    for (i=1;i<nproc;i++){
      nullv (condvect,nbdofmas[domproc[i]]);
      
      gndofe = gtop->give_ndofe (domproc[i]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[i],cn);
      globloc (lhs,condvect,cn,gndofe);
      delete [] cn;

      for (long ij=0;ij<nbdofmas[i];ij++){
	buff[ij]=condvect[ij];
      }
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
      
    }
    
    nullv (condvect,nbdofmas[ndom]);

    gndofe = gtop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    gtop->give_code_numbers (domproc[0],cn);
    globloc (lhs,condvect,cn,gndofe);
    delete [] cn;
    
    delete [] lhs;  delete [] rhs;
    
  }
  else{
    if (decomp == 0)
      MPI_Send(buff,maxnbdof*maxnbdof+maxnbdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    else
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);

    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
    for (long ij=0;ij<nbdof;ij++){
      condvect[ij]=buff[ij];
    }
  }
  
  delete [] buff;
  
}

/**
   function solves reduced system of linear algebraic equations
   
   @param domproc - domain-processor correspondence
   @param condmat - array containing reduced matrices (Schur complements)
   @param condvect - array containing reduced vectors
   @param decomp - indicator for valid condensed matrix of reduced system used for elimination of reduced system only.
                  (0=condensation have to be performed/1=only backward substitution is necessary) 
   @param out - output file (for auxiliary output)
   
   JK
*/
void schurcompl::solve_red_sys (long *domproc,double *condmat,double *condvect,long decomp,FILE *out)
{
  switch (trssol){
  case master_sol:{
    solve_red_sys_fin (domproc,condmat,condvect,decomp,out);
    break;
  }
  case paral_sol:{
    solve_red_sys_iter (domproc,condmat,condvect,out);
    break;
  }
  default:{
    par_print_err(myrank,"unknown type of solver of reduced problem is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function solves system of linear algebraic equations
   by the Schur complment method
   
   @param top - pointer to the sequential general topology
   @param gm - pointer to the subdomain %matrix stored in the form gmatrix
   @param domproc - domain-processor correspondence
   @param lhs - array containing solution
   @param rhs - array containing right hand side
   @param out - output file (for auxiliary output)
   
   JK
*/
void schurcompl::solve_system (gtopology *top,gmatrix *gm,
			       long *domproc,double *lhs,double *rhs,FILE *out)
{
  long i,j;
  double *condmat,*condvect;
  time_t t1,t2,t3,t4;
  long decomp;
  
  //  allocation of reduced matrix
  if (myrank==0){
    condmat = new double [maxnbdof*maxnbdof];
    memset (condmat,0,maxnbdof*maxnbdof*sizeof(double));
  }
  else{
    condmat = new double [nbdof*nbdof];
    memset (condmat,0,nbdof*nbdof*sizeof(double));
  }
  
  //  allocation of reduced vector
  if (myrank==0){
    condvect = new double [maxnbdof];
    memset (condvect,0,maxnbdof*sizeof(double));
  }
  else{
    condvect = new double [nbdof];
    memset (condvect,0,nbdof*sizeof(double));
  }
  
  t1 = time (NULL);

  decomp = gm->decomp();
  //  elimination of inner DOFs = static condensation of inner DOFs
  //  the Schur complement matrix is stored in array condmat
  gm->condense (top,condmat,condvect,lhs,rhs,nbdof,1,out);
  
  /*
  if (gm->dsky){
    fprintf (out,"\n\n\n\n Kontrola pole dsky");
    for (i=0;i<gm->dsky->n;i++){
      for (j=gm->dsky->adr[i];j<gm->dsky->adr[i+1];j++){
	fprintf (out,"\n %13.10le",gm->dsky->a[j]);
      }
    }
    fprintf (out,"\n\n\n\n polovina");
    long ijkolp=gm->dsky->n;
    for (i=0;i<gm->dsky->n;i++){
      for (j=gm->dsky->adr[i]+gm->dsky->adr[ijkolp];j<gm->dsky->adr[i+1]+gm->dsky->adr[ijkolp];j++){
	fprintf (out,"\n %13.10le",gm->dsky->a[j]);
      }
    }
  }
  if (gm->sky){
    fprintf (out,"\n\n\n\n Kontrola pole sky");
    for (i=0;i<gm->sky->n;i++){
      for (j=gm->sky->adr[i];j<gm->sky->adr[i+1];j++){
	fprintf (out,"\n %13.10le",gm->sky->a[j]);
      }
    }
  }
  */
  
  /*
  fprintf (out,"\n\n\n kontrola kondenzace \n");
  //fprintf (stdout,"\n\n\n kontrola kondenzace \n");
  for (i=0;i<nbdof;i++){
    //fprintf (stdout,"\n condvect %le   ",condvect[i]);
    fprintf (out,"\n condvect %le   ",condvect[i]);
    for (j=0;j<nbdof;j++){
      //fprintf (stdout," %le",condmat[i*nbdof+j]);
      fprintf (out," %le",condmat[i*nbdof+j]);
    }
  }
  //fprintf (stdout,"\n\n konec kontroly \n\n");
  fprintf (out,"\n\n konec kontroly \n\n");
  */

  t2 = time (NULL);
  

  //  solution of reduced problem
  solve_red_sys (domproc,condmat,condvect,decomp,out);


  fprintf (out,"\n\n\n kontrola po hrubem problemu \n");
  //fprintf (stdout,"\n\n\n kontrola kondenzace \n");
  for (i=0;i<nbdof;i++){
    //fprintf (stdout,"\n condvect %le   ",condvect[i]);
    fprintf (out,"\n condvect %le   ",condvect[i]);
  }
  //fprintf (out,"\n\n konec kontroly \n\n");
  //fprintf (stdout,"\n\n konec kontroly \n\n");


  t3 = time (NULL);
  
  //  back substitution on subdomains
  gm->condense (top,condmat,condvect,lhs,rhs,nbdof,2,out);
  
  //for (long ijk=0;ijk<ndof;ijk++){
  //fprintf (stdout,"\n lhs %ld   %e",ijk,lhs[ijk]);
  //}

  t4 = time (NULL);
  
  // ****************************************************
  //  gathering of elapsed time on the master processor
  // ****************************************************
  MPI_Status stat;
  long *buff;
  buff = new long[3];
  buff[0]=t2-t1;
  buff[1]=t3-t2;
  buff[2]=t4-t3;
  

  if (myrank==0){
    
    long *fact,*cp,*tottime;
    fact=new long [nproc];
    cp=new long [nproc];
    tottime = new long [nproc];
    
    fprintf (out,"\n\n\n\n\n");
    j=domproc[myrank];
    fact[j]=buff[0];
    cp[j]=buff[1];
    tottime[j]=buff[0]+buff[1]+buff[2];
    
    /*
    fprintf (out,"\n\n\n Domain %ld",j);
    fprintf (out,"\n time of condensation             %ld",buff[0]);
    fprintf (out,"\n time of solution of red. sys.    %ld",buff[1]);
    fprintf (out,"\n time of back substitution        %ld",buff[2]);
    */

    j=domproc[myrank];
//    fprintf (stdout,"\n\n\n Domain %ld",j);
//    fprintf (stdout,"\n time of condensation             %ld",buff[0]);
//    fprintf (stdout,"\n time of solution of red. sys.    %ld",buff[1]);
//    fprintf (stdout,"\n time of back substitution        %ld",buff[2]);

    for (i=1;i<nproc;i++){
      MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      fact[j]=buff[0];
      cp[j]=buff[1];
      tottime[j]=buff[0]+buff[1]+buff[2];
      
      /*
      fprintf (out,"\n\n\n Domain %ld",j);
      fprintf (out,"\n time of condensation             %ld",buff[0]);
      fprintf (out,"\n time of solution of red. sys.    %ld",buff[1]);
      fprintf (out,"\n time of back substitution        %ld",buff[2]);
      */
//      fprintf (stdout,"\n\n\n Domain %ld",j);
//      fprintf (stdout,"\n time of condensation             %ld",buff[0]);
//      fprintf (stdout,"\n time of solution of red. sys.    %ld",buff[1]);
//      fprintf (stdout,"\n time of back substitution        %ld",buff[2]);
    }
    
    long maxfact,mincp,mintot,maxtot;
    maxfact=0;
    maxtot=0;
    mincp=LONG_MAX;
    mintot=LONG_MAX;
    for (i=0;i<nproc;i++){
      if (maxfact<fact[i])
	maxfact=fact[i];

      if (mincp>cp[i])
	mincp=cp[i];

      if (mintot>tottime[i])
	mintot=tottime[i];
      if (maxtot<tottime[i])
	maxtot=tottime[i];
    }
    /*
    fprintf (out,"\n\n maximum factorization time  %ld",maxfact);
    fprintf (out,"\n coarse problem time         %ld",mincp);
    fprintf (out,"\n minimum total time          %ld",mintot);
    fprintf (out,"\n maximum total time          %ld",maxtot);
    */
    delete [] fact;
    delete [] cp;
    delete [] tottime;
  }
  else{
    MPI_Send(buff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  delete [] buff;
  
  
  delete [] condmat;
  delete [] condvect;

}

/**
   function gathers boundary/interface contributions from subdomains into one global %vector
   the global %vector must be deleted outside of the subroutine
   (the global %vector contains all components defined on the boundaries/interfaces)
   (the local %vector contains all components defined on a subdomain)
   
   @param lv - local %vector
   @param gv - global %vector, it is allocated only on the master processor
               the array must be deleted outside of the subroutine
   @param domproc - domain-processor correspondence
   
   JK, 6.11.2004
*/
void schurcompl::gather_bound_vect (double *lv,double *gv,long *domproc)
{
  long i,j,l,gndofe,*cn;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxnbdof];
  
  //  boundary/interface values from subdomain are selected
  j=0;
  for (i=nidof;i<ndof;i++){
    buff[j]=lv[i];  j++;
  }
  
  if (myrank==0){
    // ******************************
    //  contribution from the master
    // ******************************
    l=domproc[0];
    //  number of DOFs on appropriate subdomain
    gndofe = gtop->give_ndofe (l);
    cn = new long [gndofe];
    //  code numbers of subdomain
    gtop->give_code_numbers (l,cn);
    //  localization of values from local vector to the global vector
    locglob (gv,buff,cn,gndofe);
    delete [] cn;
    
    // ***************************
    //  contributions from slaves
    // ***************************
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=domproc[stat.MPI_TAG];
      
      //  number of DOFs on appropriate subdomain
      gndofe = gtop->give_ndofe (l);
      cn = new long [gndofe];
      //  code numbers of subdomain
      gtop->give_code_numbers (l,cn);
      //  localization of values from local vector to the global vector
      locglob (gv,buff,cn,gndofe);
      delete [] cn;
    }
  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
}

/**
   function computes dot product of two vectors which are split on several processors
   it gathers boundary/interface contributions from subdomains into global vectors
   it gathers results of dot product of internal components of the vectors
   the global vectors must be deleted outside of the subroutine
   
   @param lv1,lv2 - local vectors
   @param gv1,gv2 - global vectors, they are allocated only on the master processor
                    the arrays must be deleted outside of the subroutine
   
   JK, 6.11.2004, revised 2.6.2011
*/
double schurcompl::pss_gather_bound_vect (double *lv1,double *gv1,double *lv2,double *gv2,long *domproc)
{
  long i,j,l,gndofe,*cn,buffsize;
  double dp,idp,*buff;
  MPI_Status stat;
  
  //  contributions from the vector lv1 are stored on positions with id between 0 and maxnbdof,
  //  contributions from the vector lv2 are stored on positions with id between maxnbdof and 2.maxnbdof
  buffsize=2*maxnbdof+1;
  buff = new double [buffsize];
  
  idp=0.0;
  
  //  dot product of internal DOFs
  dp=0.0;
  for (i=0;i<nidof;i++){
    dp+=lv1[i]*lv2[i];
  }
  buff[buffsize-1]=dp;
  
  //  interface values from subdomain are selected
  j=0;
  for (i=nidof;i<ndof;i++){
    buff[j]=lv1[i];
    buff[j+maxnbdof]=lv2[i];  j++;
  }
  
  if (myrank==0){
    //  sum of contributions from internal DOFs
    idp=buff[buffsize-1];
    
    // ******************************
    //  contribution from the master
    // ******************************
    l=domproc[0];
    //  number of DOFs on appropriate subdomain
    gndofe = gtop->give_ndofe (l);
    cn = new long [gndofe];
    //  code numbers of subdomain
    gtop->give_code_numbers (l,cn);
    //  localization of values from local vectors to the global vectors
    locglob (gv1,buff,cn,gndofe);
    locglob (gv2,buff+maxnbdof,cn,gndofe);
    delete [] cn;
    
    // ***************************
    //  contributions from slaves
    // ***************************
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=domproc[stat.MPI_TAG];
      
      //  number of DOFs on appropriate subdomain
      gndofe = gtop->give_ndofe (l);
      cn = new long [gndofe];
      //  code numbers of subdomain
      gtop->give_code_numbers (l,cn);
      //  localization of values from local vector to the global vector
      locglob (gv1,buff,cn,gndofe);
      locglob (gv2,buff+maxnbdof,cn,gndofe);
      delete [] cn;
      
      //  sum of contributions from internal DOFs
      idp+=buff[buffsize-1];
    }
  }
  else{
    MPI_Send(buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;

  return idp;
}


/**
   function computes norm of a %vector where only selected components contribute
   the vector is split on several processors
   
   function assumes that no boundary/interface node is equipped with the Dirichlet boundary condition
   function assumes that no boundary/interface node is coupled with any other node
   
   @param cid[in] - component id (id of required component)
   @param err[in] - required residual
   @param thresh[in] - zero norm threshold for the given cid-th component (quantity)
   @param sc[in] - flag for scaling norm of vectors lv1 by the norm of rhs (sc=1)
   @param ptop[in] - pointer to the parallel topology
   @param top[in] - pointer to the subdomain topology
   @param lv1,rhs[in] - local vectors (e.g. residual and rhs) on subdomains
   @param gv1,grhs[in/out] - global vectors, they are allocated only on the master processor
                             the arrays must be deleted outside of the subroutine, 
                             dimension must be ndofcp, i.e. number of DOFs of coarse problem
   @param domproc - domain-processor correspondence
   @param n[in] - dimension of vectors lv1 and rhs
   @param out[in] - pointer to the opened text file for logging messages
   @param norfv[out] - resulting norm of lv %vector of whole domain
   
   if the vectors gv1 and grhs are not used outside of this function, they will be removed ???!!!
   
   JK, 13.7.2011
*/
long schurcompl::selected_norm_calculation (long cid, double err, double thresh, long sc,
                                            partopjk *ptop, gtopology *top,
                                            double *lv1, double *gv1, double *rhs, double *grhs,
                                            long *domproc, long /*n*/, FILE *out, double &norfv)
{
  long i,j,k,l,nid,dof,gndofe,*cn,buffsize,stop,ndofn;
  double dp,totdp,norrhs,totnorrhs,*buff;
  MPI_Status stat;
  
  //for (i=0;i<n;i++){
  //fprintf (stdout,"\n lv1 %ld   %le",i+1,lv1[i]);
  //}
  //for (i=0;i<n;i++){
  //fprintf (stdout,"\n rhs %ld   %le",i+1,rhs[i]);
  //}
  
  
  //  buff[0--maxnbdof-1] contains contributions of the vector lv1
  //  buff[maxnbdof--2maxnbdof-1] contains contributions of the vector rhs
  //  buff[2maxnbdof] contains partial dot product of lv1*lv1
  //  buff[2maxnbdof+1] contains partial dot product of rhs*rhs
  buffsize=2*maxnbdof+2;
  buff = new double [buffsize];
  nullv (buff,buffsize);
  
  //  dot product of selected internal DOFs
  dp=0.0;
  //  norm of the right hand side vector
  norrhs=0.0;
  //  loop over the internal nodes
  //printf ("ptop->nin %ld",ptop->nin);
  for (i=0;i<ptop->nin;i++){
    //  node id
    nid=ptop->inid[i];
    //printf ("\n nid(int) %ld",nid);
    //  id of selected DOF
    dof = top->give_dof (nid,cid);
    //printf ("dof %ld",dof);
    if (dof>0){
      dof--;
      dp+=lv1[dof]*lv1[dof];
      norrhs+=rhs[dof]*rhs[dof];
    }

  }

  buff[buffsize-2]=dp;
  buff[buffsize-1]=norrhs;

  //  boundary/interface values from subdomain are selected
  k=0;
  //  loop over the boundary/interface nodes
  for (i=0;i<ptop->nbn;i++){
    //  node id
    nid=ptop->bnid[i];
    //fprintf (stdout,"\n nid(bound) %ld",nid);
    //  id of selected DOF
    
    ndofn = top->give_ndofn (nid);
    for (j=0;j<ndofn;j++){
      dof = top->give_dof (nid,j);
      if (dof>0){
	if (cid==j){
	  dof--;
	  buff[k]=lv1[dof];
	  buff[maxnbdof+k]=rhs[dof];
	}
	k++;
      }
    }
  }

  //for (i=0;i<maxnbdof;i++){
  //fprintf (stdout,"\n buff %ld   %le",i+1,buff[i]);
  //}
  //for (i=maxnbdof;i<2*maxnbdof;i++){
  //fprintf (stdout,"\n buff %ld   %le",i+1,buff[i]);
  //}
  //fprintf (stdout,"\n\n\n\n\n\n\n\n\n\n\n\n\n\n");

  if (myrank==0){
    //  sum of contributions from partial dot products of internal DOFs
    totdp=buff[buffsize-2];
    totnorrhs=buff[buffsize-1];
 
    // ******************************
    //  contribution from the master
    // ******************************
    l=domproc[0];
    //  number of DOFs on appropriate subdomain
    gndofe = gtop->give_ndofe (l);
    cn = new long [gndofe];
    //  code numbers of subdomain
    gtop->give_code_numbers (l,cn);

    fprintf (out,"\n\n domain %d\n kontrola kodovych cisel na hranici\n",myrank);
    for (long ijk=0;ijk<gndofe;ijk++){
      fprintf (out," %ld",cn[ijk]);
    }
    fprintf (out,"\n\n");
    
    
    //  localization of values from local vectors to the global vectors
    locglob (gv1,buff,cn,gndofe);
    locglob (grhs,buff+maxnbdof,cn,gndofe);
    delete [] cn;

    //for (j=0;j<ndofcp;j++){
    //printf ("\n gv1 %6ld   %le    grhs %6ld   %le",j,gv1[j],j,grhs[j]);
    //}
    
    // ***************************
    //  contributions from slaves
    // ***************************
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //  id of subdomain
      l=domproc[stat.MPI_TAG];
      
      //  number of DOFs on appropriate subdomain
      gndofe = gtop->give_ndofe (l);
      cn = new long [gndofe];
      //  code numbers of subdomain
      gtop->give_code_numbers (l,cn);
      
      fprintf (out,"\n\n domain %ld\n kontrola kodovych cisel na hranici\n",l);
      for (long ijk=0;ijk<gndofe;ijk++){
	fprintf (out," %ld",cn[ijk]);
      }
      fprintf (out,"\n\n");
      
      //  localization of values from local vectors to the global vectors
      locglob (gv1,buff,cn,gndofe);
      locglob (grhs,buff+maxnbdof,cn,gndofe);
      delete [] cn;

      //for (j=0;j<ndofcp;j++){
      //printf ("\n gv1 %6ld   %le    grhs %6ld   %le",j,gv1[j],j,grhs[j]);
      //}
      
      //  sum of contributions from partial dot products of internal DOFs
      totdp+=buff[buffsize-2];
      totnorrhs+=buff[buffsize-1];
    }
    
  }
  else{
    MPI_Send(buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  delete [] buff;
  
  //  equilibriated state is not reached
  stop=0;
  if (myrank==0){
    
    fprintf (out,"\n Kontrola gv na masterovi");
    for (i=0;i<ndofcp;i++){
      fprintf (out,"\n gv1 %6ld   % 12.9le    grhs  % 12.9le",i,gv1[i],grhs[i]);
      totdp+=gv1[i]*gv1[i];
      totnorrhs+=grhs[i]*grhs[i];
    }
    

    totdp=sqrt(totdp);
    totnorrhs=sqrt(totnorrhs);

    /*
    fprintf (stdout,"\n\n\n\n component ID %ld",cid);
    fprintf (stdout,"\n totdp %le",totdp);
    fprintf (stdout,"\n totnorrhs %le",totnorrhs);
    fprintf (stdout,"\n thresh %le",thresh);
    fprintf (stdout,"\n err %le",err);
    fprintf (stdout,"\n\n\n\n\n");
    */    
    
    norfv = totdp;
    if ((totnorrhs<thresh) || (sc == 0)){
      if (sc == 0)
        fprintf (stdout,"\n phase %ld,  norm of residuum  %e,    required error  %e", cid+1, norfv, err);
      else
        fprintf (stdout,"\n phase %ld,  norm of residuum  %e,    norm of rhs  %e,  required error  %e", cid+1, totdp, totnorrhs, err);
      if (totdp<err){
	//  equilibriated state is reached
	stop = 1;        
      }
    }
    else{
      norfv = totdp/totnorrhs;
      fprintf (stdout,"\n phase %ld,   norm of residuum  %e,    norm of rhs  %e,   error  %e,   required error  %e",cid+1, totdp, totnorrhs, norfv, err);
      if (norfv<err){
	//  equilibriated state is reached
	stop=1;
      }
    }

    for (i=1;i<nproc;i++){ 
      MPI_Send(&stop,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      MPI_Send(&norfv,1,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
    
  }else{
    MPI_Recv(&stop,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    MPI_Recv(&norfv,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  return stop;
}


/**
   function gathers contributions from subdomains into one global %vector
   the global %vector must be deleted outside of the subroutine
   
   @param lv - local %vector
   @param gv - global %vector, it is allocated only on the master processor
               the array must be deleted outside of the subroutine
   
   JK, 6.11.2004
*/
double schurcompl::pss_gather_bound_vect_old (double *lv,double *gv,long *domproc)
{
  long i,j,l,gndofe,*cn,buffsize;
  double dp,idp,*buff;
  MPI_Status stat;
  
  buffsize=maxnbdof+1;
  buff = new double [buffsize];
  
  idp=0.0;
  
  //  dot product of internal DOFs
  dp=0.0;
  for (i=0;i<nidof;i++){
    dp+=lv[i]*lv[i];
  }
  buff[maxnbdof]=dp;
  
  //  interface values from subdomain are selected
  j=0;
  for (i=nidof;i<ndof;i++){
    buff[j]=lv[i];  j++;
  }
  
  if (myrank==0){
    //  sum of contributions from internal DOFs
    idp=buff[maxnbdof];
    
    // ******************************
    //  contribution from the master
    // ******************************
    l=domproc[0];
    //  number of DOFs on appropriate subdomain
    gndofe = gtop->give_ndofe (l);
    cn = new long [gndofe];
    //  code numbers of subdomain
    gtop->give_code_numbers (l,cn);
    //  localization of values from local vector to the global vector
    locglob (gv,buff,cn,gndofe);
    delete [] cn;
    
    // ***************************
    //  contributions from slaves
    // ***************************
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=domproc[stat.MPI_TAG];
      
      //  number of DOFs on appropriate subdomain
      gndofe = gtop->give_ndofe (l);
      cn = new long [gndofe];
      //  code numbers of subdomain
      gtop->give_code_numbers (l,cn);
      //  localization of values from local vector to the global vector
      locglob (gv,buff,cn,gndofe);
      delete [] cn;
      
      //  sum of contributions from internal DOFs
      idp+=buff[maxnbdof];
    }
  }
  else{
    MPI_Send(buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;

  return idp;
}



/**
   function scatters contributions from one global vector into local vectors
   function rewrites part of the local vector by contributions from boundaries
   
   @param lv - local vector
   @param gv - global vector, it is allocated only on the master processor
   
   JK, 6.11.2004
*/
void schurcompl::scatter_bound_vect (double *lv,double *gv,long *domproc)
{
  long i,j,gndofe,*cn;
  double *buff;
  MPI_Status stat;

  buff = new double [maxnbdof];

  if (myrank==0){
    for (i=1;i<nproc;i++){
      nullv (buff,nbdofmas[domproc[i]]);
      
      gndofe = gtop->give_ndofe (domproc[i]);
      cn = new long [gndofe];
      gtop->give_code_numbers (domproc[i],cn);
      globloc (gv,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
    
    nullv (buff,nbdofmas[domproc[0]]);
    
    gndofe = gtop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    gtop->give_code_numbers (domproc[0],cn);
    globloc (gv,buff,cn,gndofe);
    delete [] cn;
    
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  j=0;
  for (i=nidof;i<ndof;i++){
    lv[i]=buff[j]/2.0;  j++; // ???!!! proc 1/2, 09.2021 TKo
  }
  
  delete [] buff;
  
}

/**
   function computes norm of vector of unbalanced values
   function gathers vector of unbalanced values on boundaries
   function scatters correct values into subdomains
   function returns norm of vector of unbalanced values
   
   @param lv - local vector
   @param domproc - domain-processor correspondence
   @param out - output file
   
   JK, 6.11.2004
*/
double schurcompl::unbalanced_values (double *lv,long *domproc,FILE */*out*/)
{
  long i;
  double aux,norm,parlocnorm;
  double *gv=NULL;
  MPI_Status stat;
  
  if (myrank==0){
    gv = new double [ndofcp];
    nullv (gv,ndofcp);
  }
  
  //  norm of vector of unbalanced values at internal nodes
  parlocnorm =0;
  for (i=0;i<nidof;i++){
    parlocnorm+=lv[i]*lv[i];
  }
  
  //  gathering of partial norms from subdomains
  if (myrank==0){
    norm=parlocnorm;
    for (i=1;i<nproc;i++){
      MPI_Recv (&aux,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      norm+=aux;
    }
  }
  else{
    MPI_Send(&parlocnorm,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  gathering of unbalanced values at boundary nodes
  gather_bound_vect (lv,gv,domproc);
  
  //  contribution to the norm from boundary nodes
  if (myrank==0){
    for (i=0;i<ndofcp;i++){
      norm+=gv[i]*gv[i];
      //fprintf (out,"\n%6ld  %e",i,gv[i]);
    }
  }
  
  //  scattering of unbalanced values at boundary nodes into subdomains
  // ???!!! proc scatter, pokazi to vektor lv, 09.2021
  scatter_bound_vect (lv,gv,domproc);
  
  //  scattering of norm of vector of unbalanced values
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Send(&norm,1,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&norm,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (myrank==0)
    delete [] gv;
  
  return norm;
}

