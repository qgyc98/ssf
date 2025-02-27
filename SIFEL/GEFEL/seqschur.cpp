#include "seqschur.h"
#include "slesolv.h"
#include "gmatrix.h"
#include "skyline.h"
#include <limits.h>
#include <string.h>
#include <math.h>
#include <time.h>

/**
   constructor
   
   JK, 13.5.2009
*/
seqschur::seqschur ()
{
  //  number of subdomains
  ns=0;
  
  //  type of solver of reduced/coarse problem
  trssol = (redsystsolver) 0;
  //  type of storage of matrix of reduced/coarse problem
  rsmstor = (storagetype) 0;
  //  maximum number of iterations in the conjugate gradient method used for solution of reduced problem
  nicg=0;
  //  required/prescribed norm of residuum
  errcg=0.0;
  //  actual number of iterations (number of iterations performed during the solution)
  anicg=0;
  //  attained norm of residuum
  aerrcg=0.0;

  //  computer zero
  zero=0.0;







  //  number of DOFs on subdomain
  ndof=0;
  //  number of internal DOFs on subdomain
  nidof=0;
  //  number of boundary/interface DOFs on subdomain
  nbdof=0;
  //  maximum number of boundary/interface DOFs on subdomain
  maxnbdof=0;


  
  
  //  array containing numbers of boundary/interface DOFs on subdomains
  nbdofmas=NULL;
  
  
  gtop=NULL;
  arr=NULL;
  
  //  solver of system of linear eqautions
  ssle = NULL;
  
}

seqschur::~seqschur()
{
  delete ssle;
}


/**
   function reads basic data
   
   @param i - type of solver of system of equations
   @param in - input file
   
   JK, 13.5.2009
*/
void seqschur::read (gtopology *top,long mespr,XFILE *in)
{
  //  number of subdomains
  xfscanf (in,"%ld",&ns);
  
  //  type of reduced solver
  //  trssol = master_sol = 1 - solution by a direct method
  //  trssol = paral_sol = 2 - solution by an iterative method
  xfscanf (in,"%d",(int*)&trssol);
  switch (trssol){
  case master_sol:{
    if (mespr==1)  fprintf (stdout,"\n reduced system of equations is solved by direct method");
    xfscanf (in,"%d",(int*)&rsmstor);
    
    if (ssle==NULL)
      ssle = new slesolv ();

    ssle->read (top,in,mespr);
    
    break;	
  }
  case paral_sol:{
    if (mespr==1)  fprintf (stdout,"\n reduced system of equations is solved by the iterative method");
    xfscanf (in,"%ld %le",&nicg,&errcg);
    break;
  }
  default:{
    print_err("unknown reduced system solver is required",__FILE__,__LINE__,__func__);
  }
  }

}

/**
   function prints basic data
   
   @param out - output file
   
   JK, 13.5.2009
*/
void seqschur::print (FILE *out)
{
  //  number of subdomains
  fprintf (out,"%ld\n",ns);
  //  type of reduced solver
  fprintf (out,"%d\n",trssol);
  
  switch (trssol){
  case master_sol:{
    fprintf (out,"%d",rsmstor);
    ssle->print (out);
    break;	
  }
  case paral_sol:{
    fprintf (out,"%ld %le\n",nicg,errcg);
    break;
  }
  default:{
    print_err("unknown reduced system solver is required",__FILE__,__LINE__,__func__);
  }
  }

}




/**
   function generates ordering of unknowns on subdomains
   
   @param top - pointer to the sequential general topology
   @param ltg - local to global array
   @param out - output stream
   
   JK, 13.5.2009
*/
/*
void seqschur::subdomain_ordering (gtopology *top,long *ltg,FILE *out)
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
  // array with indices of boundary nodes
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

  fprintf (out,"\n\n\n seqschur::subdomain_ordering \n\n");
  for (i=0;i<top->nn;i++){
    fprintf (out,"\n node %6ld ",i);
    ndofn=top->give_ndofn (i);
    for (j=0;j<ndofn;j++){
      fprintf (out," %ld",top->give_dof(i,j));
    }
  }
  fprintf (out,"\n\n");
}
*/



/**
   function orders unknowns on subdomains with respect of Schur complement method
   inner variables are ordered at the beginning
   boundary variables are ordered at he end
   
   @param top - pointer to domain topology
   @param ptop - pointer to topology for paralle computation
   
   JK, 13.5.2009
*/
/*
void seqschur::subdomain_ordering (gtopology *top,partop *ptop)
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
   function initializes variables and arrays

   namely:
   
   ndofcp - number of coarse DOFs
   maxnbdof - maximum number of boundary/interface DOFs on subdomain
   nbdofmas - array of numbers of boundary/interface DOFs on subdomains
   
   @param selnodschur - pointer to the object of the class seqselnodes
   
   JK, 13.5.2009
*/
void seqschur::initiate (seqselnodes *selnodschur)
{
  long i;

  //  number of coarse DOFs = total number of boundary/interface DOFs
  //  number of DOFs of the coarse problem
  ndofcp = selnodschur->tndofsn;
  
  //  number of boundary/interface DOFs on subdomains
  nbdofmas = new long [ns];

  //  array containing numbers of boundary/interface DOFs on subdomains
  //  maximum number of reduced DOFs on subdomain
  maxnbdof=0;
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    nbdofmas[i] = selnodschur->snndofmas[i];
    if (maxnbdof<nbdofmas[i]){
      maxnbdof=nbdofmas[i];
    }
  }
  
  //  topology describing subdomains as superelements
  gtop = new gtopology;
  gtop->initiate (selnodschur->cndofmas,selnodschur->snndofmas,ns);
  gtop->cnstate=1;
  
}

/**
   function orders unknowns in coarse problem
   
   function initializates object gtop of the class gtopology
   it represents subdomains as superelements
   
   @param out - output file
   
   JK, 8.7.2005
*/
  /*
void seqschur::coarse_problem_ordering (FILE *out)
{

  long i,j,k,m,tnbn;
  long *nbnd,*ndofncn;
  long **cnbn,**nbdofnd,**cncn,**bcn;
  long ***lbcn;
  
  //  array of numbers of boundary nodes on subdomains
  nbnd = gtop->stop->nbnd;
  //  array containing coarse node numbers of boundary nodes
  cnbn = gtop->stop->lgnbn;
  //  array containing numbers of DOFs on boundary nodes
  nbdofnd = gtop->stop->lbndofn;
  //  total number of boundary nodes
  tnbn = gtop->stop->tnbn;
  //  array containing local boundary code numbers
  lbcn = gtop->stop->lbcn;
  
  
  //  number of DOFs at coarse nodes
  ndofncn = new long [tnbn];
  
  for (i=0;i<ns;i++){
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
  nbdofmas = new long [ns];
  
  maxnbdof=0;
  for (i=0;i<ns;i++){
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
  
  
  bcn = new long* [ns];
  for (i=0;i<ns;i++){
    bcn[i] = new long [nbdofmas[i]];
    nbdofmas[i]=0;
  }
  
  //  assembling of lists of boundary code numbers on subdomains
  for (i=0;i<ns;i++){
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
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %ld",i+1);
    for (j=0;j<nbdofmas[i];j++){
      fprintf (out,"\n %ld  %ld",j+1,bcn[i][j]);
    }
  }
  
  // nasledujici cast jiz byla jednou volana
  //gtop = new gtopology;
  //gtop->initiate (bcn,nbdofmas,ns);
  //gtop->cnstate=1;


  
  for (i=0;i<ns;i++){
    delete [] bcn[i];
  }
  delete [] bcn;
  
  for (i=0;i<tnbn;i++){
    delete [] cncn[i];
  }
  delete [] cncn;
  
  delete [] ndofncn;

}
  */

/**
   function assembles list of unknowns belonging to the subdomains
   
   ndofdom - numbers of DOFs on subdomains
   cndom - list of DOFs on subdomains
   
   @param top - pointer to topology
   @param out - output file
   
   JK, 21.5.2009
*/
void seqschur::assemble_subdom_unknowns (gtopology *top,FILE *out)
{
  long i,j,l,k,g,ndofn,ii;
  long maxdof;
  long **aux;
  
  aux = new long* [ns];
  
  //  number of DOFs on subdomains
  ndofdom = new long [ns];
  for (i=0;i<ns;i++){
    ndofdom[i]=0;
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  maximum DOF id
    maxdof=0;
    
    //  loop over the number of internal nodes on subdomains
    for (j=0;j<top->stop->nind[i];j++){
      //  global glued number of the node
      g=top->stop->ggnin[i][j];
      ndofn=top->give_ndofn (g);
      //  loop over the number of DOFs
      for (k=0;k<ndofn;k++){
	l=top->give_dof (g,k);
	if (l>maxdof)
	  maxdof=l;
      }
    }
    //  loop over the number of boundary/interface nodes on subdomains
    for (j=0;j<top->stop->nbnd[i];j++){
      //  global glued number of the node
      g=top->stop->ggnbn[i][j];
      ndofn=top->give_ndofn (g);
      //  loop over the number of DOFs
      for (k=0;k<ndofn;k++){
	l=top->give_dof (g,k);
	if (l>maxdof)
	  maxdof=l;
      }
    }
    

    aux[i] = new long [maxdof];
    for (j=0;j<maxdof;j++){
      aux[i][j]=-1;
    }
    
    ii=0;
    //  loop over the number of internal nodes on subdomains
    for (j=0;j<top->stop->nind[i];j++){
      //  global glued number of the node
      g=top->stop->ggnin[i][j];
      ndofn=top->give_ndofn (g);
      //  loop over the number of DOFs
      for (k=0;k<ndofn;k++){
	l=top->give_dof (g,k);
	if (l>0){
	  if (aux[i][l-1]==-1){
	    aux[i][l-1]=ii;
	    ii++;
	  }
	}
      }
    }
    //  loop over the number of boundary/interface nodes on subdomains
    for (j=0;j<top->stop->nbnd[i];j++){
      //  global glued number of the node
      g=top->stop->ggnbn[i][j];
      ndofn=top->give_ndofn (g);
      //  loop over the number of DOFs
      for (k=0;k<ndofn;k++){
	l=top->give_dof (g,k);
	if (l>0){
	  if (aux[i][l-1]==-1){
	    aux[i][l-1]=ii;
	    ii++;
	  }
	}
      }
    }
    
    ndofdom[i]=0;
    for (j=0;j<maxdof;j++){
      if (aux[i][j]>-1){
	ndofdom[i]++;
      }
    }
  }
  
  //  numbers of DOFs on subdomains are known now
  
  //  array of code numbers on subdomains
  cndom = new long* [ns];
  for (i=0;i<ns;i++){
    cndom[i] = new long [ndofdom[i]];
    //ndofdom[i]=0;
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of internal nodes on subdomains
    for (j=0;j<top->stop->nind[i];j++){
      //  global glued number of the node
      g=top->stop->ggnin[i][j];
      ndofn=top->give_ndofn (g);
      //  loop over the number of DOFs
      for (k=0;k<ndofn;k++){
	l=top->give_dof (g,k);
	if (l>0){
	  cndom[i][aux[i][l-1]]=l;
	  //ndofdom[i]++;
	}
      }
    }
    //  loop over the number of boundary/interface nodes on subdomains
    for (j=0;j<top->stop->nbnd[i];j++){
      //  global glued number of the node
      g=top->stop->ggnbn[i][j];
      ndofn=top->give_ndofn (g);
      //  loop over the number of DOFs
      for (k=0;k<ndofn;k++){
	l=top->give_dof (g,k);
	if (l>0){
	  cndom[i][aux[i][l-1]]=l;
	  //ndofdom[i]++;
	}
      }
    }
  }

  
  for (i=0;i<ns;i++){
    delete [] aux[i];
  }
  delete [] aux;


  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n array ndofdom - number of DOFs on subdomains \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain %3ld    %6ld",i,ndofdom[i]);
  }
  
  fprintf (out,"\n\n\n array cndom - array of code numbers of DOFs on subdomains \n\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n subdomain %ld",i);
    for (j=0;j<ndofdom[i];j++){
      fprintf (out,"\n DOF %6ld    %ld",j,cndom[i][j]);
    }
  }

}


/**
   function assembles subdomain matrices from the %matrix of the whole system
   
   @param top - general topology
   @param gm - pointer to the %matrix of the system
   @param out - output file
   
   JK, 13.5.2009
*/
void seqschur::subdomain_matrices (gmatrix *gm,FILE *out)
{
  long i,j;
  /*
  smsky = new skyline [ns];
  for (i=0;i<ns;i++){
    gndofe = gtop->give_ndofe (i);
    cn = new long [gndofe];
    gtop->give_code_numbers (i,cn);
    smsky[i].assemble_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,nbdofmas[i],cn);
    delete [] cn;
  }
  */
  smsky = new skyline [ns];
  for (i=0;i<ns;i++){
    smsky[i].assemble_from_scr (gm->scr->adr,gm->scr->ci,gm->scr->a,ndofdom[i],cndom[i]);
  }
  
  // *******************
  //  auxiliary output
  // *******************
  
  
  fprintf (out,"\n\n kontrola adres ve skyline \n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n podoblast %ld",i);
    for (j=0;j<smsky[i].n+1;j++){
      fprintf (out,"\n adr %4ld   %ld",j,smsky[i].adr[j]);
    }
  }
  
  fprintf (out,"\n\n kontrola ndofdom\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n pocet DOFs na domene %ld   %ld\n",i,ndofdom[i]);
    for (j=0;j<ndofdom[i];j++){
      fprintf (out,"  %ld",cndom[i][j]);
    }
  }
  
  fprintf (out,"\n\n\n");
  for (i=0;i<smsky[1].n;i++){
    fprintf (out,"\n\n%ld\n",i);
    for (j=0;j<i;j++){
      if (smsky[1].adr[i]+i-j>=smsky[1].adr[i+1])  fprintf (out," %10.5f",0.0);
      else fprintf (out," %10.5f",smsky[1].a[smsky[1].adr[i]+i-j]);
    }
    for (j=i;j<smsky[1].n;j++){
      if (smsky[1].adr[j]+j-i>=smsky[1].adr[j+1])  fprintf (out," %10.5f",0.0);
      else fprintf (out," %10.5f",smsky[1].a[smsky[1].adr[j]+j-i]);
    }
  }
  
}





/**
   function solves reduced system of equations from primal
   domain decomposition by conjugate gradient method
   %matrix of reduced system is not assembled
   iteration process is directed by the master processor
   
   @param condmat - array containing reduced matrices (Schur complements)
   @param condvect - array containing reduced vectors
   @param out - output file (for auxiliary output)
   
   JK, 13.5.2009
*/
void seqschur::solve_red_sys_iter (double **/*condmat*/,double **/*condvect*/,FILE */*out*/)
{
  /*
  int ind;
  long i,j,k,sizebuff,gndofe,*cn;
  double alpha,beta,nom,denom,norrhs;
  double *lhs,*rhs,*buff,*p,*d,*r;

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
      
      if (i%100==0)  fprintf (stdout,"\n iteration %ld    norres/norrhs %e",i,nom/norrhs);
      //fprintf (out,"\n iteration %ld    norres/norrhs %e",i,nom/norrhs);
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
    
  
    //fprintf (out,"\n\n\n REDUCED SYSTEM SOLUTION (on master)\n");
    //for (i=0;i<ndofcp;i++){
      //fprintf (out,"\n lhs %ld     %e",i,lhs[i]);
    //}

    
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
  */
}


/**
   function solves reduced system of linear algebraic equations on the master
   processor by a direct method
   
   @param condmat - array containing reduced matrices (Schur complements)
   @param condvect - array containing reduced vectors
   @param out - output file (for auxiliary output)
   
   JK, 13.5.2009
*/
void seqschur::solve_red_sys_fin (double **condmat,double **condvect,FILE *out)
{
  long i,gndofe,*cn;
  double *lhs,*rhs;
  precond prec;
 
  arr = new gmatrix;
  
  //arr->ts = rsmstor;
  arr->setval (ssle);
  //arr->tlinsol=ldl;


  //arr->alloc ();
  arr->initiate (gtop,ndofcp,rsmstor,1,out);
  
  lhs = new double [ndofcp];
  memset (lhs,0,ndofcp*sizeof(double));
  rhs = new double [ndofcp];
  memset (rhs,0,ndofcp*sizeof(double));
  
  for (i=0;i<ns;i++){
    gndofe = gtop->give_ndofe (i);
    cn = new long [gndofe];
    gtop->give_code_numbers (i,cn);
    arr->localized (condmat[i],cn,gndofe,i);
    locglob (rhs,condvect[i],cn,gndofe);
    delete [] cn;
  }

  arr->prepmat (0.0,1);
  
  
  //  solution of the reduced problem
  arr->solve_system (gtop,prec,lhs,rhs,out);
  
  
  fprintf (out,"\n\n\n REDUCED SYSTEM SOLUTION (on master)\n");
  for (i=0;i<ndofcp;i++){
    fprintf (out,"\n lhs %ld     %e",i,lhs[i]);
  }
  
  for (i=0;i<ns;i++){
      nullv (condvect[i],nbdofmas[i]);
      
      gndofe = gtop->give_ndofe (i);
      cn = new long [gndofe];
      gtop->give_code_numbers (i,cn);
      globloc (lhs,condvect[i],cn,gndofe);
      delete [] cn;
  }
  
  
  delete [] lhs;  delete [] rhs;
  delete arr;
  
}

/**
   function solves reduced system of linear algebraic equations
   
   @param condmat - array containing reduced matrices (Schur complements)
   @param condvect - array containing reduced vectors
   @param out - output file (for auxiliary output)
   
   JK, 13.5.2009
*/
void seqschur::solve_red_sys (double **condmat,double **condvect,FILE *out)
{
  switch (trssol){
  case master_sol:{
    solve_red_sys_fin (condmat,condvect,out);
    break;
  }
  case paral_sol:{
    solve_red_sys_iter (condmat,condvect,out);
    break;
  }
  default:{
    print_err ("unknown type of solver of reduced problem is required",__FILE__,__LINE__,__func__);
  }
  }
}


/**
   function solves system of linear algebraic equations
   by the Schur complment method
   
   @param top - pointer to the sequential general topology
   @param gm - pointer to the subdomain %matrix stored in the form gmatrix
   @param lhs - array containing solution
   @param rhs - array containing right hand side
   @param out - output file (for auxiliary output)
   
   JK, 13.5.2009
*/
void seqschur::solve_system (gtopology */*top*/,gmatrix *gm,
			     double *lhs,double *rhs,FILE *out)
{
  long i,j;
  double **condmat,**condvect;
  time_t t1,t2,t3,t4;
  
  //  allocation of reduced matrix
  condmat = new double* [ns];
  for (i=0;i<ns;i++){
    condmat[i] = new double [nbdofmas[i]*nbdofmas[i]];
    for (j=0;j<nbdofmas[i]*nbdofmas[i];j++){
      condmat[i][j]=0.0;
    }
  }
  
  //  allocation of reduced vector
  condvect = new double* [ns];
  for (i=0;i<ns;i++){
    condvect[i] = new double [nbdofmas[i]];
    for (j=0;j<nbdofmas[i];j++){
      condvect[i][j]=0.0;
    }
  }
  

  //  assembling of subdomain matrices
  subdomain_matrices (gm,out);

  //  allocation of vectors of right hand sides on subdomains
  double **ff;
  ff = new double* [ns];
  for (i=0;i<ns;i++){
    ff[i] = new double [ndofdom[i]];
    for (j=0;j<ndofdom[i];j++){
      ff[i][j]=0.0;
    }
  }
  
  
  //  allocation of vectors of right hand sides on subdomains
  for (i=0;i<ns;i++){
    globloc (rhs,ff[i],cndom[i],ndofdom[i]);
  }

  double **dd;
  dd = new double* [ns];
  for (i=0;i<ns;i++){
    dd[i] = new double [ndofdom[i]];
    for (j=0;j<ndofdom[i];j++){
      dd[i][j]=0.0;
    }
  }

  t1 = time (NULL);
  
  //  elimination of inner DOFs = static condensation of inner DOFs
  //  the Schur complement matrix is stored in array condmat
  for (i=0;i<ns;i++){
    smsky[i].ldlkon_sky (condmat[i],condvect[i],dd[i],ff[i],nbdofmas[i],1);
  }
  
  t2 = time (NULL);
  
  //  solution of reduced problem
  solve_red_sys (condmat,condvect,out);
  
  t3 = time (NULL);
  

  //  back substitution on subdomains
  for (i=0;i<ns;i++){
    smsky[i].ldlkon_sky (condmat[i],condvect[i],dd[i],ff[i],nbdofmas[i],2);
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<ndofdom[i];j++){
      lhs[cndom[i][j]-1]=dd[i][j];
    }
  }
  /*
  //  allocation of vectors of right hand sides on subdomains
  for (i=0;i<ns;i++){
    gndofe = gtop->give_ndofe (i);
    cn = new long [gndofe];
    gtop->give_code_numbers (i,cn);
    locglob (lhs,dd[i],cn,gndofe);
    delete [] cn;
  }
  */

  t4 = time (NULL);
  
  // ****************************************************
  //  gathering of elapsed time on the master processor
  // ****************************************************
  long *buff;
  buff = new long[3];
  buff[0]=long(t2-t1);
  buff[1]=long(t3-t2);
  buff[2]=long(t4-t3);
  
  /*
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

    fprintf (out,"\n\n\n Domain %ld",j);
    fprintf (out,"\n time of condensation             %ld",buff[0]);
    fprintf (out,"\n time of solution of red. sys.    %ld",buff[1]);
    fprintf (out,"\n time of back substitution        %ld",buff[2]);
    
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

      fprintf (out,"\n\n\n Domain %ld",j);
      fprintf (out,"\n time of condensation             %ld",buff[0]);
      fprintf (out,"\n time of solution of red. sys.    %ld",buff[1]);
      fprintf (out,"\n time of back substitution        %ld",buff[2]);

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
    
    fprintf (out,"\n\n maximum factorization time  %ld",maxfact);
    fprintf (out,"\n coarse problem time         %ld",mincp);
    fprintf (out,"\n minimum total time          %ld",mintot);
    fprintf (out,"\n maximum total time          %ld",maxtot);

    delete [] fact;
    delete [] cp;
    
  }
  else{
    MPI_Send(buff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  */
  delete [] buff;
  
  
  delete [] condmat;
  delete [] condvect;
  
}

/**
   function gathers contributions from subdomains into one global %vector
   the global %vector must be deleted outside of the subroutine
   
   @param lv - local %vector
   @param gv - global %vector, it is allocated only on the master processor
               the array must be deleted outside of the subroutine
   
   JK, 6.11.2004
*/
/*
void seqschur::gather_bound_vect (double *lv,double *gv,long *domproc)
{
  long i,j,l,gndofe,*cn;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxnbdof];
  
  //  values from subdomain are selected
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
*/

/**
   function gathers contributions from subdomains into one global %vector
   the global %vector must be deleted outside of the subroutine
   
   @param lv - local %vector
   @param gv - global %vector, it is allocated only on the master processor
               the array must be deleted outside of the subroutine
   
   JK, 6.11.2004
*/
 /*
double seqschur::pss_gather_bound_vect (double *lv,double *gv,long *domproc)
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
 */




/**
   function scatters contributions from one global vector into local vectors
   function rewrites part of the local vector by contributions from boundaries
   
   @param lv - local vector
   @param gv - global vector, it is allocated only on the master processor
   
   JK, 6.11.2004
*/
  /*
void seqschur::scatter_bound_vect (double *lv,double *gv,long *domproc)
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
    lv[i]=buff[j]/2.0;  j++;
  }
  
  delete [] buff;
  
}
  */

/**
   function computes norm of vector of unbalanced values
   function gathers vector of unbalanced values on boundaries
   function scatters correct values into subdomains
   function returns norm of vector of unbalanced values
   
   @param lv - local vector
   @param domproc - domain-processor correspondence
   
   JK, 6.11.2004
*/
   /*
double seqschur::unbalanced_values (double *lv,long *domproc,FILE *out)
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
   */





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
void seqschur::search_comcn(gtopology *top,long *ltg, long *idgnn, long nbn, long **gnnc, long id, long ccn)
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
