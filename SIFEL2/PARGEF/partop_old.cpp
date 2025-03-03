#include "partop.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include "../GEFEL/galias.h"
/**
   constructor
   
   @param np - number of processors
   @param mr - my rank
   @param nd - number of subdomain
   
   JK, 3.5.2004
*/
partop::partop(int np,int mr,long nd,meshdescription meshd)
{
  //  number of processors
  nproc=np;
  //  my rank
  myrank=mr;
  //  number of subdomain
  ndom=nd;
  
  //  mesh description
  md = meshd;
  
  //  number of nodes on subdomain
  nn=0;
  // number of elements on subdomain
  ne=0;
  //  total number of nodes in the whole problem
  tnnp=0;
  //  total number of boundary nodes
  tnbn=0;
  //  number of nodes on subdomain
  nn=0;
  //  maximum number of nodes on one subdomain
  maxnn=0;
  //  number of boundary nodes on subdomain
  nbn=0;
  //  maximum number of boundary nodes on one subdomain
  maxnbn = 0;
  //  maximum number of all DOFs on one subdomain
  maxndof=0;
  
  
  //  array containing nodal multiplicity
  nodmultip = NULL;
  //  local numbers of boundary nodes
  lnbn = NULL;
  //
  lnin = NULL;
  //
  lgnbn = NULL;

  dofmultip = NULL;
  
  nodeidentif=NULL;

  if (myrank==0){
    nnsd = NULL;
    nbnd = NULL;
    multip = NULL;
    allnodes = NULL;
    cnbn = NULL;
    nbdofnd = NULL;
    nbdofd = NULL;
    lbcn = NULL;
    bnmultip = NULL;
    llnbn = NULL;
    ldn = NULL;
    
    //  array containing numbers of all degrees of freedom on subdomains
    nalldof=NULL;
    
    gnn = NULL;
    pgcn = NULL;
    gcnbn = NULL;
  }
  
}

/**
   destructor
   
   JK, 3.5.2004
*/
partop::~partop()
{
  long i;

  delete [] nodmultip;
  delete [] lnbn;
  
  if (myrank==0){
    delete [] nnsd;
    delete [] nbnd;
    delete [] multip;
    delete [] nbdofd;
    delete [] bnmultip;
    
    
    for (i=0;i<nproc;i++){
      delete [] allnodes[i];
    }
    delete [] allnodes;
    
    for (i=0;i<nproc;i++){
      delete [] cnbn[i];
    }
    delete [] cnbn;
    
    for (i=0;i<nproc;i++){
      delete [] nbdofnd[i];
    }
    delete [] nbdofnd;
    
    for (i=0;i<tnbn;i++){
      delete [] llnbn[i];
    }
    delete [] llnbn;
    
    for (i=0;i<tnbn;i++){
      delete [] ldn[i];
    }
    delete [] ldn;
    
    //  missing deallocation of array lbcn
  }
}


/**
   function assignes node numbers to auxiliary indicators at nodes
   
   @param top - pointer to general topology
   @param ltg - array of local to global correcpondence
   
   JK, 1.7.2005
*/
void partop::initiation (gtopology *top,long *ltg)
{
  long i;
  
  nn = top->nn;
  ne = top->ne;
  

  switch (md){
  case all_nodes:{
    for (i=0;i<nn;i++){
      top->gnodes[i].ai = ltg[i];
    }
    break;
  }
  case bound_nodes:{
    for (i=0;i<nn;i++){
      top->gnodes[i].ai = ltg[i];
    }
    break;
  }
  case neg_bound_nodes:{
    for (i=0;i<nn;i++){
      top->gnodes[i].ai = abs(ltg[i]);
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of mesh description is required");
    fprintf (stderr,"\n in function initiation (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }


}

/**
   function collects numbers of all nodes on subdomains
   
   the following function has to be called before this function:
   void initiation (gtopology *top,long *ltg);
   
   
   function establishes:

   maxnn - maximum number of nodes on one subdomain

   nnsd - array containing numbers of nodes on subdomains
   
   
   @param domproc - array containing domain-processor correspondence
   
   JK, 30.6.2005
*/
void partop::numbers_of_all_nodes_on_subdomains (long *domproc,FILE *out)
{
  long i,j,k;
  MPI_Status stat;
  
  //  maximum number of nodes defined on subdomain (will be sent to the master)
  maxnn=nn;
  
  // *********
  //  master
  // *********
  if (myrank==0){
    //  array containing number of nodes on subdomains
    if (nnsd!=NULL)
      delete [] nnsd;
    nnsd = new long [nproc];
    
    //  master contribution
    j=domproc[0];
    nnsd[j]=maxnn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      if (maxnn<k)  maxnn=k;
      nnsd[j]=k;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  // ********
  // slaves
  // ********
  else{
    MPI_Send (&nn,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxnn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  if (myrank==0){
    fprintf (out,"\n\n\n numbers of nodes on each subdomain, maxnn %ld\n\n",maxnn);
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nnsd[i]);
    }
  }

}


/**
   function collects numbers of all degrees of freedom on subdomains
   supports and constraints are not taken into account
   
   function establishes:
   
   maxndof - maximum number of degrees of freedom on one subdomain
             it summarizes number of DOFs at nodes belonging to one subdomain
	     it differs from the actual degrees of freedom of subdomains in constraints
	     
   nalldof - array of numbers of degrees of freedom on subdomains
   
   JK, 20.3.2007
*/
void partop::numbers_of_all_dofs_on_subdomains (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k;
  MPI_Status stat;
  
  //  maximum number of nodes defined on subdomain (will be sent to the master)
  maxndof=0;
  for (i=0;i<nn;i++){
    maxndof+=top->give_ndofn (i);
  }
  
  // *********
  //  master
  // *********
  if (myrank==0){
    //  array containing number of nodes on subdomains
    if (nalldof!=NULL)
      delete [] nalldof;
    nalldof = new long [nproc];
    
    //  master contribution
    j=domproc[0];
    nalldof[j]=maxndof;
    
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      if (maxndof<k)  maxndof=k;
      nalldof[j]=k;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxndof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  // ********
  // slaves
  // ********
  else{
    MPI_Send (&maxndof,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxndof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  if (myrank==0){
    fprintf (out,"\n\n\n numbers of all DOFs on each subdomain, maxndof %ld\n\n",maxndof);
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nalldof[i]);
    }
  }

}



/**
   function computes node multiplicity
   function is used only for mesh description = all_nodes or neg_bound_nodes
   
   the following functions have to be called before this function:
   void initiation (gtopology *top,long *ltg);
   void numbers_of_all_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out);

   
   function establishes:
   
   tnnp - total number of nodes in problem
   
   multip - array of boundary node multiplicity (on the master)
   
   nodmultip - array of node multiplicity (on all processors)
   
   
   @param top - pointer to subdomain topology
   @param domproc - array containing domain-processor correspondence
   
   JK, 30.6.2005
*/
void partop::compute_multiplicity (long *ltg,long *domproc,FILE *out)
{
  long i,j,k,a,buffsize;
  long *buff;
  MPI_Status stat;
  
  if (maxnn==0){
    fprintf (stderr,"\n\n maximum number of nodes on one subdomain is equal to zero,\n");
    fprintf (stderr,"\n function numbers_of_all_nodes_on_subdomains has to be called first,\n");
    abort ();
  }
  
  switch (md){
  case bound_nodes:{
    
    buffsize=2;
    buff = new long [buffsize];
    
    //  determination of total number of boundary nodes
    maxnbn=0;
    tnbn=0;
    for (i=0;i<nn;i++){
      if (tnbn<ltg[i])
	tnbn=ltg[i];
      if (ltg[i]>-1)
	maxnbn++;;
    }
    
    buff[0]=tnbn;
    buff[1]=maxnbn;
    
    if (myrank==0){
      for (i=1;i<nproc;i++){
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	if (tnbn<buff[0])
	  tnbn=buff[0];
	if (maxnbn<buff[1])
	  maxnbn=buff[1];
      }
      tnbn++;
      
      buff[0]=tnbn;
      buff[1]=maxnbn;
      
      for (i=1;i<nproc;i++){
	MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    tnbn=buff[0];
    maxnbn=buff[1];
    
    delete [] buff;
    
    buffsize = maxnbn;
    buff = new long [buffsize];
    for (i=0;i<maxnbn;i++){
      buff[i]=-1;
    }
    
    j=0;
    for(i=0;i<nn;i++){
      if (ltg[i]>-1){
	buff[j] = ltg[i];
	j++;
      }
    }
    
    
    if(myrank == 0){
      if (multip!=NULL){
	delete [] multip;
      }
      multip = new long [tnbn];
      for(i=0;i<tnbn;i++){
	multip[i]=0;
      }
      
      // master contribution
      for(j=0;j<maxnbn;j++){
	if (buff[j]>-1){
	  multip[buff[j]]++;
	}
      }
      
      //  slave contributions
      for (i=1;i<nproc;i++){
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	for(j=0;j<maxnbn;j++){
	  if (buff[j]>-1){
	    multip[buff[j]]++;
	  }
	}
      }
      
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    delete [] buff;
    
    buffsize=tnbn;
    buff = new long [buffsize];
    
    if (myrank==0){
      for (i=0;i<tnbn;i++){
	buff[i]=multip[i];
      }
      
      for (i=1;i<nproc;i++){
     	MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (nodmultip!=NULL)
      delete [] nodmultip;
    nodmultip = new long [nn];
    
    for (i=0;i<nn;i++){
      if (ltg[i]>-1){
	nodmultip[i]=buff[ltg[i]];
      }
      else{
	nodmultip[i]=1;
      }
    }
    
    delete [] buff;
    
    break;
  }




  case all_nodes:
  case neg_bound_nodes:{
    
    buffsize=maxnn+1;
    buff = new long [buffsize];
    tnnp=0;
    for (i=0;i<nn;i++){
      buff[i]=ltg[i];
      if (tnnp<ltg[i])
	tnnp=ltg[i];
    }
    buff[maxnn]=tnnp;
    
    // *********
    //  master
    // *********
    if (myrank==0){
      if (allnodes != NULL){
	for (i=0;i<nproc;i++){
	  delete [] allnodes[i];
	}
	delete [] allnodes;
      }
      allnodes = new long* [nproc];
      for (i=0;i<nproc;i++){
	allnodes[i] = new long [nnsd[i]];
      }
      
      
      k=domproc[0];
      for (j=0;j<nnsd[k];j++){
	allnodes[k][j]=buff[j];
      }
      if (tnnp<buff[maxnn])
	tnnp=buff[maxnn];
      

      for (i=1;i<nproc;i++){
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	k=domproc[stat.MPI_TAG];
	for (j=0;j<nnsd[k];j++){
	  allnodes[k][j]=buff[j];
	}
	if (tnnp<buff[maxnn])
	  tnnp=buff[maxnn];
      }
    }
    
    // ********
    // slaves
    // ********
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    
    if (myrank==0){
      //  array containing number of subdomains which share node
      //  this is called node multiplicity
      
      //  must be increased due to indices from 0 instead of 1
      tnnp++;
      
      fprintf (out,"\n\n\n total number of nodes on whole problem %ld\n",tnnp);
      if (multip != NULL)
	delete [] multip;
      multip = new long [tnnp];
      for (i=0;i<tnnp;i++){
	multip[i]=0;
      }
      
      //  computation of node incidences
      for (i=0;i<nproc;i++){
	for (j=0;j<nnsd[i];j++){
	  multip[allnodes[i][j]]++;
	}
      }
    }
    
    
    if (myrank==0){
      for (i=1;i<nproc;i++){
	k=domproc[i];
	for (j=0;j<nnsd[k];j++){
	  buff[j]=multip[allnodes[k][j]];
	}
	MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
      
      k=domproc[0];
      for (j=0;j<nnsd[k];j++){
	buff[j]=multip[allnodes[k][j]];
      }
    }
    else{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    
    //  nodal multiplicity
    if (nodmultip != NULL)
      delete [] nodmultip;
    nodmultip = new long [nn];
    for (i=0;i<nn;i++){
      nodmultip[i] = buff[i];
    }
    
    delete [] buff;
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of mesh description is required in");
    fprintf (stderr,"\n function compute_multiplicity (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  /*
  fprintf(out,"\n\n\n Multiplicity of nodes on subdomain\n\n");  
  for (i=0;i<nn;i++)
    fprintf(out,"%ld   %ld\n",i+1,nodmultip[i]);  
  
  // kontrolni tisk allnodes
  if (myrank == 0){
    fprintf(out,"\n\n\n Global node numbers on each subdomain\n\n");
    fprintf(out,"the j-th node on the i-th subdomain has global number k\n");
    for (i=0;i<nproc;i++){
      k=domproc[i];
      fprintf(out,"Subdomain %ld\n",k+1);
      for (j=0;j<nnsd[k];j++){
	fprintf(out,"%ld   %ld\n",j+1,allnodes[i][j]+1);
      }
    }
  }
  */
  
}

/**
   function assmebles multiplicity of degrees of freedom
   it is different from the function compute_multiplicity, which
   computes multiplicity of nodes
   
   JK, 14.6.2007
*/
void partop::dof_multiplicity (gtopology *top,FILE *out)
{
  long i,j,k,nm;
  long ndofn;
  
  ndof = 0;
  for (i=0;i<nn;i++){
    ndofn=top->give_ndofn (i);
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>0){
	ndof++;
      }
    }
  }
  
  //printf (out,"\n\n ndof   %ld\n",ndof);
  
  
  dofmultip = new long [ndof];
  
  for (i=0;i<nn;i++){
    ndofn=top->give_ndofn (i);
    nm=nodmultip[i];
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>0){
	dofmultip[k-1]=nm;
      }
    }
  }
  /*
  fprintf (out,"\n\n\n kontrola dofmultip \n");
  for (i=0;i<ndof;i++){
    fprintf (out,"%ld   %ld\n",i+1,dofmultip[i]);
  }
  */

}


/**
   function selects boundary nodes
   
   the following functions have to be called before this function:
   void initiation (gtopology *top,long *ltg);
   void numbers_of_all_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out);
   void compute_multiplicity (gtopology *top,long *domproc,FILE *out);

   
   function establishes:
   
   maxnbn - maximum number of boundary nodes
   
   nbn - number of boundary nodes on each subdomain
   
   tnbn - total number of boundary nodes
   
   nbnd - array containing numbers of boundary nodes on subdomains
   
   lnbn - local numbers of boundary nodes
   
   @param ltg - local to global correspondence
   @param domproc - array containing domain-processor correspondence
   
   JK, 30.6.2005
   
*/
void partop::find_boundary_nodes (long *ltg,long *domproc,FILE *out)
{
  long i,j,k,l,buffsize;
  long *buff;
  MPI_Status stat;
  
  
  // *******************************************************
  //  searching for number of boundary nodes on subdomains
  // *******************************************************
  
  switch (md){
  case bound_nodes:{
    //  number of boundary nodes on subdomain
    nbn=0;  maxnbn=0;
    for (i=0;i<nn;i++){
      if (ltg[i]>-1)  nbn++;
      if (maxnbn<ltg[i])  maxnbn=ltg[i];
    }
    
    buffsize=2;
    buff = new long [buffsize];
    buff[0]=nbn;
    buff[1]=maxnbn;
    
    
    // *********
    //  master
    // *********
    if (myrank==0){
      //  list of number of boundary nodes on subdomains
      if (nbnd != NULL)
	delete [] nbnd;
      nbnd = new long [nproc];
      
      //  master contribution
      j=domproc[0];
      if (maxnbn<buff[0])  maxnbn=buff[0];
      if (tnbn<buff[1])  tnbn=buff[1];
      nbnd[j]=buff[0];
      
      for (i=1;i<nproc;i++){
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	//  slave contributions
	j=domproc[stat.MPI_TAG];
	if (maxnbn<buff[0])  maxnbn=buff[0];
	if (tnbn<buff[1])  tnbn=buff[1];
	nbnd[j]=buff[0];
	
      }
      
      //  node numbers start from 0 with respect to C language notation
      tnbn++;
      
      for (i=1;i<nproc;i++){
	MPI_Send (&maxnbn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
    }
    // ********
    // slaves
    // ********
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
      MPI_Recv (&maxnbn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    
    delete [] buff;
    fprintf (out,"\n\n maxnbn is %ld",maxnbn);
    break;
  }

  case all_nodes:
  case neg_bound_nodes:{
    
    if (myrank==0){
      if (multip==NULL){
	fprintf (stderr,"\n\n array multip is not allocated in function find_boundary_nodes,\n");
	fprintf (stderr,"\n function compute_multiplicity has to be called first (file %s, line %d),\n",__FILE__,__LINE__);
      }
    }
    
    if (myrank==0){
      //  number of all boundary nodes
      tnbn=0;
      for (i=0;i<tnnp;i++){
	if (multip[i]!=1)  tnbn++;
      }
      
      //  list of number of boundary nodes on subdomains
      if (nbnd != NULL)
	delete [] nbnd;
      nbnd = new long [nproc];

      maxnbn=0;
      for (i=0;i<nproc;i++){
	k=0;
	for (j=0;j<nnsd[i];j++){
	  if (multip[allnodes[i][j]]>1)  k++;
	}
	if (maxnbn<k)  maxnbn=k;
	nbnd[i]=k;
      }
      
      for (i=1;i<nproc;i++){
	MPI_Send (&maxnbn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (&maxnbn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    
    
    //  number of boundary nodes on subdomain
    nbn=0;
    for (i=0;i<nn;i++){
      if (nodmultip[i]>1)
	nbn++;
    }

    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of mesh description is required in");
    fprintf (stderr,"\n function find_boundary_nodes (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  
  
  
  buffsize=maxnbn;
  buff = new long [buffsize];
  
  switch (md){
  case bound_nodes:{

    //  list of coarse numbers of boundary nodes
    if (lgnbn!=NULL)
      delete [] lgnbn;
    lgnbn = new long [nbn];

    j=0;
    for (i=0;i<nn;i++){
      if (ltg[i]>-1){
	lgnbn[j]=ltg[i];
	j++;
      }
    }
    break;
  }
  case all_nodes:
  case neg_bound_nodes:{
    
    if (myrank==0){
      for (i=1;i<nproc;i++){
	k=domproc[i];
	l=0;
	for (j=0;j<nnsd[k];j++){
	  if (multip[allnodes[k][j]]>1){
	    buff[l]=allnodes[k][j];
	    l++;
	  }
	}
	
	MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
      
      k=domproc[0];
      l=0;
      for (j=0;j<nnsd[k];j++){
	if (multip[allnodes[k][j]]>1){
	  buff[l]=allnodes[k][j];
	  l++;
	}
      }
      
    }
    else{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }

    //  list of coarse numbers of boundary nodes
    if (lgnbn!=NULL)
      delete [] lgnbn;
    lgnbn = new long [nbn];
    
    for (i=0;i<nbn;i++){
      lgnbn[i]=buff[i];
    }
    
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of mesh description is required in");
    fprintf (stderr,"\n function find_boundary_nodes (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  
  
  delete [] buff;
  
  
  // *******************************
  //  indication of boundary nodes 
  // *******************************
  //  local numbers of boundary nodes
  if (lnbn != NULL)
    delete [] lnbn;
  lnbn = new long [nbn];
  
  switch (md){
  case all_nodes:
  case neg_bound_nodes:{
    j=0;
    for (i=0;i<nn;i++){
      if (nodmultip[i]>1){
	lnbn[j]=i;
	j++;
      }
    }
    
    break;
  }
    
  case bound_nodes:{
    j=0;
    for (i=0;i<nn;i++){
      if (ltg[i]>-1){
	lnbn[j]=i;
	j++;
      }
    }
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of mesh description is required in");
    fprintf (stderr,"\n function find_boundary_nodes (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  

  //  number of internal nodes
  nin = nn-nbn;
  
  //  array containing local numbers of internal nodes
  if (lnin != NULL)
    delete [] lnin;
  lnin = new long [nin];
  
  buff = new long [nn];
  for (i=0;i<nn;i++){
    buff[i]=0;
  }
  
  for (i=0;i<nbn;i++){
    j=lnbn[i];
    buff[j]=1;
  }
  
  j=0;
  for (i=0;i<nn;i++){
    if (buff[i]==0){
      lnin[j]=i;
      j++;
    }
  }

  delete [] buff;

  MPI_Barrier (MPI_COMM_WORLD);
  

  if (myrank==0){
    fprintf (out,"\n\n\n numbers of boundary nodes on each subdomain \n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nbnd[i]);
    }
  }
  
  fprintf (out,"\n\n\n numbers of boundary nodes on subdomain \n\n");
  fprintf (out,"\n %ld",nbn);
  fprintf (out,"\n\n\n local numbers of boundary nodes on subdomain");
  for (i=0;i<nbn;i++){
    fprintf (out,"\n %ld   %ld",i,lnbn[i]+1);
  }

  fprintf (out,"\n\n\n numbers of internal nodes on subdomain \n\n");
  fprintf (out,"\n %ld",nin);
  fprintf (out,"\n\n\n local numbers of internal nodes on subdomain");
  for (i=0;i<nin;i++){
    fprintf (out,"\n %ld   %ld",i,lnin[i]+1);
  }
  

}

/**
   function rewrites array ltg
   
   JK, 8.8.2007
*/
void partop::rewrite_ltg (long *ltg)
{
  long i;
  
  for (i=0;i<nn;i++){
    ltg[i]=-1;
  }
  for (i=0;i<nbn;i++){
    ltg[lnbn[i]]=lgnbn[i];
  }
}


/**
   function assembles boundary nodes on the master processor
   
   the following functions have to be called before this function:
   void initiation (gtopology *top,long *ltg);
   void numbers_of_all_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out);
   void compute_multiplicity (gtopology *top,long *domproc,FILE *out);
   void find_boundary_nodes (gtopology *top,long *domproc,FILE *out);
   
   
   function establishes:
   
   cnbn - array of coarse numbers of boundary nodes
   
   @param top - pointer to subdomain topology
   @param domproc - array containing domain-processor correspondence
   
   JK, 30.6.2005
*/
void partop::boundary_nodes_on_master (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,m,buffsize;
  long *buff;
  MPI_Status stat;
  
  if (myrank==0){
    if (nbnd==NULL){
      fprintf (stderr,"\n\n array containing numbers of boundary nodes on subdomains is not allocated,\n");
      fprintf (stderr,"\n function find_boundary_nodes has to be called first (file %s, line %d),\n",__FILE__,__LINE__);
    }
  }
  
  
  if (myrank==0){
    if (cnbn != NULL){
      for (i=0;i<nproc;i++){
	delete [] cnbn[i];
      }
      delete [] cnbn;
    }
    cnbn = new long* [nproc];
    for (i=0;i<nproc;i++){
      cnbn[i] = new long [nbnd[i]];
    }
  }
  
  switch (md){
  case all_nodes:
  case neg_bound_nodes:{
    
    if (myrank==0){
      if (gcnbn != NULL)
	delete [] gcnbn;
      gcnbn = new long [tnnp];

      j=0;
      for (i=0;i<tnnp;i++){
	if (multip[i]>1){
	  gcnbn[i]=j;
	  j++;
	}
	else{
	  gcnbn[i]=-1;
	}
      }
      
      /*
      // kontrolni tisk pole gcnbn - pro allnodes je to co gnodes.ai pro bounadry nodes
      fprintf (out,"\n\n gcnbn\n");
      for (i=0;i<tnnp;i++){
	fprintf (out,"%ld  %ld\n",i+1,gcnbn[i]); 
      }
      */
      
      if (j!=tnbn){
	fprintf (stderr,"\n\n total number of boundary nodes computed in function boundary_nodes_on_master");
	fprintf (stderr,"\n differs from the number computed in function find_boundary_nodes (file %s, line %d),\n",__FILE__,__LINE__);
      }
      
      for (i=0;i<nproc;i++){
	m=0;
	for (j=0;j<nnsd[i];j++){
	  k=allnodes[i][j];
	  if (gcnbn[k]>-1){
	    cnbn[i][m]=gcnbn[k];
	    m++;
	  }
	}
      }
      delete [] gcnbn;
      gcnbn = NULL;
    }
    
    break;
  }
  
  case bound_nodes:{
    
    if (maxnbn==0){
      fprintf (stderr,"\n\n maximum number of boundary nodes on one subdomain is equal to zero,");
      fprintf (stderr,"\n function find_boundary_nodes should be called first (file %s, line %d),\n",__FILE__,__LINE__);
    }
    
    buffsize=maxnbn;
    buff = new long [buffsize];
    
    
    j=0;
    for (i=0;i<nn;i++){
      if (top->gnodes[i].ai>-1){
	buff[j]=top->gnodes[i].ai;
	j++;
      }
    }
    
    if (myrank==0){
      
      //  master contribution
      k=domproc[0];
      for (j=0;j<nbnd[k];j++){
	cnbn[k][j]=buff[j];
      }
      
      for (i=1;i<nproc;i++){
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	//  slave contributions
	k=domproc[stat.MPI_TAG];
	for (j=0;j<nbnd[k];j++){
	  cnbn[k][j]=buff[j];
	}
	
      }
 
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    
    delete [] buff;
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown type of mesh description is required in");
    fprintf (stderr,"\n function boundary_nodes_on_master (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }

  
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  if (myrank==0){
    fprintf (out,"\n\n\n\n coarse numbering of boundary nodes  - cnbn\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain number %ld\n",i+1);
      for (j=0;j<nbnd[i];j++){
	fprintf (out,"%ld %ld\n",j,cnbn[i][j]+1);
      }
    }
  }

  
}


/**
   function assembles correspondence among coarse nodes and local boundary nodes
   
   function assembles:
   
   llnbn - list of local numbers of boundary nodes belonging to coarse node
   ldn - list of subdomains which contain boundary nodes belonging to coarse node
   
   not checked

   JK, 6.7.2005
*/
void partop::coarse_local_nodes ()
{
  long i,j,k,m;
  
  if (myrank==0){
    //  multiplicity of boundary nodes
    bnmultip = new long [tnbn];
    for (i=0;i<tnbn;i++){
      bnmultip[i]=0;
    }
    
    if (md == all_nodes){
      j=0;
      for (i=0;i<tnnp;i++){
	if (multip[i]>1){
	  bnmultip[j]=multip[i];
	  j++;
	}
      }
    }
    
    if (md == bound_nodes){
      for (i=0;i<nproc;i++){
	for (j=0;j<nbnd[i];j++){
	  bnmultip[cnbn[i][j]]++;
	}
      }
    }
    
    llnbn = new long* [tnbn];
    ldn = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      llnbn[i] = new long [bnmultip[i]];
      ldn[i] = new long [bnmultip[i]];
      bnmultip[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nbnd[i];j++){
	k=cnbn[i][j];
	m=bnmultip[k];
	llnbn[k][m]=j;
	ldn[k][m]=i;
	bnmultip[k]++;
      }
    }
  }
}

/**
   function sorts node numbers of nodes shared at coarse node increasingly
   
   not checked

   JK, 11.7.2005
*/
void partop::sort_nodes (FILE *out)
{
  long i,j,k,m,n,min;
  
  if (myrank==0){
    for (i=0;i<tnbn;i++){
      for (j=0;j<bnmultip[i];j++){
	min=nproc;
	for (k=j;k<bnmultip[i];k++){
	  if (ldn[i][k]<min){
	    min=ldn[i][k];
	    m=k;
	  }
	}
	n=ldn[i][j];
	ldn[i][j]=ldn[i][m];
	ldn[i][m]=n;
	
	n=llnbn[i][j];
	llnbn[i][j]=llnbn[i][m];
	llnbn[i][m]=n;
      }
    }
  }
  
  if (myrank==0){
    if(md == all_nodes){
      fprintf (out,"\n\n Multiplicity of all nodes on whole problem\n");
      for (i=0;i<tnnp;i++){
	fprintf (out,"%ld  %ld\n",i+1,multip[i]); 
      }
    }
    fprintf (out,"\n\n Multiplicity of all boundary nodes on whole problem\n");
    for (i=0;i<tnbn;i++){
      fprintf (out,"%ld  %ld\n",i+1,bnmultip[i]); 
    }
    
    
    /*
    fprintf (out,"\n\n  ldn and llnbn array\n");
    // llnbn ukazuje na pozici v poli lnbn
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n boundary node %4ld  ",i+1);
      for (j=0;j<bnmultip[i];j++){
	fprintf (out," ldn =  %ld  llnbn =  %ld  lnbn = %ld,",ldn[i][j],llnbn[i][j],lnbn[llnbn[i][j]]);
      }
    }*/
  }
  
}
  


/**
   function detects numbers of degrees of freedom on boundary nodes
   

   the following functions have to be called before this function:
   void initiation (gtopology *top,long *ltg);
   void numbers_of_all_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out);
   void compute_multiplicity (gtopology *top,long *domproc,FILE *out);
   void find_boundary_nodes (gtopology *top,long *domproc,FILE *out);
   void boundary_nodes_on_master (gtopology *top,long *domproc,FILE *out);

   function assembles:
   
   nbdofnd - array containing numbers of DOFs on boundary nodes
   
   @param top - pointer to general topology
   @param domproc - array containing domain-processor correspondence

   JK, 6.7.2005
*/
void partop::number_of_bdofs_on_nodes (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;
  
  if (lnbn==NULL){
    fprintf (stderr,"\n\n array containing local numbers of boundary nodes is not allocated,\n");
    fprintf (stderr,"\n function find_boundary_nodes has to be called first (file %s, line %d),\n",__FILE__,__LINE__);
  }

  buffsize = maxnbn;
  buff = new long [buffsize];
  
  //  detection of numbers of DOFs on boundary nodes
  for (i=0;i<nbn;i++){
    j=lnbn[i];
    buff[i]=top->give_ndofn (j);
  }
  
  if (myrank==0){
    if (nbdofnd != NULL){
      for (i=0;i<nproc;i++){
	delete [] nbdofnd[i];
      }
      delete [] nbdofnd;
    }
    nbdofnd = new long* [nproc];
    for (i=0;i<nproc;i++){
      nbdofnd[i] = new long [nbnd[i]];
    }
    
    //  contribution from master
    k=domproc[0];
    for (j=0;j<nbnd[k];j++){
      nbdofnd[k][j]=buff[j];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      k=domproc[stat.MPI_TAG];
      for (j=0;j<nbnd[k];j++){
	nbdofnd[k][j]=buff[j];
      }
      
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  
  if (myrank==0){
    fprintf (out,"\n\n\n numbers of DOFs at boundary nodes");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain %ld",i);
      for (j=0;j<nbnd[i];j++){
	fprintf (out,"\n %ld %ld",j+1,nbdofnd[i][j]);
      }
    }
  }
  
}

/**
   function assembles code numbers on master

   in fact, it assembles indicators 0 or 1 which will be
   used later for code number generation
   
   the following functions have to be called before this function:
   void initiation (gtopology *top,long *ltg);
   void numbers_of_all_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out);
   void compute_multiplicity (gtopology *top,long *domproc,FILE *out);
   void find_boundary_nodes (gtopology *top,long *domproc,FILE *out);
   void boundary_nodes_on_master (gtopology *top,long *domproc,FILE *out);
   void number_of_bdofs_on_nodes (gtopology *top,long *domproc,FILE *out);


   function assembles:
   
   nbdofd - array containing number of boundary DOFs on subdomains
   lbcn - array containing local boundary code numbers
   
   @param top - pointer to general topology
   @param domproc - array containing domain-processor correspondence
   
   JK, 6.7.2005
*/
void partop::code_numbers_on_master (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,m,l,ndofn,buffsize;
  long *buff;
  MPI_Status stat;
  
  if (myrank==0){
    if (nbdofnd==NULL){
      fprintf (stderr,"\n\n array ndofnd is not allocated in function code_numbers_on_master,\n");
      fprintf (stderr,"\n function number_of_dofs_on_nodes has to be called first (file %s, line %d),\n",__FILE__,__LINE__);
    }
  }
  
  //  detection of numbers of DOFs on particular subdomains
  if (myrank==0){
    if (nbdofd != NULL)
      delete [] nbdofd;
    nbdofd = new long [nproc];
    
    maxnbdof=0;
    for (i=0;i<nproc;i++){
      nbdofd[i]=0;
      for (j=0;j<nbnd[i];j++){
	nbdofd[i]+=nbdofnd[i][j];
      }
      if (maxnbdof<nbdofd[i])
	maxnbdof=nbdofd[i];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnbdof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&maxnbdof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  
  buffsize = maxnbdof;
  buff = new long [buffsize];
  
  m=0;
  for (i=0;i<nbn;i++){
    k=lnbn[i];
    ndofn = top->give_ndofn (k);
    for (j=0;j<ndofn;j++){
      buff[m]=top->give_dof (k,j);
      m++;
    }
  }
  
  if (myrank==0){
    if (lbcn != NULL){
      for (i=0;i<nproc;i++){
	for (j=0;j<nbnd[i];j++){
	  delete [] lbcn[i][j];
	}
	delete [] lbcn[i];
      }
      delete [] lbcn;
    }
    lbcn = new long** [nproc];
    for (i=0;i<nproc;i++){
      lbcn[i] = new long* [nbnd[i]];
      for (j=0;j<nbnd[i];j++){
	lbcn[i][j] = new long [nbdofnd[i][j]];
      }
    }
    
    //  master contribution
    l=0;
    m=domproc[0];
    for (j=0;j<nbnd[m];j++){
      for (k=0;k<nbdofnd[m][j];k++){
	lbcn[m][j][k]=buff[l];
	l++;
      }
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      l=0;
      m=domproc[stat.MPI_TAG];
      for (j=0;j<nbnd[m];j++){
	for (k=0;k<nbdofnd[m][j];k++){
	  lbcn[m][j][k]=buff[l];
	  l++;
	}
      }

    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  
  
  if (myrank==0){
    fprintf (out,"\n\n\n kontrola pole lbcn \n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domena  %ld",i);
      for (j=0;j<nbnd[i];j++){
	fprintf (out,"\n boundary uzel  %ld",j);
	for (k=0;k<nbdofnd[i][j];k++){
	  fprintf (out,"  %ld",lbcn[i][j][k]);
	}
	fprintf (out,"\n");
      }
    }
  }
  
  
  
  MPI_Barrier (MPI_COMM_WORLD);
  
}

/**
   function recognizes homogeneous or nonhomogeneous meshes
   
   @param top - general topology
   
   29.5.2007, JB
*/
long partop::give_whole_elemtype(gtopology *top)
{
  long i,t;
  gelemtype get;
  long *elem;
  long type = 0;
  elem = new long[ne];
  
  for(i = 0; i < ne; i++){
    get =  top->give_elem_type (i);
    switch(get){
    case linbar:{     t= 1;    break;  }
    case quadbar:{    t= 2;    break;  }
    case lintriag:{   t=21;    break;  }
    case quadtriag:{  t=22;    break;  }
    case linquad:{    t=25;    break;  }
    case quadquad:{   t=26;    break;  }
    case lintetra:{   t=41;    break;  }
    case quadtetra:{  t=42;    break;  }
    case linhexa:{    t=45;    break;  }
    case quadhexa:{   t=46;    break;  }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function give_whole_elemtype (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
    elem[i] = t;
  }
  
  for(i = 1; i < ne; i++){
    if(elem[i] != elem[i-1]){
      fprintf (stderr,"\n\n inhomogeneous mesh, function identification_node can not be used\n"); 
    }
    else{
      type = elem[i];
    }
  }
  return type; 
  delete []elem;
  
}

/**
   function recognizes spatial dimension of the mesh
   
   @param top - general topology
   
   29.5.2007, JB
*/
long partop::give_whole_dim(gtopology *top)
{
  
  gelemtype get;
  long *elem;
  long dim,i,t;
  dim = -1;
  
  elem = new long[ne];
  
  for(i = 0; i < ne; i++){
    get =  top->give_elem_type (i);
    switch(get){
    case linbar:{     t = 1;    break;  }
    case quadbar:{    t = 1;    break;  }
    case lintriag:{   t = 2;    break;  }
    case quadtriag:{  t = 2;    break;  }
    case linquad:{    t = 2;    break;  }
    case quadquad:{   t = 2;    break;  }
    case lintetra:{   t = 3;    break;  }
    case quadtetra:{  t = 3;    break;  }
    case linhexa:{    t = 3;    break;  }
    case quadhexa:{   t = 3;    break;  }
    default:{
      fprintf (stderr,"\n\n unknown element type is required in function give_whole_dim (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
    elem[i] = t;
  }
  
  
  for(i = 1; i < ne; i++){
    if(elem[i] != elem[i-1]){
      fprintf (stderr,"\n\n different spatial dimension is found in the mesh, function identification_node can not be used\n");  
    }
    else{
      dim = elem[i];
    }
  }
  
  return(dim);
  
  delete []elem;

}


void partop::identification_node(gtopology *top,long *ltg,long *domproc, FILE *out)
{
  long dim;
  dim = give_whole_dim(top);
  //dim = 3;
  fprintf(out,"\n\n\nDimension of problem is %ld\n",dim);
  switch(dim){
  case -1:
    fprintf (stderr,"\n\n Spatial dimension is not detected (file %s, line %d).\n",__FILE__,__LINE__);   
    break;
  case 1:
    break;
  case 2: identification_node_2d (top,ltg,domproc,out);
    break;
  case 3: identification_node_3d (top,ltg,domproc,out);
    break;
  default:
   fprintf (stderr,"\n unknown dimension is used in function identification_node (file %s, line %d)\n",__FILE__,__LINE__);
   break;
  }
}


/**
   function searches for corner nodes in 2D problems
   
   array nodeidentif contains node identification on subdomain
   
   1 - internal node
   2 - edge node
   3 - corner node

   @param out - pointer to output file stream
   @param top - pointer to general topology

   08.03.2007 JB
*/
void partop::identification_node_2d (gtopology *top, long *ltg,long *domproc, FILE *out)
{
  long i,j,k,l,m,n,a,b;
  
  if (nodeidentif != NULL)
    delete []nodeidentif;
  nodeidentif = new long [nn];
  
  for (i=0;i<nn;i++){
    nodeidentif[i]=1; 
  }
  
  //  assembling of list of adjacent nodes to each node
  //top->adjacnodes (out);
  top->adjacnodes_edge (out);
  
  // ********
  // trideni
  // ********
  
  for (i = 0; i < nn; i++){
    if (nodmultip[i] == 1)
      //  internal node
      nodeidentif[i] = 1;
    if(nodmultip[i] > 2){
      //  vertex node
      nodeidentif[i] = 3;
    }
    // stare a funkcni
    /*
      if(nodmultip[i] == 2){
      l=0;
      for (j = 0; j < top->nadjnodnod[i]; j++){
      a = top->adjnodnod[i][j];
      if(nodmultip[a] >= 2){
      l++;
      }
      }
      if(l == 2){
      //  vertex node
      nodeidentif[i] = 3;
      }
      else{
      //  edge node (not vertex node)
      nodeidentif[i] = 2;
      }
      }*/
    // nove
    if(nodmultip[i] == 2){
      m = 0;
      b = 0;
      // fprintf(out,"\nnode %ld multip = %ld\n",i+1,nodmultip[i]);
      for (j = 0;j < top->nadjacnodesedge[i]; j++){
	l = 0;
	for(k = 0; k < top->n_nned[i][j]; k++){
	  a = top->adjacnodesedge[i][j][k];
	  if(a != i){
	    // fprintf(out,"a = %ld  multip = %ld\n",a+1,nodmultip[a]);
	    if(nodmultip[a] == 2){
	      l++;
	    }
	    if(nodmultip[a] > 2){
	      b++;
	    }
	  }
	}
	// fprintf(out,"l = %ld\n",l);
	if(l == 1){
	  m++;
	}
      }
      // fprintf(out," m = %ld\n",m);
      if(m == 1 && b == 0){
	//vertex node
	nodeidentif[i] = 3;
      }
      else{
	//edge node
	nodeidentif[i] = 2;
      }
    }
    if(nodmultip[i] > 2){
      m = 0;
      b = 0;
      n = 0;
      // fprintf(out,"\nnode %ld multip = %ld\n",i+1,nodmultip[i]);
      for (j = 0;j < top->nadjacnodesedge[i]; j++){
	l = 0;
	for(k = 0; k < top->n_nned[i][j]; k++){
	  a = top->adjacnodesedge[i][j][k];
	  if(a != i){
	      // fprintf(out,"a = %ld  multip = %ld\n",a+1,nodmultip[a]);
	    if(nodmultip[i] == nodmultip[a]){
	      l++;
	    }
	    if( nodmultip[i] < nodmultip[a] ){
	      b++;
	    }
	    if(nodmultip[i] > nodmultip[a]){
	      n++;
	    }
	  }
	}
	// fprintf(out,"l = %ld b = %ld n = %ld\n",l,b,n);
	if( l == 1){
	  m++;
	}
      }
      // lezi na krivce - m == 1 && b == 0 bod dotyku n ==  nadjacnodesedge[i]
      // fprintf(out," m = %ld\n",m);
      if((m == 1  && b == 0) || n == top -> nadjacnodesedge[i] ){
	//vertex node
	nodeidentif[i] = 3;
      }
      else{
	//edge node
	nodeidentif[i] = 2;
      }
      
      
    }
  }
 
  
  //tisk do out
  for (i = 0; i < nn; i++){
    if(nodeidentif[i] == 0) fprintf(out,"\n node number %ld is not identification",i+1);
    if(nodeidentif[i] == 1) fprintf(out,"\n node number %ld is internal node",i+1);
    if(nodeidentif[i] == 2) fprintf(out,"\n node number %ld is edge node",i+1);
    if(nodeidentif[i] == 3) fprintf(out,"\n node number %ld is vertex node",i+1);
  }
  
  // deleting
  for(i = 0; i < nn; i++){
    for(j = 0; j < top->nadjacnodesedge[i]; j++){
      delete []top->adjacnodesedge[i][j];
      }
    delete []top->adjacnodesedge[i];
    delete []top->n_nned[i];
  }
  delete []top->adjacnodesedge;
  delete []top->n_nned;
  delete []top->nadjacnodesedge;
  
}



/**
   function searches for corner nodes in 3D problems
   
   array nodeidentif contains node identification on subdomain
   
   1 - internal node
   2 - edge and surface node
   3 - corner node
   
   @param out - pointer to output file stream
   @param top - pointer to general topology

   08.03.2007 JB
*/
void partop::identification_node_3d(gtopology *top,long *ltg,long *domproc, FILE *out)
{
  long type;
  
  type = give_whole_elemtype(top);
  switch(type){
  case  45:
    fprintf(out,"\n\n\nUsed element into problem is hexahedron\n\n",type);
    identification_node_hexa (top,ltg,domproc,out);
    break;
  case 46:
    fprintf(out,"\n\n\nUsed element into problem is hexahedron\n\n",type);
    identification_node_hexa (top,ltg,domproc,out);
    break;
  case 41:
    fprintf(out,"\n\n\nUsed element into problem is tetrahedron\n\n",type);
    identification_node_tetra (top, domproc,out);
    break;
   case 42:
     fprintf(out,"\n\n\nUsed element into problem is tetrahedron\n\n",type);
     identification_node_tetra (top, domproc,out);
     break; 
  default:
    fprintf (stderr,"\n\n unknown element type is required in function identification_node_3d (file %s, line %d).\n",__FILE__,__LINE__);
    break;
  }
}





/*
   array nodeidentif contaiting node identification on subdomain

   function only for brick elements
   

   1 - internal node
   2 - edge node
   3 - vertex node
   4 - face node

   @param out - pointer to output file stream
   @param top - pointer to general topology
JB
*/
void partop::identification_node_hexa (gtopology *top,long *ltg,long *domproc, FILE *out)
{
  long a,i,j,k,l,m,b,n;
  
  if (nodeidentif != NULL)
    delete []nodeidentif;
  nodeidentif = new long [nn];
  
  for (i = 0; i < nn; i++)
    nodeidentif[i] = 0; 

  
  top->adjacnodes_edge (out);
  
  // ********
  // trideni
  // ********
  
  
  for (i=0;i<nn;i++){
    
    if(nodmultip[i] == 1){
      //internal node
      nodeidentif[i] = 1; 
    }
    // stare a celkem funkcni
    /*
      if(nodmultip[i] == 2){
      m = 0;
      for (j = 0;j < top->nadjacnodessurf[i]; j++){
      l = 0;
      for(k = 0; k < top->n_nnsurf[i][j]; k++){
      a = top->adjacnodessurf[i][j][k];
      if(nodmultip[a] >= 2){
      l++;
      }
      }
      //fprintf(out,"   nn =  %ld l = %ld",i+1,l);
      if(l == 4){
      m++;
      }
      
      }
      //fprintf(out," m = %ld\n",m);
      if(m == 1 && top -> nadjacnodessurf[i] == 3 ){
      //vertex node
      nodeidentif[i] = 3;
      }
      else{
      //edge node
      nodeidentif[i] = 2;
      //if( top -> nadjacnodessurf[i] == 5){
      //edge node
      //nodeidentif[i] = 2;
      //}
      //else{
      //face node
      //nodeidentif[i] = 4; 
      //}
      }
      }
      if(nodmultip[i] >= 3){
      l=0;
      for (j=0;j<top->nadjnodnod[i];j++){
      a=top->adjnodnod [i][j];
      if(nodmultip[a] == nodmultip[i]){
      l++;
      }
      }
      if(l == 2){
      //vertex node
      nodeidentif[i]=3;
      }
      else{
      //edge node 
      nodeidentif[i]=2;
      }
      }*/
    
    // nove
    if(nodmultip[i] == 2){
      m = 0;
      b = 0;
      //fprintf(out,"\nnode %ld multip = %ld\n",i+1,nodmultip[i]);
      for (j = 0;j < top->nadjacnodesedge[i]; j++){
	l = 0;
	for(k = 0; k < top->n_nned[i][j]; k++){
	  a = top->adjacnodesedge[i][j][k];
	  if(a != i){
	    //fprintf(out,"a = %ld  multip = %ld\n",a+1,nodmultip[a]);
	    if(nodmultip[a] == 2){
	      l++;
	    }
	    if(nodmultip[a] > 2){
	      b++;
	    }
	  }
	}
	//fprintf(out,"l = %ld\n",l);
	if(l == 1){
	  m++;
	}
      }
      //fprintf(out," m = %ld\n",m);
      if((m == 2 && top->nadjacnodesedge[i] == 3) && b == 0){
	//vertex node
	nodeidentif[i] = 3;
      }
      else{
	//edge node
	nodeidentif[i] = 2;
      }
    }
    if(nodmultip[i] >= 3){
      m = 0;
      b = 0;
      n = 0;
      //fprintf(out,"\nnode %ld multip = %ld\n",i+1,nodmultip[i]);
      for (j = 0;j < top->nadjacnodesedge[i]; j++){
	l = 0;
	for(k = 0; k < top->n_nned[i][j]; k++){
	  a = top->adjacnodesedge[i][j][k];
	  if(a != i){
	    //fprintf(out,"a = %ld  multip = %ld\n",a+1,nodmultip[a]);
	    if(nodmultip[i] == nodmultip[a]){
	      l++;
	    }
	    if( nodmultip[i] < nodmultip[a] ){
	      b++;
	    }
	    if(nodmultip[i] > nodmultip[a]){
	      n++;
	    }
	  }
	}
	if( l == 1){
	  m++;
	}
      }
      // lezi na krivce - m == 1 && b == 0 bod dotyku n ==  nadjacnodesedge[i]
      if((m == 1  && b == 0) || n == top -> nadjacnodesedge[i] ){
	//vertex node
	nodeidentif[i] = 3;
      }
      else{
	//edge node
	nodeidentif[i] = 2;
      }
    }
    
  }
  //vypis
  for (i=0;i<nn;i++){
    if(nodeidentif[i]==0) fprintf(out,"\n node number %ld is not identification",i+1);
    if(nodeidentif[i]==1) fprintf(out,"\n node number %ld is internal node",i+1);
    if(nodeidentif[i]==2) fprintf(out,"\n node number %ld is edge node",i+1);
    if(nodeidentif[i]==4) fprintf(out,"\n node number %ld is face node",i+1);
    if(nodeidentif[i]==3) fprintf(out,"\n node number %ld is vertex node",i+1);
  }
  
  // deleting
  for(i = 0; i < nn; i++){
    for(j = 0; j < top->nadjacnodesedge[i]; j++){
      delete []top->adjacnodesedge[i][j];
    }
    delete []top->adjacnodesedge[i];
    delete []top->n_nned[i];
  }
  delete []top->adjacnodesedge;
  delete []top->n_nned;
  delete []top->nadjacnodesedge;
  
  
}


/*
   array nodeidentif contaiting node identification on subdomain

   function only for tetrahedron elements
   

   1 - internal node
   2 - edge node
   3 - vertex node
   4 - face node

   @param out - pointer to output file stream
   @param top - pointer to general topology
JB

nepracuje spolehlive. rozdil mezi oznacenim vrcholu na jednotlivych domenach
*/

void partop::identification_node_tetra (gtopology *top, long *domproc, FILE *out)
{
  long i,j,k,m,l,a,b,n;
  
  if (nodeidentif != NULL)
    delete [] nodeidentif;
  nodeidentif = new long [nn];
 
  for (i=0;i<nn;i++) 	nodeidentif[i] = 0; 
  
  top->adjacnodes_edge(out);
  
  // ********
  // trideni
  // ********
  
  for (i=0;i<nn;i++){
    if(nodmultip[i] == 1){
      //internal node
      nodeidentif[i]=1; 
    }
    
    // nove a pokusne
    if(nodmultip[i] == 2){
      m = 0;
      b = 0;
      n = 0;
      //fprintf(out,"\nnode %ld multip = %ld\n",i+1,nodmultip[i]);
      for (j = 0;j < top->nadjacnodesedge[i]; j++){
	l = 0;
	for(k = 0; k < top->n_nned[i][j]; k++){
	  a = top->adjacnodesedge[i][j][k];
	  if(a != i){
	    //fprintf(out,"a = %ld  multip = %ld\n",a+1,nodmultip[a]);
	    if(nodmultip[a] == 2){
	      l++;
	    }
	    if(nodmultip[a] > 2){
	      b++;
	    }
	    if(nodmultip[a] < 2){
	      n++;
	    }
	  }
	}
	//fprintf(out,"l = %ld b = %ld n = %ld\n",l,b,n);
	if(l == 1){
	    m++;
	}
      }
      //fprintf(out," m = %ld\n",m);
      if(m == 2 && b == 0 ){
	//vertex node
	nodeidentif[i] = 3;
      }
      else{
	//edge node
	nodeidentif[i] = 2;
      }
    }
    if(nodmultip[i] >= 3){
      m = 0;
      b = 0;
      n = 0;
      //fprintf(out,"\nnode %ld multip = %ld\n",i+1,nodmultip[i]);
      for (j = 0;j < top->nadjacnodesedge[i]; j++){
	l = 0;
	for(k = 0; k < top->n_nned[i][j]; k++){
	  a = top->adjacnodesedge[i][j][k];
	  if(a != i){
	    //fprintf(out,"a = %ld  multip = %ld\n",a+1,nodmultip[a]);
	    if(nodmultip[i] == nodmultip[a]){
	      l++;
	    }
	    if( nodmultip[i] < nodmultip[a] ){
	      b++;
	    }
	    if(nodmultip[i] > nodmultip[a]){
	      n++;
	    }
	  }
	}
	//fprintf(out,"l = %ld b = %ld n = %ld\n",l,b,n);
	if( l == 1){
	  m++;
	}
      }
      //fprintf(out," m = %ld\n",m);
      // lezi na krivce - m == 1 && b == 0 bod dotyku n ==  nadjacnodesedge[i]
      if((m == 1  && b == 0) || n == top -> nadjacnodesedge[i] ){
	//vertex node
	nodeidentif[i] = 3;
      }
      else{
	//edge node
	nodeidentif[i] = 2;
      }
    }
  }

  //vypis
  for (i=0;i<nn;i++){
    if(nodeidentif[i]==0) fprintf(out,"\n node number %ld is not identification",i+1);
    if(nodeidentif[i]==1) fprintf(out,"\n node number %ld is internal node",i+1);
    if(nodeidentif[i]==2) fprintf(out,"\n node number %ld is edge node",i+1);
    if(nodeidentif[i]==4) fprintf(out,"\n node number %ld is face node",i+1);
    if(nodeidentif[i]==3) fprintf(out,"\n node number %ld is vertex node",i+1);
  }
  
  // deleting
  for(i = 0; i < nn; i++){
    for(j = 0; j < top->nadjacnodesedge[i]; j++){
      delete []top->adjacnodesedge[i][j];
    }
    delete []top->adjacnodesedge[i];
    delete []top->n_nned[i];
  }
  delete []top->adjacnodesedge;
  delete []top->n_nned;
  delete []top->nadjacnodesedge;
  
}

void partop::control_numbers_of_vertices (gtopology *top,FILE *out)
{ 
  long i,j,k,l,m,a,b,dim;
  long n,nom,denom,nc;
  long *vertex,**vcomb,*line;
  double u1,u2,u3,v1,v2,v3,x1,x2,x3,y1,y2,y3,z1,z2,z3,scal,u,v,cosalpha;
  // ***********************
  // kontrola poctu vertexu
  // ***********************
  // n - number of vertices
  n = 0;
  for (i=0;i<nn;i++){
    if(nodeidentif[i]==3) n++;
  }
  
  fprintf(out,"\n\n\nNumber of vertices is %ld\n\n",n);
  
  dim = give_whole_dim(top);
   
  if(dim == 3){
    if( n < 3){
      fprintf (stderr,"Number of vertices is < 3\n");
      fprintf (stderr,"Problem in function control_numbers_of_vertices  %s, line %d).\n",__FILE__,__LINE__);
    }  
  }
  if(dim == 2){ 
    if( n < 2){
      fprintf (stderr,"Number of vertices is < 3\n");
      fprintf (stderr,"Problem in function control_numbers_of_vertices  %s, line %d).\n",__FILE__,__LINE__);
    } 
  }
  
  // list of vertices
  vertex = new long[n];
  m = 0;
  for (i=0;i<nn;i++){
    if(nodeidentif[i]==3){
      vertex[m] = i;
      //fprintf (out,"%ld %ld\n",m+1,vertex[m]);
      m++;
    }
  }
  
  /*
    for(i = 0; i < n; i++){
    fprintf(out,"\n node number %ld has coordinates %le %le %le\n",vertex[i]+1,top->gnodes[vertex[i]].x,top->gnodes[vertex[i]].y,top->gnodes[vertex[i]].z) ;
    }
  */
  
  if (n > 3){
    // **************************************
    // kontrola zda nelezi vertices na primce
    // **************************************
    // number of combination
    nom = 1;
    for(i = 1; i < n+1; i++){
      nom*=i;
    }
    denom = 1;
    for(i = 1; i < n-2; i++){
      denom*=i;
    }
    nc = nom/(denom*6);
    fprintf(out,"Pocet kombinaci k vyzkouseni je %ld\n",nc);
    
    // list of combination of verticies
    vcomb = new long*[nc];
    for(i = 0; i < nc; i++ ){
      vcomb[i] = new long[3];
    }
    line = new long[nc];
    
    
    b = 0;
    a = n-2;
    for(i = 0; i < n-2; i++){
      //fprintf (out,"i: %ld node %ld\n",i,vertex[i]);
      // node 1
      x1 = top->gnodes[vertex[i]].x;
      y1 = top->gnodes[vertex[i]].y;
      z1 = top->gnodes[vertex[i]].z;
      l = i+1;
      for(j = 0; j < a; j++){
	//fprintf (out,"j: %ld node %ld\n",j,vertex[l]);
	// node 2
	x2 = top->gnodes[vertex[l]].x;
	y2 = top->gnodes[vertex[l]].y;
	z2 = top->gnodes[vertex[l]].z;
	// vecotor u - nodes 1 and 2
	u1 = x1 - x2;
	//fprintf(out,"u1 = %le\n",u1);
	u2 = y1 - y2;
	//fprintf(out,"u2 = %le\n",u2);
	u3 = z1 - z2;
	//fprintf(out,"u3 = %le\n",u3);
	m = l+1;
	for(k = 0; k < a-j; k++){
	  //fprintf (out,"k : %ld node %ld\n",k,vertex[m]);
	  vcomb[b][0]= vertex[i];
	  vcomb[b][1]= vertex[l];
	  vcomb[b][2]= vertex[m];
	  // node 3
	  x3 = top->gnodes[vertex[m]].x;
	  y3 = top->gnodes[vertex[m]].y;
	  z3 = top->gnodes[vertex[m]].z;	
	  // vector v - nodes 1 and 3
	  v1 = x1 - x3;
	  //fprintf(out,"v1 = %le\n",v1);
	  v2 = y1 - y3;
	  //fprintf(out,"v2 = %le\n",v2);
	  v3 = z1 - z3;
	  //fprintf(out,"v3 = %le\n",v3);
	  // scalar product
	  scal = u1*v1+u2*v2+u3*v3;
	  //fprintf(out,"scal = %le\n",scal);
	  //  magnitude of vectoru u
	  u = u1*u1+u2*u2+u3*u3;
	  //fprintf(out,"u = %le\n",u);
	  u = sqrt(u);
	  //fprintf(out,"u = %le\n",u);
	  // magnitude of vectoru v
	  v = v1*v1+v2*v2+v3*v3;
	  //fprintf(out,"v = %le\n",v);
	  v = sqrt(v);
	  //fprintf(out,"v = %le\n",v);
	  // angle
	  cosalpha = fabs(scal/(u*v));
	  //fprintf(out,"cosalpha = %le\n",cosalpha);
	  if(cosalpha == 1.0000000){
	    // body 1 2 3 lezi na primce
	    line[b] = 1;
	    b++;
	  }
	  else{
	    line[b] = 0;
	    b++;
	  }
	  m++;
	}
	l++;
      }
      a--;
    }
    
    b = 0;
    for(i = 0; i < nc; i++){
      if(line[i] == 0){
	fprintf(out,"body %ld %ld %ld nelezi na primce\n",vcomb[i][0]+1,vcomb[i][1]+1,vcomb[i][2]+1);
      }
      else{
	b++;
	fprintf(out,"body %ld %ld %ld lezi na primce\n",vcomb[i][0]+1,vcomb[i][1]+1,vcomb[i][2]+1);
      }
    }
    fprintf(out," %ld kombinaci lezi na primce\n",b);
  }
  
  delete []vertex; 
}

void partop::control_vertices_on_master(gtopology *top, long *domproc, FILE *out)
{ 
  long i,j,k,m,n;
  long buffsize;
  long *buff,*aux,*auxnn;
  long **globnodeidentif;
  MPI_Status stat;
   
  // *************************
  // kontola vertexu na master
  // *************************
  buffsize = maxnbn;
  buff = new long[buffsize];

  j = 0;
  for(i = 0; i < nn; i++ ){
    if( nodeidentif[i] > 1){
      buff[j] = nodeidentif[i];
      j++;
    }
  }
  
  if(myrank == 0){
    
    globnodeidentif = new long*[nproc];
   
    // master contribution
    j=domproc[0]; 
    //printf(out,"nbnd[%ld] =  %ld\n",j,nbnd[j]);
    
    globnodeidentif[j] = new long[nbnd[j]];
    for(i = 0; i < nbnd[0]; i++){
      globnodeidentif[j][i] = buff[i];
    }
    
    
    // slaves contribution
    for (k=1;k<nproc;k++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG]; 
      globnodeidentif[j] = new long[nbnd[j]];
      for(i = 0; i < nbnd[j]; i++){
	globnodeidentif[j][i] = buff[i];
      }
    }
    
    // tisk kontrolni
    // for (k = 0;k < nproc; k++){ 
    // fprintf(out,"Domain number %ld\n",k+1);
    // for(i = 0; i < nbnd[k]; i++){
    // fprintf(out,"%ld %ld\n",i+1,globnodeidentif[k][i]);
    // }
    // }
    
    // kontrola
    aux = new long[nproc];
    auxnn = new long[nproc];
    for(i = 0; i < tnbn; i++){
      //fprintf(out,"node %ld X  ",i+1); 
      for(j = 0; j < nproc; j++){
	aux[j] = 0;
	auxnn[j] = -1;
	for(k = 0; k < nbnd[j]; k++){
	  if(cnbn[j][k] == i){
	    aux[j] = globnodeidentif[j][k];
	    auxnn[j] = k;
	    break;
	  }
	}
      }
      
      //for(j = 0; j < nproc; j++){	
      //fprintf(out,"%ld %ld   ",auxnn[j],aux[j]);	
      //}
      //fprintf(out,"\n");
      
      // corner node
      m = 0;
      // edge and surface node
      n = 0;
      for(j = 0; j < nproc; j++){
	if(aux[j] == 2){
	  n++;
	}
	if(aux[j] == 3){
	  m++;
	}
      }
      if( m >= n){
	//fprintf(out,"node %ld\n\n\n\n",i+1); 
	for(j = 0; j < nproc; j++){
	  if(aux[j] != 0 && auxnn[j] != -1){
	    globnodeidentif[j][auxnn[j]] = 3;
	  }
	}
      }
      else{
	for(j = 0; j < nproc; j++){
	  if(aux[j] != 0 && auxnn[j] != -1){
	    globnodeidentif[j][auxnn[j]] = 2;
	  }
	}
      }
    }
    delete []aux;
    delete []auxnn;
    
    //tisk kontrolni
    //for (k = 0;k < nproc; k++){ 
    //fprintf(out,"Domain number %ld\n",k+1);
    //for(i = 0; i < nbnd[k]; i++){
    //fprintf(out,"%ld %ld\n",i+1,globnodeidentif[k][i]);
    //}
    //}
    
    // for slaves
    for (i = 1;i < nproc; i++){
      for(k = 0; k < nbnd[i]; k++){
	buff[k] = globnodeidentif[i][k];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    // For master
    for(k = 0; k < nbnd[0]; k++){
      buff[k] = globnodeidentif[0][k];
    }
    
    // mazani
    for (k=0;k<nproc;k++){ 
      delete []globnodeidentif[k];
    }
    delete []globnodeidentif;
    
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  for(i = 0; i < nbn; i++){
    if(buff[i] == nodeidentif[lnbn[i]]){
      nodeidentif[lnbn[i]] = buff[i];
    }
    else{
      fprintf(out,"pozor zmena\n");
      nodeidentif[lnbn[i]] = buff[i];
    }
  }
  delete []buff;

}

/*
  corner detection 

08.05.07  JB
*/
void partop::corner_detection(gtopology *top, long *domproc, long *ltg, FILE *out)
{
  long i,j,k;
  long buffsize,c,b;
  long *buff,*master_corner;
  MPI_Status stat;
  
  // node identification
  identification_node(top,ltg,domproc,out);
  control_numbers_of_vertices (top,out);  
  control_vertices_on_master(top,domproc,out);
 
    
  buffsize = maxnbn;
  buff = new long[buffsize];
  
  //selection of boundary nodes on subdomain
  k = 0;
  for(i = 0; i < nn; i++){
    if(nodeidentif[i] > 1){
      buff[k] = nodeidentif[i];
      k++;
      //fprintf(out,"%ld\n",i+1);
    }
  }
  
   
  //  numbering of corner and edge nodes on master
  if(myrank==0){
    //fprintf (stderr,"tnbn = %ld\n",tnbn);
    master_corner = new long[tnbn];
    
    for(i = 0; i < tnbn; i++ ) master_corner[i] = -1;
    
    // master contributions
    j=domproc[0]; 
    for(k = 0; k < nbnd[j]; k++){
      master_corner[cnbn[j][k]]=buff[k];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG]; 
      for(k = 0; k < nbnd[j]; k++){
	if(master_corner[cnbn[j][k]] == -1){
	  master_corner[cnbn[j][k]]= buff[k];
	}
	else{
	  if(master_corner[cnbn[j][k]] != buff[k]){
	    fprintf (out,"\n %ld  %ld  %ld different nodeidentif on the %ld-th domain in the %ld-th node",cnbn[j][k],buff[k],master_corner[cnbn[j][k]],j+1,k);
	    fprintf (stderr,"different nodeidentif on the %ld-th domain in the %ld-th node\n",j,k);
	    fprintf (stderr,"\n in function corner_detection (file %s, line %d)\n",__FILE__,__LINE__); 
	  }
	}
      }
    }
    
    // fprintf(out,"\n\ncorner nodes\n");
    // for(k = 0; k < tnbn; k++){
    // fprintf(out,"%ld %ld\n",k+1,master_corner[k]);
    // }
    
    // numbering of corner and edge nodes on master
    c = 1;
    b = 1;
    for(i = 0; i < tnbn; i++){
      if(master_corner[i] == 3){
	master_corner[i] = c*(-1);
	c++;
      }
      if(master_corner[i] == 2){
	master_corner[i] = b;
	b++;
      }	
      //fprintf(out,"%ld %ld\n",i+1,master_corner[i]);
    }
    
    fprintf(out,"\n\n\n Number of corner nodes in whole problem is %ld\n",c-1);
    fprintf(out,"\n\n\n Number of boundary nodes in whole problem is %ld\n",b-1);
    
    if((c-1)+(b-1) != tnbn){
      fprintf (stderr,"Problem in function corner_detection - the number of corner and edge nodes is not equal to total number of boundary nodes (file %s, line %d).\n",__FILE__,__LINE__);
      
    }
    
    // For Slaves
    for (i = 1;i < nproc; i++){
      for(k = 0; k < nbnd[i]; k++){
	buff[k] =  master_corner[cnbn[i][k]];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    // For master
    for(k = 0; k < nbnd[0]; k++){
      buff[k] =  master_corner[cnbn[0][k]];
    }
    delete []master_corner;
  }
  
  // slave
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  // fprintf(out,"\n\ncorner nodes - buff\n");
  // for(i = 0; i < nbn; i++){
  // fprintf(out,"%ld %ld\n",i+1,buff[i]);
  // }
  
  
  fprintf(out,"\n\n\n array corneridentif\n"); 
  
  k = 0;
  for(i = 0; i < nn; i++){
    if(nodeidentif[i] == 1 || nodeidentif[i] == 0){
      ltg[i] = 0;
    }  
    else{
      ltg[i] = buff[k];
      k++;
    }
    fprintf(out,"%ld %ld\n",i+1,ltg[i]);

    ltg[i]--;
  }
  
  delete []buff;
  delete []nodeidentif;

}






/**
   function
   
   param top - pointer to gtopology
   
   10.5.2007, JB
*/
void partop::assemble_gnn(gtopology *top,long *domproc,FILE *out)
{
  long i,j,k;
  long a,b,buffsize;
  long *buff;
  MPI_Status stat;
  /*
  if (md==bound_nodes){
    // creation of array allnodes
    
    buffsize = maxnbn;
    buff = new long [buffsize];
    
    for(i = 0; i < nbn; i++){
      buff[i] = lnbn[i];
    }
    
    // fprintf(out,"\n");
    // for(i = 0; i < nbn; i++){
    // fprintf(out,"%ld buff %ld\n",i+1,buff[i]);
    // }
    
    
    if(myrank == 0){
      if (allnodes != NULL){
	for (i=0;i<nproc;i++){
	  delete [] allnodes[i];
	}
	delete [] allnodes;
      }
      allnodes = new long *[nproc];

      for(i = 0; i< nproc; i++){
	allnodes[i] = new long [nnsd[i]];
	for(j = 0; j < nnsd[i]; j++){
	  allnodes[i][j] = -1;
	}
      }
      
      // master contribution
      for(k = 0; k < nbnd[0]; k++){
	allnodes[0][buff[k]]=cnbn[0][k];
      }
      a = tnbn;
      for(k = 0; k < nnsd[0]; k++){
	if(allnodes[0][k] == -1){
	  allnodes[0][k] = a;
	  a++;
	}
      }
      
      
      //  slave contributions
      for (i=1;i<nproc;i++){
	MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	j=domproc[stat.MPI_TAG];  
	for(k = 0; k < nbnd[0]; k++){
	  allnodes[j][buff[k]]=cnbn[j][k];
	}
	for(k = 0; k < nnsd[j]; k++){
	  if(allnodes[j][k] == -1){
	    allnodes[j][k] = a;
	    a++;
	  }
	}
      }
      for (i=0;i<nproc;i++){  
	fprintf(out,"\n\nDomain number %ld\n",i+1);
	for(j = 0; j < nnsd[i]; j++){
	  fprintf(out,"%ld %ld\n",j+1,allnodes[i][j]+1);
	}
      }
      tnnp = a;
      fprintf(out,"tnnp is %ld\n",tnnp); 
      
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    delete []buff;
  }
*/

  buffsize = maxnn;
  buff = new long [buffsize];
  
  // detection of numbers of DOFs on nodes on each subdomain
  for (i=0;i<nn;i++){
    buff[i]=top -> give_ndofn(i);
  }
  
  if (myrank==0){
    
    //  allocation  of array gnn
    if (gnn == NULL)
      gnn = new long[tnnp+1];
    
    for(i = 0; i < tnnp+1; i++)
      gnn[i] = -10; 
    
    //  contribution from master
    k = domproc[0];
    for (i = 0;i < nnsd[k]; i++){
      gnn[allnodes[k][i]] = buff[i];
    }
    
    for (i= 1 ;i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      k=domproc[stat.MPI_TAG];
      for (j = 0; j < nnsd[k]; j++){
	if(gnn[allnodes[k][j]] == -10){
	  gnn[allnodes[k][j]] = buff[j];
	}
	else{
	  if (gnn[allnodes[k][j]] != buff[j]){
	    fprintf (stderr,"\n different numbers of DOFs on the %ld-th domain in the %ld-th node",k,j);
	    fprintf (stderr,"\n in function assemble_gnn (file %s, line %d)\n",__FILE__,__LINE__);
	  }
	}
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  if (myrank == 0){
    
    // kontrolni tisk
    // fprintf(out,"\n\n\n gnn array\n\n");
    // for(i = 0; i < tnnp; i++){
    // fprintf(out,"%ld %ld\n",i+1,gnn[i]);
    // }
    
    // rewrite of array gnn up to array of address for array pgcn
    a = gnn[0];
    gnn[0] = 0;
    for(i = 1; i < tnnp+1; i++){
      b = gnn[i];
      gnn[i] = gnn[i-1] + a;
      a = b;
    }
    tndof = gnn[tnnp];


    fprintf(out,"\n\n\ngnn array\n\n");
    for(i = 0; i < tnnp+1; i++){
      fprintf(out,"%ld %ld\n",i+1,gnn[i]);
    }
  }
}

/*
  function 
  
  the following functions have to be called before this function:
  void initiation (top,ltg);
  void numbers_of_all_nodes_on_subdomains (top,domproc,out);
  void compute_multiplicity (top,domproc,out);
  void find_boundary_nodes (top,domproc,out);
  void assemble_gnn(gtopology *top,long *domproc,FILE *out);
  
  10.5.2007, JB
*/

void partop::assemble_pgcn(gtopology *top,long *domproc,FILE *out)
{
  long i,ii,j,k,m;
  long adr,ndofn,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize = maxndof;
  buff = new long [buffsize];
  
  m = 0;
  for (i = 0;i < nn; i++){
    ndofn = top->give_ndofn (i);
    for (j = 0; j < ndofn; j++){
      buff[m] = top -> give_dof (i,j);
      m++;
    }
  }
  

  if (myrank==0){
        
    fprintf (stdout,"\n\n\n tndof je %ld\n\n\n",tndof);

    
    
    
    //  allocation  of array pgcn
    if (pgcn == NULL){
      pgcn = new long[tndof];
    }
    
    
    for(i = 0; i < tndof; i++) pgcn[i] = -10;
    
    //  contribution from master
    k = domproc[0];
    m = 0;
    for (i = 0;i < nnsd[k]; i++){
      adr = gnn[allnodes[k][i]];
      //fprintf(out,"adr = %ld  gnn[allnodes[k][adr]+1] = %ld\n", adr, gnn[allnodes[k][i]+1]);
      ndofn = gnn[allnodes[k][i]+1] - adr;
      //fprintf(out,"ndofn = %ld\n", ndofn); 
      for (j = 0;j < ndofn; j++){
	pgcn[adr+j] = buff[m];
	m++;
      }
    }
    
    // slave contributions
    for (i= 1 ;i < nproc; i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      k = domproc[stat.MPI_TAG];
      m = 0;
      for (ii = 0;ii< nnsd[k]; ii++){
	adr = gnn[allnodes[k][ii]];
	ndofn = gnn[allnodes[k][ii]+1] - adr;
	for (j = 0;j < ndofn; j++){
	  pgcn[adr+j] = buff[m];
	  m++;
	}
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete []buff;

  /*
  // kontrolni tisk pgcn
  if(myrank == 0){
    fprintf(out,"\n\n\n Kontrola pgcn\n\n");
    for(i = 0; i < tndof; i++){
      fprintf(out,"%ld %ld\n",i+1,pgcn[i]);
    }
  }
  */


  // global code numbers
  if(myrank == 0){
    m = 1;
    for(i = 0; i < tndof; i++){
      if(pgcn[i] > 0){
	pgcn[i] = m;
	m++;
      }
    }
    //tisk
    fprintf(out,"\n\n\n Global code numbers\n\n");
    for(i = 0; i < tndof; i++){
      fprintf(out,"%ld %ld\n",i,pgcn[i]);
    }
    tndof = m-1;
    fprintf(out,"\n\n\n number of DOFs on whole problem %ld\n\n",tndof);
  }
  

  
}

/**
   function assembles global code numbers on master processor
   
   the following functions have to be called before this function:
   void initiation (top,ltg);
   void numbers_of_all_nodes_on_subdomains (top,domproc,out);
   void compute_multiplicity (top,domproc,out);
   void find_boundary_nodes (top,domproc,out);
   void assemble_gnn(gtopology *top,long *domproc,FILE *out);
   void assemble_pgcn(gtopology *top,long *domproc,FILE *out);
   
   
   function assembles array gcnd
   
   @param top - pointer to gtopology of subdomain
   @param domproc - domain-processor correspondence
   @param out - output file
   
   10.5.2007, JB
*/
void partop::assemble_gcnd(gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,m;
  long ndofn,adr;
  
  if(myrank == 0){
    
    nud = new long[nproc];
    for(i = 0; i < nproc; i++){
      nud[i]=0;
      for(j = 0; j < nnsd[i]; j++){
	adr = gnn[allnodes[i][j]];
	ndofn = gnn[allnodes[i][j]+1] - adr;
	for (k = 0;k < ndofn; k++){
	  if(pgcn[adr+k] > 0){
	    nud[i]++;
	  }
	}
      }
    }
    
    fprintf(out,"\n\n\n Number of unknowns on subdomain\n");
    for(i = 0; i < nproc; i++){
      fprintf(out,"Domain %ld %ld\n",i+1,nud[i]);
    }
    
    gcnd = new long*[nproc];
    for(i = 0;i < nproc; i++){
      gcnd[i] = new long[nud[i]];
      m = 0;
      for(j = 0; j < nnsd[i]; j++){
	adr = gnn[allnodes[i][j]];
	ndofn = gnn[allnodes[i][j]+1] - adr;
	for (k = 0;k < ndofn; k++){
	  if(pgcn[adr+k] > 0){
	    gcnd[i][m] = pgcn[adr+k] ;
	    m++;
	  }
	}
      } 
    }
    
    // print gcnd
    fprintf(out,"\n\n\n Unknowns on subdomain gcnd\n");
    for(i = 0;i < nproc; i++){
      fprintf(out,"Domain %ld\n",i+1);
      for(j = 0; j < nud[i]; j++){
	fprintf(out,"%ld %ld\n",j,gcnd[i][j]);	
      }
    }	
    
   
    delete []gnn;
    delete []pgcn;
  }
  
}


