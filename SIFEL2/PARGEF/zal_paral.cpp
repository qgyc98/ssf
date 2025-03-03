#include "paral.h"

paral::paral(int np,int mr,int nd)
{
  nproc=np;  myrank=mr;  ndom=nd;
  domproc=NULL;

  ltg = NULL;
  gcn=NULL;
  lcngcn=NULL;
  inc=NULL;
  nodnum=NULL;
  
  nbndom=NULL;
  nrdofdom=NULL;
  masgcn=NULL;
  nrbmdom=NULL;
  rbmadr=NULL;
}

paral::~paral()
{
  delete [] domproc;

  delete [] ltg;
  delete [] gcn;
  delete [] lcngcn;
  delete [] inc;
  delete [] nodnum;

  delete [] nbndom;
  delete [] nrdofdom;
  delete [] masgcn;
  delete [] nrbmdom;
  delete [] rbmadr;
}

/**
   function reads local to global correspondence of nodes
   
   JK, 20.1.2002
*/
void paral::read (FILE *in,gtopology *top)
{
  long i;
  
  ltg = new long [top->nn];
  
  for (i=0;i<top->nn;i++){
    fscanf (in,"%ld",ltg+i);
    ltg[i]--;
  }
}

/**
   function establishes correspondence among processors and subdomains
   
   JK
*/
void paral::procdomcorr ()
{
  long i,j;
  MPI_Status stat;
  
  if (myrank==0){
    domproc = new long [nproc];
    domproc[0]=ndom;
    for (i=1;i<nproc;i++){
      MPI_Recv (&j,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      domproc[j]=stat.MPI_TAG;
    }
  }
  else{
    MPI_Send (&ndom,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
}

/**
   function orders unknowns on subdomain with respect of Schur complement method
   inner variables are ordered at the beginning
   boundary variables are ordered at he end
   
   @param top - pointer to topology
   
   JK, 26.11.2002
*/
void paral::schurordering (gtopology *top)
{
  long i,j,k,ndofn;
  
  //  contributions from nodes
  ndof=1;
  for (i=0;i<top->nn;i++){
    if (ltg[i]==-1){
      //  inner nodes
      ndofn=top->give_ndofn (i);
      for (j=0;j<ndofn;j++){
	k=top->give_dof(i,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k>0){
	  top->save_dof (i,j,ndof);  ndof++;
	}
      }
    }
  }
  indof=ndof-1;
  
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      //  boundary nodes
      ndofn=top->give_ndofn (i);
      for (j=0;j<ndofn;j++){
	k=top->give_dof (i,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k>0){
	  top->save_dof (i,j,ndof);  ndof++;
	}
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
  
}

/**
   function creates code numbers of global nodes

   function is called from primal domain decomposition
   function detects maximum number of reduced DOFs on subdomains
   function collects local to global code numbers on master processor

   JK, 30.1.2002
*/
void paral::globcnnum_pdd (gtopology *top,FILE *out)
{
  long i,j,k,l,gnn,maxndofn,maxnbn,maxgnn;
  long *buff,**gncn;
  MPI_Status stat;
  
  //***************************************************************************
  //  determination of maximum number of nodes and maximum number of nodal DOFs
  //***************************************************************************
  //  maximum number of DOFs on node
  maxndofn=0;
  //  number of boundary nodes on subdomain
  nbn=0;
  //  maximum global number of node
  gnn=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      if (maxndofn<top->give_ndofn (i))  maxndofn = top->give_ndofn (i);
      nbn++;
      if (gnn<ltg[i])  gnn=ltg[i];
    }
  }
  
  buff = new long [4];
  buff[0]=maxndofn;
  buff[1]=nbn;
  buff[2]=gnn;
  buff[3]=ndof-indof;
  
  if (myrank==0){
    nbndom = new long [nproc];
    nrdofdom = new long [nproc];
    
    //  master contributions
    nbndom[ndom]=nbn;
    nrdofdom[ndom]=ndof-indof;

    //  maximum number of degrees of freedom on node
    totmaxndofn=maxndofn;
    //  maximum number of boundary nodes
    maxnbn=nbn;
    //  maximum global number of node
    maxgnn=gnn;
    //  maximum number of boundary unknowns on subdomain
    maxnrdof=0;
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,4,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      nbndom[stat.MPI_TAG]=buff[1];
      nrdofdom[stat.MPI_TAG]=buff[3];
      if (buff[0]>totmaxndofn)  totmaxndofn=buff[0];
      if (buff[1]>maxnbn)       maxnbn=buff[1];
      if (buff[2]>maxgnn)       maxgnn=buff[2];
      if (buff[3]>maxnrdof)     maxnrdof=buff[3];
    }
    maxgnn++;

    fprintf (out,"\n\n\n Global data information \n");
    fprintf (out,"\n maxnbn %ld",maxnbn);
    fprintf (out,"\n totmaxndofn %ld",totmaxndofn);
    fprintf (out,"\n maxgnn %ld",maxgnn);
    fprintf (out,"\n maxnrdof %ld",maxnrdof);
    fprintf (out,"\n\n pocty velicin na podoblastech");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nbn %4ld  nrdof %4ld",nbndom[i],nrdofdom[i]);
    }

    buff[0]=totmaxndofn;
    buff[1]=maxnbn;
    buff[2]=maxnrdof;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,4,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,4,MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,4,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    totmaxndofn=buff[0];
    maxnbn=buff[1];
    maxnrdof=buff[2];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;



  //*******************************************
  //  gathering of indicators of DOFs to master 
  //*******************************************
  //  allocation of array buff
  buff = new long [maxnbn*(totmaxndofn+1)];
  for (i=0;i<maxnbn*(totmaxndofn+1);i++){
    buff[i]=0;
  }
  
  //  indicators are placed to array buff
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      buff[k*(totmaxndofn+1)+totmaxndofn]=ltg[i];
      for (j=0;j<top->give_ndofn (i);j++){
	buff[k*(totmaxndofn+1)+j]=top->give_dof (i,j);
      }
      k++;
    }
  }
  
  
  if (myrank==0){
    
    //  allocation of array gncn
    gncn = new long* [maxgnn];
    for (i=0;i<maxgnn;i++){
      gncn[i] = new long [totmaxndofn];
      for (j=0;j<totmaxndofn;j++){
	gncn[i][j]=0;
      }
    }
    
    //  master contribution
    for (j=0;j<nbndom[ndom];j++){
      l=buff[j*(totmaxndofn+1)+totmaxndofn];
      for (k=0;k<totmaxndofn;k++){
	gncn[l][k]=buff[j*(totmaxndofn+1)+k];
      }
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      for (j=0;j<nbndom[stat.MPI_TAG];j++){
	l=buff[j*(totmaxndofn+1)+totmaxndofn];
	for (k=0;k<totmaxndofn;k++){
	  gncn[l][k]=buff[j*(totmaxndofn+1)+k];
	}
      }
    }
    
    //  computation of total number of global unknowns
    ngdof=1;
    for (i=0;i<maxgnn;i++){
      for (j=0;j<totmaxndofn;j++){
	if (gncn[i][j]>0){
	  gncn[i][j]=ngdof;  ngdof++;
	}
      }
    }
    ngdof--;
    
  }
  else{
    MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  //********************************
  //  global code numbers scattering
  //********************************
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      for (j=0;j<nbndom[stat.MPI_TAG];j++){
	l=buff[j*(totmaxndofn+1)+totmaxndofn];
	for (k=0;k<totmaxndofn;k++){
	  buff[j*(totmaxndofn+1)+k]=gncn[l][k];
	}
      }
      l=stat.MPI_TAG;
      MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,l,ndom,MPI_COMM_WORLD);
    }

    k=0;
    for (i=0;i<top->nn;i++){
      if (ltg[i]!=-1){
	buff[k*(totmaxndofn+1)+totmaxndofn]=ltg[i];
	for (j=0;j<top->give_ndofn (i);j++){
	  buff[k*(totmaxndofn+1)+j]=top->give_dof(i,j);
	}
	k++;
      }
    }
    
    for (j=0;j<nbndom[ndom];j++){
      l=buff[j*(totmaxndofn+1)+totmaxndofn];
      for (k=0;k<totmaxndofn;k++){
	buff[j*(totmaxndofn+1)+k]=gncn[l][k];
      }
    }
    
    //  deallocation of array gncn
    for (i=0;i<maxgnn;i++){
      delete [] gncn[i];
    }
    delete [] gncn;
    
  }
  else{
    MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,maxnbn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  


  //*********************************************
  //  local to global code numbers correspondence
  //*********************************************
  //  allocation of array of local to global code numbers correspondence
  gcn = new long [ndof];
  for (i=0;i<ndof;i++){
    gcn[i]=0;
  }
  
  //  filling of local to global code numbers correspondence
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      for (j=0;j<top->give_ndofn (i);j++){
	l=top->give_dof (i,j)-1;
	if (l>-1)  gcn[l]=buff[k*(totmaxndofn+1)+j];
      }
      k++;
    }
  }
  
  delete [] buff;
  
  
  //********************************************************************
  //  filling of array of local to global code numbers of all subdomains
  //********************************************************************
  buff = new long [maxnrdof];
  j=indof;
  for (i=0;i<ndof-indof;i++){
    buff[i]=gcn[j];  j++;
  }
  
  if (myrank==0){
    //  allocation of array masgcn
    masgcn = new long* [nproc];
    for (i=0;i<nproc;i++){
      masgcn[i] = new long [nrdofdom[i]];
      for (j=0;j<nrdofdom[i];j++){
	masgcn[i][j]=0;
      }
    }
    
    //  master contribution
    for (j=0;j<nrdofdom[ndom];j++){
      masgcn[ndom][j]=buff[j];
    }
    
    //  slaves contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnrdof,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      for (j=0;j<nrdofdom[stat.MPI_TAG];j++){
	masgcn[stat.MPI_TAG][j]=buff[j];
      }
      
    }
  }
  else{
    MPI_Send (buff,maxnrdof,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
}

//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************

/**
   function creates code numbers
   
   maxnbn - maximum number of boundary nodes (each processor)
   pcor->totmaxndofn - maximum number of DOFs in one node (each processor)
   maxgnn - maximum node number in the problem (only on master)
   maxinc - maximum number of node to subdomain incidency (each processor)
   
   pcor->nbndom (nproc) - array containing numbers of boundary nodes on subdomains (only on master)
   noddom (maxnbn,nproc) - array containing node-subdomain correspondence (only on master)
   nodinc (maxgnn) - array containing number of node to subdomain incidencies (only on master)
   gncn (maxgnn,nodinc*pcor->totmaxndofn) - array containing local code numbers (only on master),
                                                       global code numbers (only on master)  
   pcor->lcngcn (nbn,inc*pcor->totmaxndofn) - array containing (each processor)
   pcor->nodnum (nbn) - array containing only numbers of boundary nodes (each processor)
   
   @param top - topology object pointer
   @param out - output stream
   
   JK, 15.1.2002
*/
void paral::globcnnum_feti (gtopology *top,FILE *out)
{
  long i,j,k,l,m,ii,jj,kk,ll,maxndofn,gnn,maxnbn,maxgnn,maxinc;
  long *buff,**noddom,*nodinc,**gncn;
  MPI_Status stat;
  
  //***************************************************************************
  //  determination of maximum number of nodes and maximum number of nodal DOFs
  //***************************************************************************
  //  maximum number of degrees of freedom on node
  maxndofn=0;
  //  number of boundary nodes
  nbn=0;
  //  maximum global number of node
  gnn=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      if (maxndofn<top->give_ndofn (i))  maxndofn = top->give_ndofn (i);
      nbn++;
      if (gnn<ltg[i])  gnn=ltg[i];
    }
  }
  
  
  buff = new long [3];
  buff[0]=maxndofn;
  buff[1]=nbn;
  buff[2]=gnn;

  if (myrank==0){
    //  array of numbers of boundary nodes
    nbndom = new long [nproc];
    nbndom[ndom]=nbn;

    //  maximum number of degrees of freedom on node
    totmaxndofn=maxndofn;
    //  maximum number of boundary nodes
    maxnbn=nbn;
    //  maximum global number of node
    maxgnn=gnn;
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      nbndom[stat.MPI_TAG]=buff[1];
      if (buff[0]>totmaxndofn)  totmaxndofn=buff[0];
      if (buff[1]>maxnbn)  maxnbn=buff[1];
      if (buff[2]>maxgnn)  maxgnn=buff[2];
      
    }
    maxgnn++;

    fprintf (out,"\n\n\n Global data information \n");
    fprintf (out,"\n maxnbn %ld",maxnbn);
    fprintf (out,"\n totmaxndofn %ld",totmaxndofn);
    fprintf (out,"\n maxgnn %ld",maxgnn);
    fprintf (out,"\n\n pocty velicin na podoblastech");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nbn %4ld",nbndom[i]);
    }
    
    buff[0]=totmaxndofn;
    buff[1]=maxnbn;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,3,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,3,MPI_LONG,0,ndom,MPI_COMM_WORLD);


    MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    totmaxndofn=buff[0];
    maxnbn=buff[1];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  //****************************************************
  //  node-subdomain correspondence assembling on master
  //****************************************************
  buff = new long [maxnbn];
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      buff[j]=ltg[i];  j++;
    }
  }
  
  if (myrank==0){
    //  array of node-subdomain correspondence
    noddom = new long* [maxgnn];
    for (i=0;i<maxgnn;i++){
      noddom[i] = new long [nproc];
      for (j=0;j<nproc;j++){
	noddom[i][j]=-1;
      }
    }
    
    k=ndom;
    for (j=0;j<nbndom[k];j++){
      noddom[buff[j]][k]=1;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=stat.MPI_TAG;
      for (j=0;j<nbndom[k];j++){
	noddom[buff[j]][k]=1;
      }
    }
    
    /*
    fprintf (Out,"\n\n\n\n kontrola pole noddom\n");
    for (i=0;i<maxgnn;i++){
      fprintf (Out,"\n");
      for (j=0;j<nproc;j++){
	fprintf (Out," %ld",noddom[i][j]);
      }
    }
    */
    
  }
  else{
    MPI_Send (buff,maxnbn,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  //*************************************
  //  evaluation of number of incidencies
  //*************************************
  if (myrank==0){
    //  array of numbers of incidences to subdomains of nodes
    nodinc = new long [maxgnn];
    //  maximum number of incidences
    maxinc=0;
    for (i=0;i<maxgnn;i++){
      k=0;
      for (j=0;j<nproc;j++){
	if (noddom[i][j]==1)  k++;
      }
      if (maxinc<k)  maxinc=k;
      nodinc[i]=k;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxinc,1,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
    
    /*
    fprintf (Out,"\n\n\n kontrola pole nodinc \n");
    for (i=0;i<maxgnn;i++){
      fprintf (Out," %ld",nodinc[i]);
    }
    */
    
  }
  else{
    MPI_Recv (&maxinc,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);


  //*************************************************
  //  code numbers gathering on master processor
  //  global code numbers generation
  //  computation of number of Lagrange multipliers
  //*************************************************
  
  //********************************************
  //  code numbers gathering on master processor
  //********************************************
  buff = new long [maxnbn*(totmaxndofn+1)];
  for (i=0;i<maxnbn*(totmaxndofn+1);i++){
    buff[i]=0;
  }
  
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      buff[k*(totmaxndofn+1)+totmaxndofn]=ltg[i];
      for (j=0;j<top->give_ndofn (i);j++){
	buff[k*(totmaxndofn+1)+j]=top->give_dof (i,j);
      }
      k++;
    }
  }

  
  if (myrank==0){
    gncn = new long* [maxgnn];
    for (i=0;i<maxgnn;i++){
      gncn[i] = new long [nodinc[i]*totmaxndofn];
      for (j=0;j<nodinc[i]*totmaxndofn;j++){
	gncn[i][j]=0;
      }
    }
    
    //  master contribution
    for (i=0;i<maxgnn;i++){
      nodinc[i]=0;
    }
    for (j=0;j<nbndom[ndom];j++){
      l=buff[j*(totmaxndofn+1)+totmaxndofn];
      m=nodinc[l]*totmaxndofn;
      for (k=0;k<totmaxndofn;k++){
	gncn[l][m]=buff[j*(totmaxndofn+1)+k];  m++;
      }
      nodinc[l]++;
    }


    // slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      for (j=0;j<nbndom[stat.MPI_TAG];j++){
	l=buff[j*(totmaxndofn+1)+totmaxndofn];
	m=nodinc[l]*totmaxndofn;
	for (k=0;k<totmaxndofn;k++){
	  gncn[l][m]=buff[j*(totmaxndofn+1)+k];  m++;
	}
	nodinc[l]++;
      }
    }
    
    /*
    fprintf (Out,"\n\n\n kontrola pole nodinc \n");
    for (i=0;i<maxgnn;i++){
      fprintf (Out," %ld",nodinc[i]);
    }
    fprintf (Out,"\n\n\n kontrola globalnich kodovych cisel  ngdof %ld \n",ngdof);
    for (i=0;i<maxgnn;i++){
      fprintf (Out,"\n");
      for (j=0;j<nodinc[i]*totmaxndofn;j++){
	fprintf (Out," %ld",gncn[i][j]);
      }
    }
    */
    
    
    //**************************************************
    //  computation of number of Lagrange multipliers
    //  global code numbers generation
    //**************************************************
    ngdof=1;
    for (i=0;i<maxgnn;i++){
      for (j=totmaxndofn;j<nodinc[i]*totmaxndofn;j++){
	if (gncn[i][j]>0){
	  gncn[i][j]=ngdof;  ngdof++;
	}
      }
    }
    ngdof--;

    
    /*
    fprintf (Out,"\n\n\n kontrola globalnich kodovych cisel  ngdof %ld \n",ngdof);
    for (i=0;i<maxgnn;i++){
      fprintf (Out,"\n");
      for (j=totmaxndofn;j<nodinc[i]*totmaxndofn;j++){
	fprintf (Out," %ld",gncn[i][j]);
      }
    }
    */
    

  }
  else{
    MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  
  //*******************************
  //  gobal code numbers scattering
  //*******************************
  
  buff = new long [maxnbn*(maxinc*totmaxndofn+1)];

  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      buff[j*(maxinc*totmaxndofn+1)]=ltg[i];  j++;
    }
  }

  //  array containing global node numbers of boundary nodes
  nodnum = new long [nbn];
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      nodnum[j]=i;
      j++;
    }
  }
  

  if (myrank==0){
    
    ii=ndom;
    
    
    for (j=0;j<nbndom[ii];j++){
      jj=buff[j*(maxinc*totmaxndofn+1)];
      for (k=0;k<nproc;k++){
	if (noddom[jj][k]==1)  break;
      }
      if (k==ii){
	kk=0;
	for (k=totmaxndofn;k<nodinc[jj]*totmaxndofn;k++){
	  buff[j*(maxinc*totmaxndofn+1)+kk]=0-gncn[jj][k];  kk++;
	}
      }
      else{
	kk=0;
	for (k=0;k<ii;k++){
	  if (noddom[jj][k]==1)  kk++;
	}
	ll=0;
	for (k=kk*totmaxndofn;k<(kk+1)*totmaxndofn;k++){
	  buff[j*(maxinc*totmaxndofn+1)+ll]=gncn[jj][k];  ll++;
	}
      }
      buff[j*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn]=nodinc[jj]-1;
    }
    
    
    //  array containing numbers of incidencies of nodes to subdomains
    inc = new long [nbn];
    //  array containing global code numbers in FETI method
    lcngcn = new long* [nbn];
    for (i=0;i<nbn;i++){
      inc[i]=buff[i*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn];
      lcngcn[i] = new long [inc[i]*totmaxndofn];
      for (j=0;j<inc[i]*totmaxndofn;j++){
	lcngcn[i][j]=buff[i*(maxinc*totmaxndofn+1)+j];
      }
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      ii=stat.MPI_TAG;
      

      for (j=0;j<nbndom[ii];j++){
	jj=buff[j*(maxinc*totmaxndofn+1)];
	for (k=0;k<nproc;k++){
	  if (noddom[jj][k]==1)  break;
	}
	if (k==ii){
	  kk=0;
	  for (k=totmaxndofn;k<nodinc[jj]*totmaxndofn;k++){
	    buff[j*(maxinc*totmaxndofn+1)+kk]=0-gncn[jj][k];  kk++;
	  }
	  buff[j*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn]=nodinc[jj]-1;
	}
	else{
	  kk=0;
	  for (k=0;k<ii;k++){
	    if (noddom[jj][k]==1)  kk++;
	  }
	  ll=0;
	  for (k=kk*totmaxndofn;k<(kk+1)*totmaxndofn;k++){
	    buff[j*(maxinc*totmaxndofn+1)+ll]=gncn[jj][k];  ll++;
	  }
	  buff[j*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn]=1;
	}
      }


      MPI_Send (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,domproc[ii],ndom,MPI_COMM_WORLD);
      
      
    }
  }
  else{
    MPI_Send (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    
    //  array containing numbers of incidencies of nodes to subdomains
    inc = new long [nbn];
    //  array containing global code numbers in FETI method
    lcngcn = new long* [nbn];
    for (i=0;i<nbn;i++){
      inc[i]=buff[i*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn];
      lcngcn[i] = new long [inc[i]*totmaxndofn];
      for (j=0;j<inc[i]*totmaxndofn;j++){
	lcngcn[i][j]=buff[i*(maxinc*totmaxndofn+1)+j];
      }
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  delete [] buff;
  
  if (myrank==0){
    delete [] nodinc;
    //*********************
    for (i=0;i<maxgnn;i++){
      delete [] gncn[i];
    }
    delete [] gncn;
    //*********************
    for (i=0;i<maxgnn;i++){
      delete [] noddom[i];
    }
    delete [] noddom;
  }

  /*
  fprintf (Out,"\n\n\n kontrola");
  for (i=0;i<nbn;i++){
    fprintf (Out,"\n uzel %ld  pocet podoblasti %ld\n",nodnum[i],inc[i]);
    fprintf (Out,"lcngcn");
    for (j=0;j<inc[i]*totmaxndofn;j++){
      fprintf (Out," %ld",lcngcn[i][j]);
    }
  }
  */
  

}





/**
   function localizes local vector components to global vector
   localization process is generated by FETI ordering
   
   @param gv - array containing global vector
   @param lv - array containing local vector

   4.12.2001
*/
void paral::locglobfeti (gtopology *top,double *gv,double *lv)
{
  long i,j,k,l,ln,lcn,llgcn;
  
  
  //fprintf (Out,"\n\n number of boundary nodes %ld",nbn);
  //fprintf (Out,"\n\n nove volani funkce locglobfeti  nbn %ld   pcor->totmaxndofn %ld\n",nbn,pcor->totmaxndofn);
  

  for (i=0;i<nbn;i++){
    
    l=0;
    ln=nodnum[i];
    
    //fprintf (Out,"\n pcor->nodnum %ld  inc %ld",nodnum[i],inc[i]);
    
    
    for (j=0;j<inc[i];j++){
      for (k=0;k<totmaxndofn;k++){
	lcn=top->give_dof (ln,k)-1;
	llgcn=lcngcn[i][l];
	
	//fprintf (Out,"    lcn %ld  gcn %ld",lcn,gcn);
	
	if (llgcn<0) gv[0-llgcn-1]-=lv[lcn];
	if (llgcn>0) gv[llgcn-1]+=lv[lcn];
	l++;
      }
    }
    
    
  
  }

 
  
}


/**
   function creates local vector from global vector
   localization process is generated by FETI ordering
   
   @param gv - array containing global vector
   @param lv - array containing local vector

   4.12.2001
*/
void paral::globlocfeti (gtopology *top,double *gv,double *lv)
{
  long i,j,k,l,ln,lcn,llgcn;
  

  for (i=0;i<nbn;i++){
    l=0;  ln=nodnum[i];
    for (j=0;j<inc[i];j++){
      for (k=0;k<totmaxndofn;k++){
	lcn=top->give_dof (ln,k)-1;
	llgcn=lcngcn[i][l];
	if (llgcn<0) lv[lcn]-=gv[0-llgcn-1];
	if (llgcn>0) lv[lcn]+=gv[llgcn-1];
	l++;
      }
    }
  }

}



//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************
//************************************************************************

/**
   function orders unknowns on subdomain with respect of DP-FETI method
   unknowns are decomposed into three parts
   first part contains inner nodes/unknowns
   second part contains boundary nodes/unknowns except corner nodes/unknowns
   third part contains corner nodes/unknowns
   function orders unknowns from two first parts at the beginning
   unknowns from the third part are ordered at the end
   nodes are marked with help of ltg array

   in the input file
   ltg[first part]=0
   ltg[second part]>0
   ltg[third part]<0
   before use
   ltg[first part]=-1
   ltg[second part]>-1
   ltg[third part]<-1

   @param top - pointer to topology
   
   JK, 26.11.2002
*/
void paral::dpfetiordering (gtopology *top)
{
  long i,j,k,ndofn;
  
  //  contributions from nodes
  ndof=1;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>=-1){
      //  inner and boundary nodes
      ndofn=top->give_ndofn (i);
      for (j=0;j<ndofn;j++){
	k=top->give_dof(i,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k>0){
	  top->save_dof (i,j,ndof);  ndof++;
	}
      }
    }
  }
  indof=ndof-1;
  
  for (i=0;i<top->nn;i++){
    if (ltg[i]<-1){
      //  corner nodes
      ndofn=top->give_ndofn (i);
      for (j=0;j<ndofn;j++){
	k=top->give_dof (i,j);
	if (k<0)  continue;
	if (k==0)  continue;
	if (k>0){
	  top->save_dof (i,j,ndof);  ndof++;
	}
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
  
}

/**
   function creates code numbers for DP-FETI method
   
   @param top - topology object pointer
   @param out - output stream
   
   JK, 15.1.2002
*/
void paral::globcnnum_dpfeti (gtopology *top,FILE *out)
{
  long i,j;
  long *buff,**noddom,*nodinc,**gncn;
  MPI_Status stat;
  
  //***************************************************************************
  //  determination of maximum number of nodes and maximum number of nodal DOFs
  //***************************************************************************
  //  maximum number of degrees of freedom on node
  maxndofn=0;
  //  number of boundary nodes
  nbn=0;
  //  maximum global number of boundary node
  gnbn=0;
  //  number of corner nodes
  ncn=0;
  //  maximum global number of corner node
  gncn=0;
  for (i=0;i<top->nn;i++){
    if (maxndofn<top->give_ndofn (i))  maxndofn = top->give_ndofn (i);
    
    if (ltg[i]>-1){
      nbn++;
      if (gnbn<ltg[i])  gnbn=ltg[i];
    }
    
    if (ltg[i]<-1){
      ncn++;
      if (gncn<0-ltg[i]-1)  gncn=0-ltg[i]-1;
    }
  }
  
  
  buff = new long [3];
  buff[0]=maxndofn;
  buff[1]=nbn;
  buff[2]=gnbn;
  buff[3]=ncn;
  buff[4]=gncn;
  buff[5]=ndof-indof;

  if (myrank==0){
    //  array of numbers of boundary nodes
    nbndom = new long [nproc];
    nbndom[ndom]=nbn;
    
    //  array of numbers of corner nodes
    ncndom = new long [nproc];
    ncndom[ndom]=ncn;
    
    //  array of numbers of corner unknowns
    ncdofdom = new long [nproc];
    ncdofdom[ndom]=ndof-indof;
    
    //  maximum number of degrees of freedom on node
    totmaxndofn=maxndofn;
    //  maximum number of boundary nodes
    maxnbn=nbn;
    //  maximum global number of boundary node
    maxgnbn=gnbn;
    //  maximum number of corner nodes
    maxncn=ncn;
    //  maximum global number of corner node
    maxgncn=gncn;
    //  maximum number of corner unknowns
    maxncdof=ndof-indof;
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,5,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      nbndom[stat.MPI_TAG]=buff[1];
      ncndom[stat.MPI_TAG]=buff[3];
      ncdofdom[stat.MPI_TAG]=buff[5];
      
      if (buff[0]>totmaxndofn)  totmaxndofn=buff[0];
      if (buff[1]>maxnbn)       maxnbn=buff[1];
      if (buff[2]>maxgnbn)      maxgnbn=buff[2];
      if (buff[3]>maxncn)       maxncn=buff[3];
      if (buff[4]>maxgncn)      maxgncn=buff[4];
      if (buff[5]>maxncdof)     maxncdof=buff[5];
    }
    maxgnbn++;  maxgncn++;
    
    //  control print
    fprintf (out,"\n\n\n Global data information \n");
    fprintf (out,"\n totmaxndofn %ld",totmaxndofn);
    fprintf (out,"\n maxnbn %ld",maxnbn);
    fprintf (out,"\n maxgnbn %ld",maxgnbn);
    fprintf (out,"\n maxncn %ld",maxncn);
    fprintf (out,"\n maxgncn %ld",maxgncn);
    fprintf (out,"\n maxncdof %ld",maxcndof);
    fprintf (out,"\n\n pocty velicin na podoblastech");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nbn %4ld  ncn %4ld  ncdof",nbndom[i],ncndom[i],ncdofdom[i]);
    }
    
    
    buff[0]=totmaxndofn;
    buff[1]=maxncn;
    buff[2]=maxncdof;
    buff[3]=maxnbn;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,5,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,5,MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,5,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    totmaxndofn=buff[0];
    maxncn=buff[1];
    maxncdof=buff[2];
    maxnbn=buff[3];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  //***************************
  //  corner nodes and unknowns
  //***************************


  //*******************************************
  //  gathering of indicators of DOFs to master 
  //*******************************************
  //  allocation of array buff
  buff = new long [maxncn*(totmaxndofn+1)];
  for (i=0;i<maxncn*(totmaxndofn+1);i++){
    buff[i]=0;
  }
  
  //  indicators are placed to array buff
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]<-1){
      buff[k*(totmaxndofn+1)+totmaxndofn]=0-ltg[i]-1;
      for (j=0;j<top->give_ndofn (i);j++){
	buff[k*(totmaxndofn+1)+j]=top->give_dof (i,j);
      }
      k++;
    }
  }
  
  
  if (myrank==0){
    
    //  allocation of array gncn
    gncn = new long* [maxgncn];
    for (i=0;i<maxgncn;i++){
      gncn[i] = new long [totmaxndofn];
      for (j=0;j<totmaxndofn;j++){
	gncn[i][j]=0;
      }
    }
    
    //  master contribution
    for (j=0;j<ncndom[ndom];j++){
      l=buff[j*(totmaxndofn+1)+totmaxndofn];
      for (k=0;k<totmaxndofn;k++){
	gncn[l][k]=buff[j*(totmaxndofn+1)+k];
      }
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxncn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      for (j=0;j<ncndom[stat.MPI_TAG];j++){
	l=buff[j*(totmaxndofn+1)+totmaxndofn];
	for (k=0;k<totmaxndofn;k++){
	  gncn[l][k]=buff[j*(totmaxndofn+1)+k];
	}
      }
    }
    
    //  computation of total number of global unknowns
    ncdof=1;
    for (i=0;i<maxgncn;i++){
      for (j=0;j<totmaxndofn;j++){
	if (gncn[i][j]>0){
	  gncn[i][j]=ncdof;  ncdof++;
	}
      }
    }
    ncdof--;
    
  }
  else{
    MPI_Send (buff,maxncn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  //********************************
  //  global code numbers scattering
  //********************************
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxncn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      for (j=0;j<ncndom[stat.MPI_TAG];j++){
	l=buff[j*(totmaxndofn+1)+totmaxndofn];
	for (k=0;k<totmaxndofn;k++){
	  buff[j*(totmaxndofn+1)+k]=gncn[l][k];
	}
      }
      l=stat.MPI_TAG;
      MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,l,ndom,MPI_COMM_WORLD);
    }

    k=0;
    for (i=0;i<top->nn;i++){
      if (ltg[i]<-1){
	buff[k*(totmaxndofn+1)+totmaxndofn]=0-ltg[i]-1;
	for (j=0;j<top->give_ndofn (i);j++){
	  buff[k*(totmaxndofn+1)+j]=top->give_dof(i,j);
	}
	k++;
      }
    }
    
    for (j=0;j<ncndom[ndom];j++){
      l=buff[j*(totmaxndofn+1)+totmaxndofn];
      for (k=0;k<totmaxndofn;k++){
	buff[j*(totmaxndofn+1)+k]=gncn[l][k];
      }
    }
    
    //  deallocation of array gncn
    for (i=0;i<maxgncn;i++){
      delete [] gncn[i];
    }
    delete [] gncn;
    
  }
  else{
    MPI_Send (buff,maxncn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,maxncn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  


  //*********************************************
  //  local to global code numbers correspondence
  //*********************************************
  //  allocation of array of local to global code numbers correspondence
  gcn = new long [ndof];
  for (i=0;i<ndof;i++){
    gcn[i]=0;
  }
  
  //  filling of local to global code numbers correspondence
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]<-1){
      for (j=0;j<top->give_ndofn (i);j++){
	l=top->give_dof (i,j)-1;
	if (l>-1)  gcn[l]=buff[k*(totmaxndofn+1)+j];
      }
      k++;
    }
  }
  
  delete [] buff;
  
  
  //********************************************************************
  //  filling of array of local to global code numbers of all subdomains
  //********************************************************************
  buff = new long [maxncdof];
  j=indof;
  for (i=0;i<ndof-indof;i++){
    buff[i]=gcn[j];  j++;
  }
  
  if (myrank==0){
    //  allocation of array masgcn
    masgcn = new long* [nproc];
    for (i=0;i<nproc;i++){
      masgcn[i] = new long [ncdofdom[i]];
      for (j=0;j<ncdofdom[i];j++){
	masgcn[i][j]=0;
      }
    }
    
    //  master contribution
    for (j=0;j<ncdofdom[ndom];j++){
      masgcn[ndom][j]=buff[j];
    }
    
    //  slaves contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxncdof,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      for (j=0;j<ncdofdom[stat.MPI_TAG];j++){
	masgcn[stat.MPI_TAG][j]=buff[j];
      }
      
    }
  }
  else{
    MPI_Send (buff,maxncdof,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  
  //  pole gcn neni jiz potreba









  
  //*******************************
  //  boundary nodes and unknowns
  //*******************************


  //****************************************************
  //  node-subdomain correspondence assembling on master
  //****************************************************
  buff = new long [maxnbn];
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1){
      buff[j]=ltg[i];  j++;
    }
  }
  
  if (myrank==0){
    //  array of node-subdomain correspondence
    noddom = new long* [maxgnbn];
    for (i=0;i<maxgnbn;i++){
      noddom[i] = new long [nproc];
      for (j=0;j<nproc;j++){
	noddom[i][j]=-1;
      }
    }
    
    k=ndom;
    for (j=0;j<nbndom[k];j++){
      noddom[buff[j]][k]=1;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=stat.MPI_TAG;
      for (j=0;j<nbndom[k];j++){
	noddom[buff[j]][k]=1;
      }
    }
    
    /*
    fprintf (Out,"\n\n\n\n kontrola pole noddom\n");
    for (i=0;i<maxgnn;i++){
      fprintf (Out,"\n");
      for (j=0;j<nproc;j++){
	fprintf (Out," %ld",noddom[i][j]);
      }
    }
    */
    
  }
  else{
    MPI_Send (buff,maxnbn,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  //*************************************
  //  evaluation of number of incidencies
  //*************************************
  if (myrank==0){
    //  array of numbers of incidences to subdomains of nodes
    nodinc = new long [maxgnbn];
    //  maximum number of incidences
    maxinc=0;
    for (i=0;i<maxgnbn;i++){
      k=0;
      for (j=0;j<nproc;j++){
	if (noddom[i][j]==1)  k++;
      }
      if (maxinc<k)  maxinc=k;
      nodinc[i]=k;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxinc,1,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
    
    /*
    fprintf (Out,"\n\n\n kontrola pole nodinc \n");
    for (i=0;i<maxgnn;i++){
      fprintf (Out," %ld",nodinc[i]);
    }
    */
    
  }
  else{
    MPI_Recv (&maxinc,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);


  //*************************************************
  //  code numbers gathering on master processor
  //  global code numbers generation
  //  computation of number of Lagrange multipliers
  //*************************************************
  
  //********************************************
  //  code numbers gathering on master processor
  //********************************************
  buff = new long [maxnbn*(totmaxndofn+1)];
  for (i=0;i<maxnbn*(totmaxndofn+1);i++){
    buff[i]=0;
  }
  
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>=-1){
      buff[k*(totmaxndofn+1)+totmaxndofn]=ltg[i];
      for (j=0;j<top->give_ndofn (i);j++){
	buff[k*(totmaxndofn+1)+j]=top->give_dof (i,j);
      }
      k++;
    }
  }

  
  if (myrank==0){
    gncn = new long* [maxgnbn];
    for (i=0;i<maxgnbn;i++){
      gncn[i] = new long [nodinc[i]*totmaxndofn];
      for (j=0;j<nodinc[i]*totmaxndofn;j++){
	gncn[i][j]=0;
      }
    }
    
    //  master contribution
    for (i=0;i<maxgnbn;i++){
      nodinc[i]=0;
    }
    for (j=0;j<nbndom[ndom];j++){
      l=buff[j*(totmaxndofn+1)+totmaxndofn];
      m=nodinc[l]*totmaxndofn;
      for (k=0;k<totmaxndofn;k++){
	gncn[l][m]=buff[j*(totmaxndofn+1)+k];  m++;
      }
      nodinc[l]++;
    }


    // slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn*(totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      for (j=0;j<nbndom[stat.MPI_TAG];j++){
	l=buff[j*(totmaxndofn+1)+totmaxndofn];
	m=nodinc[l]*totmaxndofn;
	for (k=0;k<totmaxndofn;k++){
	  gncn[l][m]=buff[j*(totmaxndofn+1)+k];  m++;
	}
	nodinc[l]++;
      }
    }
    
    /*
    fprintf (Out,"\n\n\n kontrola pole nodinc \n");
    for (i=0;i<maxgnn;i++){
      fprintf (Out," %ld",nodinc[i]);
    }
    fprintf (Out,"\n\n\n kontrola globalnich kodovych cisel  ngdof %ld \n",ngdof);
    for (i=0;i<maxgnn;i++){
      fprintf (Out,"\n");
      for (j=0;j<nodinc[i]*totmaxndofn;j++){
	fprintf (Out," %ld",gncn[i][j]);
      }
    }
    */
    
    
    //**************************************************
    //  computation of number of Lagrange multipliers
    //  global code numbers generation
    //**************************************************
    nmdof=1;
    for (i=0;i<maxgnbn;i++){
      for (j=totmaxndofn;j<nodinc[i]*totmaxndofn;j++){
	if (gncn[i][j]>0){
	  gncn[i][j]=nmdof;  nmdof++;
	}
      }
    }
    nmdof--;

    
    /*
    fprintf (Out,"\n\n\n kontrola globalnich kodovych cisel  ngdof %ld \n",ngdof);
    for (i=0;i<maxgnn;i++){
      fprintf (Out,"\n");
      for (j=totmaxndofn;j<nodinc[i]*totmaxndofn;j++){
	fprintf (Out," %ld",gncn[i][j]);
      }
    }
    */
    

  }
  else{
    MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  
  //*******************************
  //  gobal code numbers scattering
  //*******************************
  
  buff = new long [maxnbn*(maxinc*totmaxndofn+1)];

  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1){
      buff[j*(maxinc*totmaxndofn+1)]=ltg[i];  j++;
    }
  }

  //  array containing global node numbers of boundary nodes
  nodnum = new long [nbn];
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1){
      nodnum[j]=i;
      j++;
    }
  }
  

  if (myrank==0){
    
    ii=ndom;
    
    
    for (j=0;j<nbndom[ii];j++){
      jj=buff[j*(maxinc*totmaxndofn+1)];
      for (k=0;k<nproc;k++){
	if (noddom[jj][k]==1)  break;
      }
      if (k==ii){
	kk=0;
	for (k=totmaxndofn;k<nodinc[jj]*totmaxndofn;k++){
	  buff[j*(maxinc*totmaxndofn+1)+kk]=0-gncn[jj][k];  kk++;
	}
      }
      else{
	kk=0;
	for (k=0;k<ii;k++){
	  if (noddom[jj][k]==1)  kk++;
	}
	ll=0;
	for (k=kk*totmaxndofn;k<(kk+1)*totmaxndofn;k++){
	  buff[j*(maxinc*totmaxndofn+1)+ll]=gncn[jj][k];  ll++;
	}
      }
      buff[j*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn]=nodinc[jj]-1;
    }
    
    
    //  array containing numbers of incidencies of nodes to subdomains
    inc = new long [nbn];
    //  array containing global code numbers in FETI method
    lcngcn = new long* [nbn];
    for (i=0;i<nbn;i++){
      inc[i]=buff[i*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn];
      lcngcn[i] = new long [inc[i]*totmaxndofn];
      for (j=0;j<inc[i]*totmaxndofn;j++){
	lcngcn[i][j]=buff[i*(maxinc*totmaxndofn+1)+j];
      }
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      ii=stat.MPI_TAG;
      

      for (j=0;j<nbndom[ii];j++){
	jj=buff[j*(maxinc*totmaxndofn+1)];
	for (k=0;k<nproc;k++){
	  if (noddom[jj][k]==1)  break;
	}
	if (k==ii){
	  kk=0;
	  for (k=totmaxndofn;k<nodinc[jj]*totmaxndofn;k++){
	    buff[j*(maxinc*totmaxndofn+1)+kk]=0-gncn[jj][k];  kk++;
	  }
	  buff[j*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn]=nodinc[jj]-1;
	}
	else{
	  kk=0;
	  for (k=0;k<ii;k++){
	    if (noddom[jj][k]==1)  kk++;
	  }
	  ll=0;
	  for (k=kk*totmaxndofn;k<(kk+1)*totmaxndofn;k++){
	    buff[j*(maxinc*totmaxndofn+1)+ll]=gncn[jj][k];  ll++;
	  }
	  buff[j*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn]=1;
	}
      }


      MPI_Send (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,domproc[ii],ndom,MPI_COMM_WORLD);
      
      
    }
  }
  else{
    MPI_Send (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    
    //  array containing numbers of incidencies of nodes to subdomains
    inc = new long [nbn];
    //  array containing global code numbers in FETI method
    lcngcn = new long* [nbn];
    for (i=0;i<nbn;i++){
      inc[i]=buff[i*(maxinc*totmaxndofn+1)+maxinc*totmaxndofn];
      lcngcn[i] = new long [inc[i]*totmaxndofn];
      for (j=0;j<inc[i]*totmaxndofn;j++){
	lcngcn[i][j]=buff[i*(maxinc*totmaxndofn+1)+j];
      }
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  delete [] buff;
  
  if (myrank==0){
    delete [] nodinc;
    //*********************
    for (i=0;i<maxgnn;i++){
      delete [] gncn[i];
    }
    delete [] gncn;
    //*********************
    for (i=0;i<maxgnn;i++){
      delete [] noddom[i];
    }
    delete [] noddom;
  }

  /*
  fprintf (Out,"\n\n\n kontrola");
  for (i=0;i<nbn;i++){
    fprintf (Out,"\n uzel %ld  pocet podoblasti %ld\n",nodnum[i],inc[i]);
    fprintf (Out,"lcngcn");
    for (j=0;j<inc[i]*totmaxndofn;j++){
      fprintf (Out," %ld",lcngcn[i][j]);
    }
  }
  */
  

}


