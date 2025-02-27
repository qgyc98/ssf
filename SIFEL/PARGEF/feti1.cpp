#include <mpi.h>
#include "feti1.h"
#include <string.h>
#include <math.h>
#include <time.h>

feti1::feti1(int np,int mr,long nd)
{
  nproc=np;  myrank=mr;  ndom=nd;
  
  ndof=0;  maxlggl=0;  nclggl=0;  maxnrbm=0;
  enrbm=0;  lithr=0.0;  nrbm=0;
  nicg=0;  anicg=0;  errcg=0.0;  aerrcg=0.0;
  zero=1.0e-15;
  
  ndofcp=0;  hsize=0;
  
  lcn=NULL;

  lggl=NULL;
  glcor=NULL;

  nrbmdom=NULL;
  rbmadr=NULL;
  
  // ***********
  //  sobota 30.7.2005
  // ************
  edofs = NULL;
  
  ncnd = NULL;

  // ***********
  //  nedele 31.7.2005
  // ************
  ndofncn = NULL;
  cncn = NULL;

  // diagonal scaling matrix
  wscalmat = NULL;
  lwscalmat = NULL;
  
  smsky = NULL;
  
  fetiprecond = -1;
}

feti1::~feti1()
{
  long i;
  
  delete [] lcn;
  
  
  for (i=0;i<nproc;i++){
    delete [] glcor[i];
  }
  delete [] glcor;

  delete [] lggl;
  delete [] nrbmdom;
  delete [] rbmadr;
  
  // ***********
  //  sobota 30.7.2005
  // ***********
  delete [] edofs;
  
  if (myrank==0){
    delete [] ncnd;
  }
  // ***********
  //  nedele 31.7.
  // ****************
  
}

/**
   JK, 7.8.2007
*/
void feti1::initiate (selnodes *selnodfeti,FILE *out)
{
  long i;
  //  maximum number of components in local-global arrays
  maxlggl = selnodfeti->maxndof;
  //  number of components in local-global arrays
  nclggl = selnodfeti->ndof;
  
  //  maximum number of DOFs contributing to the coarse problem
  maxncdofd = selnodfeti->maxndof;

  //  number of contributions (unknowns) from particular subdomains to the coarse problem
  ncdof = selnodfeti->ndof;

  //  array containing local code numbers of boundary unknonws
  lcn = selnodfeti->ldof;

  //  array containing code numbers contributing to the coarse problem
  edofs = selnodfeti->ldof;
  
  fprintf (out,"\n\n\n osel \n");
  for (long i=0;i<ncdof;i++){
    fprintf (out,"\n edofs %ld",edofs[i]);
  }
  fprintf (out,"\n\n\n osel \n");
  
  
  if (myrank==0){
    //  number of DOFs (unknowns) in coarse problem = total number of boundary DOFs
    ndofcp =selnodfeti->tndof;
    
    ncdofd = selnodfeti->ndofdom;
    //lggl = selnodfeti->ndofdom;
    
    //glcor = selnodfeti->cndom;
    ccn = selnodfeti->cndom;
    
    wscalmat = new double [ndofcp]; 
    for (i = 0; i < ndofcp; i++){
      wscalmat[i] = double(selnodfeti->dofmultip[i]);
    }
    
  }
  
}

/**
   function generates coarse degrees of freedom

   @param ptop - pointer to parallel topology
   @param out - pointer to output file
   
   JK, 31.7.2005
*/
void feti1::coarse_dofs (partop *ptop,FILE *out)
{
  long i,j,k,m,n;
  long tnbn;
  long *bnmultip,**nbdofnd,**llnbn,**ldn,***lbcn;
  
  if (myrank==0){
    
    //  total number of boundary nodes = number of coarse nodes
    tnbn = ptop->tnbn;
    //  array containing multiplicity of boundary nodes
    bnmultip = ptop->bnmultip;
    //  array containing numbers of DOFs on boundary nodes
    nbdofnd = ptop->nbdofnd;
    //  array containing local numbers of boundary nodes belonging to coarse nodes
    llnbn = ptop->llnbn;
    //  array containing numbers of subdomain to which coarse nodes belong to
    ldn = ptop->ldn;
    //  array containing local boundary code numbers
    lbcn = ptop->lbcn;

    
    //  number of DOFs at coarse nodes
    //  ndofncn[i][j]=k - the j-th node shared by the i-th coarse node contains k DOFs
    ndofncn = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      ndofncn[i] = new long [bnmultip[i]];
    }
    
    for (i=0;i<tnbn;i++){
      for (j=0;j<bnmultip[i];j++){
	k=llnbn[i][j];
	m=ldn[i][j];
	ndofncn[i][j]=nbdofnd[m][k];
      }
    }


    //  array of code numbers at coarse nodes
    //  cncn[i][j][k]=m - the k-th DOF at the j-th node shared by the i-th coarse node has number m
    cncn = new long** [tnbn];
    for (i=0;i<tnbn;i++){
      cncn[i] = new long* [bnmultip[i]];
      for (j=0;j<bnmultip[i];j++){
	cncn[i][j] = new long [ndofncn[i][j]];
      }
    }
    
    //  assignment of constraint indicators
    for (i=0;i<tnbn;i++){
      for (j=0;j<bnmultip[i];j++){
	m=llnbn[i][j];
	n=ldn[i][j];
	for (k=0;k<ndofncn[i][j];k++){
	  cncn[i][j][k]=lbcn[n][m][k];
	}
      }
    }

    //  code numbers generation
    ndofcp=1;
    for (i=0;i<tnbn;i++){
      for (j=0;j<bnmultip[i]-1;j++){
	for (k=0;k<ndofncn[i][j];k++){
	  if (cncn[i][j][k]>0){
	    cncn[i][j][k]=ndofcp;
	    ndofcp++;
	  }
	}
      }
    }
    ndofcp--;
    
    fprintf (out,"\n\n number of unknowns (DOFs) in the coarse problem is  %ld",ndofcp);
    fprintf (stdout,"\n\n number of unknowns (DOFs) in the coarse problem is  %ld",ndofcp);
  }
}

/**
   function determines number of nodes and unknowns (DOFs) contributing to the coarse problem
   
   function determines

   maxncnd - maximum number of nodes contributing to the coarse problem
   maxncdofd - maximum number of DOFs contributing to the coarse problem
   ncn - number of nodes contributing to the coarse problem

   ncnd - array of numbers of nodes contributing to coarse problem
   ncdofd - array of numbers of unknowns (DOFs) contributing to the coarse problem
   
   @param top - pointer to generalized topology of one subdomain
   @param ptop - pointer to parallel topology
   @param domproc - array containing correspondence between processors and subdomains
   @param out - pointer to output file

   JK, 31.7.2005
*/
void feti1::number_contributing_nodes_dofs (gtopology */*top*/,partop *ptop,long *domproc,FILE *out)
{
  long i,j,k,m,ii,jj,buffsize;
  long *buff;
  long *lnbn,*bnmultip,*nbnd,**ldn,**cnbn;
  MPI_Status stat;
  
  buffsize=3;
  buff = new long [buffsize];
  
  //  list of local numbers of boundary nodes
  lnbn = ptop->lnbn;

  if (myrank==0){
    //  array containing multiplicity of boundary nodes
    bnmultip = ptop->bnmultip;
    //  array containing list of numbers of boundary nodes on subdomains
    nbnd = ptop->nbnd;
    //  array containing numbers of subdomain to which coarse nodes belong to
    ldn = ptop->ldn;
    //  array containing coarse node numbers of boundary nodes
    cnbn = ptop->cnbn;

    //  array of numbers of nodes contributing to the coarse problem
    ncnd = new long [nproc];
    //  array of numbers of DOFs contributing to the coarse problem
    ncdofd = new long [nproc];
    
    //  maximum number of nodes contributing from subdomains
    maxncnd = 0;
    //  maximum number of DOFs contributing from subdomains
    maxncdofd = 0;
    
    
    for (i=0;i<nproc;i++){
      ncnd[i]=0;
      ncdofd[i]=0;
      for (j=0;j<nbnd[i];j++){
	m=cnbn[i][j];
	for (k=0;k<bnmultip[m];k++){
	  if (ldn[m][k]==i){
	    if (k==0){
	      for (ii=0;ii<bnmultip[m]-1;ii++){
		ncnd[i]++;
		for (jj=0;jj<ndofncn[m][ii];jj++){
		  if (cncn[m][ii][jj]>0)
		    ncdofd[i]++; 
		}
	      }
	    }
	    else{
	      ncnd[i]++;
	      for (jj=0;jj<ndofncn[m][k-1];jj++){
		if (cncn[m][k-1][jj]>0)  ncdofd[i]++;
	      }
	      
	    }
	    break;
	  }
	}
      }
      
      if (maxncnd<ncnd[i])
	maxncnd=ncnd[i];
      if (maxncdofd<ncdofd[i])
	maxncdofd=ncdofd[i];
    }
    
    for (i=1;i<nproc;i++){
      j=domproc[i];
      buff[0]=maxncnd;
      buff[1]=ncnd[j];
      buff[2]=maxncdofd;
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    j=domproc[0];
    buff[0]=maxncnd;
    buff[1]=ncnd[j];
    buff[2]=maxncdofd;
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  maxncnd=buff[0];
  ncn=buff[1];
  maxncdofd=buff[2];

  delete [] buff;
  
  fprintf (out,"\n\n KONTROLA 1.  maxncnd je %ld   ncn je %ld      maxncdofd je %ld",maxncnd,ncn,maxncdofd);
  
  if (myrank==0){
    fprintf (out,"\n\n kontrola poctu prispivajicich uzlu    ncnd\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"  %ld",ncnd[i]);
    }
    fprintf (out,"\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"  %ld",ncdofd[i]);
    }
    fprintf (out,"\n");
  }
}

/**
   function assembles lists of contributing nodes and unknowns (DOFs) to the coarse problem

   @param top - pointer to generalized topology of one subdomain
   @param ptop - pointer to parallel topology
   @param domproc - array containing correspondence between processors and subdomains
   @param out - pointer to output file
   
   
   JK, 31.7.2005
*/
void feti1::contributing_nodes_dofs (gtopology *top,partop *ptop,long *domproc,FILE *out)
{
  long i,j,k,l,m,ii,jj,ndofn,buffsize;
  long tnbn;
  long *buff;
  long *lnbn,*bnmultip,**ldn,*nbnd,**llnbn,**cnbn;
  MPI_Status stat;

  //  list of local numbers of boundary nodes
  lnbn = ptop->lnbn;

  if (myrank==0){
    //  total number of boundary nodes = number of coarse nodes
    tnbn = ptop->tnbn;
    //  array containing multiplicity of boundary nodes
    bnmultip = ptop->bnmultip;
    //  array containing numbers of subdomain to which coarse nodes belong to
    ldn = ptop->ldn;
    //  array containing list of numbers of boundary nodes on subdomains
    nbnd = ptop->nbnd;
    //  array containing local numbers of boundary nodes belonging to coarse nodes
    llnbn = ptop->llnbn;
    //  array containing coarse node numbers of boundary nodes
    cnbn = ptop->cnbn;
  }

  //  node numbers contributing to the coarse problem

  buffsize = maxncnd;
  buff = new long [buffsize];
  
  if (myrank==0){
    //  array of local numbers of nodes contributing to the coarse problem
    for (i=1;i<nproc;i++){
      
      for (j=0;j<buffsize;j++){
	buff[j]=0;
      }
      
      m=0;
      for (j=0;j<tnbn;j++){
	for (k=0;k<bnmultip[j];k++){
	  if (ldn[j][k]==domproc[i]){
	    if (k==0){
	      for (l=0;l<bnmultip[j]-1;l++){
		buff[m]=llnbn[j][k];
		buff[m]=0-buff[m]-1;
		m++;
	      }
	    }
	    else{
	      buff[m]=llnbn[j][k];
	      m++;
	    }
	    break;
	  }
	}
      }
      
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    i=0;
    for (j=0;j<buffsize;j++){
      buff[j]=0;
    }
    
    m=0;
    for (j=0;j<tnbn;j++){
      for (k=0;k<bnmultip[j];k++){
	if (ldn[j][k]==domproc[i]){
	  if (k==0){
	    for (l=0;l<bnmultip[j]-1;l++){
	      buff[m]=llnbn[j][k];
	      buff[m]=0-buff[m]-1;
	      m++;
	    }
	  }
	  else{
	    buff[m]=llnbn[j][k];
	    m++;
	  }
	  break;
	}
      }
    }
    
    
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  fprintf (out,"\n\n KONTROLA 2. kontrola hranicnich uzlu \n");
  for (i=0;i<ncn;i++){
    if (buff[i]<0){
      j=0-buff[i]-1;
    }
    else
      j=buff[i];
    
    fprintf (out,"      %ld %ld",buff[i],lnbn[j]);
  }
  fprintf (out,"\n");
  
  //  number of contributions to the coarse problem
  ncdof=0;
  for (i=0;i<ncn;i++){
    if (buff[i]<0){
      j=0-buff[i]-1;
      k=lnbn[j];
      ndofn=top->give_ndofn (k);
      for (j=0;j<ndofn;j++){
	l=top->give_dof (k,j);
	if (l>0)  ncdof++;
      }
    }
    else{
      k=lnbn[buff[i]];
      ndofn=top->give_ndofn (k);
      for (j=0;j<ndofn;j++){
	l=top->give_dof (k,j);
	if (l>0)  ncdof++;
      }
    }
  }
  
  fprintf (out,"\n\n KONTROLA 3.  pocet stupnu volnosti prispivajicich do hrubeho problemu  %ld",ncdof);
  
  
  //  list of code numbers which are extracted from domain in the FETI method
  //  some code numbers are positive and some of them are negative
  //  it depends on number of subdomain
  //  these code numbers are used on subdomains, master processor contains corresponding array with coarse code numbers
  edofs = new long [ncdof];
  m=0;
  for (i=0;i<ncn;i++){
    if (buff[i]<0){
      j=0-buff[i]-1;
      k=lnbn[j];
      ndofn=top->give_ndofn (k);
      for (j=0;j<ndofn;j++){
	l=top->give_dof (k,j);
	if (l>0){
	  edofs[m]=0-l;
	  m++;
	}
      }
    }
    else{
      k=lnbn[buff[i]];
      ndofn=top->give_ndofn (k);
      for (j=0;j<ndofn;j++){
	l=top->give_dof (k,j);
	if (l>0){
	  edofs[m]=l;
	  m++;
	}
      }
    }
  }
  
  delete [] buff;

  
  fprintf (out,"\n\n\n KONTROLA 4.   extract code numbers     edofs\n");
  for (i=0;i<ncdof;i++){
    fprintf (out,"\n%ld",edofs[i]);
  }
  fprintf (out,"\n");
  
  

  //  coarse code numbers for subdomain contributions
  if (myrank==0){
    
    //  coarse code numbers
    //  ccn[i][j]=k - the j-th contribution from the i-th subdomain goes to the k-th coarse component
    ccn = new long* [nproc];
    for (i=0;i<nproc;i++){
      ccn[i] = new long [ncdofd[i]];
    }
    
    for (i=0;i<nproc;i++){
      l=0;
      for (j=0;j<nbnd[i];j++){
	m=cnbn[i][j];
	for (k=0;k<bnmultip[m];k++){
	  if (ldn[m][k]==i){
	    if (k==0){
	      for (ii=0;ii<bnmultip[m]-1;ii++){
		for (jj=0;jj<ndofncn[m][ii];jj++){
		  if (cncn[m][ii][jj]>0){
		    //ccn[i][l]=0-cncn[m][ii][jj];
		    ccn[i][l]=cncn[m][ii][jj];
		    l++;
		  }
		}
	      }
	    }
	    else{
	      for (jj=0;jj<ndofncn[m][k-1];jj++){
		if (cncn[m][k-1][jj]>0){
		  ccn[i][l]=cncn[m][k-1][jj];
		  l++;
		}
	      }
	      
	    }
	    break;
	  }
	}
      }
      
    }
    
    fprintf (out,"\n\n\n kontrola coarse code numbers    ccn\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain %ld     ncdofd %ld",i,ncdofd[i]);
      for (j=0;j<ncdofd[i];j++){
	fprintf (out,"\n%ld",ccn[i][j]);
      }
    }
    fprintf (out,"\n");

  }  

}

/**
   function orders unknowns in coarse problem and assembles connection among local and coarse unknowns
   
   @param top - pointer to generalized topology of one subdomain
   @param ptop - pointer to parallel topology
   @param domproc - array containing correspondence between processors and subdomains
   @param out - pointer to output file

   JK, 11.7.2005
*/
void feti1::coarse_problem_ordering (gtopology *top,partop *ptop,long *domproc,FILE *out)
{
  //  assembling constraints and code numbers at coarse nodes
  //  generation of unknown numbers in the coarse problem
  coarse_dofs (ptop,out);
  
  //  determination of number of contributing nodes and unknowns to the coarse problem
  number_contributing_nodes_dofs (top,ptop,domproc,out);
  
  //  assembling of contributing unknowns to the coarse problem
  //  assembling arrays for subdomains and array for coarse problem on the master processor
  contributing_nodes_dofs (top,ptop,domproc,out);
}

/**
   function orders unknowns on subdomains with respect to the FETI method
   
   @param top - pointer to domain topology
   
   JK, 31.7.2005
*/
void feti1::subdomain_ordering (gtopology *top,long *ltg,FILE *out)
{
  if(fetiprecond == dirichlet){
    // schur odrdering - see subdomain_matrix
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
    
    for (i=0;i<top->nn;i++){
      if (ltg[i]>-1){
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
    
  }
  else{
    
    ndof = top->codenum_generation (out);
  }
}









/**
   function extracts contributions from local %vector to buffer
   
   @param lv - local %vector
   @param buff - buffer
   
   JK, 31.7.2005
*/
void feti1::local_buff (double *lv,double *buff)
{
  long i,j;
  
  for (i=0;i<ncdof;i++){
    j=edofs[i];
    if (j<0){
      j=0-j-1;
      buff[i]=0.0-lv[j];
    }
    else{
      j=j-1;
      buff[i]=lv[j];
    }
  }
}

/**
   function assembles contributions to the coarse vector from buffer
   
   @param cv - coarse %vector
   @param buff - buffer
   @param nd - number of subdomain
   
   JK, 31.7.2005
*/
void feti1::buff_coarse (double *cv,double *buff,long nd)
{
  long i,j;
  
  for (i=0;i<ncdofd[nd];i++){
    j=ccn[nd][i]-1;
    cv[j]+=buff[i];
  }
}

/**
   function extracts contributions from coarse %vector to buffer
   
   @param cv - coarse %vector
   @param buff - buffer
   @param nd - number of subdomain

   JK, 31.7.2005
*/
void feti1::coarse_buff (double *cv,double *buff,long nd)
{
  long i,j;
  
  for (i=0;i<ncdofd[nd];i++){
    j=ccn[nd][i]-1;
    buff[i]=cv[j];
  }
}

/**
   function assembles contributions from buffer to local %vector
   
   @param lv - local %vector
   @param buff - buffer

   JK, 31.7.2005
*/
void feti1::buff_local (double *lv,double *buff)
{
  long i,j;
  
  for (i=0;i<ncdof;i++){
    j=edofs[i];
    if (j<0){
      j=0-j-1;
      lv[j]-=buff[i];
    }
    else{
      j=j-1;
      lv[j]+=buff[i];
    }
  }
}












/**
   JK, 10.6.2003
*/
/*
void feti1::globcnnum_feti (gtopology *top,long *ltg,long *domproc,FILE *out)
{
  long i,j,k,l,m,n,p,q;
  long nbdof,nbn,tnbn,maxnbn,totmaxnbn,totmaxnbdof;
  long *buff,**lgcor,*nbndom;
  long *nbdofdom,**lbndom,**llcndom,*nninc,**lbndofndom,*ndofnarr,**agcn,**nodomcor;
  MPI_Status stat;
  
  //  number of boundary nodes on one subdomain
  nbn=0;
  //  number of boundary unknowns on one subdomain
  nbdof=0;
  //  maximum number of global node on subdomain
  maxnbn=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      nbdof+=top->give_ndofn (i);  nbn++;
      if (ltg[i]>maxnbn)  maxnbn=ltg[i];
    }
  }
  
  //  data gathering
  buff = new long [3];
  buff[0]=nbn;
  buff[1]=nbdof;
  buff[2]=maxnbn;
  
  if (myrank==0){
    //  array of numbers of boundary nodes [nproc]
    nbndom = new long [nproc];
    //  array of numbers of boundary DOFs [nproc]
    nbdofdom = new long [nproc];
    //  maximum number of boundary nodes on one subdomain
    totmaxnbn=0;
    //  maximum number of boundary DOFs on one subdomain
    totmaxnbdof=0;
    //  total number of boundary nodes
    tnbn=0;
    
    j=domproc[0];
    nbndom[j]=buff[0];
    nbdofdom[j]=buff[1];
    totmaxnbn=buff[0];
    totmaxnbdof=buff[1];
    tnbn=buff[2];
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      nbndom[j]=buff[0];
      nbdofdom[j]=buff[1];
      if (buff[0]>totmaxnbn)    totmaxnbn=buff[0];
      if (buff[1]>totmaxnbdof)  totmaxnbdof=buff[1];
      if (buff[2]>tnbn)         tnbn=buff[2];
    }

    tnbn++;

    buff[0]=totmaxnbn;
    buff[1]=totmaxnbdof;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,3,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }

  }
  else{
    MPI_Send (buff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  totmaxnbn=buff[0];
  totmaxnbdof=buff[1];
  delete [] buff;


  // *********************************
  //  boundary node numbers gathering
  // *********************************
  buff = new long [totmaxnbn];
  for (i=0;i<totmaxnbn;i++){
    buff[i]=0;
  }
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      buff[j]=ltg[i];  j++;
    }
  }
  
  if (myrank==0){
    
    //  lbndom - list of boundary nodes on subdomains [nproc,nbndom]
    lbndom = new long* [nproc];
    for (i=0;i<nproc;i++){
      lbndom[i] = new long [nbndom[i]];
      for (j=0;j<nbndom[i];j++){
	lbndom[i][j]=0;
      }
    }
    
    k=domproc[0];
    for (j=0;j<nbndom[k];j++){
      lbndom[k][j]=buff[j];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,totmaxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      for (j=0;j<nbndom[k];j++){
	lbndom[k][j]=buff[j];
      }
    }
    

    //  control printing
    fprintf (out,"\n\n\n tnbn %ld \n\n",tnbn);
    for (i=0;i<nproc;i++){
      fprintf (out,"\n domain n. %ld, number of boundary nodes (nbndom) %ld\n(lbndom)",i,nbndom[i]);
      for (j=0;j<nbndom[i];j++){
	fprintf (out," %ld",lbndom[i][j]);
      }
    }

    
    
  }
  else{
    MPI_Send (buff,totmaxnbn,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // ****************************************
  //  number of boundary nodes DOFs gathering
  // ****************************************
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      buff[j]=top->give_ndofn (i);  j++;
    }
  }
  
  if (myrank==0){

    //  lbndofndom - list of numbers of DOFs on boundary nodes [nproc,nbndom]
    lbndofndom = new long* [nproc];
    for (i=0;i<nproc;i++){
      lbndofndom[i] = new long [nbndom[i]];
      for (j=0;j<nbndom[i];j++){
	lbndofndom[i][j]=0;
      }
    }
    
    k=domproc[0];
    for (j=0;j<nbndom[k];j++){
      lbndofndom[k][j]=buff[j];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,totmaxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      for (j=0;j<nbndom[k];j++){
	lbndofndom[k][j]=buff[j];
      }
    }
    

    //  control printing
    fprintf (out,"\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n domain %ld (ndofn)\n(lbndofndom)",i);
      for (j=0;j<nbndom[i];j++){
	fprintf (out," %ld",lbndofndom[i][j]);
      }
    }


  }
  else{
    MPI_Send (buff,totmaxnbn,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  
  // **********************************
  //  boundary nodes local code numbers
  // **********************************
  buff = new long [totmaxnbdof];
  for (i=0;i<totmaxnbdof;i++){
    buff[i]=0;
  }

  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]!=-1){
      for (j=0;j<top->give_ndofn (i);j++){
	buff[k]=top->give_dof(i,j);  k++;
      }
    }
  }
  
  if (myrank==0){
    
    //  llcndom - list of local code numbers of boundary nodes [nproc,nbdofdom]
    llcndom= new long* [nproc];
    for (i=0;i<nproc;i++){
      llcndom[i] = new long [nbdofdom[i]];
      for (j=0;j<nbdofdom[i];j++){
	llcndom[i][j]=0;
      }
    }
    
    k=domproc[0];
    for (j=0;j<nbdofdom[k];j++){
      llcndom[k][j]=buff[j];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,totmaxnbdof,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      for (j=0;j<nbdofdom[k];j++){
	llcndom[k][j]=buff[j];
      }
    }
    

    // control printing
    fprintf (out,"\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n domain n. %ld, number of boundary DOFs (nbdofdom)  %ld\n(llcndom)",i,nbdofdom[i]);
      for (j=0;j<nbdofdom[i];j++){
	fprintf (out," %ld",llcndom[i][j]);
      }
    }

    
  }
  else{
    MPI_Send (buff,totmaxnbdof,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  
  //  block on master

  // *************************
  //  number of incidences
  // *************************
  

    
  if (myrank==0){
    


    //  array containing number of incidences [tnbn]
    nninc = new long [tnbn];
    //  array containing number of DOFs on boundary nodes [tnbn]
    ndofnarr = new long [tnbn];
    for (i=0;i<tnbn;i++){
      nninc[i]=0;
      ndofnarr[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	nninc[lbndom[i][j]]++;
	if (ndofnarr[lbndom[i][j]]==0){
	  ndofnarr[lbndom[i][j]]=lbndofndom[i][j];
	}
	else{
	  if (ndofnarr[lbndom[i][j]]!=lbndofndom[i][j]){
	    fprintf (stderr,"\n\n wrong number of DOFs on node, program terminates, (%s, line%d).\n\n",__FILE__,__LINE__);
	    abort ();
	  }
	}
      }
    }
    

    //  control printing
    fprintf (out,"\n\n array of ndofn for boundary nodes (ndofnarr)\n");
    for (i=0;i<tnbn;i++){
      fprintf (out," %ld",ndofnarr[i]);
    }
    fprintf (out,"\n\n number of incidences of boundary nodes (nninc)\n");
    for (i=0;i<tnbn;i++){
      fprintf (out," %ld",nninc[i]);
    }

    
    
    //  auxiliary array of global code numbers

    agcn = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      if (nninc[i]<2){
	fprintf (stderr,"\n\n wrong number of node incidences, program terminates, (%s, line%d).\n\n",__FILE__,__LINE__);
	abort ();
      }

      //fprintf (out,"\n blbec %ld %ld",ndofnarr[i]*(nninc[i]-1),ndofnarr[i]*(nninc[i]-1));

      agcn[i] = new long [ndofnarr[i]*(nninc[i]-1)];

      for (j=0;j<ndofnarr[i]*(nninc[i]-1);j++){
	agcn[i][j]=0;
      }

    }
    
    
    fprintf (out,"\n\n");
    for (i=0;i<nproc;i++){
      l=0;
      for (j=0;j<nbndom[i];j++){
	for (k=0;k<lbndofndom[i][j];k++){
	  //l=llcndom[i][l];
	  //fprintf (out," %ld",lbndom[i][j]);
	  agcn[lbndom[i][j]][k]=llcndom[i][l];  l++;
	}
      }
      //fprintf (out,"\n l po skonceni %ld\n",l);
    }
    
    for (i=0;i<tnbn;i++){
      l=ndofnarr[i];
      for (j=0;j<nninc[i]-2;j++){
	for (k=0;k<ndofnarr[i];k++){
	  agcn[i][l]=agcn[i][k];  l++;
	}
      }
    }
    

    // *******************************
    //  global code numbers generation
    // *******************************
    ndofcp=1;
    for (i=0;i<tnbn;i++){
      l=0;
      for (j=0;j<nninc[i]-1;j++){
	for (k=0;k<ndofnarr[i];k++){
	  if (agcn[i][l]!=0){
	    agcn[i][l]=ndofcp;  ndofcp++;
	  }
	  l++;
	}
      }
    }
    ndofcp--;
    

    // control printing
    fprintf (out,"\n\n\n global code numbers \n");
    for (i=0;i<tnbn;i++){
      l=0;
      fprintf (out,"\n global node number %ld",i);
      for (j=0;j<nninc[i]-1;j++){
	for (k=0;k<ndofnarr[i];k++){
	  fprintf (out," %ld",agcn[i][l]);  l++;
	}
      }
    }

    
    //  node-domain correspondence [tnbn,nninc]
    nodomcor = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      nodomcor[i] = new long [nninc[i]];
      nninc[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	nodomcor[lbndom[i][j]][nninc[lbndom[i][j]]]=i;
	nninc[lbndom[i][j]]++;
      }
    }
    

    //  control printing
    fprintf (out,"\n\n number of incidences of boundary nodes (nninc)\n");
    for (i=0;i<tnbn;i++){
      fprintf (out," %ld",nninc[i]);
    }
    fprintf (out,"\n\n list of node-domain correspondence");
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n node number %ld   ",i);
      for (j=0;j<nninc[i];j++){
	fprintf (out," %ld",nodomcor[i][j]);
      }
    }

    
    //  local-global code numbers correspondence
    lgcor = new long* [nproc];
    glcor = new long* [nproc];
    lggl = new long [nproc];
    maxlggl = 0;
    for (i=0;i<nproc;i++){
      k=0;
      for (j=0;j<nbndom[i];j++){
	l=lbndom[i][j];
	if (nodomcor[l][0]==i){
	  k+=(nninc[l]-1)*ndofnarr[l];
	}
	else{
	  k+=ndofnarr[l];
	}
      }
      lggl[i] = k;
      if (maxlggl<k)  maxlggl=k;
      lgcor[i] = new long [lggl[i]];
      glcor[i] = new long [lggl[i]];
    }
    
    for (i=0;i<nproc;i++){
      m=0;  p=0;
      for (j=0;j<nbndom[i];j++){
	l=lbndom[i][j];
	
	for (k=0;k<nninc[l];k++){
	  if (nodomcor[l][k]==i){
	    if (k==0)  q=0;
	    else       q=ndofnarr[l]*(k-1);
	  }
	}
	
	if (nodomcor[l][0]==i){
	  for (n=0;n<nninc[l]-1;n++){
	    for (k=0;k<ndofnarr[l];k++){
	      lgcor[i][m]=llcndom[i][p];  p++;
	      glcor[i][m]=0-agcn[l][q];  m++;  q++;
	    }
	    p-=ndofnarr[l];
	  }
	  p+=ndofnarr[l];
	}
	else{
	  for (k=0;k<ndofnarr[l];k++){
	    lgcor[i][m]=llcndom[i][p];  p++;
	    glcor[i][m]=agcn[l][q];  m++;  q++;
	  }
	}
      }
    }
    

    // control printing
    fprintf (out,"\n\n\n loc-glob corresp");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain %ld",i);
      for (j=0;j<lggl[i];j++){
	fprintf (out,"\n %ld %ld",lgcor[i][j],glcor[i][j]);
      }
    }

    
  }
  


  buff = new long [2];
  if (myrank==0){
    buff[0]=maxlggl;
    for (i=1;i<nproc;i++){
      buff[1]=lggl[domproc[i]];
      MPI_Send (buff,2,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    nclggl=lggl[domproc[0]];
  }
  else{
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    maxlggl=buff[0];
    nclggl=buff[1];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  buff = new long [maxlggl];
  lcn = new long [nclggl];
  if (myrank==0){
    for (i=1;i<nproc;i++){
      k=domproc[i];
      for (j=0;j<lggl[k];j++){
	buff[j]=lgcor[k][j];
      }
      MPI_Send (buff,maxlggl,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    k=domproc[0];
    for (j=0;j<lggl[k];j++){
      buff[j]=lgcor[k][j];
    }
  }
  else{
    MPI_Recv (buff,maxlggl,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  for (i=0;i<nclggl;i++){
    lcn[i]=buff[i];
  }
  
  delete [] buff;
  

  //  control printing
  fprintf (out,"\n\n number of components of array lcn %ld\n",nclggl);
  for (i=0;i<nclggl;i++){
    fprintf (out," %ld",lcn[i]);
  }


  if (myrank==0){
    // ***********************
    for (i=0;i<tnbn;i++){
      delete [] agcn[i];
    }
    delete [] agcn;
    // ***********************
    for (i=0;i<tnbn;i++){
      delete [] nodomcor[i];
    }
    delete [] nodomcor;
    // ***********************
    for (i=0;i<nproc;i++){
      delete [] llcndom[i];
    }
    delete [] llcndom;
    // ***********************
    for (i=0;i<nproc;i++){
      delete [] lbndofndom[i];
    }
    delete [] lbndofndom;
    // ***********************
    for (i=0;i<nproc;i++){
      delete [] lbndom[i];
    }
    delete [] lbndom;
    // ***********************
    for (i=0;i<nproc;i++){
      delete [] lgcor[i];
    }
    delete [] lgcor;
    // ***********************
    
    delete [] nbndom;
    delete [] nbdofdom;
    delete [] nninc;
    delete [] ndofnarr;
  }

}
*/



/**
   function localizes components from local %vector (belonging to one subdomain)
   to the buffer (which will be usualy sent to master)
   
   @param buff - buffer
   @param lv - local %vector
   
   JK, 10.6.2003
*/
/*
void feti1::locbuff (double *buff,double *lv)
{
  long i,j;
  

//  for (i=0;i<nclggl;i++){
//  j=lcn[i]-1;
//  if (j>-1){
//    buff[i]=lv[j];
//  }
// }


  for (i=0;i<nclggl;i++){
    j=lcn[i];
    if (j>0){
      buff[i]=lv[j-1];
    }
    if (j<0){
      buff[i]=0.0-lv[0-j-1];
    }
  }
}
*/
/**
   function localizes components from buffer to local %vector (belonging to one subdomain)
   components of local array must be set to zero before application ot this function
   
   @param buff - buffer
   @param lv - local %vector
   
   JK, 10.6.2003
*/
/*
void feti1::buffloc (double *buff,double *lv)
{
  long i,j;
  

//  for (i=0;i<nclggl;i++){
//  j=lcn[i]-1;
//  if (j>-1){
//    lv[j]+=buff[i];
//  }
// }


  for (i=0;i<nclggl;i++){
    j=lcn[i];
    if (j>0){
      lv[j-1]+=buff[i];
    }
    if (j<0){
      lv[0-j-1]-=buff[i];
    }
  }
}
*/
/**
   function localizes components from global %vector (belonging to reduced problem)
   to the buffer
   
   @param buff - buffer
   @param gv - global %vector
   @param nd - number of subdomain
   
   JK, 10.6.2003
*/
/*
void feti1::globbuff (double *buff,double *gv,long nd)
{
  long i,j;
  

//for (i=0;i<lggl[nd];i++){
//  j=glcor[nd][i];
//  if (j>0){
//    buff[i]=gv[j-1];
//  }
//  if (j<0){
//    buff[i]=0.0-gv[0-j-1];
//  }
// }
  
  for (i=0;i<lggl[nd];i++){
    j=glcor[nd][i];
    if (j>0){
      buff[i]=gv[j-1];
    }
  }

}
*/
/**
   function localizes components from buffer to global %vector (belonging to reduced problem)
   components of global array must be set to zero before application of this function

   @param buff - buffer
   @param gv - global %vector
   @param nd - number of subdomain
   
   JK, 10.6.2003
*/
/*
void feti1::buffglob (double *buff,double *gv,long nd)
{
  long i,j;
  

//  for (i=0;i<lggl[nd];i++){
//  j=glcor[nd][i];
//  if (j>0){
//    gv[j-1]+=buff[i];
//  }
//  if (j<0){
//    gv[0-j-1]-=buff[i];
//  }
// }

  
  for (i=0;i<lggl[nd];i++){
    j=glcor[nd][i];
    if (j>0){
      gv[j-1]+=buff[i];
    }
  }
}
*/









//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************

/**
   function computes size of matrix H
   
   @param rbm - array containing rigid body motions
   @param maxnrbm - maximum number of rigid body motions
   
   2.12.2001 - origin
   JK, 10.6.2003 - modification
*/
void feti1::hmatrixsize (double */*rbm*/,long *domproc,FILE */*out*/)
{
  long i,j,k;
  MPI_Status stat;
  
  if (myrank==0){
    //  array containing numbers of RBMs
    nrbmdom = new long [nproc];
    //  addresses computed from nrbmdom array
    rbmadr = new long [nproc+1];
    for (i=0;i<nproc+1;i++){
      rbmadr[i]=0;
    }
    
    j=domproc[0];
    maxnrbm=nrbm;  nrbmdom[j]=nrbm;
    for (i=1;i<nproc;i++){
      MPI_Recv(&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (maxnrbm<k)  maxnrbm=k;
      j=domproc[stat.MPI_TAG];
      nrbmdom[j]=k;
    }
    
    //  computation of addresses
    rbmadr[0]=0;
    for (i=0;i<nproc;i++){
      rbmadr[i+1]=rbmadr[i]+nrbmdom[i];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send(&maxnrbm,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    //  size of the matrix H
    hsize=rbmadr[nproc];
  }
  else{
    MPI_Send(&nrbm,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv(&maxnrbm,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  
}


/**
   function assembles %matrix H
   columns are RBM after localization
   
   @param h - %matrix H
   @param rbm - array containing rigid body motions
   @param domproc - domain-processor correspondence
   
   2.12.2001
   JK, 10.6.2003; 31.7.2005 - modification
*/
void feti1::hmatrix (double *h,double *rbm,long *domproc)
{
  long i,j,k,l,buffsize;
  double *buff,*av;
  MPI_Status stat;
  
  buffsize = maxnrbm*maxncdofd;
  buff = new double [buffsize];
  for (i=0;i<buffsize;i++){
    buff[i]=0.0;
  }
  
  //  extraction of RBM from local vectors to buffer
  for (i=0;i<nrbm;i++){
    local_buff (rbm+i*ndof,buff+i*maxncdofd);
  }
  
  fprintf (stdout,"\n myrank %d  ndof %ld     nrbm %ld    maxnrbm %ld     ndofcp %ld   maxncdofd %ld",myrank,ndof,nrbm,maxnrbm,ndofcp,maxncdofd);
  

  if (myrank==0){
    
    av = new double [ndofcp];
    
    //  master contribution
    j=domproc[0];
    for (k=0;k<nrbmdom[j];k++){
      nullv (av,ndofcp);
      buff_coarse (av,buff+k*maxncdofd,j);
      for (l=0;l<ndofcp;l++){
	h[l*hsize+rbmadr[j]+k]=0.0-av[l];
      }
    }
    
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      fprintf (stdout,"\n j %ld  hsize %ld  nrbmdom %ld %ld %ld %ld",j,hsize,nrbmdom[0],nrbmdom[1],nrbmdom[2],nrbmdom[3]);
      fprintf (stdout,"\n j %ld  rbmadr %ld %ld %ld %ld",j,rbmadr[0],rbmadr[1],rbmadr[2],rbmadr[3]);

      fprintf (stdout,"\n j %ld  ncdofd %ld %ld %ld %ld",j,ncdofd[0],ncdofd[1],ncdofd[2],ncdofd[3]);

      for (k=0;k<nrbmdom[j];k++){
	nullv (av,ndofcp);
	//buff_coarse (av,buff+k*maxncdofd,j);
	for (l=0;l<ndofcp;l++){
	  //h[l*hsize+rbmadr[j]+k]=0.0-av[l];
	}
      }

    }
    
    delete [] av;
  }
  else{
    MPI_Send(buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  delete [] buff;
}




/**
   function assembles vector q
   
   @param q - vector q
   @param rbm - array containing rigid body motions
   @param f - right hand side
   
   6.3.2002 - origin
   JK, 10.6.2003 - modification
*/
void feti1::qvector (double *q,double *rbm,double *f,long *domproc)
{
  long i,j,k,l;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxnrbm];
  
  for (i=0;i<nrbm;i++){
    buff[i]=0.0-ss(rbm+i*ndof,f,ndof);
  }
  
  if (myrank==0){
    
    k=domproc[0];
    l=rbmadr[k];
    for (j=0;j<nrbmdom[k];j++){
      q[l]=buff[j];  l++;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnrbm,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      l=rbmadr[k];
      for (j=0;j<nrbmdom[k];j++){
	q[l]=buff[j];  l++;
      }
    }
    
  }
  else{
    MPI_Send(buff,maxnrbm,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  
  delete [] buff;
}

/**
   function computes projection in FETI method
   v_new = v_old - H . (H^T . H)^{-1} . H^T . v_old
   function overwrites vector v by projected vector
   
   @param v - %vector
   @param h - matrix H
   @param h1 - matrix (H^T.H)^{-1}
   @param ndofcp,hsize - number of rows and columns of matrix H
   
   ndofcp - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 7.3.2002
*/
void feti1::feti_projection (double *v,double *h,double *h1)
{
  long i;
  double *p,*q,*w;
  
  p = new double [hsize];
  q = new double [hsize];
  w = new double [ndofcp];

  //  p = H^T.v
  //mtvc (h,v,p,ndofcp,hsize);
  mtxv (h,v,p,ndofcp,hsize);
  
  
  /*
  fprintf (out,"\n\n\n kontrola H^T.v\n");
  for (i=0;i<hsize;i++){
    fprintf (out," %lf",p[i]);
  }
  */

  
  //  q = (H^T.H)^{-1}.p
  mxv (h1,p,q,hsize,hsize);


  /*
  fprintf (out,"\n\n\n kontrola (H^T.H)^{-1}.H^T.v\n");
  for (i=0;i<hsize;i++){
    fprintf (out," %lf",q[i]);
  }
  */


  //  w = H.q
  mxv (h,q,w,ndofcp,hsize);
  
  
  /*
  fprintf (out,"\n\n\n kontrola H.(H^T.H)^{-1}.H^T.v\n");
  for (i=0;i<ndofcp;i++){
    fprintf (out," %lf",w[i]);
  }
  */
  
  //  v_new = v_old - w
  for (i=0;i<ndofcp;i++){
    v[i]-=w[i];
  }
  
  delete [] w;  delete [] q;  delete [] p;
}




void feti1::subdomain_matrix (gmatrix *gm,long */*domproc*/,FILE *out)
{
  long i,j,k;
  long *aux;
  MPI_Status stat;
  long buffsize;
  double *buff,*auxscal;
  
  fprintf(out,"\n\n\nkontrola subdomain_matrix\n\n\n\n kontrola ndof %ld\n\n\n",ndof);
  aux = new long [ndof];
  for(i = 0; i < ndof; i++){
    aux[i]=0;
  }
  
  fprintf(out,"kontrola ncdof %ld\n\n\n",ncdof);
  
  for (i = 0; i < ncdof; i++){
    k = edofs[i];
    //fprintf (out,"edofs[%ld]=%ld\n",i,k);
    if (k<0)
      k=0-k-1;
    else
      k--;
    if (k>-1)
      aux[k]++;
  }
  
  
  fprintf (out,"aux  domene\n");
  for (i = 0; i < ndof; i++){
    fprintf (out," %ld\n",aux[i]);
  }
  
  ndofprec = 0;
  for (i  = 0;i < ndof; i++){
    if (aux[i] > 0)
      ndofprec++;
  }
  fprintf (out,"ndofprec %ld\n",ndofprec);
  
  
  cnprec = new long [ndofprec];
  //cpreccn = new long [ndofprec];
  
  k=0;
  for (i = 0; i < ndof;i++){
    if (aux[i]>0){
      cnprec[k]=i;
      //cpreccn[k]=i+1;
      k++;
    }
  }
  delete [] aux;
  
  /*
  fprintf (out,"Kontrola cpreccn\n");
  for (i = 0; i < ndofprec; i++){
    fprintf (out,"  cpreccn[%ld] = %ld\n",i,cpreccn[i]);
  }
  */
  
  fprintf (out,"Kontrola cnprec\n");
  for (i = 0; i < ndofprec; i++){
    fprintf (out,"  cnprec[%ld] = %ld\n",i,cnprec[i]);
  }
  
  switch (fetiprecond){
  case lumped:{
    smsky = new skyline;
    gm->sky->select_submatrix(smsky,ndofprec,cnprec,out);
    break;
  }
  case dirichlet:{
    // schur ordering - kondenzace vnitrnich stupnu volnosti
    // tvorba pomocne matice
    // kopirovani matice
    smsky = new skyline;
    smsky->n=gm->sky->n;
    //fprintf(out,"smsky->n = %ld gm->sky->n %ld\n",smsky->n,gm->sky->n);
    
    smsky->adr = new long [smsky->n+1];
    for(i = 0; i < smsky->n+1;i++){
      smsky->adr[i]=gm->sky->adr[i];
    }
    /*
    for(i = 0; i <  smsky->n+1; i++){
      fprintf(out,"%ld\n",smsky->adr[i]); 
    }
    */
    
    smsky->adr[smsky->n]= gm->sky->negm;
    smsky->negm =  gm->sky->negm;
    //   fprintf(out,"smsky->n = %ld\n",smsky->negm);
    smsky->a = new double [smsky->negm];
    
    for(i = 0; i < smsky->negm;i++){
      smsky->a[i]=gm->sky->a[i];
    }
    /*
    for(i = 0; i <  smsky->negm; i++){
      fprintf(out,"%le\n",smsky->a[i]); 
    }
    */
    // kondenzace matice
    smdm = new densemat;
    double *x,*y,*av;
    x = new double[ndof];
    y = new double[ndof];
    smdm->a = new double [ndofprec*ndofprec];
    av = new double [ndofprec];
    smdm->negm=ndofprec*ndofprec;
    smdm->n=ndofprec;
    
    smsky->ldlkon_sky (smdm->a,av,x,y,ndofprec,1);
    // mazani pomocne matice
    smsky->~skyline ();
    
    for(i = 0; i < smdm->n*smdm->n; i++){
      fprintf(out,"%le\n",smdm->a[i]); 
    }
    
    delete []x;
    delete []y;
    delete [] av;

    break;
  }
  }
  
//   buffsize=maxncdofd;
//   buff = new double [buffsize];
  
  
//   if(myrank == 0){
    
//     for ( i = 1; i  < nproc; i++){
//       nullv (buff,buffsize);
//       j=domproc[i];
//       coarse_buff (wscalmat,buff,j);
//       MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
      
//       //fprintf(out,"kontrolni tisk buff %ld\n",j);
//       //for (k = 0; k < buffsize; k++){
//       //	fprintf(out,"buff[%ld] %le\n",k,buff[k]);
//       //     }
 
//     }
    
//     nullv (buff,buffsize);
//     j=domproc[0];
//     coarse_buff (wscalmat,buff,j);
//   }
//   else{
//     MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
//   }
  
//   //for(i = 0; i < buffsize; i++){
//   // fprintf(out,"buff[%ld] = %lf\n",i,buff[i]); 
//   //}
  
//   auxscal = new double [ndof];
//   nullv (auxscal,ndof);
  
//   for (i=0;i<ncdof;i++){
//     j=edofs[i];
//     if (j<0){
//       j=0-j-1;
//       if(auxscal[j] == 0){
// 	auxscal[j] = buff[i];
//       }
//       else{
// 	if(auxscal[j] != buff[i]){
// 	  fprintf(out,"pozor chyba v plneni auxscal auxscal[%ld]= %lf buff[%ld] = %lf\n",j,auxscal[j],i,buff[i]);
// 	  auxscal[j] = buff[i];
// 	}
//       }
//     }
//     else{
//       j=j-1;
//       if(auxscal[j] == 0){
// 	auxscal[j] = buff[i];
//       }
//       else{
// 	if(auxscal[j] != buff[i]){
// 	  fprintf(out,"pozor chyba v plneni auxscal auxscal[%ld]= %lf buff[%ld] = %lf\n",j,auxscal[j],i,buff[i]);
// 	  auxscal[j] = buff[i];
// 	}
//       }
//     }
//   }
  
//   //for(i = 0; i < ndof; i++){
//    // fprintf(out,"auxscal[%ld] = %lf\n",i,auxscal[i]); 
//   //}
  
//   // lokalni skalovaci matice
//   lwscalmat = new double[ndofprec];
  
//   for (i = 0; i < ndofprec; i++){
//     lwscalmat[i] = auxscal[cnprec[i]];
//   }
  
  
//   for(i = 0; i < ndofprec; i++){
//     fprintf(out,"lwscalmat[%ld] = %lf\n",i,lwscalmat[i]); 
//   }
  
//   delete []auxscal;
//   delete []buff;
}





void feti1::scaling(double *invect,double *outvect,long n,FILE *out)
{
  double *aux;
  long i;
  aux = new double[n];
  
  for(i = 0; i < n; i++){
    aux[i] = invect[i]/wscalmat[i];
  }
  
  
  fprintf (out,"\n\n\n kontrola skalovani \n");
  for(i = 0; i < n; i++){
    fprintf (out,"%ld  pred %le po %le wscalmat %lf\n",i+1,invect[i],aux[i],wscalmat[i]); 
  }
  
  
  for(i = 0; i < n; i++){
    outvect[i] = aux[i];
  }
  
  delete []aux;
}

void feti1::locscaling (double *invect,double *outvect,FILE *out)
{
  long i;
  double *aux;
  
  aux = new double[ndofprec];
  
  for(i = 0; i < ndofprec; i++){
    aux[i] = invect[i]/lwscalmat[i];
  }
  
  
  fprintf (out,"\n\n\n kontrola skalovani \n");
  for(i = 0; i < ndofprec; i++){
    fprintf (out,"%ld  pred %le po %le lwscalmat %lf\n",i+1,invect[i],aux[i],lwscalmat[i]); 
  }
  
  
  for(i = 0; i < ndofprec; i++){
    outvect[i] = aux[i];
  }
  
  delete []aux;
}

/**

*/
void feti1::lumpedprec (gmatrix */*gm*/,double *dd,double *pp,FILE */*out*/)
{
  
  long i;
  double *d,*p;
  
  d = new double [ndofprec];
  p = new double [ndofprec];
  
  for (i = 0; i < ndofprec; i++){
    d[i]=0.0;
    p[i]=0.0;
  }
  
  for (i=0;i<ndofprec;i++){
    d[i]=dd[cnprec[i]];
  }
  
  //locscaling (d,d,out);
  /*
  for (i = 0;i < ndofprec; i++){
    p[i]=d[i];
  }
  */
  /*
  long j,acb,aci,aci1;
  double s,g;
  for (i = 0; i < smsky->n; i++){
    aci=smsky->adr[i];  aci1=smsky->adr[i+1]-1;
    g=d[i];  s=0.0;  acb=i-aci1+aci;
    printf("aci1 %ld aci %ld\n",aci1,aci);
    for (j=aci1;j>aci;j--){
      s+=smsky->a[j]*d[acb];
      p[acb]+=smsky->a[j]*g;
      printf("acb %ld\n",acb);
      acb++;
    }
    p[i]=s+smsky->a[aci]*g;
  }
  */
  smsky->mxv_sky (d,p);
    
  //gm->gmxv (dd,pp);
  
  //locscaling (p,p,out);
  for (i = 0;i < ndofprec; i++){
    pp[cnprec[i]]=p[i];
  }
  
  delete [] p;
  delete [] d;
  
}


void feti1::dirichletprec (double *dd,double *pp,FILE */*out*/)
{
  long i;
  double *d,*p;
  
  d = new double [ndofprec];
  p = new double [ndofprec];
  
  for (i=0;i<ndofprec;i++){
    d[i]=0.0;
    p[i]=0.0;
  }
  
  for (i=0;i<ndofprec;i++){
    d[i]=dd[cnprec[i]];
  }
  //locscaling (d,d,out);
  /*
  for (i = 0;i < ndofprec; i++){
    p[i]=d[i];
  }
  */
  
  smdm->mxv_dm (d,p);
  
  //locscaling (p,p,out);
  
  for (i=0;i<ndofprec;i++){
    pp[cnprec[i]]=p[i];
  }
  
  delete [] p;
  delete [] d;
}



/**
   function performs modified conjugate gradient method
   
   @param w - vector of Lagrange multipliers
   @param h - matrix H
   @param h1 - matrix (H^T.H)^{-1}
   @param ndofcp,hsize - number of rows and columns of matrix H
   
   ndofcp - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 7.3.2002
   10.6.2003 - modification
*/
void feti1::mpcg (gtopology */*top*/,gmatrix *gm,long *domproc,
		  double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out)
{
  long i,j,k,buffsize;
  long l;
  double nom,denom,alpha,beta;
  double *d,*dd,*g,*p,*pp,*hq,*buff;
  MPI_Status stat;

  //  vectors defined on each subdomain
  dd = new double [ndof];
  pp = new double [ndof];

  
  buffsize=maxncdofd+1;
  buff = new double [buffsize];
  buff[maxncdofd]=0.0;
  
  if (myrank==0){

    //  direction vector allocation
    d = new double [ndofcp];
    //  residuum vector allocation
    g = new double [ndofcp];
    //  auxiliary vector allocation
    p = new double [ndofcp];
    
    //fprintf (stdout,"\n\n\n ndofcp %ld       hsize  %ld",ndofcp,hsize);

    //  initiation of conjugate gradients

    hq = new double [hsize];
    //  (H^T.H)^{-1} q
    mxv (h1,q,hq,hsize,hsize);
    //  H [(H^T.H)^{-1} q]
    mxv (h,hq,w,ndofcp,hsize);
    delete [] hq;

    for (j=1;j<nproc;j++){
      k=domproc[j];
      coarse_buff (w,buff,k);
      MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }

    k=domproc[0];
    coarse_buff (w,buff,k);
    nullv (dd,ndof);
    buff_local (dd,buff);
    // f - \lambda_0
    subv (dd,rhs,ndof);
    
    fprintf (out,"\n");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n zoufalec %4ld    %20.15e",i,dd[i]);
    }

    //  K^+ (f - \lambda_0) = r
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    fprintf (out,"\n");
    for (i=0;i<nrbm;i++){
      fprintf (out,"\n zoufalec    se %4ld",rbmi[i]);
    }
    fprintf (out,"\n");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n zoufalec %4ld    %20.15e",i,pp[i]);
    }

    //nullvr (buff,maxlggl+1);
    local_buff (pp,buff);
    
    
    nullv (g,ndofcp);
    k=domproc[0];
    buff_coarse (g,buff,k);

    for (j=1;j<nproc;j++){
      MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      buff_coarse (g,buff,k);
    }

    
    fprintf (out,"\n\n\n");
    for (i=0;i<ndofcp;i++){
      fprintf (out,"\n zoufalec %4ld    %20.15e",i,g[i]);
    }


    //  P r
    feti_projection (g,h,h1);
    
    //  initialization of direction vector
    for (i=0;i<ndofcp;i++){
      d[i]=0.0-g[i];
    }
    
    //  nominator evaluation
    nom = ss (g,g,ndofcp);
    

  }
  else{

    //  receiving of part of initial approximation of Lagrange multipliers
    MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullv (dd,ndof);
    //  localization from buffer to local vector
    buff_local (dd,buff);
    // f - \lambda_0
    subv (dd,rhs,ndof);
    
    fprintf (out,"\n");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n zoufalec %4ld    %20.15e",i,dd[i]);
    }

    //  K^+ (f - \lambda_0) = r
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    fprintf (out,"\n");
    for (i=0;i<nrbm;i++){
      fprintf (out,"\n zoufalec    se %4ld",rbmi[i]);
    }
    fprintf (out,"\n");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n zoufalec %4ld    %20.15e",i,pp[i]);
    }
    
    //nullvr (buff,maxlggl+1);
    local_buff (pp,buff);
    MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  
  if (myrank==0){

    // *************************
    //  main iteration loop   **
    // *************************
    for (i=0;i<nicg;i++){
      
      // ******************************************
      //   auxiliary vector computation K.d = p  **
      // ******************************************
      for (j=1;j<nproc;j++){
	//nullvr (buff,maxlggl+1);
	k=domproc[j];
	coarse_buff (d,buff,k);
	MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }
      
      
      //  computation of K^+ d on master
      //nullvr (buff,maxlggl+1);
      k=domproc[0];
      coarse_buff (d,buff,k);
      nullv (dd,ndof);
      buff_local (dd,buff);
      
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      //nullvr (buff,maxlggl+1);
      local_buff (pp,buff);
      
      nullv (p,ndofcp);
      k=domproc[0];
      buff_coarse (p,buff,k);
      
      for (j=1;j<nproc;j++){
	MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	k=domproc[stat.MPI_TAG];
	buff_coarse (p,buff,k);
      }
      
      
      //  denominator of alpha
      denom = ss (d,p,ndofcp);
      
      //fprintf (stdout,"\n denominator u alpha   %le",denom);
      
      fprintf (out,"\n\n kontrola citatele a jmenovatele pred alpha %e / %e",nom,denom);
      //fprintf (stderr,"\n\n kontrola citatele a jmenovatele pred alpha %lf / %lf",nom,denom);
      
      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator in modified conjugate gradient method (%s, line %d).\n",__FILE__,__LINE__);
	
	//if (i!=nicg-1){
	buff[maxncdofd]=1.0;
	for (j=1;j<nproc;j++){
	  MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	}
	//}
	break;
      }
      
      //fprintf (stderr," SEM SE NELZE DOSTAT\n");
      
      // *******************
      //   vypocet alpha  **
      // *******************
      alpha = nom/denom;
      
      // **************************************************************
      //  vypocet noveho gradientu g a nove aproximace neznamych x   **
      // **************************************************************
      for (j=0;j<ndofcp;j++){
	w[j]+=alpha*d[j];
	g[j]+=alpha*p[j];
      }
      
      feti_projection (g,h,h1);
      
      
      // *******************************************************
      // *******************************************************
      //  vlozka s predpodminenim
      // *******************************************************
      // *******************************************************
      
      

      
      nullv (p,ndofcp);
      switch (fetiprecond){
      	//without preconditioning
      case nofetiprecond:{
	
	copyv(g,p,ndofcp);
	
	/*
	for (j=1;j<nproc;j++){
	  //nullvr (buff,maxlggl+1);
	  k=domproc[j];
	  coarse_buff (g,buff,k);
	  MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	}
	
	//nullvr (buff,maxlggl+1);
	k=domproc[0];
	coarse_buff (g,buff,k);
	nullv (dd,ndof);
	buff_local (dd,buff);
	
	
	nullv (pp,ndof);
	for (j=0;j<ndof;j++){
	  pp[j]=dd[j];
	}
	
	//nullvr (buff,maxlggl+1);
	local_buff (pp,buff);
	
	nullv (p,ndofcp);
	k=domproc[0];
	buff_coarse (p,buff,k);
	
	for (j=1;j<nproc;j++){
	  MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	  k=domproc[stat.MPI_TAG];
	  buff_coarse (p,buff,k);
	}
	*/
      	break;
      }
      case lumped:{
	// scaling
	scaling(g,g,ndofcp,out);
	
	/*
	fprintf(out,"\n\n\n kontrolni tisk g ndofcp = %ld\n",ndofcp);
	for (j = 0; j < ndofcp; j++){
	  fprintf(out,"g[%ld] %le\n",j,g[j]);
	}
	*/
	// test
	//copyv(g,p,ndofcp);
	
	for (j=1;j<nproc;j++){
	  //nullvr (buff,maxlggl+1);
	  k=domproc[j];
	  coarse_buff (g,buff,k);
	  MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  
	  
	 // fprintf(out,"kontrolni tisk buff %ld\n",j);
	 // for (l = 0; l < buffsize; l++){
	 //   fprintf(out,"buff[%ld] %le\n",l,g[l]);
	 // }
	  
	
	}
	
	//nullvr (buff,maxlggl+1);
	k=domproc[0];
	coarse_buff (g,buff,k);
	nullv (dd,ndof);
	buff_local (dd,buff);
	
	//for (l = 0; l < ndof; l++){
	 // fprintf(out,"dd[%ld] %le\n",l,dd[l]);
	//}
	
	
	nullv (pp,ndof);
	lumpedprec (gm,dd,pp,out);
	
	
	//for (j=0;j<ndof;j++){
	// pp[j]=dd[j];
	//}
	
	
	//for (l = 0; l < ndof; l++){
	 // fprintf(out,"pp[%ld] %le\n",l,pp[l]);
	//}
	
	
	//nullvr (buff,maxlggl+1);
	local_buff (pp,buff);
	
	nullv (p,ndofcp);
	k=domproc[0];
	buff_coarse (p,buff,k);
	
	for (j=1;j<nproc;j++){
	  MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	  k=domproc[stat.MPI_TAG];
	  buff_coarse (p,buff,k);
	  //fprintf(out,"kontrolni tisk buff %ld\n",j);
	 // for (l = 0; l < buffsize; l++){
	  //  fprintf(out,"buff[%ld] %le\n",l,g[l]);
	    //}
	}

	//for (l = 0; l < ndofcp; l++){
	 // fprintf(out,"p[%ld] %le\n",l,p[l]);
	//}
	  
	// scaling
	scaling(p,p,ndofcp,out);
	
	break;
      }
      case dirichlet:{
	// scaling
	scaling(g,g,ndofcp,out);
	
	/*
	fprintf(out,"\n\n\n kontrolni tisk g ndofcp = %ld\n",ndofcp);
	for (j = 0; j < ndofcp; j++){
	  fprintf(out,"g[%ld] %le\n",j,g[j]);
	}
	*/
	
	for (j=1;j<nproc;j++){
	  //nullvr (buff,maxlggl+1);
	  k=domproc[j];
	  coarse_buff (g,buff,k);
	  MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  
	  /*
	  fprintf(out,"kontrolni tisk buff %ld\n",j);
	  for (l = 0; l < buffsize; l++){
	    fprintf(out,"buff[%ld] %le\n",l,g[l]);
	  }
	  */
	
	}
	
	//nullvr (buff,maxlggl+1);
	k=domproc[0];
	coarse_buff (g,buff,k);
	nullv (dd,ndof);
	buff_local (dd,buff);
	/*
	for (l = 0; l < ndof; l++){
	  fprintf(out,"dd[%ld] %le\n",l,dd[l]);
	}
	*/
	
	nullv (pp,ndof);
	dirichletprec (dd,pp,out);
	
	/*
	for (j=0;j<ndof;j++){
	  pp[j]=dd[j];
	}
	*/
	
	for (l = 0; l < ndof; l++){
	  fprintf(out,"pp[%ld] %le\n",l,pp[l]);
	}
	
	
	//nullvr (buff,maxlggl+1);
	local_buff (pp,buff);
	
	nullv (p,ndofcp);
	k=domproc[0];
	buff_coarse (p,buff,k);
	
	for (j=1;j<nproc;j++){
	  MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	  k=domproc[stat.MPI_TAG];
	  buff_coarse (p,buff,k);
	  /*fprintf(out,"kontrolni tisk buff %ld\n",j);
	  for (l = 0; l < buffsize; l++){
	    fprintf(out,"buff[%ld] %le\n",l,g[l]);
	    }*/
	}

	for (l = 0; l < ndofcp; l++){
	  fprintf(out,"p[%ld] %le\n",l,p[l]);
	}
	  
	// scaling
	scaling(p,p,ndofcp,out);
	break;
      }
      }
      
    
    
	
      feti_projection (p,h,h1);
      
      // *******************************************************
      // *******************************************************
      //   konec vlozky s predpodminenim
      // *******************************************************
      // *******************************************************
      
      
      denom = nom;
      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator of beta in function (%s, line %d).\n",__FILE__,__LINE__);
	
	if (i!=nicg-1){
	  buff[maxncdofd]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      //nom = ss (g,g,ndofcp);
      nom=ss(p,g,ndofcp);
      
      fprintf (stdout,"\n iteration  %4ld        norm g   %e",i,sqrt(nom));
      fprintf (out,"\n iteration  %4ld        norm g   %e",i,sqrt(nom));
      
      fprintf (out,"\n kontrola citatele a jmenovatele pred betou  %e / %e",nom,denom);
      //fprintf (stderr,"\n kontrola citatele a jmenovatele pred betou  %lf / %lf",nom,denom);
      
      if (sqrt(nom)<errcg){
	if (i!=nicg-1){
	  buff[maxncdofd]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      
      //  computation of beta coefficient
      beta = nom/denom;
      
      //  new direction vector
      for (j=0;j<ndofcp;j++){
	//d[j]=beta*d[j]-g[j];
	d[j]=beta*d[j]-p[j];
      }
      
    }
    
    anicg=i;  aerrcg=nom;
    
    fprintf (out,"\n\n\n\n kontrola Lagrangeovych multiplikatoru \n");
    for (i=0;i<ndofcp;i++){
      fprintf (out,"\n lagr. mult %4ld    %e",i,w[i]);
    }
    
    delete [] d;
    delete [] p;
    delete [] g;
  }
  else{
    // SLAVES
    for (i=0;i<nicg;i++){
      
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (buff[maxncdofd]>0.5){  break;  }
      
      nullv (dd,ndof);
      //  localization from buffer to local vector
      buff_local (dd,buff);
      
      //  computation of K^+ d on master
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      //nullvr (buff,maxlggl+1);
      local_buff (pp,buff);
      
      MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
      // *******************************************************
      // *******************************************************
      //  vlozka s predpodminenim
      // *******************************************************
      // *******************************************************
      
      switch (fetiprecond){
      	// without preconditioning
      case nofetiprecond:{

	/*
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	if (buff[maxncdofd]>0.5){  break;  }
	
	nullv (dd,ndof);
	buff_local (dd,buff);

	nullv (pp,ndof);
	for (j=0;j<ndof;j++){
	  pp[j]=dd[j];
	}
	//nullvr (buff,maxlggl+1);
	local_buff (pp,buff);
	MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	*/
	break;
      }
	// lumped preconditioning
      case lumped:{
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	if (buff[maxncdofd]>0.5){  break;  }
	
	nullv (dd,ndof);
	buff_local (dd,buff);
	
	//for (l = 0; l < ndof; l++){
	 // fprintf(out,"dd[%ld] %le\n",l,dd[l]);
	//}
	
	nullv (pp,ndof);
	lumpedprec (gm,dd,pp,out);
	//for (l = 0; l < ndof; l++){
	 // fprintf(out,"pp[%ld] %le\n",l,pp[l]);
	//}

	
	//for (j=0;j<ndof;j++){
	 // pp[j]=dd[j];
	//}
	
	//nullvr (buff,maxlggl+1);
	local_buff (pp,buff);
	MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	break;
      }
	// dirichlet preconditioning
      case dirichlet:{
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	if (buff[maxncdofd]>0.5){  break;  }
	
	nullv (dd,ndof);
	buff_local (dd,buff);
	/*
	for (l = 0; l < ndof; l++){
	  fprintf(out,"dd[%ld] %le\n",l,dd[l]);
	}
	*/
	nullv (pp,ndof);
	dirichletprec (dd,pp,out);
	
	for (l = 0; l < ndof; l++){
	  fprintf(out,"pp[%ld] %le\n",l,pp[l]);
	}

	/*
	for (j=0;j<ndof;j++){
	  pp[j]=dd[j];
	}
	*/
	//nullvr (buff,maxlggl+1);
	local_buff (pp,buff);
	MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	break;
      }
      }
      // *******************************************************
      // *******************************************************
      //   konec vlozky s predpodminenim
      // *******************************************************
      // *******************************************************
      
    }
    
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] dd;
  delete [] pp;
  delete [] buff;
}

/**
   function computes displacements from Lagrange multipliers
   function is used in FETI method
   
   @param w - array containing Lagrange multipliers
   @param d - array containing displacements
   @param f - array containing load forces
   @param rbm - array containing rigid body motions
   @param h - matrix H (columns are rigid body motions)
   @param h1 - matrix (H^T . H)^{-1}
   @param ndofcp - number of Lagrange multipliers
   @param hsize - number of columns of matrix H
   
   13.3.2002 - origin
   JK, 14.6.2003 - modification
*/
void feti1::lagrmultdispl (gtopology */*top*/,gmatrix *gm,long *domproc,
			   double *w,double *d,double *f,double *rbm,long *rbmi,
			   double *h,double *h1,FILE *out)
{
  long i,j,k,l,buffsize;
  double *g,*pp,*av,*gamma,*buff,*u;
  MPI_Status stat;
  
  pp = new double [ndof];
  av = new double [maxnrbm];
  
  buffsize=maxncdofd;
  buff = new double [buffsize];
  
  if (myrank==0){
    g = new double [ndofcp];
    
    for (j=1;j<nproc;j++){
      coarse_buff (w,buff,domproc[j]);
      MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }
    
    coarse_buff (w,buff,domproc[0]);
    nullv (pp,ndof);
    buff_local (pp,buff);
    subv (pp,f,ndof);
    
    gm->ldl_feti (d,pp,nrbm,rbmi,zero);
    
    local_buff (d,buff);
    
    nullv (g,ndofcp);
    buff_coarse (g,buff,domproc[0]);
    
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      buff_coarse (g,buff,domproc[stat.MPI_TAG]);
    }
  }
  else{
    MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullv (pp,ndof);
    buff_local (pp,buff);
    
    subv (pp,f,ndof);
    
    gm->ldl_feti (d,pp,nrbm,rbmi,zero);
    
    local_buff (d,buff);

    MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  if (myrank==0){
    
    fprintf (out,"\n");
    for (i=0;i<ndofcp;i++){
      fprintf (out,"\n g  %20.15le",g[i]);
    }

    gamma = new double [hsize];
    u = new double [hsize];
    
    mtxv (h,g,u,ndofcp,hsize);
    mxv (h1,u,gamma,hsize,hsize);
    cmulv (-1.0,gamma,hsize);
    
    for (j=1;j<nproc;j++){
      k=domproc[j];
      l=rbmadr[k];
      for (i=0;i<nrbmdom[k];i++){
	av[i] = gamma[l];  l++;
      }

      fprintf (out,"\n");
      for (i=0;i<nrbmdom[k];i++){
	fprintf (out,"\n av %20.15le",av[i]);
      }
      

      MPI_Send (av,maxnrbm,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }
    
    
    
    
    k=domproc[0];
    l=rbmadr[k];
    for (i=0;i<nrbmdom[k];i++){
      av[i] = gamma[l];  l++;
    }
    
    fprintf (out,"\n");
    for (i=0;i<nrbmdom[k];i++){
      fprintf (out,"\n av %20.15le",av[i]);
    }
    
    delete [] u;
    delete [] gamma;
    delete [] g;
  }
  else{
    MPI_Recv (av,maxnrbm,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  
  cmulv (-1.0,d,ndof);
  
  for (i=0;i<ndof;i++){
    for (j=0;j<nrbm;j++){
      d[i]+=rbm[j*ndof+i]*av[j];
    }
  }
  
  
  delete [] buff;
  delete [] av;  delete [] pp;
}






void feti1::solve_system (gtopology */*top*/,gmatrix *gm,
			  long *domproc,double */*lhs*/,double */*rhs*/,FILE *out)
{
  long i,j;
  long *rbmi;
  double *rbm,*h,*q,*hh,*ih,*ee,*lm;
  time_t t1,t2,t3,t4,t5,t6;
  

  //  rigid body motion array allocation
  rbm = new double [enrbm*ndof];
  memset (rbm,0,enrbm*ndof*sizeof(double));
  //  array of RBM indices
  rbmi = new long [enrbm];
  for (i=0;i<enrbm;i++){
    rbmi[i]=-1;
  }
  
  // assemble matrix for preconditioning
  switch (fetiprecond){
  case nofetiprecond:{
    fprintf(out,"\n\n\n No preconditioning is used for FETI method\n\n");
    break;
  }
  case lumped:{
    fprintf(out,"\n\n\n Lumped preconditioning is used for FETI method\n\n");
    subdomain_matrix (gm,domproc,out);
    break;
  }
  case dirichlet:{
    fprintf(out,"\n\n\n Dirichlet preconditioning is used for FETI method\n\n");
    subdomain_matrix (gm,domproc,out);
    break;
  }
  }
  
  
  t1 = time (NULL);
  
  //gm->printmat (out);

  //  computation of kernel of system matrix
  gm->kernel (rbm,nrbm,rbmi,enrbm,lithr,3);

  //gm->printmat (out);

  t2 = time (NULL);

  fprintf (stdout,"\n proc %d  nrbm   %ld",myrank,nrbm);
  fprintf (stdout,"\n proc %d  enrbm  %ld",myrank,enrbm);
  fprintf (stdout,"\n proc %d  lithr  %e",myrank,lithr);
  fprintf (stdout,"\n proc %d  ndof   %ld",myrank,ndof);
  
  fprintf (out,"\n proc %d  nrbm   %ld",myrank,nrbm);
  fprintf (out,"\n proc %d  enrbm  %ld",myrank,enrbm);
  fprintf (out,"\n proc %d  lithr  %e",myrank,lithr);
  fprintf (out,"\n proc %d  ndof   %ld",myrank,ndof);
  

  //  investigation of H matrix sizes
  hmatrixsize (rbm,domproc,out);
  


  fprintf (out,"\n\n\n kontrola RBM");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n");
    for (j=0;j<nrbm;j++){
      fprintf (out," %f",rbm[j*ndof+i]);
    }
  }
  

  
  //if (myrank==0){
  //hsize=rbmadr[nproc];
  //
  //fprintf (stdout,"\n\n hsize = 0 \n");
  //fprintf (stderr,"\n\n hsize = 0 \n");
  //}
  
  if (myrank==0){
    hsize=rbmadr[nproc];
    
    if (hsize==0){
      fprintf (stdout,"\n\n hsize = 0 \n");
      fprintf (stderr,"\n\n hsize = 0 \n");
      abort ();
    }
    
    h = new double [ndofcp*hsize];
    q = new double [hsize];
    lm = new double [ndofcp];
    
    fprintf (out,"\n hsize   %ld",hsize);
    fprintf (out,"\n ndofcp   %ld",ndofcp);
    
    for (i=0;i<=nproc;i++){
      fprintf (out,"\n rbmadr %4ld   %4ld",i,rbmadr[i]);
    }
    
  }
  

  //  assembling of H matrix
  hmatrix (h,rbm,domproc);
    
  /*
  if (myrank==0){
    fprintf (out,"\n\n\n kontrola matice H \n");
    for (i=0;i<ndofcp;i++){
      fprintf (out,"\n %4ld",i);
      for (j=0;j<hsize;j++){
	fprintf (out,"  %f",h[i*hsize+j]);
      }
    }
  }

  //  assembling of q vector
  qvector (q,rbm,rhs,domproc);
  
  
  if (myrank==0){
    fprintf (out,"\n\n\n kontrola vektoru q \n");
    for (j=0;j<hsize;j++){
      fprintf (out,"\n %4ld   %f",j,q[j]);
    }
  }
  
  t3 = time (NULL);
  
  //  H^T . H and (H^T . H)^{-1.0} evaluation
  if (myrank==0){
    ih = new double [hsize*hsize];
    hh = new double [hsize*hsize];
    ee = new double [hsize*hsize];
    for (i=0;i<hsize;i++){
      for (j=0;j<hsize;j++){
	ee[i*hsize+j]=0.0;
      }
      ee[i*hsize+i]=1.0;
    }
    
    mtxm (h,h,hh,ndofcp,hsize,hsize);
    
    fprintf (out,"\n\n\n kontrola matice HH \n");
    for (i=0;i<hsize;i++){
      fprintf (out,"\n %4ld ",i);
      for (j=0;j<hsize;j++){
	fprintf (out," %lf",hh[i*hsize+j]);
      }
    }

    gemp (hh,ih,ee,hsize,hsize,zero,1);
    
    //mtm (h,h,hh,ndofcp,hsize,hsize);
    //mm (hh,ih,ee,hsize,hsize,hsize);
    
    
    fprintf (out,"\n\n\n kontrola matice HHinv \n");
    for (i=0;i<hsize;i++){
      fprintf (out,"\n %4ld ",i);
      for (j=0;j<hsize;j++){
	fprintf (out," %lf",ih[i*hsize+j]);
      }
    }
    
    delete [] ee;
  }
  
  //  modified conjugate gradient method
  
  t4 = time (NULL);
  

  mpcg (top,gm,domproc,lm,rhs,q,h,ih,rbmi,out);
  
  t5 = time (NULL);
  
  lagrmultdispl (top,gm,domproc,lm,lhs,rhs,rbm,rbmi,h,ih,out);

  t6 = time (NULL);


  fprintf (out,"\n\n evaluation of kernel             %ld",t2-t1);
  fprintf (out,"\n\n matrix H computation             %ld",t3-t2);
  fprintf (out,"\n\n decomposition of H               %ld",t4-t3);
  fprintf (out,"\n\n modified conjug. gradient m.     %ld",t5-t4);
  fprintf (out,"\n\n displacements computation        %ld",t6-t5);

  fprintf (stdout,"\n\n evaluation of kernel             %ld",t2-t1);
  fprintf (stdout,"\n\n matrix H computation             %ld",t3-t2);
  fprintf (stdout,"\n\n decomposition of H               %ld",t4-t3);
  fprintf (stdout,"\n\n modified conjug. gradient m.     %ld",t5-t4);
  fprintf (stdout,"\n\n displacements computation        %ld",t6-t5);


  MPI_Status stat;
  long *buff;
  buff = new long[5];
  
  buff[0]=t2-t1;
  buff[1]=t3-t2;
  buff[2]=t4-t3;
  buff[3]=t5-t4;
  buff[4]=t6-t5;

  if (myrank==0){

    fprintf (out,"\n\n\n\n\n");
    j=domproc[0];
    fprintf (out,"\n\n\n Domain %ld",j);
    fprintf (out,"\n evaluation of kernel             %ld",buff[0]);
    fprintf (out,"\n matrix H computation             %ld",buff[1]);
    fprintf (out,"\n decomposition of H               %ld",buff[2]);
    fprintf (out,"\n modified conjug. gradient m.     %ld",buff[3]);
    fprintf (out,"\n displacements computation        %ld",buff[4]);
    
    
    j=domproc[0];
    fprintf (stdout,"\n\n\n Domain %ld",j);
    fprintf (stdout,"\n evaluation of kernel             %ld",buff[0]);
    fprintf (stdout,"\n matrix H computation             %ld",buff[1]);
    fprintf (stdout,"\n decomposition of H               %ld",buff[2]);
    fprintf (stdout,"\n modified conjug. gradient m.     %ld",buff[3]);
    fprintf (stdout,"\n displacements computation        %ld",buff[4]);

    for (i=1;i<nproc;i++){
      MPI_Recv (buff,5,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];

      fprintf (out,"\n\n\n Domain %ld",j);
      fprintf (out,"\n evaluation of kernel             %ld",buff[0]);
      fprintf (out,"\n matrix H computation             %ld",buff[1]);
      fprintf (out,"\n decomposition of H               %ld",buff[2]);
      fprintf (out,"\n modified conjug. gradient m.     %ld",buff[3]);
      fprintf (out,"\n displacements computation        %ld",buff[4]);
      
      fprintf (stdout,"\n\n\n Domain %ld",j);
      fprintf (stdout,"\n evaluation of kernel             %ld",buff[0]);
      fprintf (stdout,"\n matrix H computation             %ld",buff[1]);
      fprintf (stdout,"\n decomposition of H               %ld",buff[2]);
      fprintf (stdout,"\n modified conjug. gradient m.     %ld",buff[3]);
      fprintf (stdout,"\n displacements computation        %ld",buff[4]);

    }
  }
  else{
    MPI_Send(buff,5,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  */


}
