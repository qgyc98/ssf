#include "feti1.h"
#include <string.h>
#include <math.h>
#include <time.h>

feti1::feti1(int np,int mr,int nd)
{
  nproc=np;  myrank=mr;  ndom=nd;
  
  ndof=0;  nbn=0;  totmaxndofn=0;
  enrbm=0;  lithr=0.0;  nrbm=0;  hsize=0;
  nicg=0;  anicg=0;  errcg=0.0;  aerrcg=0.0;
  zero=0.0;
  
  ngdof=0;
  
  lcngcn=NULL;
  inc=NULL;
  nodnum=NULL;
  
  nbndom=NULL;
  nrbmdom=NULL;
  rbmadr=NULL;
}

feti1::~feti1()
{
  delete [] lcngcn;
  delete [] inc;
  delete [] nodnum;

  delete [] nbndom;
  delete [] nrbmdom;
  delete [] rbmadr;
}

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
void feti1::globcnnum_feti (gtopology *top,long *ltg,long *domproc,FILE *out)
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
    nbndom[domproc[0]]=nbn;

    //  maximum number of degrees of freedom on node
    totmaxndofn=maxndofn;
    //  maximum number of boundary nodes
    maxnbn=nbn;
    //  maximum global number of node
    maxgnn=gnn;
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      nbndom[domproc[stat.MPI_TAG]]=buff[1];
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
      MPI_Send (buff,3,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);


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
    
    k=domproc[0];
    for (j=0;j<nbndom[k];j++){
      noddom[buff[j]][k]=1;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
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
    MPI_Send (buff,maxnbn,MPI_LONG,0,myrank,MPI_COMM_WORLD);
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
      MPI_Send (&maxinc,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
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
    for (j=0;j<nbndom[domproc[0]];j++){
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
      for (j=0;j<nbndom[domproc[stat.MPI_TAG]];j++){
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
    MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,0,myrank,MPI_COMM_WORLD);
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
    
    ii=domproc[0];
    
    
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
      ii=domproc[stat.MPI_TAG];
      

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


      MPI_Send (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,domproc[ii],myrank,MPI_COMM_WORLD);
      
      
    }
  }
  else{
    MPI_Send (buff,maxnbn*(maxinc*totmaxndofn+1),MPI_LONG,0,myrank,MPI_COMM_WORLD);
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
void feti1::locglobfeti (gtopology *top,double *gv,double *lv)
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
void feti1::globlocfeti (gtopology *top,double *gv,double *lv)
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
   @param ngdof - number of global DOFs = number of Lagrange multipliers
   
   2.12.2001
*/
void feti1::hmatrixsize (double *rbm,long &maxnrbm,long *domproc)
{
  long i;
  long *ibuff;
  MPI_Status stat;
  
  ibuff = new long [3];
  ibuff[0]=nrbm;
  ibuff[1]=0;
  ibuff[2]=0;
  
  if (myrank==0){
    nrbmdom = new long [nproc];
    rbmadr = new long [nproc+1];
    for (i=0;i<nproc+1;i++){
      rbmadr[i]=0;
    }
    maxnrbm=nrbm;  nrbmdom[domproc[0]]=nrbm;
    for (i=1;i<nproc;i++){
      MPI_Recv(ibuff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (maxnrbm<ibuff[0])  maxnrbm=ibuff[0];
      nrbmdom[domproc[stat.MPI_TAG]]=ibuff[0];
    }
    
    rbmadr[0]=0;
    //  chyba, melo by byt i=0
    for (i=1;i<nproc;i++){
      rbmadr[i+1]=rbmadr[i]+nrbmdom[i];
      
      //fprintf (stderr,"\n i %ld    nrbmdom %ld      rbmadr %ld",i,nrbmdom[i],rbmadr[i+1]);


    }
    
    ibuff[0]=maxnrbm;
    ibuff[1]=ngdof;
    ibuff[2]=rbmadr[nproc];
    for (i=1;i<nproc;i++){
      MPI_Send(ibuff,3,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    ngdof=ngdof;
    hsize=rbmadr[nproc];
  }
  else{
    MPI_Send(ibuff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv(ibuff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    maxnrbm=ibuff[0];  ngdof=ibuff[1];  hsize=ibuff[2];
  }
  delete [] ibuff;
  
}


/**
   function assembles matrix H
   columns are RBM after localization
   
   @param h - matrix H
   @param rbm - array containing rigid body motions
   @param maxnrbm - maximum number of rigid body motions
   @param ngdof - number of global DOFs = number of Lagrange multipliers
   
   2.12.2001
*/
void feti1::hmatrix (gtopology *top,double *h,double *rbm,long maxnrbm,long *domproc)
{
  long i,j,k,l,m;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxnrbm*ngdof];
  for (i=0;i<maxnrbm*ngdof;i++){
    buff[i]=0.0;
  }


  /*
  fprintf (out,"\n\n\n kontrola RBM PRED LOCGLOBFETI");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n");
    for (j=0;j<ppd->nrbm;j++){
      fprintf (out," %lf",rbm[j*ndof+i]);
    }
  }
  */


  for (i=0;i<nrbm;i++){
    locglobfeti (top,buff+i*ngdof,rbm+i*ndof);
  }
  
  
  /*
  fprintf (out,"\n\n\n KONTROLA RBM PO LOCGLOBFETI");
  for (i=0;i<ngdof;i++){
    fprintf (out,"\n");
    for (j=0;j<ppd->nrbm;j++){
      fprintf (out," %lf",buff[j*ngdof+i]);
    }
  }
  */


  
  if (myrank==0){
    hsize = rbmadr[nproc];
    

    m=0;
    for (j=0;j<nrbmdom[domproc[0]];j++){
      l=rbmadr[domproc[0]]+j;
      for (k=0;k<ngdof;k++){
	h[l]=0.0-buff[m];
	l+=hsize;  m++;
      }

      //fprintf (stderr,"\n konecne l je rovno %ld",l);

    }



    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnrbm*ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      

      m=0;
      for (j=0;j<nrbmdom[domproc[stat.MPI_TAG]];j++){
	l=rbmadr[domproc[stat.MPI_TAG]]+j;
	for (k=0;k<ngdof;k++){
	  h[l]=0.0-buff[m];
	  l+=hsize;  m++;
	}

	//fprintf (stderr,"\n j %ld    konecne l je rovno %ld",j,l);

      }

    }
  }
  else{
    MPI_Send(buff,maxnrbm*ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);

}



/**
   function assembles vector q
   
   @param q - vector q
   @param rbm - array containing rigid body motions
   @param f - right hand side
   @param maxnrbm - maximum number of rigid body motions
   
   6.3.2002
*/
void feti1::qvector (double *q,double *rbm,double *f,long maxnrbm,long *domproc)
{
  long i,j,k,l,m;
  double *buff;
  MPI_Status stat;
  

  buff = new double [maxnrbm];
  
  for (i=0;i<nrbm;i++){
    buff[i]=0.0-ss(rbm+i*ndof,f,ndof);
  }
  

  if (myrank==0){
    

    l=rbmadr[domproc[0]];  m=0;
    for (j=0;j<nrbmdom[domproc[0]];j++){
      q[l]=buff[m];  l++;  m++;
    }

    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnrbm,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      

      l=rbmadr[k];  m=0;
      for (j=0;j<nrbmdom[k];j++){
	q[l]=buff[m];  l++;  m++;
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
   @param ngdof,hsize - number of rows and columns of matrix H
   
   ngdof - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 7.3.2002
*/
void feti1::feti_projection (double *v,double *h,double *h1)
{
  long i;
  double *p,*q,*w;
  
  p = new double [hsize];
  q = new double [hsize];
  w = new double [ngdof];

  //  p = H^T.v
  //mtvc (h,v,p,ngdof,hsize);
  mtv (h,v,p,ngdof,hsize);
  
  
  /*
  fprintf (out,"\n\n\n kontrola H^T.v\n");
  for (i=0;i<hsize;i++){
    fprintf (out," %lf",p[i]);
  }
  */

  
  //  q = (H^T.H)^{-1}.p
  mv (h1,p,q,hsize,hsize);


  /*
  fprintf (out,"\n\n\n kontrola (H^T.H)^{-1}.H^T.v\n");
  for (i=0;i<hsize;i++){
    fprintf (out," %lf",q[i]);
  }
  */


  //  w = H.q
  mv (h,q,w,ngdof,hsize);
  
  
  /*
  fprintf (out,"\n\n\n kontrola H.(H^T.H)^{-1}.H^T.v\n");
  for (i=0;i<ngdof;i++){
    fprintf (out," %lf",w[i]);
  }
  */
  
  //  v_new = v_old - w
  for (i=0;i<ngdof;i++){
    v[i]-=w[i];
  }
  
  delete [] w;  delete [] q;  delete [] p;
}




/**
   function performs modified conjugate gradient method
   
   @param w - vector of Lagrange multipliers
   @param h - matrix H
   @param h1 - matrix (H^T.H)^{-1}
   @param ngdof,hsize - number of rows and columns of matrix H
   
   ngdof - number of Lagrange multipliers
   hsize - total number of all rigid body motions
   
   JK, 7.3.2002
*/
void feti1::mpcg (gtopology *top,gmatrix *gm,
		  double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out)
{
  long i,j;
  double nom,denom,alpha,beta;
  double *d,*dd,*g,*p,*pp,*buff;
  MPI_Status stat;
  
  //  direction vector allocation
  d = new double [ngdof+1];
  d[ngdof]=0.0;

  dd = new double [ndof];
  //  residuum vector allocation
  g = new double [ngdof];
  //  auxiliary vector allocation
  p = new double [ngdof+1];
  p[ngdof]=0.0;

  pp = new double [ndof];
  
  
  /*
  fprintf (stdout,"\n\n pocet iteraci %ld",ni);
  fprintf (out,"\n\n\n\n kontrola RHS\n");
  for (i=0;i<ndof;i++){
    fprintf (out," %lf",rhs[i]);
  }

  fprintf (out,"\n\n\n kontrola indexu zavislych rovnic");
  for (i=0;i<ppd->enrbm;i++){
    fprintf (out,"   %ld",rbmi[i]);
  }
  */
  
  /*
  fprintf (out,"\n\n\n kontrola matice tuhosti po eliminaci");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n");
    for (j=Sm_sky->adr[i];j<Sm_sky->adr[i+1];j++){
      fprintf (out," %lf",Sm_sky->a[j]);
    }
  }
  */
  
  
  if (myrank==0){
    
    buff = new double [hsize];
    mv (h1,q,buff,hsize,hsize);
    mv (h,buff,w,ngdof,hsize);
    delete [] buff;
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru q \n");
    for (i=0;i<hsize;i++){
      fprintf (out," %le",q[i]);
    }

    fprintf (out,"\n\n\n lagrangeovy multiplikatory pred ldl_feti_sky \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",w[i]);
    }
    */
    
    
    buff = new double [ngdof+1];
    
    for (j=1;j<nproc;j++){
      MPI_Send (w,ngdof,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }
    nullvr (dd,ndof);
    globlocfeti (top,w,dd);
    subv (dd,rhs,ndof);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru prave strany dd pred funkci funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",dd[i]);
    }
    */
    
    
    //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru reseni pp funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",pp[i]);
    }
    */
    
    
    nullvr (g,ngdof);
    locglobfeti (top,g,pp);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (g,buff,ngdof);
    }
    
    
    /*
    fprintf (out,"\n\n\n  inicializace    kontrola vektoru reziduii v inicializacni fazi pred projekci \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",g[i]);
    }
    */
    
    
    feti_projection (g,h,h1);
    
    
    /*
    fprintf (out,"\n\n\n  inicializace    kontrola vektoru reziduii v inicializacni fazi po projekci \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",g[i]);
    }
    */
    
    
    //  initialization of direction vector
    for (i=0;i<ngdof;i++){
      d[i]=0.0-g[i];
    }
    
    //  nominator evaluation
    nom = ss (g,g,ngdof);
    

    //fprintf (stdout,"\n\n nominator %lf",nom);

    
    /*************************/
    /*  main iteration loop  */
    /*************************/
    for (i=0;i<nicg;i++){
      
      //fprintf (stdout,"\n\n JEDEME GRADIENTY, proc %d, iterace %ld",myrank,i);

      /******************************************/
      /*  auxiliary vector computation K.d = p  */
      /******************************************/
      
      
      for (j=1;j<nproc;j++){
	MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }
      

      //  computation of K^+ d on master
      nullvr (dd,ndof);
      globlocfeti (top,d,dd);
      
      /*
      fprintf (out,"\n\n\n krok %ld        kontrola vektoru prave strany pred funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",dd[j]);
      }
      
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru reseni pp PRED  PRED funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      
      //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru reseni pp funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      
      nullvr (p,ngdof);
      locglobfeti (top,p,pp);
      
      
      for (j=1;j<nproc;j++){
	MPI_Recv(buff,ngdof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	addv (p,buff,ngdof);
      }
      
      
      /*
      fprintf (out,"\n\n\n krok %ld        kontrola pomocneho vektoru K.d\n",i);
      for (j=0;j<ngdof;j++){
	fprintf (out," %le",p[j]);
      }
      */
      /*
      fprintf (out,"\n\n\n krok %ld        kontrola pomocneho vektoru d\n",i);
      for (j=0;j<ngdof;j++){
	fprintf (out," %lf",d[j]);
      }
      */
      
      
      //  denominator of alpha
      denom = ss (d,p,ngdof);
      
      //fprintf (stdout,"\n denominator u alpha   %le",denom);

      fprintf (out,"\n\n kontrola citatele a jmenovatele pred alpha %e / %e",nom,denom);
      //fprintf (stderr,"\n\n kontrola citatele a jmenovatele pred alpha %lf / %lf",nom,denom);

      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator in modified conjugate gradient method (%s, line %d).\n",__FILE__,__LINE__);
	
	fprintf (stderr,"\n\n POKUD JSME TADY, TAK");
	
	if (i!=nicg-1){
	  d[ngdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      //fprintf (stderr," SEM SE NELZE DOSTAT\n");

      /*******************/
      /*  vypocet alpha  */
      /*******************/
      alpha = nom/denom;
      
      /**************************************************************/
      /*  vypocet noveho gradientu g a nove aproximace neznamych x  */
      /**************************************************************/
      for (j=0;j<ngdof;j++){
	w[j]+=alpha*d[j];
	g[j]+=alpha*p[j];
      }
      
      feti_projection (g,h,h1);
      
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  vlozka s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************
      
      
      for (j=1;j<nproc;j++){
	MPI_Send (g,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }

      nullvr (dd,ndof);
      globlocfeti (top,g,dd);
      nullvr (pp,ndof);
      //lumpedprec (dd,pp);
      for (j=0;j<ndof;j++){
	pp[j]=dd[j];
      }
      nullvr (p,ngdof);
      locglobfeti (top,p,pp);

      for (j=1;j<nproc;j++){
	MPI_Recv(buff,ngdof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	addv (p,buff,ngdof);
      }

      feti_projection (p,h,h1);

      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  konec vlozky s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************


      denom = nom;
      if (fabs(denom)<zero){
	fprintf (stderr,"\n\n zero denominator of beta in function (%s, line %d).\n",__FILE__,__LINE__);
	
	if (i!=nicg-1){
	  d[ngdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      //nom = ss (g,g,ngdof);
      nom=ss(p,g,ngdof);
      
      fprintf (stdout,"\n iteration  %4ld        norm g   %e",i,sqrt(nom));
      fprintf (out,"\n iteration  %4ld        norm g   %e",i,sqrt(nom));
      
      fprintf (out,"\n kontrola citatele a jmenovatele pred betou  %e / %e",nom,denom);
      //fprintf (stderr,"\n kontrola citatele a jmenovatele pred betou  %lf / %lf",nom,denom);
      
      if (sqrt(nom)<errcg){
	if (i!=nicg-1){
	  d[ngdof]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (d,ngdof+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
	  }
	}
	break;
      }
      
      
      //  computation of beta coefficient
      beta = nom/denom;
      
      //  new direction vector
      for (j=0;j<ngdof;j++){
	//d[j]=beta*d[j]-g[j];
	d[j]=beta*d[j]-p[j];
      }
      
    }
    
    anicg=i;  aerrcg=nom;
    
    fprintf (out,"\n\n\n\n kontrola Lagrangeovych multiplikatoru \n");
    for (i=0;i<ngdof;i++){
      fprintf (out,"\n lagr. mult %4ld    %e",i,w[i]);
    }
    
  }
  else{
    
    MPI_Recv (d,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullvr (dd,ndof);
    globlocfeti (top,d,dd);
    subv (dd,rhs,ndof);
    
    
    /*
      fprintf (out,"\n\n\n inicializace --- kontrola vektoru prave strany dd pred funkci funkci ldl_feti_sky \n");
      for (i=0;i<ndof;i++){
      fprintf (out," %le",dd[i]);
      }
    */
	
	
    //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru reseni pp funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",pp[i]);
    }
    */
    

    nullvr (p,ngdof);
    locglobfeti (top,p,pp);
    MPI_Send (p,ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    

    for (i=0;i<nicg;i++){

      //fprintf (stdout,"\n\n JEDEME GRADIENTY, proc %d, iterace %ld",myrank,i);
      
      MPI_Recv (d,ngdof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (d[ngdof]>0.5){  break;  }
      
      nullvr (dd,ndof);
      globlocfeti (top,d,dd);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru prave strany pred funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",dd[j]);
      }
      */
      
      
      //Sm_sky->ldl_feti_sky (pp,dd,ppd->nrbm,rbmi,ppd->zero);
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld          kontrola vektoru reseni pp funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      
      nullvr (p,ngdof);
      locglobfeti (top,p,pp);
      
      
      MPI_Send (p,ngdof+1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  vlozka s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************

      MPI_Recv (p,ngdof+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (p[ngdof]>0.5){  break;  }
      
      nullvr (dd,ndof);
      globlocfeti (top,p,dd);
      nullvr (pp,ndof);
      //lumpedprec (dd,pp);
      for (j=0;j<ndof;j++){
	pp[j]=dd[j];
      }
      nullvr (p,ngdof);
      locglobfeti (top,p,pp);
      MPI_Send (p,ngdof+1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  konec vlozky s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************
      
    }
  }
  
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
   @param ngdof - number of Lagrange multipliers
   @param hsize - number of columns of matrix H
   
   13.3.2002
*/
void feti1::lagrmultdispl (gtopology *top,gmatrix *gm,long *domproc,
			     double *w,double *d,double *f,double *rbm,long *rbmi,
			     double *h,double *h1)
{
  long i,j,k;
  double *g,*p,*pp,*av,*gamma;
  MPI_Status stat;
  
  g = new double [ngdof];
  p = new double [ngdof];
  pp = new double [ndof];
  av = new double [hsize];
  

  if (myrank==0){
    
    
    /*
    fprintf (out,"\n\n\n kontrola Lagrangeovych multiplikatoru \n");
    for (j=0;j<ngdof;j++){
      fprintf (out," %lf",w[j]);
    }
    */
    
    
    for (j=1;j<nproc;j++){
      MPI_Send (w,ngdof,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }
    nullvr (pp,ndof);
    globlocfeti (top,w,pp);

    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    

    subv (pp,f,ndof);
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    
    /*
    fprintf (out,"\n\n\n kontrola matice tuhosti po eliminaci");
    for (i=0;i<ndof;i++){
      fprintf (out,"\n");
      for (j=Sm_sky->adr[i];j<Sm_sky->adr[i+1];j++){
	fprintf (out," %lf",Sm_sky->a[j]);
      }
    }
    */
    
    //Sm_sky->ldl_feti_sky (d,pp,pdd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (d,pp,nrbm,rbmi,zero);
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }

    fprintf (out,"\n\n\n kontrola vektoru d \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    

    nullvr (g,ngdof);
    locglobfeti (top,g,d);
    for (j=1;j<nproc;j++){
      MPI_Recv(p,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (g,p,ngdof);
    }
    
    gamma = new double [hsize];
    
    mtv (h,g,av,ngdof,hsize);
    mv (h1,av,gamma,hsize,hsize);
    scalarray (gamma,-1.0,hsize);
    
    for (j=1;j<nproc;j++){
      k=rbmadr[domproc[j]];
      for (i=0;i<nrbmdom[domproc[j]];i++){
	av[i]=gamma[k];  k++;
      }
      MPI_Send (av,hsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }

  }

  else{
    MPI_Recv (p,ngdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullvr (pp,ndof);
    globlocfeti (top,p,pp);

    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */
    
    subv (pp,f,ndof);
    
    /*
    fprintf (out,"\n\n\n kontrola vektoru pp \n");
    for (j=0;j<ndof;j++){
      fprintf (out," %lf",pp[j]);
    }
    */


    //Sm_sky->ldl_feti_sky (d,pp,ppd->nrbm,rbmi,ppd->zero);
    gm->ldl_feti (d,pp,nrbm,rbmi,zero);
    nullvr (g,ngdof);
    locglobfeti (top,g,d);
    MPI_Send (g,ngdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    
    
    MPI_Recv (av,hsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
  }
  
  
  /*
  fprintf (out,"\n\n kontrola posunu bez RBM\n");
  for (i=0;i<ndof;i++){
    fprintf (out," %lf",d[i]);
  }
  */
  
  scalarray (d,-1.0,ndof);
  
  for (i=0;i<ndof;i++){
    for (j=0;j<nrbm;j++){
      d[i]+=rbm[j*ndof+i]*av[j];
    }
  }
  
  
  /*
  fprintf (out,"\n\n kontrola posunu s RBM\n");
  for (i=0;i<ndof;i++){
    fprintf (out," %lf",d[i]);
  }
  */
  
  
  delete [] av;  delete [] pp;  delete [] p;  delete [] g;
}

/**
   function makes preconditioning by lumped matrix
   
   @param v - input vector
   @param pv - preconditioned vector
   
   8.4.2002
*/
void feti1::lumpedprec (gtopology *top,double *v,double *pv)
{
/*
  long i,j,k,ndofe;
  ivector cn;
  vector av,lv;
  matrix sm;
  
  for (i=0;i<top->ne;i++){
    ndofe=top->get_ndofe (i);
    allocm (ndofe,ndofe,sm);
    allocv (ndofe,cn);
    allocv (ndofe,av);
    allocv (ndofe,lv);
    stiffmat (i,0,sm);
    top->give_code_numbers (i,cn.a);
    for (j=0;j<ndofe;j++){
      av[j]=0.0;
      k=cn[j]-1;
      if (k>-1){
	av[j]=v[k];
      }
    }
    mxv (sm,av,lv);
    for (j=0;j<ndofe;j++){
      k=cn[j]-1;
      if (k>-1){
	pv[k]+=lv[j];
      }
    }
    destrm (sm);
    destrv (lv);
    destrv (av);
    destrv (cn);
  }
*/
}









void feti1::solve_system (gtopology *top,gmatrix *gm,
			  long *domproc,double *lhs,double *rhs,FILE *out)
{
  long i,j,maxnrbm;
  long *rbmi;
  double *rbm,*h,*q,*hh,*ih,*ee,*lm;
  time_t t1,t2,t3,t4,t5,t6;
  
  //  rigid body motion array allocation
  rbm = new double [enrbm*ndof];
  memset (rbm,0,enrbm*ndof*sizeof(double));
  rbmi = new long [enrbm];
  for (i=0;i<enrbm;i++){
    rbmi[i]=-1;
  }
  
  fprintf (out,"\n enrbm   %ld",enrbm);
  fprintf (out,"\n ndof    %ld",ndof);
  
  t1 = time (NULL);
  
  //  computation of kernel of system matrix
  //Sm_sky->ker (rbm,nse,ppd->se,ppd->ense,ppd->lithr,2);
  //Sm_sky->ker (rbm,ppd->nrbm,rbmi,ppd->enrbm,ppd->lithr,3);
  gm->kernel (rbm,nrbm,rbmi,enrbm,lithr,3);
  
  t2 = time (NULL);

  fprintf (stdout,"\n proc %d  nrbm   %ld",myrank,nrbm);
  fprintf (stdout,"\n proc %d  enrbm  %ld",myrank,enrbm);
  fprintf (stdout,"\n proc %d  lithr  %e",myrank,lithr);
  
  fprintf (out,"\n proc %d  nrbm   %ld",myrank,nrbm);
  fprintf (out,"\n proc %d  enrbm  %ld",myrank,enrbm);
  fprintf (out,"\n proc %d  lithr  %e",myrank,lithr);
  
  //  investigation of H matrix sizes
  hmatrixsize (rbm,maxnrbm,domproc);
  
  
  /*
  fprintf (out,"\n\n\n kontrola RBM");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n");
    for (j=0;j<nrbm;j++){
      fprintf (out," %f",rbm[j*ndof+i]);
    }
  }
  */
  
  
  if (myrank==0){
    hsize=rbmadr[nproc];
    h = new double [ngdof*hsize];
    q = new double [hsize];
    lm = new double [ngdof];
    
    fprintf (out,"\n hsize   %ld",hsize);
    
    for (i=0;i<=nproc;i++){
      fprintf (out,"\n rbmadr %4ld   %4ld",i,rbmadr[i]);
    }
    
  }
  
  
  //  assembling of H matrix
  hmatrix (top,h,rbm,maxnrbm,domproc);
  
  
  /*
  if (myrank==0){
    fprintf (out,"\n\n\n\n kontrola matice H");
    for (i=0;i<ngdof;i++){
      fprintf (out,"\n");
      for (j=0;j<hsize;j++){
	fprintf (out," %f",h[i*hsize+j]);
      }
    }
  }
  */
  
  
  
  //  assembling of q vector
  qvector (q,rbm,rhs,maxnrbm,domproc);
  

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
    
    //  old
    //mtmccr (h,h,hh,ngdof,hsize,hsize);
    
    mtm (h,h,hh,ngdof,hsize,hsize);
    
    /*
    if (myrank==0){
      fprintf (out,"\n\n\n\n kontrola matice (H^T.H)");
      for (i=0;i<hsize;i++){
	fprintf (out,"\n");
	for (j=0;j<hsize;j++){
	  fprintf (out," %f",hh[i*hsize+j]);
	}
      }
    }
    */
    
    gemp (hh,ih,ee,hsize,hsize,zero,1);
    
    
    
    /*
    if (myrank==0){
      fprintf (out,"\n\n\n\n kontrola matice (H^T.H)^-1");
      for (i=0;i<hsize;i++){
	fprintf (out,"\n");
	for (j=0;j<hsize;j++){
	  fprintf (out," %f",ih[i*hsize+j]);
	}
      }
    }
    */
    
    
    mtm (h,h,hh,ngdof,hsize,hsize);
    mm (hh,ih,ee,hsize,hsize,hsize);
    
    /*
      if (myrank==0){
      fprintf (out,"\n\n\n\n kontrola matice jednotkove matice");
      for (i=0;i<hsize;i++){
	fprintf (out,"\n");
	for (j=0;j<hsize;j++){
	  fprintf (out," %f",ee[i*hsize+j]);
	}
      }
    }
    */
    
    
    delete [] ee;
  }
  
  //  modified conjugate gradient method
  
  t4 = time (NULL);
  
  mpcg (top,gm,lm,rhs,q,h,ih,rbmi,out);
  
  t5 = time (NULL);

  fprintf (out,"\n\n hsize %ld",hsize);
  
  lagrmultdispl (top,gm,domproc,lm,lhs,rhs,rbm,rbmi,h,ih);

  t6 = time (NULL);

  fprintf (out,"\n\n evaluation of kernel             %ld",t2-t1);
  fprintf (out,"\n\n matrix H computation             %ld",t3-t2);
  fprintf (out,"\n\n decomposition of H               %ld",t4-t3);
  fprintf (out,"\n\n modified conjug. gradient m.     %ld",t5-t4);
  fprintf (out,"\n\n displacements computation        %ld",t6-t5);
  

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
    j=domproc[myrank];
    fprintf (out,"\n\n\n Domain %ld",j);
    fprintf (out,"\n evaluation of kernel             %ld",buff[0]);
    fprintf (out,"\n matrix H computation             %ld",buff[1]);
    fprintf (out,"\n decomposition of H               %ld",buff[2]);
    fprintf (out,"\n modified conjug. gradient m.     %ld",buff[3]);
    fprintf (out,"\n displacements computation        %ld",buff[4]);
    
    
    j=domproc[myrank];
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
      
      
      j=domproc[myrank];
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

}
