#include "feti1n.h"
#include <string.h>
#include <math.h>
#include <time.h>

feti1::feti1(int np,int mr,int nd)
{
  nproc=np;  myrank=mr;  ndom=nd;
  
  ndof=0;  nbn=0;
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
   JK, 10.6.2003
*/
void feti1::globcnnum_feti (gtopology *top,long *ltg,long *domproc,FILE *out)
{
  long i,j,k,l,m,n,p,q;
  long maxnbn,gndof,totmaxnbn,totmaxnbdof;
  long *buff,**lgcor,**glcor,*lggl;
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
    gndof=1;
    for (i=0;i<tnbn;i++){
      l=0;
      for (j=0;j<nninc[i]-1;j++){
	for (k=0;k<ndofnarr[i];k++){
	  if (agcn[i][l]!=0){
	    agcn[i][l]=gndof;  gndof++;
	  }
	  l++;
	}
      }
    }
    
    
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
  
  
  
}




/**
   function localizes components from local %vector (belonging to one subdomain)
   to the buffer (which will be usualy sent to master)
   
   @param buff - buffer
   @param lv - local %vector
   
   JK, 10.6.2003
*/
void feti1::locbuff (double *buff,double *lv)
{
  long i,j;
  
  for (i=0;i<nclggl;i++){
    j=lcn[i]-1;
    if (j>-1){
      buff[i]=lv[j];
    }
  }
}
/**
   function localizes components from buffer to local %vector (belonging to one subdomain)
   
   @param buff - buffer
   @param lv - local %vector
   
   JK, 10.6.2003
*/
void feti1::buffloc (double *buff,double *lv)
{
  long i,j;
  
  for (i=0;i<nclggl;i++){
    j=lcn[i]-1;
    if (j>-1){
      lv[j]+=buff[i];
    }
  }
}
/**
   function localizes components from global %vector (belonging to reduced problem)
   to the buffer
   
   @param buff - buffer
   @param gv - global %vector
   @param nd - number of subdomain
   
   JK, 10.6.2003
*/
void feti1::globbuff (double *buff,double *gv,long nd)
{
  long i,j;
  
  for (i=0;i<lggl[nd];i++){
    j=glcor[i];
    if (j>0){
      buff[i]=gv[j-1];
    }
    if (j<0){
      buff[i]=0.0-gv[0-j-1];
    }
  }
}
/**
   function localizes components from buffer to global %vector (belonging to reduced problem)
   
   @param buff - buffer
   @param gv - global %vector
   @param nd - number of subdomain
   
   JK, 10.6.2003
*/
void feti1::buffglob (double *buff,double *gv,long nd)
{
  long i,j;
  
  for (i=0;i<lggl[nd];i++){
    j=glcor[i];
    if (j>0){
      gv[j-1]+=buff[i];
    }
    if (j<0){
      gv[0-j-1]-=buff[i];
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
   
   2.12.2001 - origin
   JK, 10.6.2003 - modification
*/
void feti1::hmatrixsize (double *rbm,long *domproc,FILE *out)
{
  long i;
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
    for (i=1;i<nproc;i++){
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
   function assembles matrix H
   columns are RBM after localization
   
   @param h - matrix H
   @param rbm - array containing rigid body motions
   
   2.12.2001
   JK, 10.6.2003 - modification
*/
void feti1::hmatrix (double *h,double *rbm)
{
  long i,j,k,l,m;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxnrbm*maxlggl];
  for (i=0;i<maxnrbm*maxlggl;i++){
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
    //locglobfeti (top,buff+i*ngdof,rbm+i*ndof);
    locbuff (buff+i*maxlggl,rbm+i*ndof);
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

    j=domproc[0];
    for (k=0;k<nrbmdom[j];k++){
      buffglob (buff+k*maxlggl,h+(rbmadr[j]+k)*ngdof,j);
    }
    
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxnrbm*maxlggl,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      j=domproc[stat.MPI_TAG];
      for (k=0;k<nrbmdom[j];k++){
	buffglob (buff+k*maxlggl,h+(rbmadr[j]+k)*ngdof,j);
      }
      
    }
  }
  else{
    MPI_Send(buff,maxnrbm*maxlggl,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
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
void feti1::qvector (double *q,double *rbm,double *f)
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
   10.6.2003 - modification
*/
void feti1::mpcg (gtopology *top,gmatrix *gm,
		  double *w,double *rhs,double *q,double *h,double *h1,long *rbmi,FILE *out)
{
  long i,j;
  double nom,denom,alpha,beta;
  double *d,*dd,*g,*p,*pp,*hq,*buff;
  MPI_Status stat;
  
  //  direction vector allocation
  d = new double [ngdof+1];
  d[ngdof]=0.0;

  //  residuum vector allocation
  g = new double [ngdof];
  //  auxiliary vector allocation
  p = new double [ngdof+1];
  p[ngdof]=0.0;



  //  vectors defined on each subdomain
  dd = new double [ndof];
  pp = new double [ndof];

  
  
  buff = new double [maxlggl+1];

  if (myrank==0){
    
    //  initiation of conjugate gradients
    
    hq = new double [hsize];
    //  (H^T.H)^{-1} q
    mv (h1,q,hq,hsize,hsize);
    //  H [(H^T.H)^{-1} q]
    mv (h,hq,w,ngdof,hsize);
    delete [] hq;
    
    for (j=1;j<nproc;j++){
      k=domproc[j];
      globbuff (buff,w,k);
      MPI_Send (buff,maxlggl+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
    }

    k=domproc[0];
    globbuff (buff,w,k);
    nullvr (dd,ndof);
    buffloc (buff,dd);
    subv (dd,rhs,ndof);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru prave strany dd pred funkci funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",dd[i]);
    }
    */
    
    
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru reseni pp funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",pp[i]);
    }
    */

    nullvr (buff,maxlggl+1);
    locbuff (buff,pp);
    
    
    nullvr (g,ngdof);
    k=domproc[0];
    buffglob (buff,g,k);

    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxlggl+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      buffglob (buff,g,k);
    }
    
    
    /*
    fprintf (out,"\n\n\n  inicializace    kontrola vektoru reziduii v inicializacni fazi pred projekci \n");
    for (i=0;i<ngdof;i++){
      fprintf (out," %le",g[i]);
    }
    */
    
    
    //  P r
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
    

    /*************************/
    /*  main iteration loop  */
    /*************************/
    for (i=0;i<nicg;i++){
      
      //fprintf (stdout,"\n\n JEDEME GRADIENTY, proc %d, iterace %ld",myrank,i);

      /******************************************/
      /*  auxiliary vector computation K.d = p  */
      /******************************************/
      
      
      for (j=1;j<nproc;j++){
	nullvr (buff,maxlggl+1);
	k=domproc[j];
	globbuff (buff,d,k);
	MPI_Send (buff,maxlggl+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }
      
      
      //  computation of K^+ d on master
      nullvr (buff,maxlggl+1);
      k=domproc[0];
      globbuff (buff,d,k);
      nullvr (dd,ndof);
      buffloc (buff,dd);
      
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
      
      
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      
      /*
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru reseni pp funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      nullvr (buff,maxlggl+1);
      locbuff (buff,pp);
      
      nullvr (p,ngdof);
      k=domproc[0];
      buffglob (buff,p,k);
      
      for (j=1;j<nproc;j++){
	MPI_Recv(buff,maxlggl+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	k=domproc[stat.MPI_TAG];
	buffglob (buff,p,k);
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
	  buff[maxlggl]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (buff,maxlggl+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
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
	nullvr (buff,maxlggl+1);
	k=domproc[j];
	globbuff (buff,g,k);
	MPI_Send (buff,maxlggl+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
      }

      nullvr (buff,maxlggl+1);
      k=domproc[0];
      globbuff (buff,g,k);
      nullvr (dd,ndof);
      buffloc (buff,dd);

      nullvr (pp,ndof);
      //lumpedprec (dd,pp);
      for (j=0;j<ndof;j++){
	pp[j]=dd[j];
      }

      nullvr (buff,maxlggl+1);
      locbuff (buff,pp);
      
      nullvr (p,ngdof);
      k=domproc[0];
      buffglob (buff,p,k);

      for (j=1;j<nproc;j++){
	MPI_Recv(buff,maxlggl+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	k=domproc[stat.MPI_TAG];
	buffglob (buff,p,k);
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
	  buff[maxlggl]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (buff,maxlggl+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
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
	  buff[maxlggl]=1.0;
	  for (j=1;j<nproc;j++){
	    MPI_Send (buff,maxlggl+1,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);
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
    //  receiving of part of initial approximation of Lagrange multipliers
    MPI_Recv (buff,maxlggl+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    nullvr (dd,ndof);
    //  localization from buffer to local vector
    buffloc (buff,dd);
    // f - \lambda_0
    subv (dd,rhs,ndof);
    
    //  K^+ (f - \lambda_0) = r
    gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
    
    
    /*
    fprintf (out,"\n\n\n inicializace --- kontrola vektoru reseni pp funkci ldl_feti_sky \n");
    for (i=0;i<ndof;i++){
      fprintf (out," %le",pp[i]);
    }
    */
    
    nullvr (buff,maxlggl+1);
    locbuff (buff,pp);
    MPI_Send (buff,maxlggl+1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    
    
        
    for (i=0;i<nicg;i++){

      //fprintf (stdout,"\n\n JEDEME GRADIENTY, proc %d, iterace %ld",myrank,i);
      
      MPI_Recv (buff,maxlggl+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (buff[maxlggl]>0.5){  break;  }
      
      nullvr (dd,ndof);
      //  localization from buffer to local vector
      buffloc (buff,dd);
      
      /*
      fprintf (out,"\n\n\n krok %ld         kontrola vektoru prave strany pred funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",dd[j]);
      }
      */
      
      
      gm->ldl_feti (pp,dd,nrbm,rbmi,zero);
      
      /*
      fprintf (out,"\n\n\n krok %ld          kontrola vektoru reseni pp funkci ldl_feti_sky \n",i);
      for (j=0;j<ndof;j++){
	fprintf (out," %le",pp[j]);
      }
      */
      
      nullvr (buff,maxlggl+1);
      locbuff (buff,pp);
      
      MPI_Send (buff,maxlggl+1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
      //*******************************************************
      //*******************************************************
      //*******************************************************
      //  vlozka s predpodminenim
      //*******************************************************
      //*******************************************************
      //*******************************************************

      MPI_Recv (buff,maxlggl+1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      if (buff[maxlggl]>0.5){  break;  }
      
      nullvr (dd,ndof);
      buffloc (buff,dd);
      nullvr (pp,ndof);
      //lumpedprec (dd,pp);
      for (j=0;j<ndof;j++){
	pp[j]=dd[j];
      }
      nullvr (buff,maxlggl+1);
      locbuff (buff,pp);
      MPI_Send (buff,maxlggl+1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
      
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
    
    for (j=0;j<nproc;j++){
      k=rbmadr[j];
      for (i=0;i<nrbmdom[j];i++){
	av[i]=gamma[k];  k++;
      }
      if (domproc[j]!=0)
	MPI_Send (av,hsize,MPI_DOUBLE,domproc[j],myrank,MPI_COMM_WORLD);
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
  //  array of RBM indices
  rbmi = new long [enrbm];
  for (i=0;i<enrbm;i++){
    rbmi[i]=-1;
  }
  
  fprintf (out,"\n enrbm   %ld",enrbm);
  fprintf (out,"\n ndof    %ld",ndof);
  
  t1 = time (NULL);
  
  //  computation of kernel of system matrix
  gm->kernel (rbm,nrbm,rbmi,enrbm,lithr,3);
  
  t2 = time (NULL);

  fprintf (stdout,"\n proc %d  nrbm   %ld",myrank,nrbm);
  fprintf (stdout,"\n proc %d  enrbm  %ld",myrank,enrbm);
  fprintf (stdout,"\n proc %d  lithr  %e",myrank,lithr);
  
  fprintf (out,"\n proc %d  nrbm   %ld",myrank,nrbm);
  fprintf (out,"\n proc %d  enrbm  %ld",myrank,enrbm);
  fprintf (out,"\n proc %d  lithr  %e",myrank,lithr);
  

  //  investigation of H matrix sizes
  hmatrixsize (rbm,domproc,out);
  
  
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
  hmatrix (h,rbm);
  
  
  
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
  qvector (q,rbm,rhs);
  


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



  
  // provizorni konec

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
  
}
