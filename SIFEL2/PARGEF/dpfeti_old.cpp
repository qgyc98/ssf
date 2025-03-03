#include "dpfeti.h"
#include <string.h>
#include <math.h>
#include <time.h>

dpfeti::dpfeti(int np,int mr,int nd)
{
  nproc=np;  myrank=mr;  ndom=nd;

  ndof=0;  nidof=0;  ncdof=0;  tncdof=0;  tnmdof=0;
  ncn=0;  nbn=0;  maxncdof=0;  totmaxndofn=0;  maxinc=0;
  
  nicg=0;  anicg=0;  errcg=0.0;  aerrcg=0.0;
  zero=0.0;

  lcngcn=NULL;
  inc=NULL;
  nodnum=NULL;
  
  nbndom=NULL;
  ncndom=NULL;
  ncdofdom=NULL;
  masgcn=NULL;
}

dpfeti::~dpfeti()
{
  delete [] lcngcn;
  delete [] inc;
  delete [] nodnum;

  delete [] nbndom;
  delete [] ncndom;
  delete [] ncdofdom;
  delete [] masgcn;

}


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
void dpfeti::dpfetiordering (gtopology *top,long *ltg)
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
  nidof=ndof-1;
  
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
  
  ncdof=ndof-nidof;
}

/**
   function creates code numbers for DP-FETI method
   
   @param top - topology object pointer
   @param out - output stream
   
   JK, 15.1.2002
*/
void dpfeti::globcnnum_dpfeti (gtopology *top,long *ltg,long *domproc,FILE *out)
{
  long i,j,k,l,m,ii,jj,kk,ll;
  long gnbn,gncnn,maxndofn,maxnbn,maxncn,maxgnbn,maxgncn;
  long *buff,*gcn,**noddom,*nodinc,**gncn;
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
  gncnn=0;
  for (i=0;i<top->nn;i++){
    if (maxndofn<top->give_ndofn (i))  maxndofn = top->give_ndofn (i);
    
    if (ltg[i]>-1){
      nbn++;
      if (gnbn<ltg[i])  gnbn=ltg[i];
    }
    
    if (ltg[i]<-1){
      ncn++;
      if (gncnn<0-ltg[i]-2)  gncnn=0-ltg[i]-2;
    }
  }
  
  
  buff = new long [6];
  buff[0]=maxndofn;
  buff[1]=nbn;
  buff[2]=gnbn;
  buff[3]=ncn;
  buff[4]=gncnn;
  buff[5]=ncdof;

  if (myrank==0){
    //  array of numbers of boundary nodes
    nbndom = new long [nproc];
    nbndom[ndom]=nbn;
    
    //  array of numbers of corner nodes
    ncndom = new long [nproc];
    ncndom[ndom]=ncn;
    
    //  array of numbers of corner unknowns
    ncdofdom = new long [nproc];
    ncdofdom[ndom]=ncdof;
    
    //  maximum number of degrees of freedom on node
    totmaxndofn=maxndofn;
    //  maximum number of boundary nodes
    maxnbn=nbn;
    //  maximum global number of boundary node
    maxgnbn=gnbn;
    //  maximum number of corner nodes
    maxncn=ncn;
    //  maximum global number of corner node
    maxgncn=gncnn;
    //  maximum number of corner unknowns
    maxncdof=ncdof;
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,6,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
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
    fprintf (out,"\n maxncdof %ld",maxncdof);
    fprintf (out,"\n\n pocty velicin na podoblastech");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nbn %4ld  ncn %4ld  ncdof %4ld",nbndom[i],ncndom[i],ncdofdom[i]);
    }
    
    
    buff[0]=totmaxndofn;
    buff[1]=maxncn;
    buff[2]=maxncdof;
    buff[3]=maxnbn;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,6,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,6,MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,6,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
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
      buff[k*(totmaxndofn+1)+totmaxndofn]=0-ltg[i]-2;
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
    
    /*
    fprintf (out,"\n\n kontrola sebranych kodovych cisel \n");
    for (i=0;i<maxgncn;i++){
      fprintf (out,"\n corner node %ld   ",i);
      for (j=0;j<totmaxndofn;j++){
	fprintf (out," %ld",gncn[i][j]);
      }
    }
    */
    
    //  computation of total number of corner unknowns
    tncdof=1;
    for (i=0;i<maxgncn;i++){
      for (j=0;j<totmaxndofn;j++){
	if (gncn[i][j]>0){
	  gncn[i][j]=tncdof;  tncdof++;
	}
      }
    }
    tncdof--;
    
    /*
    fprintf (out,"\n\n kontrola rohovych kodovych cisel \n");
    for (i=0;i<maxgncn;i++){
      fprintf (out,"\n corner node %ld   ",i);
      for (j=0;j<totmaxndofn;j++){
	fprintf (out," %ld",gncn[i][j]);
      }
    }
    */
    
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
      MPI_Send (buff,maxncn*(totmaxndofn+1),MPI_LONG,l,ndom,MPI_COMM_WORLD);
    }

    k=0;
    for (i=0;i<top->nn;i++){
      if (ltg[i]<-1){
	buff[k*(totmaxndofn+1)+totmaxndofn]=0-ltg[i]-2;
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
  
  /*
  fprintf (out,"\n\n kontrola bufferu\n");
  long ijk=0;
  for (i=0;i<maxncn;i++){
    fprintf (out,"\n");
    for (j=0;j<totmaxndofn+1;j++){
      fprintf (out,"  %ld",buff[ijk]);  ijk++;
    }
  }
  */

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
  j=nidof;
  for (i=0;i<ncdof;i++){
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
    
    ptop = new gtopology;
    ptop->initiate (masgcn,ncdofdom,nproc);
    
    /*
    fprintf (out,"\n\n\nKONTROLA TOPOLOGIE\n");
    for (i=0;i<ptop->ne;i++){
      fprintf (out,"\n prvek cislo %ld\n",i);
      for (j=0;j<ptop->gelements[i].ndofe;j++){
	fprintf (out,"  %ld",ptop->gelements[i].cn[j]);
      }
    }
    */

    
    
  }
  else{
    MPI_Send (buff,maxncdof,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  delete [] gcn;
  







  //  boundary nodes and unknowns




  //  node-subdomain correspondence assembling on master

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
    fprintf (out,"\n\n kontrola incidenci \n");
    for (i=0;i<maxgnbn;i++){
      fprintf (out,"\n");
      for (j=0;j<nproc;j++){
	fprintf (out,"  %ld",noddom[i][j]);
      }
    }
    */
    
  }
  else{
    MPI_Send (buff,maxnbn,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  

  //  evaluation of number of incidencies

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
    
    /*
    fprintf (out,"\n\n kontrola poctu incidenci\n");
    for (i=0;i<maxgnbn;i++){
      fprintf (out,"   %ld",nodinc[i]);
    }
    */
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxinc,1,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
    
  }
  else{
    MPI_Recv (&maxinc,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);



  //  code numbers gathering on master processor
  //  global code numbers generation
  //  computation of number of Lagrange multipliers

  

  //  code numbers gathering on master processor

  buff = new long [maxnbn*(totmaxndofn+1)];
  for (i=0;i<maxnbn*(totmaxndofn+1);i++){
    buff[i]=0;
  }
  
  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1){
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
    
    
    

    //  computation of number of Lagrange multipliers
    //  global code numbers generation

    tnmdof=1;
    for (i=0;i<maxgnbn;i++){
      for (j=totmaxndofn;j<nodinc[i]*totmaxndofn;j++){
	if (gncn[i][j]>0){
	  gncn[i][j]=tnmdof;  tnmdof++;
	}
      }
    }
    tnmdof--;
    
    /*
    fprintf (out,"\n\n kontrola kodovych cisel multiplikatoru \n");
    for (i=0;i<maxgnbn;i++){
      fprintf (out,"\n");
      for (j=0;j<nodinc[i]*totmaxndofn;j++){
	fprintf (out,"   %ld",gncn[i][j]);
      }
    }
    */

  }
  else{
    MPI_Send (buff,maxnbn*(totmaxndofn+1),MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  

  

  //  gobal code numbers scattering

  
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
    
    /*
    fprintf (out,"\n\n kontrola bufferu \n");
    long ijk=0;
    for (i=0;i<nbndom[ndom];i++){
      fprintf (out,"\n");
      for (j=0;j<maxinc*totmaxndofn+1;j++){
	fprintf (out,"  %ld",buff[ijk]);  ijk++;
      }
    }
    */
    
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
      
      /*
      fprintf (out,"\n\n kontrola bufferu %ld \n",ii);
      long ijk=0;
      for (long ij=0;ij<nbndom[ii];ij++){
	fprintf (out,"\n");
	for (j=0;j<maxinc*totmaxndofn+1;j++){
	  fprintf (out,"  %ld",buff[ijk]);  ijk++;
	}
      }
      */
      
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

    for (i=0;i<maxgnbn;i++){
      delete [] gncn[i];
    }
    delete [] gncn;

    for (i=0;i<maxgnbn;i++){
      delete [] noddom[i];
    }
    delete [] noddom;
  }
  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Send(&tnmdof,1,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(&tnmdof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

}


/**
   function localizes local vector components to global vector
   localization process is generated by FETI ordering
   
   @param gv - array containing global vector
   @param lv - array containing local vector

   4.12.2001
*/
void dpfeti::locglobfeti (gtopology *top,double *gv,double *lv)
{
  long i,j,k,l,ln,lcn,llgcn;

  for (i=0;i<nbn;i++){
    l=0;
    ln=nodnum[i];
    
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
   
   @param top - topology
   @param gv - array containing global vector
   @param lv - array containing local vector

   4.12.2001
*/
void dpfeti::globlocfeti (gtopology *top,double *gv,double *lv)
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


//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************
//****************************************************************

void dpfeti::arrmatrix (double *condmat)
{
  long i,j,k,l,ii,gndofe,*cn;
  double *buff;
  MPI_Status stat;
  
  buff = new double [maxncdof*maxncdof];
  ii=0;
  for (i=0;i<ncdof;i++){
    for (j=0;j<ncdof;j++){
      buff[ii]=condmat[i*ncdof+j];  ii++;
    }
  }
  
  if (myrank==0){
    arr = new gmatrix;

    arr->ts = rsmstor;
    arr->tlinsol = tlinsol;
    
    arr->zero = zero;
    arr->limit = limit;
    
    arr->initiate (ptop,tncdof,1);
    arr->prepmat (0.0,1);
    
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    arr->localized (condmat,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof*maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=stat.MPI_TAG;
      
      ii=0;
      for (j=0;j<ncdofdom[l];j++){
	for (k=0;k<ncdofdom[l];k++){
	  condmat[j*ncdofdom[l]+k]=buff[ii];  ii++;
	}
      }
      
      gndofe = ptop->give_ndofe (l);
      cn = new long [gndofe];
      ptop->give_code_numbers (l,cn);
      arr->localized (condmat,cn,gndofe);
      delete [] cn;
      
    }
    
    arr->decompose_matrix ();
    
  }
  else{
    MPI_Send(buff,maxncdof*maxncdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
}


void dpfeti::vectors_br_bm (gtopology *top,gmatrix *gm,
			    double *condmat,double *condvect,double *rhs)
{
  long i,j,l,gndofe,*cn;
  double *lhs,*buff,*zaloha;
  MPI_Status stat;
  
  //************
  //  vector b_r
  //************
  buff = new double [maxncdof];
  for (i=0;i<ncdof;i++){
    buff[i]=condvect[i];
  }
  
  if (myrank==0){
    br = new double [tncdof];
    nullvr (br,tncdof);
    
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    locglob (br,condvect,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=stat.MPI_TAG;
      
      for (j=0;j<ncdofdom[l];j++){
	condvect[j]=buff[j];
      }
      
      gndofe = ptop->give_ndofe (l);
      cn = new long [gndofe];
      ptop->give_code_numbers (l,cn);
      locglob (br,condvect,cn,gndofe);
      delete [] cn;
      
    }
    
  }
  else{
    MPI_Send(buff,maxncdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  //************
  //  vector b_m
  //************

  buff = new double [tnmdof];
  lhs = new double [nidof];
  zaloha = new double [nidof];
  for (i=0;i<nidof;i++){
    zaloha[i]=rhs[i];
  }
  
  
  //  (K_{rr}^{j})^{-1} . f_{r}^{j}
  gm->condense (condmat,lhs,zaloha,ncdof,3);
  delete [] zaloha;
  
  
  //  G^{j} . (K_{rr}^{j})^{-1} . f_{r}^{j}
  nullvr (buff,tnmdof);
  locglobfeti (top,buff,lhs);
  
  if (myrank==0){
    bm = new double [tnmdof];
    nullvr (bm,tnmdof);
    locglobfeti (top,bm,lhs);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,tnmdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (bm,buff,tnmdof);
    }
  }
  else{
    MPI_Send(buff,tnmdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  delete [] lhs;
  
}


void dpfeti::rhs_dpfeti (gtopology *top,long *domproc,gmatrix *gm)
{
  long i,j,gndofe,*cn;
  double *ur,*wr,*utc,*vtc,*buff;
  MPI_Status stat;
  
  ur = new double [nidof];
  wr = new double [nidof];
  buff = new double [maxncdof];
  
  if (myrank==0){
    utc = new double [tncdof];
    vtc = new double [tncdof];
    
    for (i=0;i<tncdof;i++){
      utc[i]=br[i];
    }
    
    //  A_{RR}^{-1} . b_{R} = \mu
    arr->back_substitution (vtc,utc);

    for (i=0;i<nproc;i++){
      j=domproc[i];
      
      //  F_c^j . \mu
      if (j==myrank)  continue;
      nullvr (buff,maxncdof);
      gndofe = ptop->give_ndofe (i);
      cn = new long [gndofe];
      ptop->give_code_numbers (i,cn);
      globloc (vtc,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
    }
    
    nullvr (buff,maxncdof);
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    globloc (vtc,buff,cn,gndofe);
    delete [] cn;
    
    delete [] utc;
    delete [] vtc;
  }
  else{
    MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  //  K_{rc}^j . [F_c^j . \mu]
  mv (krc,buff,ur,nidof,ncdof);
  delete [] buff;
  
  //  (K_{rr}^j)^{-1} . [K_{rc}^j . F_c^j . \mu]
  gm->condense (buff,wr,ur,ncdof,3);
  

  buff = new double [tnmdof];

  //  (G^j) . [(K_{rr}^j)^{-1} . K_{rc}^j . F_c^j . \mu]
  nullvr (buff,tnmdof);
  locglobfeti (top,buff,wr);
  
  if (myrank==0){
    nullvr (cgrhs,tnmdof);
    locglobfeti (top,cgrhs,wr);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,tnmdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (cgrhs,buff,tnmdof);
    }
  }
  else{
    MPI_Send(buff,tnmdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  delete [] ur;
  delete [] wr;
  
  if (myrank==0){
    for (i=0;i<tnmdof;i++){
      cgrhs[i]=bm[i]-cgrhs[i];
    }
  }
  
}

void dpfeti::matxvect (gtopology *top,long *domproc,
		       gmatrix *gm,double *input,double *output)
{
  long i,j,gndofe,*cn;
  double *utc,*vtc,*ur,*vr,*wr,*buff;
  MPI_Status stat;
  
  ur = new double [nidof];
  vr = new double [nidof];
  wr = new double [nidof];
  
  buff = new double [tnmdof];

  if (myrank==0){
    for (i=0;i<tnmdof;i++){
      buff[i]=input[i];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send(buff,tnmdof,MPI_DOUBLE,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(buff,tnmdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullvr (ur,nidof);
  globlocfeti (top,buff,ur);
  delete [] buff;
  
  //  (K_{rr}^j)^{-1} . [(G^j)^T . \lambda]
  gm->condense (buff,vr,ur,ncdof,3);
  
  buff = new double [maxncdof];

  //  (K_{rc}^j)^T . [(K_{rr}^j)^{-1} . (G^j)^T . \lambda]
  mtv (krc,vr,buff,nidof,ncdof);
  
  if (myrank==0){
    utc = new double [tncdof];
    vtc = new double [tncdof];
    nullvr (utc,tncdof);
    
    //  (F_c^j)^T . [(K_{rc}^j)^T . (K_{rr}^j)^{-1} . (G^j)^T . \lambda] = \nu
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    locglob (utc,buff,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      j=stat.MPI_TAG;
      gndofe = ptop->give_ndofe (j);
      cn = new long [gndofe];
      ptop->give_code_numbers (j,cn);
      locglob (utc,buff,cn,gndofe);
      delete [] cn;
    }

    //  A_{RR}^{-1} . A_{RM} . \lambda = \mu
    arr->back_substitution (vtc,utc);
    
    for (i=0;i<nproc;i++){
      j=domproc[i];
      
      //  F_c^j . \mu
      if (j==myrank)  continue;
      nullvr (buff,maxncdof);
      gndofe = ptop->give_ndofe (i);
      cn = new long [gndofe];
      ptop->give_code_numbers (i,cn);
      globloc (vtc,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
    }

    nullvr (buff,maxncdof);
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    globloc (vtc,buff,cn,gndofe);
    delete [] cn;
    
    delete [] utc;
    delete [] vtc;
  }
  else{
    MPI_Send(buff,maxncdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
    MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  K_{rc}^j . [F_c^j . \mu]
  mv (krc,buff,ur,nidof,ncdof);
  delete [] buff;
  
  //  (K_{rr}^j)^{-1} . [K_{rc}^j . F_c^j . \mu]
  gm->condense (buff,wr,ur,ncdof,3);
  
  //  puvodni verze
  //for (i=0;i<nidof;i++){
  // wr[i]=vr[i]-wr[i];
  //}
  for (i=0;i<nidof;i++){
    wr[i]=vr[i]+wr[i];
  }

  buff = new double [tnmdof];

  //  (G^j) . {[(K_{rr}^j)^{-1} . K_{rc}^j . F_c^j . \mu] + [(K_{rr}^{-1})^{-1} . (G^j)^T . \lambda]}
  nullvr (buff,tnmdof);
  locglobfeti (top,buff,wr);
  
  if (myrank==0){
    nullvr (output,tnmdof);
    locglobfeti (top,output,wr);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,tnmdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      addv (output,buff,tnmdof);
    }
  }
  else{
    MPI_Send(buff,tnmdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  delete [] buff;
  delete [] ur;
  delete [] vr;
  delete [] wr;
}






void dpfeti::cg (gtopology *top,long *domproc,gmatrix *gm,long iv,FILE *out)
{
  long i,j,stop;
  double nom,denom,norrhs,alpha,beta;
  double *d,*r,*p;
  MPI_Status stat;
  
  if (myrank==0){
    d = new double [tnmdof];
    r = new double [tnmdof];
    p = new double [tnmdof];
    
    //  initial values
    if (iv==0){
      for (i=0;i<tnmdof;i++)
	cglhs[i]=0.0;
    }
  }
  matxvect (top,domproc,gm,cglhs,p);
  
  if (myrank==0){
    norrhs=0.0;  nom=0.0;
    for (i=0;i<tnmdof;i++){
      norrhs+=cgrhs[i]*cgrhs[i];
      r[i]=p[i]-cgrhs[i];
      nom+=r[i]*r[i];
      d[i]=-1.0*r[i];
    }
    
    if (norrhs<zero){
      fprintf (stderr,"\n\n norm of right hand side in conjugate gradient method is smaller than %e",zero);
      fprintf (stderr,"\n see file %s, line %d.\n",__FILE__,__LINE__);
      aerrcg=norrhs;  anicg=0;
      return;
    }
  }
  
  
  //  iteration loop
  stop=0;
  for (i=0;i<nicg;i++){
    
    
    //  new coefficient alpha
    matxvect (top,domproc,gm,d,p);
    
    if (myrank==0){
      denom = ss (d,p,tnmdof);
      if (fabs(denom)<zero){
	fprintf (stdout,"\n there is zero denominator in alpha computation in conjugate gradient method (%s, line %d)\n",__FILE__,__LINE__);
	stop=1;
      }
      
      alpha = nom/denom;
      
      //  new approximation of x and r
      for (j=0;j<tnmdof;j++){
	cglhs[j]+=alpha*d[j];
	r[j]+=alpha*p[j];
      }
      
      denom=nom;
      
      nom = ss (r,r,tnmdof);
      
      fprintf (out,"\n iteration   %ld  norres/norrhs %e",i,sqrt(nom/norrhs));
      fprintf (stdout,"\n iteration   %ld  norres/norrhs %e",i,sqrt(nom/norrhs));
      //printf ("\n iteration   %ld  norres/norrhs %e",i,sqrt(nom/norrhs));
      
      if (sqrt(nom/norrhs)<errcg)  stop=1;
      //if (fabs(nom)<limit)  break;
      
      
      beta = nom/denom;
      
      //  new vector of direction
      for (j=0;j<tnmdof;j++){
	d[j]=beta*d[j]-r[j];
      }

      for (j=1;j<nproc;j++){
	MPI_Send(&stop,1,MPI_LONG,j,ndom,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv(&stop,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (stop==1)  break;
  }
  
  anicg=i;  aerrcg=nom;
  
  if (myrank==0){
    delete [] p;  delete [] r;  delete [] d;
  }
  
}


void dpfeti::corner_displ (gtopology *top,long *domproc,gmatrix *gm)
{
  long i,j,gndofe,*cn;
  double *utc,*ur,*vr,*buff;
  MPI_Status stat;
  
  ur = new double [nidof];
  vr = new double [nidof];
  
  buff = new double [tnmdof];
  
  if (myrank==0){
    for (i=0;i<tnmdof;i++){
      buff[i]=cglhs[i];
    }
    for (i=1;i<nproc;i++){
      MPI_Send(buff,tnmdof,MPI_DOUBLE,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(buff,tnmdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullvr (ur,nidof);
  globlocfeti (top,buff,ur);
  
  //  (K_{rr}^j)^{-1} . [(G^j)^T . \lambda]
  gm->condense (buff,vr,ur,ncdof,3);
  
  delete [] buff;
  buff = new double [maxncdof];

  //  (K_{rc}^j)^T . [(K_{rr}^j)^{-1} . (G^j)^T . \lambda]
  mtv (krc,vr,buff,nidof,ncdof);
  
  if (myrank==0){
    utc = new double [tncdof];
    nullvr (utc,tncdof);
    
    //  (F_c^j)^T . [(K_{rc}^j)^T . (K_{rr}^j)^{-1} . (G^j)^T . \lambda] = \nu
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    locglob (utc,buff,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      j=stat.MPI_TAG;
      gndofe = ptop->give_ndofe (j);
      cn = new long [gndofe];
      ptop->give_code_numbers (j,cn);
      locglob (utc,buff,cn,gndofe);
      delete [] cn;
    }
    
    //  puvodni verze
    //for (i=0;i<tncdof;i++){
    // utc[i]=br[i]-utc[i];
    //}
    for (i=0;i<tncdof;i++){
      utc[i]=br[i]+utc[i];
    }
    
    //  A_{RR}^{-1} . (b_R - A_{RM} . \lambda) = r_C
    arr->back_substitution (displ,utc);
    
    delete [] utc;
  }
  else{
    MPI_Send(buff,maxncdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  delete [] ur;
  delete [] vr;
}

void dpfeti::compute_displ (gtopology *top,gmatrix *gm,long *domproc,
			    double *subdispl,double *rhs)
{
  long i,j,gndofe,*cn;
  double *ur,*vr,*buff;
  MPI_Status stat;
  
  ur = new double [nidof];
  vr = new double [nidof];
  
  buff = new double [maxncdof];
  if (myrank==0){
    for (i=0;i<nproc;i++){
      j=domproc[i];
      
      //  F_c^j . \mu
      if (j==myrank)  continue;
      nullvr (buff,maxncdof);
      gndofe = ptop->give_ndofe (i);
      cn = new long [gndofe];
      ptop->give_code_numbers (i,cn);
      globloc (displ,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,j,ndom,MPI_COMM_WORLD);
    }

    nullvr (buff,maxncdof);
    gndofe = ptop->give_ndofe (ndom);
    cn = new long [gndofe];
    ptop->give_code_numbers (ndom,cn);
    globloc (displ,buff,cn,gndofe);
    delete [] cn;
    
  }
  else{
    MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  j=nidof;
  for (i=0;i<ncdof;i++){
    subdispl[j]=buff[i];  j++;
  }
  
  //  K_{rc}^j . [F_c^j . r_C]
  mv (krc,buff,ur,nidof,ncdof);
  delete [] buff;
  
  
  buff = new double [tnmdof];
  
  if (myrank==0){
    for (i=0;i<tnmdof;i++){
      buff[i]=cglhs[i];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send(buff,tnmdof,MPI_DOUBLE,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(buff,tnmdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullvr (vr,nidof);
  globlocfeti (top,buff,vr);
  delete [] buff;
  
  for (i=0;i<nidof;i++){
    rhs[i]=rhs[i]-ur[i]-vr[i];
  }
  
  gm->condense (buff,subdispl,rhs,ncdof,3);
  
  delete [] ur;
  delete [] vr;
}


/**
   function solves system of algebraic equations by DP-FETI method
   
   24.7.2002
*/
void dpfeti::solve_system (gtopology *top,gmatrix *gm,long *domproc,
			   double *lhs,double *rhs,FILE *out,long mespr)
{
  long i,j;
  double *condmat,*condvect,*rrhs,*crhs;
  time_t t1,t2,t3,t4,t5,t6,t7,t8;
  
  t1 = time (NULL);
  
  if (myrank==0){
    condmat = new double [maxncdof*maxncdof];
    memset (condmat,0,maxncdof*maxncdof*sizeof(double));
    condvect = new double [maxncdof];
    memset (condvect,0,maxncdof*sizeof(double));
  }
  else{
    condmat = new double [ncdof*ncdof];
    memset (condmat,0,ncdof*ncdof*sizeof(double));
    condvect = new double [ncdof];
    memset (condvect,0,ncdof*sizeof(double));
  }
  
  if (myrank==0){

    fprintf (out,"\n\n tncdof %ld",tncdof);

    cglhs = new double [tnmdof];
    cgrhs = new double [tnmdof];
    displ = new double [tncdof];
  }

  //  allocation of matrix K_rc
  krc = new double [nidof*ncdof];
  
  //  K_rc matrix assembling
  gm->a12block (krc,ncdof);
  
  crhs = new double [ncdof];
  rrhs = new double [nidof];
  
  //  copy of right hand side vector
  for (i=0;i<nidof;i++){
    rrhs[i]=rhs[i];
  }
  j=0;
  for (i=nidof;i<ndof;i++){
    crhs[j]=rhs[i];  j++;
  }
  
  t2 = time (NULL);
 
  //  matrix condensation
  gm->condense (condmat,lhs,rhs,ncdof,1);
  
  j=0;
  for (i=nidof;i<ndof;i++){
    condvect[j]=rhs[i];
    j++;
  }
  
  t3 = time (NULL);
  
  arrmatrix (condmat);
  
  t4 = time (NULL);

  //  computation of vectors b_R and b_M
  vectors_br_bm (top,gm,condmat,condvect,rrhs);
  
  //  computation of right hand side of final system of equations
  rhs_dpfeti (top,domproc,gm);
  
  t5 = time (NULL);
    
  //  solution of final problem by conjugate gradient method
  cg (top,domproc,gm,0,out);
  
  t6 = time (NULL);
  
  if (myrank==0){
    fprintf (out,"\n\n L A G R A N G E    M U L T I P L I E R S \n");
    for (i=0;i<tnmdof;i++){
      fprintf (out,"\n %7ld  %30.15f",i,cglhs[i]);
    }
  }
  
  //  corner displacement computation
  corner_displ (top,domproc,gm);
  
  t7 = time (NULL);

  if (myrank==0){
    fprintf (out,"\n\n C O R N E R    D I S P L A C E M E N T S \n");
    fprintf (out,"\n\n tncdof  %ld",tncdof);
    for (i=0;i<tncdof;i++){
      fprintf (out,"\n %7ld   %30.15f",i,displ[i]);
    }
    fprintf (out,"\n\n konec tisku\n");
  }
  
  //  subdomain displacement computation
  compute_displ (top,gm,domproc,lhs,rrhs);

  t8 = time (NULL);
  
  fprintf (out,"\n\n A_rc block assembling            %ld",t2-t1);
  fprintf (out,"\n\n time of condensation             %ld",t3-t2);
  fprintf (out,"\n\n A_rr assembling                  %ld",t4-t3);
  fprintf (out,"\n\n b_R and b_M assembling           %ld",t5-t4);
  fprintf (out,"\n\n solution of reduced system       %ld",t6-t5);
  fprintf (out,"\n\n corner displacements calcul.     %ld",t7-t6);
  fprintf (out,"\n\n all displacements calcul.        %ld",t8-t7);

  
  delete [] condmat;
  delete [] condvect;
  delete [] rrhs;
  delete [] crhs;
  
}
