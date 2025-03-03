#include "dpfeti_new.h"
#include <string.h>
#include <math.h>
#include <time.h>

dpfeti::dpfeti(int np,int mr,int nd)
{
  nproc=np;  myrank=mr;  ndom=nd;

}

dpfeti::~dpfeti()
{
  delete [] cgrhs;
  delete [] cglhs;
  delete [] displ;
  delete [] br;
  delete [] bm;
  
  delete [] krc;
  
  for (i=0;i<nproc;i++){
    delete [] loccn[i];
    delete [] globcn[i];
  }
  delete [] loccn;
  delete [] globcn;
  
  delete [] nlbdof;
  delete [] ncdofdom;
  
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
   
   JK, 26.11.2002, 10.2.2003
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
  long i,j,k,l,m,n,ndofn,stop,nbn,ncn,maxnbn,maxncn,tnbn,tncn;
  long *nibn,*lgnbn,*llnbn,*lgncn,*llncn,***lcnbn,***lcncn,**gnbndom,**gncndom;
  long **cncn,**cnbn,*nbndom,*ncndom,**dominc,**supelemcn;
  long buffsize,*buff;
  MPI_Status stat;
  
  //  ndofn - number of degrees of freedom of one node
  ndofn = top->give_ndofn (0);
  
  //  nbn - number of boundary nodes (except corner nodes)
  //  ncn - number of corner nodes
  nbn=0;  ncn=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1)  nbn++;
    if (ltg[i]<-1)  ncn++;
  }
  
  //  lgnbn - list of global numbers of boundary nodes
  //  lgncn - list of global numbers of corner nodes
  lgnbn = new long [nbn];
  lgncn = new long [ncn];
  //  llnbn - list of local numbers of boundary nodes
  //  llncn - list of local numbers of corner nodes
  llnbn = new long [nbn];
  llncn = new long [ncn];

  j=0;  k=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1){
      lgnbn[j]=ltg[i];
      llnbn[j]=i;
      j++;
    }
    if (ltg[i]<-1){
      lgncn[k]=0-ltg[i]-2;
      llncn[k]=i;
      k++;
    }
  }
  
  buff = new long [2];
  buff[0]=nbn;
  buff[1]=ncn;
  
  if (myrank==0){
    //  array containing numbers of boundary nodes on subdomains
    nbndom = new long [nproc];
    //  array containing numbers of corner nodes on subdomains
    ncndom = new long [nproc];
    
    nbndom[ndom]=nbn;
    ncndom[ndom]=ncn;
    
    //  maxnbn - maximum number of boundary nodes on one subdomain
    //  maxncn - maximum number of corner nodes on one subdomain
    maxnbn=nbn;  maxncn=ncn;
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      nbndom[stat.MPI_TAG]=buff[0];
      ncndom[stat.MPI_TAG]=buff[1];
      
      if (maxnbn<buff[0])  maxnbn=buff[0];
      if (maxncn<buff[1])  maxncn=buff[1];
    }
    
    buff[0]=maxnbn;
    buff[1]=maxncn;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,2,MPI_LONG,i,ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,2,MPI_LONG,0,ndom,MPI_COMM_WORLD);
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    maxnbn=buff[0];
    maxncn=buff[1];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  buffsize = (maxnbn+maxncn)*(ndofn+1);
  buff = new long [buffsize];
  k=0;
  for (i=0;i<nbn;i++){
    buff[k]=lgnbn[i];  k++;
    for (j=0;j<ndofn;j++){
      buff[k]=top->give_dof (llnbn[i],j);
      k++;
    }
  }
  for (i=0;i<ncn;i++){
    buff[k]=lgncn[i];  k++;
    for (j=0;j<ndofn;j++){
      buff[k]=top->give_dof (llncn[i],j);
      k++;
    }
  }
  
  delete [] llnbn;
  delete [] llncn;
  delete [] lgnbn;
  delete [] lgncn;
  
  if (myrank==0){
    
    //  list of local code numbers of boundary nodes from each subdomain
    lcnbn = new long** [nproc];
    //  list of local code numbers of corner nodes from each subdomain
    lcncn = new long** [nproc];
    //  list of global numbers of boundary nodes from each subdomain
    gnbndom = new long* [nproc];
    //  list of global numbers of corner nodes from each subdomain
    gncndom = new long* [nproc];
    for (i=0;i<nproc;i++){
      gnbndom[i] = new long [nbndom[i]];
      gncndom[i] = new long [ncndom[i]];
      
      lcnbn[i] = new long* [nbndom[i]];
      lcncn[i] = new long* [ncndom[i]];
      for (j=0;j<nbndom[i];j++){
	lcnbn[i][j] = new long [ndofn];
      }
      for (j=0;j<ncndom[i];j++){
	lcncn[i][j] = new long [ndofn];
      }
    }
    
    //  master contribution
    l=0;
    for (j=0;j<nbndom[ndom];j++){
      gnbndom[ndom][j]=buff[l];  l++;
      for (k=0;k<ndofn;k++){
	lcnbn[ndom][j][k]=buff[l];
	l++;
      }
    }
    for (j=0;j<ncndom[ndom];j++){
      gncndom[ndom][j]=buff[l];  l++;
      for (k=0;k<ndofn;k++){
	lcncn[ndom][j][k]=buff[l];
	l++;
      }
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      l=0;
      for (j=0;j<nbndom[stat.MPI_TAG];j++){
	gnbndom[stat.MPI_TAG][j]=buff[l];  l++;
	for (k=0;k<ndofn;k++){
	  lcnbn[stat.MPI_TAG][j][k]=buff[l];
	  l++;
	}
      }
      for (j=0;j<ncndom[stat.MPI_TAG];j++){
	gncndom[stat.MPI_TAG][j]=buff[l];  l++;
	for (k=0;k<ndofn;k++){
	  lcncn[stat.MPI_TAG][j][k]=buff[l];
	  l++;
	}
      }
      
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,ndom,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  
  
  if (myrank==0){
    

    //  tnbn - total number of all boundary nodes
    tnbn = 0;
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	if (gnbndom[i][j]>tnbn)  tnbn=gnbndom[i][j];
      }
    }
    tnbn++;
    //  tncn - total number of all corner nodes
    tncn = 0;
    for (i=0;i<nproc;i++){
      for (j=0;j<ncndom[i];j++){
	if (gncndom[i][j]>tncn)  tncn=gncndom[i][j];
      }
    }
    tncn++;
    
    //  code numbers of corner nodes
    cncn = new long* [tncn];
    for (i=0;i<tncn;i++){
      cncn[i] = new long [ndofn];
      for (j=0;j<ndofn;j++){
	cncn[i][j]=0;
      }
    }
    
    //  auxiliary information about supports of corner nodes
    for (i=0;i<nproc;i++){
      for (j=0;j<ncndom[i];j++){
	for (k=0;k<ndofn;k++){
	  cncn[gncndom[i][j]][k]=lcncn[i][j][k];
	}
      }
    }
    
    //  generation of code numbers of corner nodes
    tncdof=1;
    for (i=0;i<tncn;i++){
      for (j=0;j<ndofn;j++){
	if (cncn[i][j]>0){
	  cncn[i][j]=tncdof;
	  tncdof++;
	}
      }
    }
    tncdof--;
    
    //  array containing numbers of degrees of freedom of corner nodes on subdomains
    ncdofdom = new long [nproc];
    for (i=0;i<nproc;i++){
      ncdofdom[i]=0;
    }
    
    //  computation of numbers of superelement code numbers
    for (i=0;i<nproc;i++){
      for (j=0;j<ncndom[i];j++){
	for (k=0;k<ndofn;k++){
	  if (cncn[gncndom[i][j]][k]>0)  ncdofdom[i]++;
	}
      }
    }
    
    maxncdof=0;
    for (i=0;i<nproc;i++){
      if (maxncdof<ncdofdom[i])  maxncdof=ncdofdom[i];
    }
    
    
    //  array containing superelement code numbers
    supelemcn = new long* [nproc];
    for (i=0;i<nproc;i++){
      supelemcn[i] = new long [ncdofdom[i]];
      for (j=0;j<ncdofdom[i];j++){
	supelemcn[i][j]=0;
      }
    }
    
    //  generation of superelement code numbers
    for (i=0;i<nproc;i++){
      l=0;
      for (j=0;j<ncndom[i];j++){
	for (k=0;k<ndofn;k++){
	  m=cncn[gncndom[i][j]][k];
	  if (m>0){
	    supelemcn[i][l]=m;
	    l++;
	  }
	}
      }
    }
    
    //  supertopology
    ptop = new gtopology;
    ptop->initiate (supelemcn,ncdofdom,nproc);
    
    for (i=0;i<nproc;i++){
      delete [] gncndom[i];
    }
    delete [] gncndom;

    for (i=0;i<nproc;i++){
      delete [] supelemcn[i];
    }
    delete [] supelemcn;
    
    for (i=0;i<tncn;i++){
      delete [] cncn[i];
    }
    delete [] cncn;
    
    
    
    
    
    //  array containing number of incidences of boundary nodes
    nibn = new long [tnbn];
    for (i=0;i<tnbn;i++){
      nibn[i]=0;
    }
    
    //  computation of numbers of incidences
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	nibn[gnbndom[i][j]]++;
      }
    }
    
    //  array containing incidences of nodes to subdomains
    dominc = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      dominc[i] = new long [nibn[i]];
      for (j=0;j<nibn[i];j++){
	dominc[i][j]=-1;
      }
    }
    
    for (i=0;i<tnbn;i++){
      nibn[i]=0;
    }
    
    //  searching for incidences
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	dominc[gnbndom[i][j]][nibn[gnbndom[i][j]]]=i;
	nibn[gnbndom[i][j]]++;
      }
    }
    
    
    //  array containing code numbers of boundary nodes (FETI system)
    cnbn = new long** [tnbn];
    for (i=0;i<tnbn;i++){
      cnbn[i] = new long* [nibn[i]-1];
      for (j=0;j<nibn[i]-1;j++){
	cnbn[i][j] = new long [ndofn];
	for (k=0;k<ndofn;k++){
	  cnbn[i][j][k]=0;
	}
      }
    }
    
    //  auxiliary information about supports
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	m=gnbndom[i][j];
	for (k=0;k<nibn[m]-1;k++){
	  for (l=0;l<ndofn;l++){
	    cnbn[m][k][l]=lcnbn[i][j][l];
	  }
	}
      }
    }
    
    //  generation of code numbers of boundary nodes
    tnmdof=1;
    for (i=0;i<tnbn;i++){
      for (j=0;j<nibn[i]-1;j++){
	for (k=0;k<ndofn;k++){
	  if (cnbn[i][j][k]>0){
	    cnbn[i][j][k]=tnmdof;
	    tnmdof++;
	  }
	}
      }
    }
    tnmdof--;
    
    
    nlbdofdom = new long [nproc];
    for (i=0;i<nproc;i++){
      nlbdofdom[i]=0;
    }
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	m=gnbndom[i][j];
	if (dominc[m][0]==i){
	  for (k=0;k<nibn[m]-1;k++){
	    for (l=0;l<ndofn;l++){
	      if (cnbn[m][k][l]>0)  nlbdofdom[i]++;
	    }
	  }
	}
	else{
	  for (l=0;l<ndofn;l++){
	    if (cnbn[m][0][l]>0)  nlbdofdom[i]++;
	  }
	}
      }
    }
    
    
    //  local code numbers
    loccn = new long* [nproc];
    globcn = new long* [nproc];
    maxnbdof=0;
    for (i=0;i<nproc;i++){
      loccn[i] = new long [nlbdofdom[i]];
      globcn[i] = new long [nlbdofdom[i]];
      if (maxnbdof<nlbdofdom[i])  maxnbdof=nlbdofdom[i];
    }
    
    for (i=0;i<nproc;i++){
      n=0;
      for (j=0;j<nbndom[i];j++){
	m=gnbndom[i][j];
	if (dominc[m][0]==i){
	  for (k=0;k<nibn[m]-1;k++){
	    for (l=0;l<ndofn;l++){
	      if (cnbn[m][k][l]>0){
		globcn[i][n]=0-cnbn[m][k][l];
		loccn[i][n]=lcnbn[i][j][l];
		n++;
	      }
	    }
	  }
	}
	else{
	  stop=0;
	  for (k=0;k<nibn[m]-1;k++){
	    for (l=0;l<ndofn;l++){
	      if (cnbn[m][k][l]==0){
		cnbn[m][k][k]=-3;
		stop=1;
	      }
	      if (cnbn[m][k][l]>0){
		globcn[i][n]=cnbn[m][k][l];
		loccn[i][n]=lcnbn[i][j][l];
		n++;
		cnbn[m][k][k]=-3;
		stop=1;
	      }
	    }
	    if (stop==1)  break;
	  }
	}
      }
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	delete [] lcnbn[i][j];
      }
      delete [] lcnbn[i];
    }
    delete [] lcnbn;
    
    for (i=0;i<nproc;i++){
      for (j=0;j<ncndom[i];j++){
	delete [] lcncn[i][j];
      }
      delete [] lcncn[i];
    }
    delete [] lcncn;

    for (i=0;i<nproc;i++){
      delete [] nbndom[i];
    }
    delete [] nbndom;

    for (i=0;i<nproc;i++){
      delete [] ncndom[i];
    }
    delete [] ncndom;

    for (i=0;i<nproc;i++){
      delete [] gnbndom[i];
    }
    delete [] gnbndom;
    
    for (i=0;i<tnbn;i++){
      delete [] dominc[i];
    }
    delete [] dominc;
    
    for (i=0;i<tnbn;i++){
      for (j=0;j<nibn[i]-1;j++){
	delete [] cnbn[i][j];
      }
      delete [] cnbn[i];
    }
    delete [] cnbn;
    
    delete [] nibn;
    
  }
  
  buff = new long [3];
  if (myrank==0){
    buff[0]=maxncdof;
    buff[1]=maxnbdof;
    
    for (i=0;i<nproc;i++){
      if (domproc[i]==0)  continue;
      buff[2]=nlbdofdom[i];
      MPI_Send (buff,3,MPI_LONG,domproc[i],ndom,MPI_COMM_WORLD);
    }
    nlbdof=nlbdofdom[domproc[nproc]];
  }
  else{
    MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    maxncdof=buff[0];
    maxnbdof=buff[1];
    nlbdof=buff[2];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  localcn = new long [nlbdof];
  globalcn = new long [nlbdof];
  
  buffsize=maxnbdof;
  buff = new long [buffsize];
  if (myrank==0){
    for (i=0;i<nproc;i++){
      if (domproc[i]==0)  continue;
      for (j=0;j<nlbdofdom[i];j++){
	buff[j]=loccn[i][j];
      }
      MPI_Send (buff,buffsize,MPI_LONG,domproc[i],ndom,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    for (j=0;j<nlbdof;j++){
      localcn[j]=buff[j];
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
}


/**
   function extracts boundary components from local %vector
   localization process is generated by FETI ordering
   
   @param ev - array containing extracted %vector
   @param lv - array containing local %vector
   
   4.12.2001, 11.2.2003
*/
void dpfeti::extract_from_local_vector (double *ev,double *lv)
{
  long i,lcn;
  
  for (i=0;i<nlbdof;i++){
    lcn=localcn[i];
    if (lcn==0)  continue;
    ev[i]=lv[lcn-1];
  }
  
}

/**
   function extracts boundary components from local %vector
   localization process is generated by FETI ordering
   
   @param ev - array containing extracted %vector
   @param lv - array containing local %vector
   
   4.12.2001, 11.2.2003
*/
void dpfeti::put_into_local_vector (double *ev,double *lv)
{
  long i,lcn;
  
  for (i=0;i<nlbdof;i++){
    lcn=localcn[i];
    if (lcn==0)  continue;
    lv[lcn-1]+=ev[i];
  }
  
}


/**
   function extracts boundary components from local %vector
   localization process is generated by FETI ordering
   
   @param ev - array containing extracted %vector
   @param gv - array containing global %vector
   @param nsub - number of subdomain
   
   4.12.2001, 11.2.2003
*/
void dpfeti::extract_from_global_vector (double *ev,double *gv,long nsub)
{
  long i,gcn;
  
  for (i=0;i<nlbdofdom[nsub];i++){
    gcn=globcn[nsub][i];
    if (gcn>0)
      ev[i]=gv[gcn-1];
    if (gcn<0)
      ev[i]=0-gv[0-gcn-1];
  }
  
}

/**
   function extracts boundary components from local %vector
   localization process is generated by FETI ordering
   
   @param ev - array containing extracted %vector
   @param lv - array containing global %vector
   @param nsub - number of subdomain
   
   4.12.2001, 11.2.2003
*/
void dpfeti::put_into_global_vector (double *ev,double *gv,long nsub)
{
  long i,gcn;
  
  for (i=0;i<nlbdofdom[nsub];i++){
    gcn=globcn[nsub][i];
    if (gcn>0)
      gv[gcn-1]+=ev[i];
    if (gcn<0)
      gv[gcn-1]-=ev[i];
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

  buff = new double [maxnbdof];
  lhs = new double [nidof];
  zaloha = new double [nidof];
  for (i=0;i<nidof;i++){
    zaloha[i]=rhs[i];
  }
  
  
  //  (K_{rr}^{j})^{-1} . f_{r}^{j}
  gm->condense (condmat,lhs,zaloha,ncdof,3);
  delete [] zaloha;
  
  
  //  G^{j} . (K_{rr}^{j})^{-1} . f_{r}^{j}
  nullvr (buff,maxnbdof);
  extract_from_local_vector (buff,lhs);
  
  if (myrank==0){
    bm = new double [tnmdof];
    nullvr (bm,tnmdof);
    put_into_global_vector (buff,bm,ndom);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      put_into_global_vector (buff,bm,stat.MPI_TAG);
    }
  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
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
      //  F_c^j . \mu
      if (domproc[i]==0)  continue;
      nullvr (buff,maxncdof);
      gndofe = ptop->give_ndofe (i);
      cn = new long [gndofe];
      ptop->give_code_numbers (i,cn);
      globloc (vtc,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,domproc[i],ndom,MPI_COMM_WORLD);
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
  

  buff = new double [maxnbdof];

  //  (G^j) . [(K_{rr}^j)^{-1} . K_{rc}^j . F_c^j . \mu]
  nullvr (buff,maxnbdof);
  extract_from_local_vector (buff,wr);
  
  if (myrank==0){
    nullvr (cgrhs,tnmdof);
    put_into_global_vector (buff,cgrhs,stat.MPI_TAG);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      put_into_global_vector (buff,cgrhs,stat.MPI_TAG);
    }
  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
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
  
  buff = new double [maxnbdof];

  if (myrank==0){
    for (i=0;i<nproc;i++){
      if (domproc[i]==0)  continue;
      extract_from_global_vector (buff,input,i);
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,domproc[i],ndom,MPI_COMM_WORLD);
    }
    extract_from_global_vector (buff,input,domproc[nproc]);
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullvr (ur,nidof);
  put_into_local_vector (buff,ur);
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
      //  F_c^j . \mu
      if (domproc[i]==0)  continue;
      nullvr (buff,maxncdof);
      gndofe = ptop->give_ndofe (i);
      cn = new long [gndofe];
      ptop->give_code_numbers (i,cn);
      globloc (vtc,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,domproc[i],ndom,MPI_COMM_WORLD);
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

  buff = new double [maxnbdof];

  //  (G^j) . {[(K_{rr}^j)^{-1} . K_{rc}^j . F_c^j . \mu] + [(K_{rr}^{-1})^{-1} . (G^j)^T . \lambda]}
  extract_from_local_vector (buff,wr);
  
  if (myrank==0){
    nullvr (output,tnmdof);
    put_into_global_vector (buff,output,ndom);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      put_into_global_vector (buff,output,stat.MPI_TAG);
    }
  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,ndom,MPI_COMM_WORLD);
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
  
  buff = new double [maxnbdof];
  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      if (domproc[i]==0)  continue;
      extract_from_global_vector (buff,cglhs,i);
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,ndom,MPI_COMM_WORLD);
    }
    extract_from_global_vector (buff,cglhs,domproc[nproc]);
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  (G^j)^T . \lambda
  nullvr (ur,nidof);
  put_into_local_vector (buff,ur);
  
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
      //  F_c^j . \mu
      if (domproc[i]==0)  continue;
      nullvr (buff,maxncdof);
      gndofe = ptop->give_ndofe (i);
      cn = new long [gndofe];
      ptop->give_code_numbers (i,cn);
      globloc (displ,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,domproc[i],ndom,MPI_COMM_WORLD);
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
  
  
  buff = new double [maxnbdof];
  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      if (domproc[i]==0)  continue;
      extract_from_global_vector (buff,cglhs,i);
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,ndom,MPI_COMM_WORLD);
    }
    extract_from_global_vector (buff,cglhs,domproc[nproc]);
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullvr (vr,nidof);
  put_into_local_vector (buff,vr);
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

