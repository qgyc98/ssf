#include "dpfeti.h"
#include <string.h>
#include <math.h>
#include <time.h>

dpfeti::dpfeti(int np,int mr,long nd)
{
  nproc=np;  myrank=mr;  ndom=nd;

}

dpfeti::~dpfeti()
{
  long i;
  
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
  
  delete [] ncdofdom;
  delete [] nlbdofdom;
  
  delete [] localcn;
  delete [] globalcn;
  
}


/**
   function detects corner nodes in general problem
   
   in input file
   ltg[i]=0 - inner node
   ltg[i]>0 - boundary node
   in the code
   ltg[i]=-1 - inner node
   ltg[i]>-1 - boundary node
   
   22.9.2003, JK
   
*/
void dpfeti::cornodedetection (gtopology *top,long *ltg,long *domproc,FILE *out)
{
  long i,j,k,nbn,maxnbn,totnbn,totncn,totnrn,stop;
  long *buff,**gnn,*cnid,*nbndom;
  MPI_Status stat;
  
  nbn=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1)  nbn++;
  }
  
  if (myrank==0){
    nbndom = new long [nproc];
    j=domproc[0];
    nbndom[j]=nbn;
    
    //  maxnbn - maximum number of boundary nodes on one subdomain
    maxnbn=nbn;
    for (i=1;i<nproc;i++){
      MPI_Recv (&nbn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      nbndom[j]=nbn;
      
      if (maxnbn<nbn)  maxnbn=nbn;
    }
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnbn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    /*
    fprintf (out,"\n\n kontrola poctu hranicnich bodu  ");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %4ld  %4ld",i,nbndom[i]);
    }
    fprintf (out,"\n konec kontroly \n");
    fprintf (out,"\n\n maximalni pocet hranicnich bodu  %4ld",maxnbn);
    */
  }
  else{
    MPI_Send (&nbn,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxnbn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);


  buff = new long [maxnbn];
  j=0;
  for (i=0;i<top->nn;i++){
    if (ltg[i]>-1){
      buff[j]=ltg[i];
      j++;
    }
  }
  //fprintf (out,"\n index j po plneni bufferu  %ld",j);

  
  if (myrank==0){
    //  array containing global numbers of boundary nodes
    gnn = new long* [nproc];
    for (i=0;i<nproc;i++){
      gnn[i] = new long [nbndom[i]];
    }
    


    j=domproc[0];
    for (k=0;k<nbndom[j];k++){
      gnn[j][k]=buff[k];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,maxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];

      for (k=0;k<nbndom[j];k++){
	gnn[j][k]=buff[k];
      }

    }
    
    
    fprintf (out,"\n\n kontrola cisel uzlu");
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	fprintf (out,"\n %4ld %4ld %4ld",i,j,gnn[i][j]);
      }
    }
    fprintf (out,"\n\n konec kontroly");


    
    //  total number of boundary nodes
    totnbn=0;
    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	if (totnbn<gnn[i][j])  totnbn=gnn[i][j];
      }
    }
    
    totnbn++;

    //fprintf (out,"\n\n pocet vsech hranicnich uzlu   %ld",totnbn);
    


    cnid = new long [totnbn];
    for (i=0;i<totnbn;i++){
      cnid[i]=0;
    }
    

    for (i=0;i<nproc;i++){
      for (j=0;j<nbndom[i];j++){
	cnid[gnn[i][j]]++;
      }
    }
    


    fprintf (out,"\n\n kontrola cnid");
    for (i=0;i<totnbn;i++){
      fprintf (out,"\n %4ld %4ld",i,cnid[i]);
    }
    fprintf (out,"\n\n konec kontroly");


    
    totncn=1;  totnrn=1;
    for (i=0;i<totnbn;i++){
      if (cnid[i]>2){
	cnid[i]=0-totncn;
	totncn++;
      }
      else{
	cnid[i]=totnrn;
	totnrn++;
      }
    }
    totncn--;  totnrn--;		
    

    fprintf (out,"\n\n pocet corner points  %ld",totncn);
    fprintf (out,"\n\n pocet remaining points  %ld",totnrn);

    fprintf (out,"\n\n kontrola cnid");
    for (i=0;i<totnbn;i++){
      fprintf (out,"\n %4ld %4ld",i,cnid[i]);
    }
    fprintf (out,"\n\n konec kontroly");



    for (i=1;i<nproc;i++){
      j=domproc[i];
      for (k=0;k<nbndom[j];k++){
	buff[k]=cnid[gnn[j][k]];
      }
      MPI_Send (buff,maxnbn,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    j=domproc[0];
    for (k=0;k<nbndom[j];k++){
      buff[k]=cnid[gnn[j][k]];
    }

  }
  else{

    MPI_Send (buff,maxnbn,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (buff,maxnbn,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  //fprintf (out,"\n\n kontrola ltg \n");

  j=0;
  for(i=0;i<top->nn;i++){
    if (ltg[i]>-1){
      ltg[i]=buff[j]-1;
      j++;
    }
    //fprintf (out,"\n ltg %4ld   %ld",i,ltg[i]);
    if (ltg[i]>=-1)  top->gnodes[i].ai = 0;
    if (ltg[i]<-1)   top->gnodes[i].ai = 0-ltg[i]-1;
  }
  
  j=0;
  for(i=0;i<top->nn;i++){
    if (ltg[i]<-1){
      j++;
    }
  }
  
  fprintf (stdout,"\n pocet corner points %ld",j);
  

  stop=0;
  if (myrank==0){
    if (j<3)  stop=1;
    for (i=1;i<nproc;i++){
      MPI_Recv (&j,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (j<3)  stop=1;
    }
    for (i=1;i<nproc;i++){
      MPI_Send (&stop,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (&j,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&stop,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);


  if (stop==1){
    fprintf (stderr,"\n\n insufficient number of corner points (domain %d, number of corner points %ld).\n",ndom,j);
    fprintf (stderr,"\n\n insufficient number of corner points (domain %d, number of corner points %ld).\n",ndom,j);
    //MPI_Finalize ();
    abort ();
  }


  if (myrank==0){
    for (i=0;i<nproc;i++){
      delete [] gnn[i];
    }
    delete [] gnn;
    
    delete [] cnid;
    delete [] nbndom;
  }
  delete [] buff;


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
   JK, 6.8.2007
*/
void dpfeti::initiate (selnodes *selnodschur,selnodes *selnodfeti)
{
  //  number of Lagrange multipliers on subdomain
  nlbdof = selnodfeti->ndof;
  //  maximum number of multipliers on subdomains
  maxnbdof = selnodfeti->maxndof;
  //  local code numbers for multipliers on subdomain
  localcn = selnodfeti->ldof;
  
  
  if (myrank==0){
    //  total number of multipliers in problem
    tnmdof = selnodfeti->tndof;
    //  number of Lagrange multipliers on subdomains
    nlbdofdom = selnodfeti->ndofdom;
    //  global code numbers for multipliers on subdomains
    globcn = selnodfeti->cndom;
  }
  
  
  
  //  maximum number of coarse DOFs on subdomains
  maxncdof = selnodschur->maxndof;
  //  number of coarse DOFs on subdomain
  ncdof = selnodschur->ndof;
  
  if (myrank==0){
    //  total number of coarse DOFs in problem
    tncdof = selnodschur->tndof;
    //  number of coarse DOFs on subdomains
    ncdofdom=selnodschur->ndofdom;
    
    ptop = new gtopology;
    ptop->initiate (selnodschur->cndom,selnodschur->ndofdom,nproc);

  }
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
  long **cncn,***cnbn,*nbndom,*ncndom,**dominc,**supelemcn;
  long buffsize,*buff;
  MPI_Status stat;
  
  fprintf (out,"\n\n kontrola kodovych cisel \n");
  for (i=0;i<top->nn;i++){
    fprintf (out,"\n uzel %5ld   ",i+1);
    for (j=0;j<top->give_ndofn (i);j++){
      fprintf (out," %ld",top->give_dof (i,j));
    }
  }
  fprintf (out,"\n\n konec kontroly kodovych  cisel\n");

  
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
    
    nbndom[domproc[0]]=nbn;
    ncndom[domproc[0]]=ncn;
    
    //  maxnbn - maximum number of boundary nodes on one subdomain
    //  maxncn - maximum number of corner nodes on one subdomain
    maxnbn=nbn;  maxncn=ncn;
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];

      nbndom[j]=buff[0];
      ncndom[j]=buff[1];
      
      if (maxnbn<buff[0])  maxnbn=buff[0];
      if (maxncn<buff[1])  maxncn=buff[1];
    }
    
    buff[0]=maxnbn;
    buff[1]=maxncn;
    for (i=1;i<nproc;i++){
      MPI_Send (buff,2,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,2,MPI_LONG,0,myrank,MPI_COMM_WORLD);
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
    for (j=0;j<nbndom[domproc[0]];j++){
      gnbndom[domproc[0]][j]=buff[l];  l++;
      for (k=0;k<ndofn;k++){
	lcnbn[domproc[0]][j][k]=buff[l];
	l++;
      }
    }
    for (j=0;j<ncndom[domproc[0]];j++){
      gncndom[domproc[0]][j]=buff[l];  l++;
      for (k=0;k<ndofn;k++){
	lcncn[domproc[0]][j][k]=buff[l];
	l++;
      }
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      m=domproc[stat.MPI_TAG];
      l=0;
      for (j=0;j<nbndom[m];j++){
	gnbndom[m][j]=buff[l];  l++;
	for (k=0;k<ndofn;k++){
	  lcnbn[m][j][k]=buff[l];
	  l++;
	}
      }
      for (j=0;j<ncndom[m];j++){
	gncndom[m][j]=buff[l];  l++;
	for (k=0;k<ndofn;k++){
	  lcncn[m][j][k]=buff[l];
	  l++;
	}
      }
      
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
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
    
    /*
    fprintf (out,"\n\n tncn %ld",tncn);
    for (i=0;i<tncn;i++){
      fprintf (out,"\n corner node %ld  supports",i);
      for (j=0;j<ndofn;j++){
	fprintf (out,"  %ld",cncn[i][j]);
      }
    }
    */
    
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
    
    /*
    fprintf (out,"\n\n tncn %ld",tncn);
    for (i=0;i<tncn;i++){
      fprintf (out,"\n corner node %ld  code numbers",i);
      for (j=0;j<ndofn;j++){
	fprintf (out,"  %ld",cncn[i][j]);
      }
    }
    */
    
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
    
    /*
    fprintf (out,"\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n ncdofdom %ld",ncdofdom[i]);
    }
    */
    
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
    

    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n procesor %ld,  ncdofdom %ld",i,ncdofdom[i]);
      for (j=0;j<ncdofdom[i];j++){
	fprintf (out,"\n supelemcn %ld",supelemcn[i][j]);
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
    
    /*
    fprintf (out,"\n\n node-domain incidences");
    fprintf (out,"\n\n tnbn %ld",tnbn);
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n boundary node %ld,  nibn %ld, domains",i,nibn[i]);
      for (j=0;j<nibn[i];j++){
	fprintf (out,"  %ld",dominc[i][j]);
      }
    }
    */
    
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
    
    /*
    fprintf (out,"\n\n support information\n");
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n boundary node %ld, supp",i);
      for (j=0;j<nibn[m]-1;j++){
	fprintf (out,"   ");
	for (k=0;k<ndofn;k++){
	  fprintf (out," %ld",cnbn[i][j][k]);
	}
      }
    }
    */
    
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
    
    /*
    fprintf (out,"\n\n code numbers\n");
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n boundary node %ld, cod num",i);
      for (j=0;j<nibn[m]-1;j++){
	fprintf (out,"   ");
	for (k=0;k<ndofn;k++){
	  fprintf (out," %ld",cnbn[i][j][k]);
	}
      }
    }
    */
    
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
    

    fprintf (out,"\n\n nlbdofdom \n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nlbdofdom %ld",nlbdofdom[i]);
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
    

    fprintf (out,"\n\n kontrola loccn na masterovi");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n procesor %ld  loccn",i);
      fprintf (out,"\n pole nlbdofdom %ld",nlbdofdom[i]);
      for (j=0;j<nlbdofdom[i];j++){
	fprintf (out," %ld",loccn[i][j]);
      }
      fprintf (out,"  globcn");
      for (j=0;j<nlbdofdom[i];j++){
	fprintf (out," %ld",globcn[i][j]);
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

    delete [] nbndom;

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
    
    for (i=1;i<nproc;i++){
      buff[2]=nlbdofdom[domproc[i]];
      MPI_Send (buff,3,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    nlbdof=nlbdofdom[domproc[0]];
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
    for (i=1;i<nproc;i++){
      k=domproc[i];
      for (j=0;j<nlbdofdom[k];j++){
	buff[j]=loccn[k][j];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    k=domproc[0];
    for (j=0;j<nlbdofdom[k];j++){
      buff[j]=loccn[k][j];
    }
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  for (j=0;j<nlbdof;j++){
    localcn[j]=buff[j];
  }

  delete [] buff;
  

  fprintf (out,"\n\n\n kontrolni tisk lokalnich kodovych cisel");
  for (i=0;i<nlbdof;i++){
    fprintf (out,"\n localcn %ld",localcn[i]);
  }

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
  
  /*
  for (i=0;i<nlbdof;i++){
    lcn=localcn[i];
    if (lcn==0)  continue;
    ev[i]=lv[lcn-1];
  }
  */

  for (i=0;i<nlbdof;i++){
    lcn=localcn[i];
    if (lcn<0){
      lcn=0-lcn-1;
      ev[i]=0.0-lv[lcn];
    }
    else{
      lcn=lcn-1;
      ev[i]=lv[lcn];
    }
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
  
  /*
  for (i=0;i<nlbdof;i++){
    lcn=localcn[i];
    if (lcn==0)  continue;
    lv[lcn-1]+=ev[i];
  }
  */

  for (i=0;i<nlbdof;i++){
    lcn=localcn[i];
    if (lcn<0){
      lcn=0-lcn-1;
      lv[lcn]-=ev[i];
    }
    else{
      lcn=lcn-1;
      lv[lcn]+=ev[i];
    }
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
  
  /*
  for (i=0;i<nlbdofdom[nsub];i++){
    gcn=globcn[nsub][i];
    if (gcn>0)
      ev[i]=gv[gcn-1];
    if (gcn<0)
      ev[i]=0-gv[0-gcn-1];
  }
  */

  for (i=0;i<nlbdofdom[nsub];i++){
    gcn=globcn[nsub][i]-1;
    ev[i]=gv[gcn];
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

  /*
  for (i=0;i<nlbdofdom[nsub];i++){
    gcn=globcn[nsub][i];
    if (gcn>0)
      gv[gcn-1]+=ev[i];
    if (gcn<0)
      gv[0-gcn-1]-=ev[i];
  }
  */

  for (i=0;i<nlbdofdom[nsub];i++){
    gcn=globcn[nsub][i]-1;
    gv[gcn]+=ev[i];
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

void dpfeti::arrmatrix (double *condmat,long *domproc,FILE *out)
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
    arr->setval (ssle);
    arr->initiate (ptop,tncdof,1);
    arr->prepmat (ptop,0.0,1);
    
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
    arr->localized (condmat,cn,gndofe,domproc[0]);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof*maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=domproc[stat.MPI_TAG];
      
      ii=0;
      for (j=0;j<ncdofdom[l];j++){
	for (k=0;k<ncdofdom[l];k++){
	  condmat[j*ncdofdom[l]+k]=buff[ii];  ii++;
	}
      }
      
      //fprintf (out,"\n\n\n\n\n\n\n prispevek od domeny %ld",l);
      //for (j=0;j<ncdof;j++){
      //fprintf (out,"\n condmat %4ld  %f",j,condmat[j*ncdofdom[l]+j]);
      //}
  
  

      gndofe = ptop->give_ndofe (l);
      cn = new long [gndofe];
      ptop->give_code_numbers (l,cn);
      arr->localized (condmat,cn,gndofe,l);
      delete [] cn;
      
    }
    
    arr->decompose_matrix ();
    //fprintf (out,"\n\n kontrola matice v arr");
    //arr->printmat (out);
  }
  else{
    MPI_Send(buff,maxncdof*maxncdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
}


void dpfeti::vectors_br_bm (gtopology *top,gmatrix *gm,
			    double *condmat,double *condvect,double *rhs,long *domproc,FILE *out)
{
  long i,j,l,gndofe,*cn;
  double *lhs,*buff,*zaloha;
  MPI_Status stat;
  
  // ************
  //  vector b_r
  // ************
  buff = new double [maxncdof];
  for (i=0;i<ncdof;i++){
    buff[i]=condvect[i];
  }
  
  if (myrank==0){
    br = new double [tncdof];
    nullv (br,tncdof);
    
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
    locglob (br,condvect,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      l=domproc[stat.MPI_TAG];
      
      for (j=0;j<ncdofdom[l];j++){
	condvect[j]=buff[j];
      }
      
      gndofe = ptop->give_ndofe (l);
      cn = new long [gndofe];
      ptop->give_code_numbers (l,cn);
      locglob (br,condvect,cn,gndofe);
      delete [] cn;
      
    }
    
    /*
    fprintf (out,"\n\n kontrola br \n");
    for (i=0;i<tncdof;i++){
      fprintf (out,"\n br %ld %le",i,br[i]);
    }
    */

  }
  else{
    MPI_Send(buff,maxncdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;
  
  // ************
  //  vector b_m
  // ************
  

    
  buff = new double [maxnbdof];

  //lhs = new double [nidof];
  //zaloha = new double [nidof];
  //for (i=0;i<nidof;i++){
  //zaloha[i]=rhs[i];
  //}
  
  lhs = new double [ndof];
  zaloha = new double [ndof];
  for (i=0;i<nidof;i++){
    zaloha[i]=rhs[i];
  }
  
  //fprintf (out,"\n\n kontrola Richarda");
  //fprintf (out,"\n nidof   %ld",nidof);
  //fprintf (out,"\n ncdof   %ld",ncdof);
  //fprintf (out,"\n ndof    %ld",ndof);

  //  (K_{rr}^{j})^{-1} . f_{r}^{j}
  gm->condense (top,condmat,lhs,zaloha,ncdof,3);
  delete [] zaloha;
  

  //  G^{j} . (K_{rr}^{j})^{-1} . f_{r}^{j}
  nullv (buff,maxnbdof);
  extract_from_local_vector (buff,lhs);
  
  if (myrank==0){
    bm = new double [tnmdof];
    nullv (bm,tnmdof);
    put_into_global_vector (buff,bm,domproc[0]);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      put_into_global_vector (buff,bm,domproc[stat.MPI_TAG]);
    }

    /*
    fprintf (out,"\n\n\n vector bm \n");
    for (i=0;i<tnmdof;i++){
      fprintf (out,"\n bm %ld   %e",i,bm[i]);
    }
    */

  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);


  delete [] buff;
  delete [] lhs;


}


void dpfeti::rhs_dpfeti (gtopology *top,long *domproc,gmatrix *gm,FILE *out)
{
  long i,j,gndofe,*cn;
  double *ur,*wr,*utc,*vtc,*buff;
  MPI_Status stat;
  
  //ur = new double [nidof];
  //wr = new double [nidof];
  ur = new double [ndof];
  wr = new double [ndof];
  buff = new double [maxncdof];
  
  if (myrank==0){
    utc = new double [tncdof];
    vtc = new double [tncdof];
    
    for (i=0;i<tncdof;i++){
      utc[i]=br[i];
    }
    
    /*
    fprintf (out,"\n\n\n vector utc \n");
    for (i=0;i<tncdof;i++){
      fprintf (out,"\n utc %ld   %e",i,utc[i]);
    }
    */
    
    //  A_{RR}^{-1} . b_{R} = \mu
    arr->back_substitution (vtc,utc);
    
    /*
    fprintf (out,"\n\n\n vector vtc \n");
    for (i=0;i<tncdof;i++){
      fprintf (out,"\n vtc %ld   %e",i,vtc[i]);
    }
    */

    for (i=1;i<nproc;i++){
      //  F_c^j . \mu
      nullv (buff,maxncdof);
      gndofe = ptop->give_ndofe (domproc[i]);
      cn = new long [gndofe];
      ptop->give_code_numbers (domproc[i],cn);
      globloc (vtc,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
    
    nullv (buff,maxncdof);
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
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
  mxv (krc,buff,ur,nidof,ncdof);
  delete [] buff;
  
  /*
  fprintf (out,"\n\n\n vector ur \n");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n ur %ld   %e",i,ur[i]);
  }
  */

  //  (K_{rr}^j)^{-1} . [K_{rc}^j . F_c^j . \mu]
  gm->condense (top,buff,wr,ur,ncdof,3);
  
  /*
  fprintf (out,"\n\n\n vector wr \n");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n wr %ld   %e",i,wr[i]);
  }
  */

  buff = new double [maxnbdof];

  //  (G^j) . [(K_{rr}^j)^{-1} . K_{rc}^j . F_c^j . \mu]
  nullv (buff,maxnbdof);
  extract_from_local_vector (buff,wr);
  
  if (myrank==0){
    nullv (cgrhs,tnmdof);
    put_into_global_vector (buff,cgrhs,domproc[0]);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      put_into_global_vector (buff,cgrhs,domproc[stat.MPI_TAG]);
    }
  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  delete [] ur;
  delete [] wr;
  
  if (myrank==0){
    for (i=0;i<tnmdof;i++){
      cgrhs[i]=bm[i]-cgrhs[i];
    }
    
    /*
    fprintf (out,"\n\n\n vector cgrhs \n");
    for (i=0;i<tnmdof;i++){
      fprintf (out,"\n cgrhs %ld   %e",i,cgrhs[i]);
    }
    */

  }

}

void dpfeti::matxvect (gtopology *top,long *domproc,
		       gmatrix *gm,double *input,double *output)
{
  long i,j,gndofe,*cn;
  double *utc,*vtc,*ur,*vr,*wr,*buff;
  MPI_Status stat;
  
  ur = new double [ndof];
  vr = new double [ndof];
  wr = new double [ndof];
  //ur = new double [nidof];
  //vr = new double [nidof];
  //wr = new double [nidof];

  
  buff = new double [maxnbdof];

  if (myrank==0){
    for (i=1;i<nproc;i++){
      extract_from_global_vector (buff,input,domproc[i]);
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
    extract_from_global_vector (buff,input,domproc[0]);
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullv (ur,nidof);
  put_into_local_vector (buff,ur);
  delete [] buff;
  
  //  (K_{rr}^j)^{-1} . [(G^j)^T . \lambda]
  gm->condense (top,buff,vr,ur,ncdof,3);
  
  buff = new double [maxncdof];
  
  //  (K_{rc}^j)^T . [(K_{rr}^j)^{-1} . (G^j)^T . \lambda]
  mtxv (krc,vr,buff,nidof,ncdof);
  
  if (myrank==0){
    utc = new double [tncdof];
    vtc = new double [tncdof];
    nullv (utc,tncdof);
    
    //  (F_c^j)^T . [(K_{rc}^j)^T . (K_{rr}^j)^{-1} . (G^j)^T . \lambda] = \nu
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
    locglob (utc,buff,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      j=domproc[stat.MPI_TAG];
      gndofe = ptop->give_ndofe (j);
      cn = new long [gndofe];
      ptop->give_code_numbers (j,cn);
      locglob (utc,buff,cn,gndofe);
      delete [] cn;
    }

    //  A_{RR}^{-1} . A_{RM} . \lambda = \mu
    arr->back_substitution (vtc,utc);
    
    for (i=1;i<nproc;i++){
      //  F_c^j . \mu
      nullv (buff,maxncdof);
      gndofe = ptop->give_ndofe (domproc[i]);
      cn = new long [gndofe];
      ptop->give_code_numbers (domproc[i],cn);
      globloc (vtc,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }

    nullv (buff,maxncdof);
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
    globloc (vtc,buff,cn,gndofe);
    delete [] cn;
    
    delete [] utc;
    delete [] vtc;
  }
  else{
    MPI_Send(buff,maxncdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
    MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  K_{rc}^j . [F_c^j . \mu]
  mxv (krc,buff,ur,nidof,ncdof);
  delete [] buff;
  
  //  (K_{rr}^j)^{-1} . [K_{rc}^j . F_c^j . \mu]
  gm->condense (top,buff,wr,ur,ncdof,3);
  
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
    nullv (output,tnmdof);
    put_into_global_vector (buff,output,domproc[0]);
    for (j=1;j<nproc;j++){
      MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      put_into_global_vector (buff,output,domproc[stat.MPI_TAG]);
    }
  }
  else{
    MPI_Send(buff,maxnbdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
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
      aerrcgdpfeti=norrhs;  anicgdpfeti=0;
      return;
    }
  }
  
  
  //  iteration loop
  stop=0;
  for (i=0;i<nicgdpfeti;i++){
    
    
    //  new coefficient alpha
    matxvect (top,domproc,gm,d,p);
    
    if (myrank==0){
      denom = ss (d,p,tnmdof);
      if (fabs(denom)<1.0e-25){
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
      
      if (sqrt(nom/norrhs)<errcgdpfeti)  stop=1;
      //if (fabs(nom)<limit)  break;
      
      
      beta = nom/denom;
      
      //  new vector of direction
      for (j=0;j<tnmdof;j++){
	d[j]=beta*d[j]-r[j];
      }

      for (j=1;j<nproc;j++){
	MPI_Send(&stop,1,MPI_LONG,j,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv(&stop,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (stop==1)  break;
  }
  
  anicgdpfeti=i;  aerrcgdpfeti=nom;
  
  if (myrank==0){
    delete [] p;  delete [] r;  delete [] d;
  }
  
}


void dpfeti::corner_displ (gtopology *top,long *domproc,gmatrix *gm)
{
  long i,j,gndofe,*cn;
  double *utc,*ur,*vr,*buff;
  MPI_Status stat;
  
  //ur = new double [nidof];
  //vr = new double [nidof];
  ur = new double [ndof];
  vr = new double [ndof];
  
  buff = new double [maxnbdof];
  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      extract_from_global_vector (buff,cglhs,domproc[i]);
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
    extract_from_global_vector (buff,cglhs,domproc[0]);
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  (G^j)^T . \lambda
  nullv (ur,nidof);
  put_into_local_vector (buff,ur);
  
  //  (K_{rr}^j)^{-1} . [(G^j)^T . \lambda]
  gm->condense (top,buff,vr,ur,ncdof,3);
  
  delete [] buff;
  buff = new double [maxncdof];

  //  (K_{rc}^j)^T . [(K_{rr}^j)^{-1} . (G^j)^T . \lambda]
  mtxv (krc,vr,buff,nidof,ncdof);
  
  if (myrank==0){
    utc = new double [tncdof];
    nullv (utc,tncdof);
    
    //  (F_c^j)^T . [(K_{rc}^j)^T . (K_{rr}^j)^{-1} . (G^j)^T . \lambda] = \nu
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
    locglob (utc,buff,cn,gndofe);
    delete [] cn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv(buff,maxncdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      j=domproc[stat.MPI_TAG];
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
    MPI_Send(buff,maxncdof,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
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
    for (i=1;i<nproc;i++){
      //  F_c^j . \mu
      nullv (buff,maxncdof);
      gndofe = ptop->give_ndofe (domproc[i]);
      cn = new long [gndofe];
      ptop->give_code_numbers (domproc[i],cn);
      globloc (displ,buff,cn,gndofe);
      delete [] cn;
      
      MPI_Send(buff,maxncdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }

    nullv (buff,maxncdof);
    gndofe = ptop->give_ndofe (domproc[0]);
    cn = new long [gndofe];
    ptop->give_code_numbers (domproc[0],cn);
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
  mxv (krc,buff,ur,nidof,ncdof);
  delete [] buff;
  
  
  buff = new double [maxnbdof];
  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      extract_from_global_vector (buff,cglhs,domproc[i]);
      MPI_Send(buff,maxnbdof,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD);
    }
    extract_from_global_vector (buff,cglhs,domproc[0]);
  }
  else{
    MPI_Recv(buff,maxnbdof,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  //  (G^j)^T . \lambda
  nullv (vr,nidof);
  put_into_local_vector (buff,vr);
  delete [] buff;
  
  for (i=0;i<nidof;i++){
    rhs[i]=rhs[i]-ur[i]-vr[i];
  }
  
  gm->condense (top,buff,subdispl,rhs,ncdof,3);
  
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
  for (i=0;i<nidof*ncdof;i++){
    krc[i]=0.0;
  }

  //  K_rc matrix assembling
  gm->a12block (krc,ncdof);
  
  
  //fprintf (out,"\n\n\n KONTROLA Krc \n");
  //for (i=0;i<nidof;i++){
  //fprintf (out,"\n %4ld",i);
  //for (j=0;j<ncdof;j++){
  //fprintf (out," %f",krc[i*ncdof+j]);
  //}
  //}


  crhs = new double [ncdof];
  //rrhs = new double [nidof];
  rrhs = new double [ndof];
  
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
  gm->condense (top,condmat,lhs,rhs,ncdof,1);
  
  
  j=0;
  for (i=nidof;i<ndof;i++){
    condvect[j]=rhs[i];
    j++;
  }
  
  t3 = time (NULL);
  
  /*
  fprintf (out,"\n\n CONDMAT");
  for (i=0;i<ncdof;i++){
    fprintf (out,"\n");
    for (j=0;j<ncdof;j++){
      fprintf (out,"  %le",condmat[i*ncdof+j]);
    }
  }
  
  fprintf (out,"\n\n CONDVECT");
  for (i=0;i<ncdof;i++){
    fprintf (out,"  %le",condvect[i]);
  }
  */
  
  arrmatrix (condmat,domproc,out);


  t4 = time (NULL);

  //  computation of vectors b_R and b_M
  vectors_br_bm (top,gm,condmat,condvect,rrhs,domproc,out);
  
  
  //  computation of right hand side of final system of equations
  rhs_dpfeti (top,domproc,gm,out);

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

  
  MPI_Status stat;
  long *buff;
  buff = new long[7];
  buff[0]=t2-t1;
  buff[1]=t3-t2;
  buff[2]=t4-t3;
  buff[3]=t5-t4;
  buff[4]=t6-t5;
  buff[5]=t7-t6;
  buff[6]=t8-t7;

  if (myrank==0){
    
    fprintf (out,"\n\n\n\n\n");
    j=domproc[myrank];
    fprintf (out,"\n\n\n Domain %ld",j);
    fprintf (out,"\n A_rc block assembling            %ld",buff[0]);
    fprintf (out,"\n time of condensation             %ld",buff[1]);
    fprintf (out,"\n A_rr assembling                  %ld",buff[2]);
    fprintf (out,"\n b_R and b_M assembling           %ld",buff[3]);
    fprintf (out,"\n solution of reduced system       %ld",buff[4]);
    fprintf (out,"\n corner displacements calcul.     %ld",buff[5]);
    fprintf (out,"\n all displacements calcul.        %ld",buff[6]);

    
    j=domproc[myrank];
    fprintf (stdout,"\n\n\n Domain %ld",j);
    fprintf (stdout,"\n A_rc block assembling            %ld",buff[0]);
    fprintf (stdout,"\n time of condensation             %ld",buff[1]);
    fprintf (stdout,"\n A_rr assembling                  %ld",buff[2]);
    fprintf (stdout,"\n b_R and b_M assembling           %ld",buff[3]);
    fprintf (stdout,"\n solution of reduced system       %ld",buff[4]);
    fprintf (stdout,"\n corner displacements calcul.     %ld",buff[5]);
    fprintf (stdout,"\n all displacements calcul.        %ld",buff[6]);


    for (i=1;i<nproc;i++){
      MPI_Recv (buff,7,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      fprintf (out,"\n\n\n Domain %ld",j);
      fprintf (out,"\n A_rc block assembling            %ld",buff[0]);
      fprintf (out,"\n time of condensation             %ld",buff[1]);
      fprintf (out,"\n A_rr assembling                  %ld",buff[2]);
      fprintf (out,"\n b_R and b_M assembling           %ld",buff[3]);
      fprintf (out,"\n solution of reduced system       %ld",buff[4]);
      fprintf (out,"\n corner displacements calcul.     %ld",buff[5]);
      fprintf (out,"\n all displacements calcul.        %ld",buff[6]);
      

      fprintf (stdout,"\n\n\n Domain %ld",j);
      fprintf (stdout,"\n A_rc block assembling            %ld",buff[0]);
      fprintf (stdout,"\n time of condensation             %ld",buff[1]);
      fprintf (stdout,"\n A_rr assembling                  %ld",buff[2]);
      fprintf (stdout,"\n b_R and b_M assembling           %ld",buff[3]);
      fprintf (stdout,"\n solution of reduced system       %ld",buff[4]);
      fprintf (stdout,"\n corner displacements calcul.     %ld",buff[5]);
      fprintf (stdout,"\n all displacements calcul.        %ld",buff[6]);
      
    }
  }
  else{
    MPI_Send(buff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;


  delete [] condmat;
  delete [] condvect;
  delete [] rrhs;
  delete [] crhs;


}

