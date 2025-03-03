#include <mpi.h>
#include "partopjk.h"
#include "psolver.h"
#include <stdlib.h>
#include <string.h>

/**
   constructor
   
   @param np - the number of processors
   @param mr - my rank (process id)
   @param meshd - type of mesh description
   @param *nameproc - name of processor
   @param nl - length of processor name  
   
   JK, 22.8.2011
*/
partopjk::partopjk(int np,int mr,meshdescription meshd,char *nameproc, int nl)
{
  //  the number of processors
  nproc=np;
  //  my rank
  myrank=mr;
  //  mesh description
  md = meshd;
  // name of processor
  strcpy (procName,nameproc);
  nameLength=nl;
  
  //  the number of nodes on subdomain
  nn = 0;
  //  the maximum number of nodes on subdomain
  maxnn = 0;
  //  the total number of nodes in the whole problem
  tnnp = 0;
  //  the maximum number of DOFs on subdomains
  maxndofd = 0;
  //  the number of internal nodes
  nin = 0;
  //  the number of boundary/interface nodes on subdomain
  nbn = 0;
  //  the number of internal DOFs on subdomain
  nidof = 0;
  //  the number of boundary/interface DOFs on subdomain
  nbdof = 0;
  //  the number of DOFs on subdomain
  ndof = 0;
  //  the maximum number of boundary/interface DOFs on subdomain
  maxnbdofd = 0;
  //  maximum number of boundary/interface nodes with constraints collected from all subdomains
  maxnbnwcd = 0;
  //  maximum number of DOFs at boundary/interface nodes with constraints collected from all subdomains
  maxndofbn = 0;
  // number of boundary/interface nodes with constraints on the given subdomain, i.e. dimension of array bnwc
  nbnwc = 0;
   
  
  // **************************
  //  variables on the master
  // **************************
  //  the total number of nodes in the whole problem
  tnnp = 0;
  //  the number of DOFs in the coarse problem
  ndofc = 0;
  

  // ***************************
  //  arrays on all processors
  // ***************************
  //  array containing boundary/interface nodes on subdomains
  bnid=NULL;
  //  array containing internal nodes on subdomains
  inid=NULL;
  // array of local index numbers of boundary/interface nodes with constraints on the given subdomain
  bnwc = NULL;

  // ***********************
  //  arrays on the master
  // ***********************
  //  array containing the numbers of nodes on subdomains (M)
  nnsd = NULL;
  //  array of local to global correspondence on the master (M)
  mltg = NULL;
  //  array of nodal multiplicity, the multiplicity is the number of subdomains which share the node (M)
  nodmultip = NULL;
  //  array of node-domain correspondence
  noddom = NULL;
  //  array of numbers of DOFs on nodes (M)
  ndofnm = NULL;
  //  array of DOF indicators on nodes
  dofm = NULL;

  // array of numbers of boundary/interface nodes with constraints on individual subdomains
  nbnwcd = NULL;
  // array of number of boundary/interface nodes on individual subdomains
  bnwcd = NULL;
  //  array containing the numbers of boundary/interface nodes on subdomains
  nbnd = NULL;  
  //  array containing the numbers of internal nodes on subdomains
  nind = NULL;
  //  array containing the numbers of boundary/interface DOFs on subdomains
  nbdofd = NULL;
  //  array containing coarse code numbers of boundary/interface DOFs on subdomains
  bdofd = NULL;
  ///  array containing coarse code numbers of boundary/interface DOFs on subdomain
  bdof = NULL;
  //  array of full node-domain correspondence
  fullnoddom = NULL;
  //  array of Lagrange multipliers, it contains coarse code numbers
  lagrmultip = NULL;
}

/**
   destructor
   
   JK, 22.8.2011
*/
partopjk::~partopjk()
{
  long i;

  if (myrank==0){
    //  array containing the numbers of nodes on subdomains (M)
    delete [] nnsd;

    //  array of local to global correspondence on the master (M)
    if (mltg) {
      for (i=0;i<nproc;i++){
        delete [] mltg[i];
      }
      delete [] mltg;
    }
    
    //  array of nodal multiplicity, the multiplicity is the number of subdomains which share the node (M)
    delete [] nodmultip;
    
    //  array of node-domain correspondence
    delete [] noddom;

    //  array of numbers of DOFs on nodes (M)
    delete [] ndofnm;
    
    //  array of DOF indicators on nodes
    if (dofm){
      for (i=0;i<tnnp;i++){
        delete [] dofm[i];
      }
      delete [] dofm;
    }
    
    //  array containing the numbers of boundary/interface nodes on subdomains
    delete [] nbnd;
    //  array containing the numbers of internal nodes on subdomains
    delete [] nind;
    //  array containing the numbers of boundary/interface DOFs on subdomains
    delete [] nbdofd;
    
    //  array containing coarse code numbers of boundary/interface DOFs on subdomains
    if (bdofd){
      for (i=0;i<nproc;i++){
        delete [] bdofd[i];
      }
      delete [] bdofd;
    }
    
    
    //  array of full node-domain correspondence
    if (fullnoddom){
      for (i=0;i<tnnp;i++){
        delete [] fullnoddom[i];
      }
      delete [] fullnoddom;
    }

    //  array of Lagrange multipliers, it contains coarse code numbers
    if (lagrmultip){
      for (i=0;i<tnnp;i++){
        delete [] lagrmultip[i];
      }
      delete [] lagrmultip;
    }
    
    delete [] bdof;

    delete [] nbnwcd;
    for (i=0; i<nproc; i++)
      delete [] bnwcd[i];
    delete [] bnwcd;
  }
  
  //  array containing boundary/interface nodes on subdomains
  delete [] bnid;
  //  array containing internal nodes on subdomains
  delete [] inid;
  delete [] bnwc;
}

/**
   function assembles the array nnsd on the master processor
   variable nn (the number of particular subdomain) is assigned
   
   @param top - pointer to general topology
   @param domproc - domain-processor correspondence
   @param out - output file
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::numbers_of_nodes_on_subdomains (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k;
  MPI_Status stat;
  
  //  the number of nodes on subdomain
  nn = top->nn;
  
  if (myrank==0){
    if (nnsd!=NULL){
      delete [] nnsd;
    }
    nnsd = new long [nproc];
    
    //  master contribution
    j=domproc[0];
    nnsd[j]=nn;
    if (maxnn<nn){
      maxnn=nn;
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      nnsd[j]=k;
      if (maxnn<k){
	maxnn=k;
      }
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
  }
  else{
    MPI_Send (&nn,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxnn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  auxiliary output
  fprintf (out,"\n\n maxnn (maximum number of nodes on subdomain) %ld",maxnn);
  if (myrank==0){
    fprintf (out,"\n\n array nnsd (number of nodes on subdomains");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nnsd %6ld  %ld",i+1,nnsd[i]);
    }
  }//  end of myrank
  //  end of auxiliary output
}

/**
   function assembles the array mltg on the master processor
   
   @param ltg - local to global correspondence
   @param domproc - domain-processor correspondence
   @param out - output file
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::ltg_on_master (long *ltg,long *domproc,FILE *out)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize=maxnn;
  buff = new long [buffsize];
  
  for (i=0;i<nn;i++){
    buff[i]=ltg[i];
  }
  
  if (myrank==0){
    //  array on the master for ltg arrays
    mltg = new long* [nproc];
    for (i=0;i<nproc;i++){
      mltg[i] = new long [nnsd[i]];
    }
    
    //  master contribution
    j=domproc[0];
    for (k=0;k<nnsd[j];k++){
      mltg[j][k]=buff[k];
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      for (k=0;k<nnsd[j];k++){
	mltg[j][k]=buff[k];
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n mltg (arrays ltg collected on the master)");
    for (i=0;i<nproc;i++){
      for (j=0;j<nnsd[i];j++){
	fprintf (out,"\n mltg %6ld %6ld   %ld",i+1,j+1,mltg[i][j]);
      }
    }
  }//  end of myrank
  //  end of auxiliary output
}

/**
   function converts any mesh description to the all_nodes mesh description
   
   @param ltg - local to global correspondence
   @param domproc - domain-processor correspondence
   @param out - output file
   @param proc_name - processor name
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::ltg_conversion (long */*ltg*/,long */*domproc*/,FILE *out,char *proc_name)
{
  long i,j,max;
  
  if (myrank==0){
    
    // ******************************************
    //  the maximum node number is searched for
    // ******************************************
    switch (md){
    case all_nodes:{
      max=0;
      for (i=0;i<nproc;i++){
	for (j=0;j<nnsd[i];j++){
	  if (max<mltg[i][j]){
	    max = mltg[i][j];
	  }
	}
      }
      break;
    }
    case bound_nodes:{
      max=0;
      for (i=0;i<nproc;i++){
	for (j=0;j<nnsd[i];j++){
	  if (max<mltg[i][j]){
	    max = mltg[i][j];
	  }
	}
      }
      break;
    }
    case neg_bound_nodes:{
      max=0;
      for (i=0;i<nproc;i++){
	for (j=0;j<nnsd[i];j++){
	  if (max<abs(mltg[i][j])){
	    max = abs(mltg[i][j]);
	  }
	}
      }
      break;
    }
    default:{
      par_print_err(myrank,proc_name,"wrong type of mesh description is required",__FILE__, __LINE__, __func__);
    }
    }//  end of the switch
    max++;

    // ******************************************
    //  conversion to all_nodes
    // ******************************************
    switch (md){
    case all_nodes:{
      //  no conversion is needed
      break;
    }
    case bound_nodes:{
      for (i=0;i<nproc;i++){
	for (j=0;j<nnsd[i];j++){
	  if (mltg[i][j]==-1){
	    mltg[i][j]=max;
	    max++;
	  }
	}
      }
      break;
    }
    case neg_bound_nodes:{
      for (i=0;i<nproc;i++){
	for (j=0;j<nnsd[i];j++){
	  if (mltg[i][j]<-1){
	    mltg[i][j]=0-mltg[i][j];
	  }
	}
      }
      break;
    }
    default:{
      par_print_err(myrank,proc_name,"wrong type of mesh description is required",__FILE__, __LINE__, __func__);
    }
    }//  end of the switch
    
    //  the total number of nodes in the whole problem
    tnnp=max;
    
  }//  end of myrank

  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n tnnp (the total number of nodes in the whole problem) %ld",tnnp);
    fprintf (out,"\n\n\n mltg (arrays ltg collected on the master) after conversion");
    for (i=0;i<nproc;i++){
      for (j=0;j<nnsd[i];j++){
	fprintf (out,"\n mltg %6ld %6ld   %ld",i+1,j+1,mltg[i][j]);
      }
    }
  }//  end of myrank
  //  end of auxiliary output
}

/**
   function assembles array nodmultip
   
   @param out - output file
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::node_multiplicity (FILE *out)
{
  long i,j;
  
  if (myrank==0){
    nodmultip = new long [tnnp];
    for (i=0;i<tnnp;i++){
      nodmultip[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nnsd[i];j++){
	nodmultip[mltg[i][j]]++;
      }
    }
  }
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n nodmultip (array of nodal multiplicity for all nodes in the problem)");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n nodmultip %6ld   %ld",i+1,nodmultip[i]);
    }
  }//  end of myrank
  //  end of auxiliary output
}

/**
   function assembles array noddom
   
   @param out - output file
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::node_domain (FILE *out)
{
  long i,j,nid;
  
  if (myrank==0){
    noddom = new long [tnnp];
    
    //  loop over the number of subdomains/processors
    for (i=0;i<nproc;i++){
      //  loop over the number of nodes on actual subdomains
      for (j=0;j<nnsd[i];j++){
	//  node id
	nid = mltg[i][j];
	noddom[nid]=i;
      }
    }
  }
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n noddom (node-domain correspondence)");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n noddom %6ld   %ld",i+1,noddom[i]);
    }
  }//  end of myrank
  //  end of auxiliary output
}


/**
   function assembles array ndofnm
   function determines the variable maxndofd
   
   @param top - pointer to general topology
   @param domproc - domain-processor correspondence
   @param out - output file
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::ndofn_on_master (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize = maxnn+1;
  buff = new long [buffsize];
  
  maxndofd=0;
  for (i=0;i<nn;i++){
    buff[i] = top->give_ndofn (i);
    maxndofd+=buff[i];
  }
  buff[nn]=maxndofd;
  
  if (myrank==0){
    ndofnm = new long [tnnp];
    
    //  master contribution
    j=domproc[0];
    for (k=0;k<nnsd[j];k++){
      ndofnm[mltg[j][k]]=buff[k];
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      for (k=0;k<nnsd[j];k++){
	ndofnm[mltg[j][k]]=buff[k];
      }
      if (maxndofd<buff[nnsd[j]]){
	maxndofd=buff[nnsd[j]];
      }
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxndofd,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxndofd,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  // neposila se maxndofd nize podruhe ?? - TKo
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Send (&maxndofd,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }  
  else{
    MPI_Recv (&maxndofd,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  //  auxiliary output
  fprintf (out,"\n\n maxndofd (maximum number of DOFs on subdomain)  %ld",maxndofd);
  if (myrank==0){
    fprintf (out,"\n\n\n ndofnm (array of nodal ndofn)");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n ndofnm %6ld   %ld",i+1,ndofnm[i]);
    }
  }//  end of myrank
  //  end of auxiliary output
  
}

/**
   function assembles array dofm
   
   @param top - pointer to general topology
   @param domproc - domain-processor correspondence
   @param out - output file
   
   JK, 1.7.2005, revision 22.8.2011
*/
void partopjk::dof_on_master (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,l,m,nid,ndofn,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize=maxndofd;
  buff = new long [buffsize];
  
  k=0;
  for (i=0;i<nn;i++){
    ndofn = top->give_ndofn (i);
    for (j=0;j<ndofn;j++){
      buff[k]=top->give_dof (i,j);
      k++;
    }
  }
  
  if (myrank==0){
    dofm = new long* [tnnp];
    for (i=0;i<tnnp;i++){
      dofm[i] = new long [ndofnm[i]];
    }
    
    //  master contribution
    j=domproc[0];
    l=0;
    for (k=0;k<nnsd[j];k++){
      //  node id
      nid=mltg[j][k];
      //  the number of DOFs on the actual node
      ndofn=ndofnm[nid];
      for (m=0;m<ndofn;m++){
	dofm[nid][m]=buff[l];
	l++;
      }
    }
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      l=0;
      for (k=0;k<nnsd[j];k++){
	//  node id
	nid=mltg[j][k];
	//  the number of DOFs on the actual node
	ndofn=ndofnm[nid];
	for (m=0;m<ndofn;m++){
	  dofm[nid][m]=buff[l];
	  l++;
	}
      }
    }//  end of the loop
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n dofm (array of nodal DOFs on the master)");
    for (i=0;i<tnnp;i++){
      for (j=0;j<ndofnm[i];j++){
	fprintf (out,"\n dofm %6ld %6ld   %ld",i+1,j+1,dofm[i][j]);
      }
    }
  }//  end of myrank
  //  end of auxiliary output
   
}

/**
   function detects coupled DOFs in the problem
   coupled DOFs within one subdomain are not important
   the coupled DOFs on two or more subdomains have to be detected
   also DOFs coupled with any boundary/interface DOF have to be found
   such coupled DOFs have to be denoted as the boundary/interface DOFs
   
   DOF indicators
   dofm[i][j]<0 - nonzero prescribed value of the DOF
   dofm[i][j]=0 - zero prescribed value of the DOF
   dofm[i][j]=1 - unknown degree of freedom
   dofm[i][j]>1 - coupled DOF
   
   @param out - output file
   
   JK, 1.7.2005, revision 23.8.2011
*/
void partopjk::coupled_dofs_detection (FILE *out)
{
  long i,j,k,max,dofid,stop,nid,nid2,sdid,sdid2;
  long *ncdof,**cdof;
  
  if (myrank==0){
    //  the maximum DOF indicator is determined
    max=0;
    for (i=0;i<tnnp;i++){
      for (j=0;j<ndofnm[i];j++){
	if (max<dofm[i][j]){
	  max=dofm[i][j];
	}
      }
    }
    
    if (max>1){
      //  there are coupled DOFs in the problem

      ncdof = new long [max];
      for (i=0;i<max;i++){
	ncdof[i]=0;
      }
      for (i=0;i<tnnp;i++){
	for (j=0;j<ndofnm[i];j++){
	  if (dofm[i][j]>1){
	    ncdof[dofm[i][j]-2]++;
	  }
	}
      }
      
      cdof = new long* [max-1];
      for (i=0;i<max-1;i++){
	cdof[i] = new long [ncdof[i]];
	ncdof[i]=0;
      }
      

      //  loop over the total number of nodes in the problem
      for (i=0;i<tnnp;i++){
	//  loop over the number of DOFs in nodes
	for (j=0;j<ndofnm[i];j++){
	  if (dofm[i][j]>1){
	    //  the actual DOF is coupled
	    
	    //  DOF id
	    dofid=dofm[i][j]-2;
	    //  node id is stored in the array cdof
	    cdof[dofid][ncdof[dofid]]=i;
	    //  actual position in the array cdof is incremented
	    ncdof[dofid]++;
	  }
	}
      }//  end of the loop over the total number of nodes in the problem
      
      fprintf (out,"\n\n\n jouda \n\n");
      for (i=0;i<max-1;i++){
	fprintf (out,"\n ncdof %ld  %ld",i,ncdof[i]);
	for (j=0;j<ncdof[i];j++){
	  fprintf (out,"   cdof %ld %ld %ld",i,j,cdof[i][j]);
	}
      }
      fprintf (out,"\n jouda \n\n\n\n\n");

      //  checking of the nodes with coupled DOFs
      //  loop over the suspicious DOFs
      for (i=0;i<max-1;i++){

	stop=0;
	//  loop over the suspicious DOFs with the same DOF indicator
	for (j=0;j<ncdof[i];j++){
	  //  node id
	  nid=cdof[i][j];
	  if (nodmultip[nid]>1){
	    //  there is a node with multiplicity greater than 1
	    //  it means, some coupled DOFs are connected with DOF on the boundary/interface
	    //  all DOFs with the actual DOF indicator are denoted as the boundary/interface DOFs
	    
	    //  loop over the DOFs with actual DOF indicator
	    for (k=0;k<ncdof[i];k++){
	      nid2=cdof[i][k];
	      if (nodmultip[nid2]==1){
		nodmultip[nid2]=-1;
	      }
	    }//  end of the loop over the DOFs with actual DOF indicator
	    stop=1;
	    break;
	  }//  end of if
	}//  end of the loop over the suspicious DOFs with the same DOF indicator
	
	
	if (stop==0){
	  //  the actual DOF is not connected with any boundary/interface node
	  //  the subdomain id's are checked now
	  
	  //  node id
	  nid=cdof[i][0];
	  //  subdomain id
	  sdid=noddom[nid];
	  //  loop over the suspicious DOFs with the actual DOF indicator
	  for (j=0;j<ncdof[i];j++){
	    //  node id
	    nid2=cdof[i][j];
	    //  subdomain id
	    sdid2=noddom[nid2];
	    if (sdid!=sdid2){
	      //  there is coupled DOF on two different subdomains
	      
	      //  loop over the DOFs with actual DOF indicator
	      for (k=0;k<ncdof[i];k++){
		nid2=cdof[i][k];
		if (nodmultip[nid2]==1){
		  nodmultip[nid2]=-1;
		}
	      }//  end of the loop over the DOFs with actual DOF indicator
	      
	      stop=1;
	      break;
	    }//  end of if statement
	    
	  }//  end of the loop over the suspicious DOFs with the actual DOF indicator
	}//  end of if statement
      }//  end of the loop over the suspicious DOFs with the same DOF indicator
      

      for (i=0;i<max-1;i++){
	delete [] cdof[i];
      }
      delete [] cdof;
      
      delete [] ncdof;
      
    }//  end of the if statement
  }//  end of myrank
  
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n nodmultip (array of nodal multiplicity for all nodes in the problem)");
    fprintf (out,"\n the array was modified due to coupled DOFs");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n nodmultip %6ld   %ld",i+1,nodmultip[i]);
    }
  }//  end of myrank
  //  end of auxiliary output
  
}



/**
   function detects internal and boundary/interface nodes
   
   @param out - output file
   
   JK, 24.8.2011
*/
void partopjk::internal_boundary_nodes_detection (FILE *out)
{
  long i,j,k,nid,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize = maxnn+2;
  buff = new long [buffsize];
  
  if (myrank==0){
    fprintf (out,"\n nbnd = %p\n",(void*)nbnd);
    if (nbnd!=NULL){ // proc test na realokaci ? priprava na rostouci kce ? TKo
      delete [] nbnd;
    }
    nbnd = new long [nproc];
    fprintf (out,"\n nbnd = %p\n",(void*)nbnd);
    if (nind!=NULL){
      delete [] nind;
    }
    nind = new long [nproc];
    
    //  slave contributions

    //  loop over the number of subdomains
    for (i=1;i<nproc;i++){
      nbnd[i]=0;
      nind[i]=0;

      //  loop over the number of nodes on the actual subdomain
      for (j=0;j<nnsd[i];j++){
	//  node id
	nid = mltg[i][j];
	
	if (nodmultip[nid]>1 || nodmultip[nid]<0){
	  buff[j]=2;
	  nbnd[i]++;
	}
	else{
	  buff[j]=1;
	  nind[i]++;
	}//  end of if statement (nodmultip[nid]>1 || nodmultip[nid]<0)
      }//  end of the loop over the number of nodes on the actual subdomain
      buff[buffsize-2]=nind[i];
      buff[buffsize-1]=nbnd[i];
      
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }//  end of the loop over the number of subdomains
    
    //  master contribution
    i=0;
    nbnd[i]=0;
    nind[i]=0;
    
    //  loop over the number of nodes on the actual subdomain
    for (j=0;j<nnsd[i];j++){
      //  node id
      nid = mltg[i][j];
      
      if (nodmultip[nid]>1 || nodmultip[nid]<0){
	buff[j]=2;
	nbnd[i]++;
      }
      else{
	buff[j]=1;
	nind[i]++;
      }//  end of if statement (nodmultip[nid]>1 || nodmultip[nid]<0)
    }//  end of the loop over the number of nodes on the actual subdomain
    buff[buffsize-2]=nind[i];
    buff[buffsize-1]=nbnd[i];
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }//  end of if statement (myrank)
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  the number of internal node on actual subdomain
  nin = buff[buffsize-2];
  //  the number of boundary/interface nodes on actual subdomain
  nbn = buff[buffsize-1];
  //  array of internal node numbers
  inid = new long [nin];
  //  array of boundary/interface node numbers
  bnid = new long [nbn];
  
  //  loop over the number of nodes on actual subdomain
  j=0;  k=0;
  for (i=0;i<nn;i++){
    if (buff[i]==1){
      inid[j]=i;
      j++;
    }
    if (buff[i]==2){
      bnid[k]=i;
      k++;
    }
  }//  end of the loop over the number of nodes on actual subdomain
  
  delete [] buff;
  
  
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n nbnd (the number of boundary/interface nodes on subdomains)");
    fprintf (out,"\n nind (the number of internal nodes on subdomains)");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld nbnd  %6ld   nind  %6ld",i+1,nbnd[i],nind[i]);
    }
  }

  fprintf (out,"\n\n\n nin (the number of internal nodes  %ld",nin);
  fprintf (out,"\n nbn (the number of boundary/interface nodes  %ld",nbn);
  fprintf (out,"\n\n\n inid (array of internal node numbers)");
  for (i=0;i<nin;i++){
    fprintf (out,"\n inid %6ld   %ld",i+1,inid[i]);
  }
  fprintf (out,"\n\n\n bnid (array of boundary/interface node numbers)");
  for (i=0;i<nbn;i++){
    fprintf (out,"\n bnid %6ld   %ld",i+1,bnid[i]);
  }
  //  end of auxiliary output
   
}



/**
   function detects boundary/interface nodes with some constraints, i.e. some dof of boundary node is non-positive.
   
   @param out - output file
   
   Created by Tomas Koudelka, 8.2021
*/
void partopjk::boundary_nodes_with_constraints_detection (FILE *out)
{
  long i, j, k, l, nid, constr;
  ivector buff(3);
  MPI_Status stat;
  
  
  
  if (myrank==0){    

    nbnwcd = new long[nproc];
    bnwcd  = new long*[nproc];
    
    maxndofbn = 0;    
    maxnbnwcd = 0;
    //  loop over the number of all subdomains
    for (i=0; i<nproc; i++){
      nbnwcd[i] = 0;      
      //  loop over the number of nodes on the actual subdomain
      for (j=0;j<nnsd[i];j++){
	//  node id
	nid = mltg[i][j];
        constr = 0;
	if (nodmultip[nid]>1 || nodmultip[nid]<0){
          for(k=0; k<ndofnm[nid]; k++){
            if (dofm[nid][k] <= 0){
              constr = 1;
              break;
            }
          }
          if (constr){ // boundary node with constraint(s)            
            nbnwcd[i]++; // increase the number of boundary nodes with constraints on the given domain
            if (maxndofbn < ndofnm[nid]) // check maximum number of DOFs at boundary nodes
              maxndofbn = ndofnm[nid];
          }
	}
      }
      if (maxnbnwcd < nbnwcd[i])
        maxnbnwcd = nbnwcd[i];
    }
    
    //  loop over the number of all subdomains
    for (i=0; i<nproc; i++){
      bnwcd[i] = new long[nbnwcd[i]];
      l = 0;
      //  loop over the number of nodes on the actual subdomain
      for (j=0;j<nnsd[i];j++){
	//  node id
	nid = mltg[i][j];
        constr = 0;
	if (nodmultip[nid]>1 || nodmultip[nid]<0){
          for(k=0; k<ndofnm[nid]; k++){
            if (dofm[nid][k] <= 0){
              constr = 1;
              break;
            }
          }
          if (constr){ // boundary node with constraint(s)
            bnwcd[i][l] = j;
            l++;
          }
	}
      }
    }

    nbnwc = nbnwcd[0];
    buff(1) = maxnbnwcd;
    buff(2) = maxndofbn;
    //  loop over the number of slave subdomains - scatter number of boundary nodes with constraints and maximum values
    for (i=1; i<nproc; i++){
      buff(0) = nbnwcd[i];
      MPI_Send (buff.a, 3, MPI_LONG, i, myrank, MPI_COMM_WORLD);
    }//  end of the loop over the number of subdomains
  }    
  else{
    MPI_Recv (buff.a, 3, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    nbnwc     = buff(0);
    maxnbnwcd = buff(1);
    maxndofbn = buff(2);
  }//  end of if statement (myrank)
  MPI_Barrier (MPI_COMM_WORLD);

  // scatter arrays of local index numbers of boundary nodes with constraints to individual subdomains
  reallocv(maxnbnwcd, buff);
  bnwc = new long[nbnwc];
  if (myrank==0){
    // scatter on master
    i = 0;
    copyv(bnwcd[i], bnwc, nbnwcd[i]);
    nbnwc = nbnwcd[i];
    
    //  loop over the number of slave subdomains - scatter 
    for (i=1; i<nproc; i++){
      copyv(bnwcd[i], buff.a, nbnwcd[i]);
      MPI_Send (buff.a, maxnbnwcd, MPI_LONG, i, myrank, MPI_COMM_WORLD);
    }//  end of the loop over the number of subdomains
  }
  else{
    MPI_Recv (buff.a, maxnbnwcd, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    copyv(buff.a, bnwc, nbnwc);
  }//  end of if statement (myrank)
  MPI_Barrier (MPI_COMM_WORLD);

  
  //  auxiliary output
  if(myrank == 0){
    fprintf (out,"\n\n\n maximum number of boundary/interface nodes with constraints from all subdomains: maxnbnwcd = %ld", maxnbnwcd);
    fprintf (out,"\n maximum number of DOFs at boundary/interface nodes with constraints from all subdomains: maxndofbn = %ld", maxndofbn);
    fprintf(out, "\n Master: nbnwc = %ld", nbnwc);
    for(i=0; i<nproc; i++){
      fprintf(out, "\n nbnwcd[%ld] = %ld", i+1, nbnwcd[i]);
    }
    fprintf(out, "\n");
    for(i=0; i<nproc; i++){
      for(j=0; j<nbnwcd[i]; j++)
        fprintf(out, "\nbnwcd[%ld][%ld] = %ld (gn=%ld)", i+1, j+1, bnwcd[i][j]+1, mltg[i][bnwcd[i][j]]+1);
    }
  }
  else{
    fprintf (out,"\n\n\n maximum number of boundary/interface nodes with constraints from all subdomains: maxnbnwcd = %ld", maxnbnwcd);
    fprintf (out,"\n maximum number of DOFs at boundary/interface nodes with constraints from all subdomains: maxndofbn = %ld", maxndofbn);
    fprintf(out, "\n Slave: nbnwc = %ld", nbnwc);
    fprintf(out, "\n");
    for(i=0; i<nbnwc; i++)
      fprintf(out, "\nbnwc[%ld] = %ld", i+1, bnwc[i]);
  }
  //  end of auxiliary output
}



/**
   function assembles ordering of unknowns on subdomains
   
   @param top - pointer to general topology
   @param out - output file
   
   JK, 25.8.2011
*/
long partopjk::schur_local_ordering (gtopology *top,FILE *out)
{
  long i,j,k,nid,ndofn;
  long *aux;
  
  // ********************************
  //  detection of coupled DOFs
  // ********************************
  ndof=0;
  //  loop over the all nodes on subdomains
  for (i=0;i<nn;i++){
    //  number of DOFs in required node
    ndofn = top->give_ndofn (i);
    //  loop over the number of DOFs in actual node
    for (j=0;j<ndofn;j++){
      //  DOF indicator
      k=top->give_dof (i,j);
      if (ndof<k){
	ndof=k;
      }
    }//  end of the loop over the number of DOFs in actual node
  }//  end of the loop over the all nodes on subdomains
  
  //  conversion to C notation
  ndof--;
  if (ndof<0){
    ndof=0;
  }
  aux = new long [ndof];
  for (i=0;i<ndof;i++){
    aux[i]=-1;
  }
  
  
  // ********************************
  //  ordering of internal unknowns
  // ********************************
  ndof=1;
  //  loop over the internal nodes
  for (i=0;i<nin;i++){
    //  node id (in local ordering)
    nid=inid[i];
    //  number of DOFs at required node
    ndofn = top->give_ndofn (nid);
    //  loop over the number of DOFs in actual node
    for (j=0;j<ndofn;j++){
      //  DOF indicator
      k=top->give_dof (nid,j);
      if (k<0){
	//  non-zero prescribed value
	continue;
      }
      if (k==0){
	//  prescribed value equal to zero
	continue;
      }
      if (k==1){
	//  un-coupled DOF
	top->save_dof (nid,j,ndof);
	ndof++;
      }
      if (k>1){
	//  coupled DOF
	if (aux[k-2]==-1){
	  top->save_dof (nid,j,ndof);
	  aux[k-2]=ndof;
	  ndof++;
	}
	else{
	  top->save_dof (nid,j,aux[k-2]);
	}
      }
    }//  end of the loop over the number of DOFs in actual node
  }//  end of the loop over the internal nodes

  //  number of internal DOFs (unknowns)
  nidof = ndof-1;
  
  // ******************************************
  //  ordering of boundary/interface unknowns
  // ******************************************
  //  loop over the boundary/interface nodes
  for (i=0;i<nbn;i++){
    //  node id (in local ordering)
    nid=bnid[i];
    //  number of DOFs at required node
    ndofn = top->give_ndofn (nid);
    //  loop over the number of DOFs in actual node
    for (j=0;j<ndofn;j++){
      //  DOF indicator
      k=top->give_dof (nid,j);
      if (k<0){
	//  non-zero prescribed value
	continue;
      }
      if (k==0){
	//  prescribed value equal to zero
	continue;
      }
      if (k==1){
	//  un-coupled DOF
	top->save_dof (nid,j,ndof);
	ndof++;
      }
      if (k>1){
	//  coupled DOF
	if (aux[k-2]==-1){
	  top->save_dof (nid,j,ndof);
	  aux[k-2]=ndof;
	  ndof++;
	}
	else{
	  top->save_dof (nid,j,aux[k-2]);
	}
      }
    }//  end of the loop over the number of DOFs in actual node
  }//  end of the loop over the boundary/interface nodes
  
  delete [] aux;
  
  //  the number of DOFs on actual subdomain
  ndof--;
  
  //  number of boundary/interface DOFs (unknowns)
  nbdof = ndof-nidof;

  //  state of code numbers is changed
  //  code numbers are generated
  top->cnstate=1;
  
  //  auxiliary output
  fprintf (out,"\n\n\n nidof (the number of internal DOFs   %ld",nidof);
  fprintf (out,"\n nbdof (the number of boundary/interface DOFs   %ld",nbdof);
  fprintf (out,"\n ndof (the number of all DOFs   %ld",ndof);
  fprintf (out,"\n\n\n code number");
  //  loop over the all nodes on subdomains
  for (i=0;i<nn;i++){
    fprintf (out,"\n node  %6ld   ",i+1);
    //  number of DOFs in required node
    ndofn = top->give_ndofn (i);
    //  loop over the number of DOFs in actual node
    for (j=0;j<ndofn;j++){
      fprintf (out,"  %6ld",top->give_dof (i,j));
    }//  end of the loop over the number of DOFs in actual node
  }//  end of the loop over the all nodes on subdomains
  fprintf (out,"\n\n");
  //  end of auxiliary output
  
  return ndof;
}

/**
   function assembles ordering of the coarse problem with respect to the Schur complement method
   
   only nodes on boundary/interface or nodes with coupled DOFs across
   subdomains are taken into account
   
   @param out - output file
   
   JK, 24.8.2011
*/
void partopjk::schur_coarse_ordering (FILE *out)
{
  long i,j,k;
  long *aux;
  
  if (myrank==0){
    //  the number of DOFs in the coarse problem
    ndofc=0;
    
    //  loop over the number of all nodes in the problem
    for (i=0;i<tnnp;i++){
      if (nodmultip[i]>1 || nodmultip[i]<0){
	//  node is shared by more than two subdomains
	//  or the node contains coupled DOF which is connected with other domain
	
	//  loop over the number of DOFs defined in the actual node
	for (j=0;j<ndofnm[i];j++){
	  if (ndofc<dofm[i][j]){
	    ndofc=dofm[i][j];
	  }
	}//  end of the loop over the number of DOFs defined in the actual node
      }//  end of if statement (nodmultip[i]>1 || nodmultip[i]<0)
    }//  end of the loop over the number of all nodes in the problem
    
    //  conversion to C notation
    ndofc--;
    if (ndofc<0){
      ndofc=0;
    }
    aux = new long [ndofc];
    for (i=0;i<ndofc;i++){
      aux[i]=-1;
    }
    
    ndofc=1;
    //  loop over the number of all nodes in the problem
    for (i=0;i<tnnp;i++){
      if (nodmultip[i]>1 || nodmultip[i]<0){
	//  node is shared by more than two subdomains
	//  or the node contains coupled DOF which is connected with other domain
	
	//  loop over the number of DOFs in the actual node
	for (j=0;j<ndofnm[i];j++){
	  //  DOF indicator
	  k=dofm[i][j];
	  if (k<0){
	    //  non-zero prescribed value
	    continue;
	  }
	  if (k==0){
	    //  prescribed value equal to zero
	    continue;
	  }
	  if (k==1){
	    //  un-coupled DOF
	    dofm[i][j]=ndofc;
	    ndofc++;
	  }
	  if (k>1){
	    //  coupled DOF
	    if (aux[k-2]==-1){
	      dofm[i][j]=ndofc;
	      aux[k-2]=ndofc;
	      ndofc++;
	    }
	    else{
	      dofm[i][j]=aux[k-2];
	    }
	  }
	}//  end of the loop over the number of DOFs in the actual node
      }//  end of if statement (nodmultip[i]>1 || nodmultip[i]<0)
    }//  end of the loop over the number of all nodes in the problem
    ndofc--;
    
    delete [] aux;
    
  }//  end of if statement (myrank)
  


  
  //  DOF indicators of internal nodes are deleted because of check of coarse code numbers
  for (i=0;i<tnnp;i++){
    if (nodmultip[i]==1){
      for (j=0;j<ndofnm[i];j++){
	dofm[i][j]=0;
      }
    }
  }
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n ndofc   %ld",ndofc);
    fprintf (out,"\n\n\n dofm after code numbers generation ");
    fprintf (out,"\n dofm (array of nodal DOFs on the master)");
    for (i=0;i<tnnp;i++){
      for (j=0;j<ndofnm[i];j++){
	fprintf (out,"\n dofm %6ld %6ld   %ld",i+1,j+1,dofm[i][j]);
      }
    }
  }//  end of myrank
  //  end of auxiliary output
  
}

/**
   function assembles coarse code numbers for particular subdomains
   
   only nodes on boundary/interface or nodes with coupled DOFs across
   subdomains are taken into account
   
   @param out - output file
   
   JK, 24.8.2011
*/
void partopjk::coarse_code_numbers (FILE *out)
{
  long i,j,k,nid,dofid,buffsize;
  long *aux,*buff;
  MPI_Status stat;
  
  if (myrank==0){
    if (nbdofd!=NULL){
      delete [] nbdofd;
    }
    nbdofd = new long [nproc];
    
    aux = new long [ndofc];
    
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      for (j=0;j<ndofc;j++){
	aux[j]=0;
      }
      
      nbdofd[i]=0;
      
      //  loop over the number of nodes on actual subdomain
      for (j=0;j<nnsd[i];j++){
	//  node id
	nid = mltg[i][j];
	if (nodmultip[nid]>1 || nodmultip[nid]<0){
	  //  node is shared by more than two subdomains
	  //  or the node contains coupled DOF which is connected with other domain
	  for (k=0;k<ndofnm[nid];k++){
	    //  DOF id
	    dofid=dofm[nid][k];
	    if (dofid>0){
	      aux[dofid-1]++;
	      if (aux[dofid-1]==1){
		nbdofd[i]++;
	      }
	    }
	  }
	}//  end of if statement (nodmultip[nid]>1 || nodmultip[nid]<0)
      }//  end of the loop over the number of nodes on actual subdomain
    }//  end of the loop over the number of subdomains
    
    if (bdofd!=NULL){
      for (i=0;i<nproc;i++){
	delete [] bdofd[i];
      }
      delete [] bdofd;
    }
    bdofd = new long* [nproc];
    for (i=0;i<nproc;i++){
      fprintf (stdout,"\n nbdofd %ld \n",nbdofd[i]);
      bdofd[i] = new long [nbdofd[i]];
    }
    

    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      for (j=0;j<ndofc;j++){
	aux[j]=0;
      }

      nbdofd[i]=0;
      
      //  loop over the number of nodes on actual subdomain
      for (j=0;j<nnsd[i];j++){
	//  node id
	nid = mltg[i][j];
	if (nodmultip[nid]>1 || nodmultip[nid]<0){
	  //  node is shared by more than two subdomains
	  //  or the node contains coupled DOF which is connected with other domain
	  for (k=0;k<ndofnm[nid];k++){
	    //  DOF id
	    dofid=dofm[nid][k];
	    if (dofid>0){
	      aux[dofid-1]++;
	      if (aux[dofid-1]==1){
		bdofd[i][nbdofd[i]]=dofm[nid][k];
		nbdofd[i]++;
	      }
	    }
	  }
	}//  end of if statement (nodmultip[nid]>1 || nodmultip[nid]<0)
      }//  end of the loop over the number of nodes on actual subdomain
    }//  end of the loop over the number of subdomains

    delete [] aux;
  }//  end of if statement (myrank)
  

  if (myrank==0){
    maxnbdofd = 0;
    for (i=0;i<nproc;i++){
      fprintf (stdout,"\n nbdofd %ld \n",nbdofd[i]);
      if (maxnbdofd<nbdofd[i]){
	maxnbdofd=nbdofd[i];
      }
    }
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnbdofd,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&maxnbdofd,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  
  //  auxiliary output
  fprintf (out,"\n\n\n maxnbdofd (numbers of boundary/interface DOFs on subdomains)  %ld",maxnbdofd);
  if (myrank==0){
    fprintf (out,"\n\n\n nbdofd (numbers of boundary/interface DOFs on subdomains)");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n nbdofd %6ld   %ld",i+1,nbdofd[i]);
    }
    fprintf (out,"\n\n\n nbdofd (numbers of boundary/interface DOFs on subdomains)");
    for (i=0;i<nproc;i++){
      for (j=0;j<nbdofd[i];j++){
	fprintf (out,"\n bdofd %6ld %6ld   %ld",i+1,j+1,bdofd[i][j]);
      }
    }
  }//  end of myrank
  //  end of auxiliary output

  
  buffsize=maxnbdofd;
  buff = new long [buffsize];
  
  if (myrank==0){
    
    fprintf (out,"\n nbnd = %(void*)p\n",nbnd);
    for (i=1;i<nproc;i++){
      fprintf (out," smycka  i %ld",i);
      for (j=0;j<nbdofd[i];j++){
	buff[j]=bdofd[i][j];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    i=0;
    for (j=0;j<nbdofd[i];j++){
      buff[j]=bdofd[i][j];
    }
    
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  bdof = new long [nbdof];
  for (i=0;i<nbdof;i++){
    bdof[i]=buff[i];
  }

  delete [] buff;
  

  //  auxiliary output
  fprintf (out,"\n\n\n bdof (boundary/interface DOFs on subdomain");
  for (j=0;j<nbdof;j++){
    fprintf (out,"\n bdof %6ld    %ld",j+1,bdof[j]);
  }
  //  end of auxiliary output

  
}



long partopjk::vse (long *ltg,gtopology *top,long *domproc,FILE *out,char *proc_name)
{
  long ndof;

  fprintf (out,"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n ZACATEK NOVEHO PARTOPU \n\n\n\n\n\n\n\n\n");

  //  function assembles the array nnsd which contains the numbers of nodes on subdomains
  numbers_of_nodes_on_subdomains (top,domproc,out);

  //  function assembles the array mltg on the master which stores all ltg arrays from processors
  ltg_on_master (ltg,domproc,out);
  
  //  function converts any mesh description to the all_nodes mesh description
  ltg_conversion (ltg,domproc,out,proc_name);
  
  //  function assembles array nodmultip which contains nodal multiplicity
  node_multiplicity (out);
  
  //  function assembles array noddom
  node_domain (out);

  //  function assembles array ndofnm which contains the number of DOFs for each node
  ndofn_on_master (top,domproc,out);
  
  //  function assembles array dofm which contains DOF indicators on each node
  dof_on_master (top,domproc,out);

  //  function detects coupled DOFs in the problem
  coupled_dofs_detection (out);

  //  function detects internal and boundary/interface nodes
  internal_boundary_nodes_detection (out);
  
  //  function detects internal and boundary/interface nodes
  boundary_nodes_with_constraints_detection (out);

  //  function assembles ordering of unknowns on subdomains
  ndof = schur_local_ordering (top,out);

  //  function assembles ordering of the coarse problem with respect to the Schur complement method
  schur_coarse_ordering (out);

  //  function assembles coarse code numbers for particular subdomains
  coarse_code_numbers (out);

  fprintf (out,"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n KONEC NOVEHO PARTOPU \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");

  return ndof;

}





/**
   function assembles array fullnoddom
   FETI
   
   @param out - output file

   JK, 15.10.2011
*/
void partopjk::full_node_domain (FILE *out)
{
  long i,j,k,l,min;
  
  if (myrank==0){
    if (fullnoddom!=NULL){
      for (i=0;i<tnnp;i++){
	delete [] fullnoddom[i];
      }
      delete [] fullnoddom;
    }
    fullnoddom = new long* [tnnp];
    for (i=0;i<tnnp;i++){
      fullnoddom[i] = new long [nodmultip[i]];
      //  array nodmultip will be used for definition of actual position in array fullnoddom
      //  at the end of this subroutine, nodmultip will be restored
      nodmultip[i]=0;
    }
    
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      //  loop over the number of nodes on particular subdomain
      for (j=0;j<nnsd[i];j++){
	//  the j-th node on the i-th subdomain has global number k
	k=mltg[i][j];
	//  actual position for the k-th node (in global ordering)
	l=nodmultip[k];
	fullnoddom[k][l]=i;
	nodmultip[k]++;
      }//  end of the loop over the number of nodes on particular subdomain
    }//  end of the loop over the number of subdomains
    
    //  sorting of subdomain numbers
    //  loop over the number of nodes
    for (i=0;i<tnnp;i++){
      //  loop over the number of domains which share the node
      for (j=0;j<nodmultip[i];j++){
	min=nproc;  l=-1;
	for (k=j;k<nodmultip[i];k++){
	  if (min>fullnoddom[i][k]){
	    min=fullnoddom[i][k];
	    l=k;
	  }
	}
	min=fullnoddom[i][j];
	fullnoddom[i][j]=fullnoddom[i][l];
	fullnoddom[i][l]=min;
      }//  end of the loop over the number of domains which share the node
    }//  end of the loop over the number of nodes
    
  }//  end of myrank=0
  
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n fullnoddom (array of full node-domain correspondence)");
    for (i=0;i<tnnp;i++){
      for (j=0;j<nodmultip[i];j++){
	fprintf (out,"\n fullnoddom  %6ld %6ld   %ld",i+1,j+1,fullnoddom[i][j]);
      }
    }
  }//  end of myrank
  //  end of auxiliary output
  
}


/**
   function generates Lagrange multipliers
   FETI
   
   @param ps - parallel solver
   @param out - output file

   JK, 16.10.2011
*/
void partopjk::lagrange_multipliers (psolver ps,FILE *out)
{
  long i,j,k,l,m,ndofn;
  
  if (myrank==0){
    lagrmultip = new long* [tnnp];
    for (i=0;i<tnnp;i++){
      if (ps.f1->fetiimpl == nonredundant){
	lagrmultip[i] = new long [(nodmultip[i]-1)*ndofnm[i]];
	for (j=0;j<(nodmultip[i]-1)*ndofnm[i];j++){
	  lagrmultip[i][j]=0;
	}
      }
      if (ps.f1->fetiimpl == redundant){
	lagrmultip[i] = new long [nodmultip[i]*(nodmultip[i]-1)*ndofnm[i]/2];
	for (j=0;j<nodmultip[i]*(nodmultip[i]-1)*ndofnm[i]/2;j++){
	  lagrmultip[i][j]=0;
	}
      }
    }

    //  the number of DOFs in the coarse problem
    ndofc=1;
    
    //  loop over the number of all nodes
    for (i=0;i<tnnp;i++){
      if (nodmultip[i]>1){
	//  only nodes shared by at least two subdomains are taken into account
	
	if (ps.f1->fetiimpl == nonredundant){
	  //  loop over the number of DOFs
	  ndofn=ndofnm[i];
	  for (j=0;j<ndofn;j++){
	    if (dofm[i][j]>0){
	      //  loop over the number of neighbours
	      for (k=0;k<nodmultip[i]-1;k++){
		lagrmultip[i][k*ndofn+j]=ndofc;
		ndofc++;
	      }//  end of the loop over the number of neighbours
	    }//  end of the if statement dofm[i][j]>0
	  }//  end of the loop over the number of DOFs
	}//  end of the if statement ps.f1->fetiimpl == nonredundant

	if (ps.f1->fetiimpl == redundant){
	  //  loop over the number of DOFs
	  ndofn=ndofnm[i];
	  for (j=0;j<ndofn;j++){
	    if (dofm[i][j]>0){
	      m=0;
	      //  loop over the number of neighbours
	      //  the following loop generates nodmultip[i]*(nodmultip[i]-1)/2 cycles
	      for (k=0;k<nodmultip[i]-1;k++){
		for (l=k;l<nodmultip[i]-1;l++){
		  lagrmultip[i][m*ndofn+j]=ndofc;
		  ndofc++;
		  m++;
		}
	      }//  end of the loop over the number of neighbours
	    }//  end of the if statement dofm[i][j]>0
	  }//  end of the loop over the number of DOFs
	}//  end of the if statement ps.f1->fetiimpl == redundant
	
	
      }//  end of the if statement nodmultip[i]>1
    }//  end of the loop over the number of all nodes
    
  }//  end of myrank=0
  
  
  //  auxiliary output
  if (myrank==0){
    fprintf (out,"\n\n\n lagrmultip ()");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n node %8ld (nodmultip %2ld)",i+1,nodmultip[i]);
      
      if (ps.f1->fetiimpl == nonredundant){
	for (j=0;j<(nodmultip[i]-1)*ndofnm[i];j++){
	  fprintf (out," %ld",lagrmultip[i][j]);
	}
      }
      
      if (ps.f1->fetiimpl == redundant){
	for (j=0;j<nodmultip[i]*(nodmultip[i]-1)*ndofnm[i]/2;j++){
	  fprintf (out," %ld",lagrmultip[i][j]);
	}
      }
      
    }
    fprintf (out,"\n\n");
  }//  end of myrank
  //  end of auxiliary output

}


/**
   function assembles coarse code numbers for particular subdomains
   
   only nodes on boundary/interface or nodes with coupled DOFs across
   subdomains are taken into account
   
   @param ps - parallel solver
   @param out - output file
   
   JK, 16.10.2011
*/
void partopjk::coarse_code_numbers_feti (psolver ps,FILE */*out*/)
{
  long i,j,k,l,nid,aux,buffsize,ndofn;
  long *buff;
  MPI_Status stat;
  
  if (myrank==0){
    maxnbdofd=0;
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      aux=0;
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
	//  global node number
	nid=mltg[i][j];
	if (nodmultip[nid]>1){
	  
	  if (ps.f1->fetiimpl == nonredundant){
	    if (fullnoddom[nid][0]==i){
	      for (k=0;k<ndofnm[nid];k++){
		if (dofm[nid][k]>0){
		  aux+=nodmultip[i]-1;
		}
	      }
	    }
	    else{
	      for (k=0;k<ndofnm[nid];k++){
		if (dofm[nid][k]>0){
		  aux++;
		}
	      }
	    }
	  }//  end of the if statement ps.f1->fetiimpl == nonredundant
	  
	  if (ps.f1->fetiimpl == redundant){
	    for (k=0;k<ndofnm[nid];k++){
	      if (dofm[nid][k]>0){
		aux+=nodmultip[i]-1;
	      }
	    }
	  }
	  
	}//  end of the if statement nodmultip[nid]>1
      }//  end of the loop over the number of nodes on subdomains
      
      if (maxnbdofd<aux){
	maxnbdofd=aux;
      }
    }//  end of the loop over the number of subdomains
  }//  end of the if statement myrank==0
  
  
  if (myrank==0){
    //  loop over the number of subdomains
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnbdofd,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }//  end of the loop over the number of subdomains
  }//  end of the if statement myrank==0
  else{
    MPI_Recv (&maxnbdofd,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }//  end of the else statement
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  buffsize=maxnbdofd;
  buff=new long [buffsize];
  
  
  if (myrank==0){
    //  loop over the number of subdomains
    for (i=1;i<nproc;i++){
      aux=0;
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
	//  global node number
	nid=mltg[i][j];
	if (nodmultip[nid]>1){
	  
	  if (ps.f1->fetiimpl == nonredundant){
	    if (fullnoddom[nid][0]==i){
	      ndofn=ndofnm[nid];
	      for (k=0;k<ndofn;k++){
		if (dofm[nid][k]>0){
		  for (l=0;l<nodmultip[nid];l++){
		    buff[aux]=lagrmultip[nid][l*ndofn+k];
		    aux++;
		  }
		}
	      }
	    }
	    else{
	      for (k=0;k<ndofnm[nid];k++){
		if (dofm[nid][k]>0){
		  for (l=0;l<nodmultip[nid];l++){
		    if (fullnoddom[nid][l]==i){
		      break;
		    }
		  }
		  buff[aux]=lagrmultip[nid][(l-1)*ndofn+k];
		  aux++;
		}
	      }
	    }
	  }//  end of the if statement ps.f1->fetiimpl == nonredundant
	  
	  if (ps.f1->fetiimpl == redundant){
	    for (k=0;k<ndofnm[nid];k++){
	      if (dofm[nid][k]>0){
		aux+=nodmultip[i]-1;
	      }
	    }
	  }//  end of the if statement ps.f1->fetiimpl == redundant
	  
	}//  end of the if statement nodmultip[nid]>1
      }//  end of the loop over the number of nodes on subdomains
      
      if (maxnbdofd<aux){
	maxnbdofd=aux;
      }
    }//  end of the loop over the number of subdomains
  }//  end of the if statement myrank==0
  

  
  
}
