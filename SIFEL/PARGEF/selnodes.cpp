#include "mpi.h"
#include "selnodes.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "galias.h"
#include "iotools.h"



/**
   constructor
   
   @param np - number of processors
   @param mr - my rank
   @param nd - number of subdomain
   @param k - number of nodes on subdomain
   @param j - array containing list of selected nodes
              it contains k components
   @param d - type of mesh description
   
   JK, 30.7.2007
*/
selnodes::selnodes(int np,int mr,long nd,long kk,long *j,meshdescription d,FILE *out,long mespr)
{
  long i,k;
  
  //  number of processors
  nproc=np;
  //  my rank
  myrank=mr;
  //  mesh description
  md = d;
  //  number of nodes
  nn=kk;

  //  number of selected nodes
  nsn=0;
  for (i=0;i<nn;i++){
    if (j[i]>-1)
      nsn++;
  }
  
  if (nsn==0){
    print_err("wrong number of selected nodes in constructor", __FILE__, __LINE__, __func__);
    abort();
  }
  
  //  list of selected nodes - local numbers
  lnsn=new long [nsn];
  //  list of selected nodes - global numbers
  lsng=new long [nsn];
  
  k=0;
  for (i=0;i<nn;i++){
    if (j[i]>-1){
      lnsn[k]=i;
      lsng[k]=j[i];
      k++;
    }
  }
  
  if (mespr==0){
    fprintf (out,"\n\n\n list of local and global/group numbers of selected nodes (lnsn, lsng)\n");
    fprintf (out,"\n number of selected nodes is %ld",nsn);
    for (i=0;i<nsn;i++){
      fprintf (out,"\n local node number %6ld       global/group node number %6ld",lnsn[i],lsng[i]);
    }
  }

  //  number of all DOFs including prescribed DOFs
  ndof=0;
  //  maximum number of selected nodes on subdomain
  maxnsn=0;
  //  maximum number of DOFs on subdomain
  maxndof=0;
  //  maximum number of contributing nodes on subdomain in the FETI method
  maxncn=0;
  
  ncn=0;

  ldof = NULL;
  
  ldofmultip = NULL;  // not used

  //  numbers of DOFs on subdomains
  ndofdom = NULL;  // not used
  //  code numbers at selected nodes on master
  cnm = NULL;  // not used

  //  coarse numbers of selected nodes on one subdomain
  cnsn = NULL;
  //  global numbers of selected nodes on one subdomain
  gnsn = NULL;

  if (myrank==0){
    //  total number of selected nodes
    tnsn=0;
    
    //  group node numbers
    gnn = NULL;
    //  node multiplicty
    nodmultip = NULL;
    // dof multiplicity
    dofmultip = NULL; 
    //  code numbers on master
    cndom = NULL;
    //  joint nodes to selected nodes assumed as coarse nodes
    ljn = NULL;
    //  subdomain numbers which contain connected nodes to coarse nodes
    lsn = NULL;
    //  code numbers / indicators for FETI method
    doffeti = NULL;
    //  number of contributing nodes in the FETI method
    ncndom = NULL;

    //  number of multiplicity of all boundary/interface nodes
    icmultip = NULL;
    //  numbers of selected nodes
    nsnmas = NULL;
    //  coarse numbers of selected nodes on the master
    cnsnmas = NULL;
    //  global numbers of selected nodes on the master
    gnsnmas = NULL;
    //  number of multiplicity of selected boundary/interface nodes
    snicmultip = NULL;
    //  numbers of DOFs
    snndofmas = NULL;
    //  number of DOFs on nodes
    snndofnmas = NULL;
    //  DOFs or indicators on master
    sndofmas = NULL;
    //  global numbers of boundary/interface nodes appropriate to all coarse nodes
    gnbncn = NULL;
    //  subdomain id of interface/boundary nodes appropriate to all coarse nodes
    sid = NULL;
    //  global numbers of selected boundary/interface nodes appropriate to coarse node
    sngnbncn = NULL;
    //  subdomain id of selected interface/boundary nodes of the coarse node
    snsid = NULL;
    //  numbers of DOFs for selected nodes
    ndofnsn = NULL;
    //  code numbers at selected nodes on master
    codensn = NULL;
    //  code numbers on master
    cnmas = NULL;
  }
  
  // ****************************************************
  //  to co je nize neni zkontrolovano
  // ****************************************************
  //  number of subdomain
  ndom=nd;
}

/**
   destructor
   
   JK, 30.7.2007
*/
selnodes::~selnodes()
{
  long i, j;

  delete [] lsng;
  delete [] ldof;
//  delete [] ldofmultip; // not used (commented out)
  delete [] lnsn;
  delete [] cnsn;
  delete [] gnsn;
  
  if (myrank==0){
    if (gnn)
    {
      for (i=0;i<nproc;i++)
        delete [] gnn[i];
      delete [] gnn;
    }

    delete [] dofmultip; 
/*  not used

    delete [] ndofdom;
    if (cnm)
    {
      for(i=0;)
        delete [] cnm[i];
      delete [] cnm;
    }
*/
    if (cndom)
    {
      for(i=0; i<nproc; i++)
        delete [] cndom[i];
      delete [] cndom;
    }

    for(i=0; i<tnsn; i++)
    {
      if (ljn)
        delete [] ljn[i];
      if (lsn)
        delete [] lsn[i];
      if (doffeti)
      {
        for (j=0; j<nodmultip[i]; j++)
          delete [] doffeti[i][j];
        delete [] doffeti[i];
      }
    }
    delete [] ljn;
    delete [] lsn;
    delete [] doffeti;
    delete [] nodmultip; // must not be deleted before here due to doffeti array
    delete [] ncndom;
    delete [] icmultip;

    if (cnsnmas)
    {
      for(i=0; i<nproc; i++)
        delete [] cnsnmas[i];
      delete [] cnsnmas;
    }

    if (gnsnmas)
    {
      for(i=0; i<nproc; i++)
        delete [] gnsnmas[i];
      delete [] gnsnmas;
    }

    delete [] snicmultip;

    delete [] snndofmas;

    if (snndofnmas)
    {
      for (i=0;i<nproc;i++)
        delete [] snndofnmas[i];
      delete [] snndofnmas;
    }

    if (sndofmas)
    {
      for (i=0;i<nproc;i++)
      {
        for (j=0;j<nsnmas[i];j++)
          delete [] sndofmas[i][j];
        delete [] sndofmas[i];
      }
      delete [] sndofmas;
    }
    delete [] nsnmas; // must not be deleted befor here due to sndofmas array

    if (gnbncn)
    {
      for (i=0;i<tnbn;i++)
        delete [] gnbncn[i];
      delete [] gnbncn;
    }
    
    if(sid)
    {
      for (i=0;i<tnbn;i++)
        delete [] sid[i];
      delete [] sid;
    }

    if (sngnbncn)
    {
      for (i=0;i<tnsn;i++)
        delete [] sngnbncn[i];
      delete [] sngnbncn;
    }

    if (snsid!=NULL)
    {
      for (i=0;i<tnsn;i++)
        delete [] snsid[i];
      delete [] snsid;
    }
    
    delete [] ndofnsn;

    if (codensn)
    {
      for (i=0;i<tnsn;i++)
        delete [] codensn[i];
      delete [] codensn;
    }

    if (cnmas)
    {
      for (i=0;i<nproc;i++)
        delete [] cnmas[i];
      delete [] cnmas;
    }
  }
}



/**
   function computes the number of all selected nodes and assembles
   the array nsnmas on the master
   maxnsn is obtained

   @param domproc - domain-processor correpondence
   @param out - output file
   
   30.7.2007, 29.7.2009, JK
*/
void selnodes::number_of_selected_nodes (long *domproc,FILE *out)
{
  long i,j,k;
  MPI_Status stat;
  
  if (myrank==0){
    //  array containing numbers of selected nodes on subdomains
    if (nsnmas!=NULL)
      delete [] nsnmas;
    nsnmas = new long [nproc];
    //  maximum number of selected nodes on subdomain
    maxnsn=nsn;
    //  master contributions
    j=domproc[0];
    nsnmas[j]=nsn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      nsnmas[j]=k;
      if (maxnsn<k)
        maxnsn=k;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnsn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (&nsn,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxnsn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n the maximum number of selected nodes on subdomain (maxnsn) %ld\n",maxnsn);
  if (myrank==0){
    fprintf (out,"\n\n the numbers of selected nodes on subdomains (nsnmas)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain %6ld     %ld",i+1,nsnmas[i]);
    }
  }
  
}



/**
   function collects node numbers on master
   
   function assembles array gnn
   
   @param tnnp - total number of nodes in problem
   
   30.7.2007, JK
*/
void selnodes::nodes_on_master (long *domproc,FILE *out)
{
  long i,j,k,max,tnnp,buffsize;
  long *buff,*nl;
  MPI_Status stat;
  
  tnnp=0;
  buffsize=maxnsn;
  buff = new long [buffsize];
  for (i=0;i<nsn;i++){
    buff[i]=lsng[i];
  }
  
  switch (md){
  case all_nodes:
  case neg_bound_nodes:{
    
    max=0;
    for (i=0;i<nsn;i++){
      if (max<lsng[i])
        max=lsng[i];
    }
    
    if (myrank==0){
      for (i=1;i<nproc;i++){
        MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        if (max<k)
          max=k;
      }
      tnnp=max;
      tnnp++;
    }
    else{
      MPI_Send (&max,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);

    if (myrank==0){
      //  group node numbers on master
      gnn = new long* [nproc];
      for (i=0;i<nproc;i++){
        gnn[i] = new long [nsnmas[i]];
      }
      
      nl=new long [tnnp];
      for (i=0;i<tnnp;i++){
        nl[i]=-1;
      }
      
      //  master contribution
      k=domproc[0];
      for (j=0;j<nsnmas[k];j++){
        nl[buff[j]]++;
        gnn[k][j]=buff[j];
      }
      
      //  slave contribution
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        k=domproc[stat.MPI_TAG];
        for (j=0;j<nsnmas[k];j++){
          nl[buff[j]]++;
          gnn[k][j]=buff[j];
        }
      }
      //  total number of selected nodes
      //  generation of new node numbers
      tnsn=0;
      for (i=0;i<tnnp;i++){
        if (nl[i]>-1){
          nl[i]=tnsn;
          tnsn++;
        }
      }

      for (i=0;i<nproc;i++){
        for (j=0;j<nsnmas[i];j++){
          gnn[i][j]=nl[gnn[i][j]];
        }
      }
      
      delete [] nl;
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    break;
  }
    
  case bound_nodes:{
    
    if (myrank==0){
      //  group node numbers on master
      gnn = new long* [nproc];
      for (i=0;i<nproc;i++){
        gnn[i] = new long [nsnmas[i]];
      }
      //  total number of selected nodes
      tnsn=0;
      //  master contribution
      k=domproc[0];
      for (j=0;j<nsnmas[k];j++){
        gnn[k][j]=buff[j];
        if (tnsn<buff[j])
          tnsn=buff[j];
      }
      //  slave contribution
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        k=domproc[stat.MPI_TAG];
        for (j=0;j<nsnmas[k];j++){
          gnn[k][j]=buff[j];
          if (tnsn<buff[j])
            tnsn=buff[j];
        }
      }
      tnsn++;
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    break;
  }
  default:{
    fprintf (stderr,"\n\n wrong type of mesh description is required in functio");
    fprintf (stderr,"\n node_on_master (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }

  if (myrank==0){
    fprintf (out,"\n\n total number of selected nodes is %ld",tnsn);
    fprintf (out,"\n\n group node numbers (gnn)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n subdomain %6ld",i+1);
      for (j=0;j<nsnmas[i];j++){
        fprintf (out,"\n %6ld   %ld",j+1,gnn[i][j]);
      }
    }
  }

  delete [] buff;
}



/**
   function computes multiplicity of selected nodes
   
   function assembles array nodmultip
   
   JK, 31.7.2007
*/
void selnodes::node_multiplicity (FILE *out)
{
  long i,j;
  
  if (myrank==0){
    if (nodmultip!=NULL)
      delete [] nodmultip;
    nodmultip = new long [tnsn];
    for (i=0;i<tnsn;i++){
      nodmultip[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nsnmas[i];j++){
        nodmultip[gnn[i][j]]++;
      }
    }
    fprintf (out,"\n\n node multiplicity (nodmultip)\n");
    for (i=0;i<tnsn;i++){
      fprintf (out,"\n selected node %6ld    %ld",i+1,nodmultip[i]);
    }
  }
}



void selnodes::dof_multiplicity (FILE *out)
{
  long i,j,k,n,nm,ndofn;
  MPI_Status stat;
  long buffsize;
  long *buff;
  
  buffsize = maxndof+1;
  buff = new long [buffsize];

  if (myrank == 0){
    if(dofmultip != NULL){
      delete []dofmultip;
    }
    
    
    dofmultip = new long[tndof];
    n = 0;
    for(i = 0; i < tnsn; i++){
      nm = nodmultip[i];
      for(j = 0; j < nm-1; j++){
        ndofn = ndofnsn[i];
        for(k = 0; k < ndofn; k++){
          if(doffeti[i][j][k] > 0){
            dofmultip[n] = nm;
            n++;
          }
        }
      }
    }
    // pro ladeni
    fprintf (out,"\n\n\n Control of array dofmultip\n\n\nMultiplicity of coarse DOF n = %ld\n\n",n); 
    for(i = 0; i < tndof; i++){
      fprintf (out,"coarse DOF %ld has multiplicity %ld\n",i+1,dofmultip[i]); 
    }
    /*
    for(i = 1; i < nproc; i++){
      
      for(j = 0; j < snndofmas[i];j++){
	//fprintf (out,"cndom = %ld \n",cndom[i][j]);
	buff[j] = dofmultip[cndom[i][j]-1];
      }
      buff[maxndof] = snndofmas[i];
      fprintf (out,"snndofmas[%ld] = %ld %ld\n",i,snndofmas[i],buff[maxndof]);
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD); 
    }
    
    ldofmultip = new long [snndofmas[0]];
    for(j = 0; j < snndofmas[0];j++){
      ldofmultip[j] = dofmultip[cndom[0][j]-1];
    }
    
    // pro ladeni
    fprintf (out,"ndof = %ld\n",snndofmas[0]);
    for(i = 0; i < snndofmas[0]; i++){
      fprintf (out,"ldofmultip %ld\n",ldofmultip[i]); 
    }
    */
  }
  delete [] buff;

  /*
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    if(ldofmultip != NULL){
      delete []dofmultip;
    }
    ldofmultip = new long [ndof];
    for(i = 0; i < ndof; i++){
      ldofmultip[i] = buff[i];
    }
       // pro ladeni
    fprintf (out,"ndof = %ld\n",ndof);
    for(i = 0; i < ndof; i++){
      fprintf (out," ldofmultipl %ld\n",ldofmultip[i]); 
    }
   
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete []buff;
  */
}



/**
   function assembles lists of connected nodes to coarse nodes
   and numbers of subdomains which contain connected nodes
   
   arrays ljn and lsn are assembled
   
   JK, 1.8.2007
*/
void selnodes::group_local_nodes (FILE *out)
{
  long i,j,k,m,n,min;
  
  if (myrank==0){
    ljn = new long* [tnsn];
    lsn = new long* [tnsn];
    for (i=0;i<tnsn;i++){
      ljn[i] = new long [nodmultip[i]];
      lsn[i] = new long [nodmultip[i]];
      nodmultip[i]=0;
    }
    
    //  list of local node numbers of connected nodes
    //  list of numbers of subdomains which contains connected nodes
    for (i=0;i<nproc;i++){
      for (j=0;j<nsnmas[i];j++){
        k=gnn[i][j];
        ljn[k][nodmultip[k]]=j;
        lsn[k][nodmultip[k]]=i;
        nodmultip[k]++;
      }
    }
    
    //  node sorting
    for (i=0;i<tnsn;i++){
      for (j=0;j<nodmultip[i];j++){
        min=nproc;
        for (k=j;k<nodmultip[i];k++){
          if (lsn[i][k]<min){
            min=lsn[i][k];
            m=k;
          }
        }
        n=lsn[i][j];
        lsn[i][j]=lsn[i][m];
        lsn[i][m]=n;
        n=ljn[i][j];
        ljn[i][j]=ljn[i][m];
        ljn[i][m]=n;
      }
    }
  }
  

  if (myrank==0){
    fprintf (out,"\n\n list of local node numbers of connected nodes (ljn)\n");
    fprintf (out,"\n\n list of numbers of subdomains which contains connected nodes (lsn)\n");
    for (i=0;i<tnsn;i++){
      fprintf (out,"\n coarse node %6ld   multip %6ld   ",i+1,nodmultip[i]);
      for (j=0;j<nodmultip[i];j++){
        fprintf (out,"    %ld %ld",ljn[i][j],lsn[i][j]);
      }
    }
  }
}



/**
   function assembles indicators of code numbers and then generates
   code numbers
   
   array doffeti is assembled

   JK, 1.8.2007
*/
void selnodes::dof_feti (FILE *out)
{
  long i,j,k,l,m;
  
  if (myrank==0){
    //  determination of numbers of DOFs at selected nodes
    if (ndofnsn!=NULL)
      delete [] ndofnsn;
    ndofnsn = new long [tnsn];
    for (i=0;i<tnsn;i++){
      ndofnsn[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nsnmas[i];j++){
        k=snndofnmas[i][j];
        l=gnn[i][j];
        if (ndofnsn[l]>0){
          if (ndofnsn[l]!=k){
            fprintf (stderr,"\n\n incompatible number of DOFs at selected node number %ld (file %s, line %d)\n",l,__FILE__,__LINE__);
          }
        }
        else{
          ndofnsn[l]=k;
        }
      }
    }
    
    //fprintf (out,"\n\n kontrola ndofnsn");
    //for (i=0;i<tnsn;i++){
    //fprintf (out,"\n %ld   %ld",i+1,ndofnsn[i]);
    //}

    //  array allocation
    doffeti = new long** [tnsn];
    for (i=0;i<tnsn;i++){
      doffeti[i] = new long* [nodmultip[i]];
      for (j=0;j<nodmultip[i];j++){
        doffeti[i][j] = new long [ndofnsn[i]];
      }
    }
    
    //  assembling of code numbers indicators
    for (i=0;i<tnsn;i++){
      for (j=0;j<nodmultip[i];j++){
        l=ljn[i][j];
        m=lsn[i][j];
        for (k=0;k<ndofnsn[i];k++){
          doffeti[i][j][k]=sndofmas[m][l][k];
        }
      }
    }
    
    //  code numbers generation
    tndof=1;
    for (i=0;i<tnsn;i++){
      for (j=0;j<nodmultip[i]-1;j++){
        for (k=0;k<ndofnsn[i];k++){
          if (doffeti[i][j][k]>0){
            doffeti[i][j][k]=tndof;
            tndof++;
          }
        }
      }
    }
    tndof--;
  }

  if (myrank==0){
    fprintf (out,"\n\n code numbers for FETI method (doffeti)\n");
    for (i=0;i<tnsn;i++){
      fprintf (out,"\n coarse node %6ld",i+1);
      for (j=0;j<nodmultip[i]-1;j++){
        fprintf (out,"  con.n. %6ld  sub.n. %6ld   ",ljn[i][j],lsn[i][j]);
        for (k=0;k<ndofnsn[i];k++){
          fprintf (out," %ld",doffeti[i][j][k]);
        }
      }
    }
  }
}



/**
   function determines number of contributions to coarse problem
   from particular subdomains
   
   maxndof is rewritten
   snndofmas is rewritten
   
   array ncn is assembled
   
   JK, 1.8.2007
*/
void selnodes::number_contrib (long *domproc,FILE *out)
{
  long i,j,k,l,m,ii,buffsize;
  long *buff;
  MPI_Status stat;
  
  buffsize=3;
  buff = new long [buffsize];
  
  if (myrank==0){
    //  allocation of array
    if (snndofmas!=NULL)
      delete [] snndofmas;
    snndofmas = new long [nproc];
    for (i=0;i<nproc;i++){
      snndofmas[i]=0;
    }
    
    if (ncndom!=NULL)
      delete [] ncndom;
    ncndom = new long [nproc];
    for (i=0;i<nproc;i++){
      ncndom[i]=0;
    }
     
    //  number of contributing DOFs from subdomains to coarse problem
    for (i=0;i<nproc;i++){
      for (j=0;j<nsnmas[i];j++){
        ii=gnn[i][j];
        for (k=0;k<nodmultip[ii];k++){
          if (lsn[ii][k]==i){
            if (k==0){
              for (l=0;l<nodmultip[ii]-1;l++){
                ncndom[i]++;
                for (m=0;m<ndofnsn[ii];m++){
                  if (doffeti[ii][l][m]>0)
                    snndofmas[i]++;
                }
              }
            }
            else{
              ncndom[i]++;
              for (m=0;m<ndofnsn[ii];m++){
                if (doffeti[ii][k-1][m]>0)
                  snndofmas[i]++;
              }
            }
            break;
          }
        }
      }
    }
    
    //  maximum number of contributions from subdomain
    maxndof=0;
    //  maximum number of contributing nodes from one subdomain
    maxncn=0;
    for (i=0;i<nproc;i++){
      if (maxndof<snndofmas[i])
        maxndof=snndofmas[i];
      if (maxncn<ncndom[i])
        maxncn=ncndom[i];
    }
    
    buff[0]=maxndof;
    buff[1]=maxncn;
    
    for (i=1;i<nproc;i++){
      j=domproc[i];
      buff[2]=ncndom[j];
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    j=domproc[0];
    buff[2]=ncndom[j];
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  maxndof=buff[0];
  maxncn=buff[1];
  ncn=buff[2];

  if (myrank==0){
    fprintf (out,"\n\n number of contributing DOFs to coarse problem (snndofmas)\n");
    fprintf (out,"\n maximum number of DOFs on subdomain %ld",maxndof);
    fprintf (out,"\n maximum number of contributing nodes on subdomain %ld",maxncn);
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %6ld   %ld",i+1,snndofmas[i]);
    }
  }
  fprintf (out,"\n\n number of nodes which contribute to coarse FETI problem (ncn) %ld",ncn);
  delete [] buff;
}



/**
   function assmebles nodes contributing to the coarse proble in the FETI method
      
   JK, 1.8.2007
*/
void selnodes::contrib_dofs (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,l,m,ii,jj,ndofn,buffsize,min;
  long *buff;
  MPI_Status stat;
  
  buffsize=maxncn;
  buff = new long [buffsize];
  
  if (myrank==0){
    
    //  array of local numbers of nodes contributing to the coarse problem
    for (i=1;i<nproc;i++){
      
      for (j=0;j<buffsize;j++){
        buff[j]=0;
      }
      
      m=0;
      for (j=0;j<tnsn;j++){
        for (k=0;k<nodmultip[j];k++){
          if (lsn[j][k]==domproc[i]){
            if (k==0){
              for (l=0;l<nodmultip[j]-1;l++){
                buff[m]=ljn[j][k];
                buff[m]=0-buff[m]-1;
                m++;
              }
            }
            else{
              buff[m]=ljn[j][k];
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
    for (j=0;j<tnsn;j++){
      for (k=0;k<nodmultip[j];k++){
        if (lsn[j][k]==domproc[i]){
          if (k==0){
            for (l=0;l<nodmultip[j]-1;l++){
              buff[m]=ljn[j][k];
              buff[m]=0-buff[m]-1;
              m++;
            }
          }
          else{
            buff[m]=ljn[j][k];
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
  

  fprintf (out,"\n\n kontrola bufferu na vsech uzlech");
  for (i=0;i<ncn;i++){
    fprintf (out,"\n %ld   %ld",i+1,buff[i]);
  }
  for (i=0;i<ncn;i++){
    if (buff[i]<0)
      min=0-buff[i]-1;
    else
      min=buff[i];
      
    l=i;
    for (j=i+1;j<ncn;j++){
      if (buff[j]<0)
        k=0-buff[j]-1;
      else
        k=buff[j];
      if (min>k){
        l=j;
        min=k;
      }
    }
    k=buff[i];
    buff[i]=buff[l];
    buff[l]=k;
  }

  fprintf (out,"\n\n kontrola bufferu na vsech uzlech");
  for (i=0;i<ncn;i++){
    fprintf (out,"\n %ld   %ld",i+1,buff[i]);
  }
  fprintf (out,"\n konec kontroly\n");
  /*
  fprintf (out,"\n\n kontrola lnsn na vsech uzlech");
  for (i=0;i<nsn;i++){
    fprintf (out,"\n %ld   %ld",i+1,lnsn[i]);
  }
  fprintf (out,"\n konec kontroly\n");
  */
  
  
  //  number of contributions to the coarse problem
  //  number of DOFs which contribute to the coarse problem
  ndof=0;
  for (i=0;i<ncn;i++){
    j=buff[i];
    if (j<0){
      j=0-j-1;
    }
    l=lnsn[j];
    ndofn=top->give_ndofn (l);
    for (k=0;k<ndofn;k++){
      if (top->give_dof (l,k)>0)
        ndof++;
    }
  }
  
  //  list of code numbers which are extracted from domain in the FETI method
  //  some code numbers are positive and some of them are negative
  //  it depends on number of subdomain
  //  these code numbers are used on subdomains, master processor contains corresponding array with coarse code numbers

  if (ldof!=NULL)
    delete [] ldof;
  ldof = new long [ndof];
  
  m=0;
  for (i=0;i<ncn;i++){
    ii=buff[i];
    if (ii<0){
      ii=0-ii-1;
      jj=lnsn[ii];
      ndofn=top->give_ndofn (jj);
      for (j=0;j<ndofn;j++){
        l=top->give_dof (jj,j);
        if (l>0){
          ldof[m]=0-l;
          m++;
        }
      }
    }
    else{
      jj=lnsn[ii];
      ndofn=top->give_ndofn (jj);
      for (j=0;j<ndofn;j++){
        l=top->give_dof (jj,j);
        if (l>0){
          ldof[m]=l;
          m++;
        }
      }
    }
  }
  delete [] buff;
  
  if (myrank==0){
    cndom = new long* [nproc];
    for (i=0;i<nproc;i++){
      cndom[i] = new long [snndofmas[i]];
    }
    
    for (i=0;i<nproc;i++){
      l=0;
      for (j=0;j<nsnmas[i];j++){
        m=gnn[i][j];
        for (k=0;k<nodmultip[m];k++){
          if (lsn[m][k]==i){
            if (k==0){
              for (ii=0;ii<nodmultip[m]-1;ii++){
                for (jj=0;jj<ndofnsn[m];jj++){
                  if (doffeti[m][ii][jj]>0){
                    cndom[i][l]=doffeti[m][ii][jj];
                    l++;
                  }
                }
              }
            }
            else{
              for (jj=0;jj<ndofnsn[m];jj++){
                if (doffeti[m][k-1][jj]>0){
                  cndom[i][l]=doffeti[m][k-1][jj];
                  l++;
                }
              }
            }
            break;
          }
        }
      }
    }
  }

  fprintf (out,"\n\n number of DOFs which contribute to coarse FETI problem %ld",ndof);
  fprintf (out,"\n array ldof");
  for (i=0;i<ndof;i++){
    fprintf (out,"\n %6ld   %ld",i+1,ldof[i]);
  }
  
  if (myrank==0){
    fprintf (out,"\n\n code numbers on master (cndom)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %6ld (number of DOFs %ld)",i+1,snndofmas[i]);
      for (j=0;j<snndofmas[i];j++){
        fprintf (out,"\n %6ld   %ld",j+1,cndom[i][j]);
      }
    }
  }
}




































































































/*
  ZACATEK ROZSAHLYCH UPRAV
  27.7.2009
*/

/**
   constructor
   
   @param nd - the number of subdomains/processors
   @param mr - my rank
   @param d - type of mesh description
   @param domproc - domain-processor correspondence
   @param nnsd - array containing the numbers of all nodes on subdomains
   @param jj - array containing list of selected nodes
              it contains nd components
   @param itnbn - total number of boundary/interface nodes
   @param iicmultip - array containing node multiplicity of boundary/interface nodes
   @param anodmultip - array containing node multiplicity of nodes on subdomain
   @param icnbn - array containing coarse numbers of boundary/interface nodes on one subdomain
   @param gnbndom - array containing global numbers of boundary/interface nodes on one subdomain
   @param out - output stream
   @param proc_name - processor name
   @param mespr - message printing indicator
   
   jj[i][k]>-1 - the k-th node on the i-th subdomain is selected
   
   JK, 29.7.2009
*/
selnodes::selnodes(int nd,int mr,meshdescription d,long *domproc,long *nnsd,long *jj,
                   long itnbn,long *iicmultip,long *anodmultip,long *icnbn,long *gnbndom,
                   long **ignbncn,long **isid,
                   FILE *out,char *proc_name,long /*mespr*/)
{
  long i,j,k,kk;
  long *buff,buffsize;
  MPI_Status stat;

  //  the number of processors = the number of subdomains
  nproc=nd;
  //  my rank
  myrank=mr;
  //  mesh description
  md = d;

  //  the number of all nodes on subdomain
  if (myrank==0){
    for (i=1;i<nproc;i++){
      nn=nnsd[domproc[i]];
      MPI_Send (&nn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    nn=nnsd[domproc[0]];
  }
  else{
    MPI_Recv (&nn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);


  //  the number of selected nodes on subdomain
  nsn=0;
  for (i=0;i<nn;i++){
    if (jj[i]>-1)
      nsn++;
  }
  
  if (nsn==0){
    par_print_err(myrank,proc_name,"wrong number of selected nodes in constructor",__FILE__, __LINE__, __func__);
    abort();
  }
  
  // *****************************
  //  assembling of array nsnmas
  // *****************************
  if (myrank==0){
    //  array containing numbers of selected nodes on subdomains
    nsnmas = new long [nproc];
    
    //  maximum number of selected nodes on subdomain
    maxnsn=nsn;
    
    //  master contributions
    j=domproc[0];
    nsnmas[j]=nsn;
    
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      nsnmas[j]=k;
      if (maxnsn<k)
        maxnsn=k;
    }
    
    for (i=1;i<nproc;i++){
      MPI_Send (&maxnsn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Send (&nsn,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxnsn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  // ************************************
  //  end of assembling of array nsnmas
  // ************************************
  

  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n the maximum number of selected nodes on subdomain (maxnsn) %ld\n",maxnsn);
  if (myrank==0){
    fprintf (out,"\n\n the numbers of selected nodes on subdomains (nsnmas)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain %6ld     %ld",i,nsnmas[i]);
    }
  }
  
  // *********************************************
  //  local and coarse numbers of selected nodes
  //  these numbers are always known
  // *********************************************
  buffsize=maxnsn;
  buff = new long [buffsize];
  
  //  local numbers of selected nodes
  lnsn = new long [nsn];
  //  coarse numbers of selected nodes on subdomain
  cnsn = new long [nsn];
  
  //  index
  k=0;
  //  index
  kk=0;
  //  loop over all nodes on subdomain
  for (i=0;i<nn;i++){
    if (anodmultip[i]>1){
      //  boundary/interface node
      if (jj[i]>-1){
        //  the node is selected
        //  local numbers of selected nodes
        lnsn[k]=i;
        //  coarse numbers of selected nodes
        cnsn[k]=icnbn[kk];
        buff[k]=cnsn[k];
        k++;
      }
      kk++;
    }
    else{
      //  internal node
      if (jj[i]>-1){
        //  the node is selected
        par_print_err(myrank,proc_name,"internal node has been selected",__FILE__, __LINE__, __func__);
      }
    }
  }
  
  
  if (myrank==0){
    //  coarse numbers of selected nodes
    cnsnmas = new long* [nproc];
    for (i=0;i<nproc;i++){
      cnsnmas[i]=new long [nsnmas[i]];
    }
    
    //  subdomain id
    k=domproc[0];
    //  loop over the number of nodes
    for (j=0;j<nsnmas[k];j++){
      cnsnmas[k][j]=buff[j];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //  slave contributions
      //  subdomain id
      k=domproc[stat.MPI_TAG];
      //  loop over the number of nodes
      for (j=0;j<nsnmas[k];j++){
        cnsnmas[k][j]=buff[j];
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // ****************************************
  //  global numbers of selected nodes
  //  these numbers are known only for
  //  mesh description denoted by all_nodes
  // ****************************************
  
  if (md==all_nodes){
    //  global numbers of selected nodes
    gnsn = new long [nsn];
    //  index
    k=0;
    //  index
    kk=0;
    //  loop over all nodes on subdomain
    for (i=0;i<nn;i++){
      if (anodmultip[i]>1){
        //  boundary/interface node
        if (jj[i]>-1){
          //  the node is selected
          //  global numbers of selected nodes
          gnsn[k]=gnbndom[kk];
          buff[k]=gnsn[k];
          k++;
        }
        kk++;
      }
      else{
        //  internal node
        if (jj[i]>-1){
          //  the node is selected
          par_print_err(myrank,proc_name,"internal node has been selected",__FILE__, __LINE__, __func__);
        }
      }
    }
    if (myrank==0){
      //  global numbers of selected nodes
      gnsnmas = new long* [nproc];
      for (i=0;i<nproc;i++){
        gnsnmas[i]=new long [nsnmas[i]];
      }
      //  subdomain id
      k=domproc[0];
      //  loop over the number of nodes
      for (j=0;j<nsnmas[k];j++){
        gnsnmas[k][j]=buff[j];
      }
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        //  slave contributions
        //  subdomain id
        k=domproc[stat.MPI_TAG];
        //  loop over the number of nodes
        for (j=0;j<nsnmas[k];j++){
          gnsnmas[k][j]=buff[j];
        }
      }
    }
    else{
      MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
  }
  if (myrank==0){
    //  total number of boundary/interafce nodes
    tnbn=itnbn;
    //  number of multiplicity of boundary/interface nodes
    icmultip = new long [tnbn];
    for (i=0;i<tnbn;i++){
      icmultip[i]=iicmultip[i];
    }
  }
  
  if (md==all_nodes){
    if (myrank==0){
      //  global numbers of boundary/interface nodes appropriate to coarse node
      gnbncn = new long* [tnbn];
      for (i=0;i<tnbn;i++){
        gnbncn[i]=new long [icmultip[i]];
      }
      //  subdomain id of interface/boundary nodes of the coarse node
      sid = new long* [tnbn];
      for (i=0;i<tnbn;i++){
        sid[i] = new long [icmultip[i]];
      }
      for (i=0;i<tnbn;i++){
        for (j=0;j<icmultip[i];j++){
          gnbncn[i][j]=ignbncn[i][j];
          sid[i][j]=isid[i][j];
        }
      }
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete [] buff;
  //  the total number of selected nodes
  tnsn=0;
  //  the number of all selected DOFs including prescribed DOFs on one subdomain
  snndof=0;
  //  the maximum number of selected DOFs on subdomain
  maxsnndof=0;
  //  total number of DOFs on selected nodes
  tndofsn=0;

  lsng = NULL;
  ldof = NULL;

  ldofmultip = NULL;  // not used
  //  numbers of DOFs on subdomains
  ndofdom = NULL;  // not used
  //  code numbers at selected nodes on master
  cnm = NULL;  // not used

  if (myrank==0){
    //  group node numbers
    gnn = NULL;
    //  node multiplicty
    nodmultip = NULL;
    // dof multiplicity
    dofmultip = NULL; 
    //  code numbers on master
    cndom = NULL;
    //  joint nodes to selected nodes assumed as coarse nodes
    ljn = NULL;
    //  subdomain numbers which contain connected nodes to coarse nodes
    lsn = NULL;
    //  code numbers / indicators for FETI method
    doffeti = NULL;
    //  number of contributing nodes in the FETI method
    ncndom = NULL;

    //  number of multiplicity of selected boundary/interface nodes
    snicmultip = NULL;
    //  numbers of DOFs
    snndofmas = NULL;
    //  number of DOFs on nodes
    snndofnmas = NULL;
    //  DOFs or indicators on master
    sndofmas = NULL;
    //  global numbers of selected boundary/interface nodes appropriate to coarse node
    sngnbncn = NULL;
    //  subdomain id of selected interface/boundary nodes of the coarse node
    snsid = NULL;
    //  numbers of DOFs for selected nodes
    ndofnsn = NULL;
    //  code numbers at selected nodes on master
    codensn = NULL;
    //  code numbers on master
    cnmas = NULL;
  }
}



/**
   function collects coarse numbers of nodes on the master
   
   function assembles arrays:

   snicmultip - number of multiplicity of selected boundary/interface nodes
   
   @param out - output file (for auxiliary output)
   
   29.7.2009, JK
*/
void selnodes::node_coarse_numbers (FILE *out)
{
  long i,j,k;
  long *av;

  if (myrank==0){
    //  auxiliary array
    av = new long [tnbn];
    for (i=0;i<tnbn;i++){
      av[i]=0;
    }
    
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      //  loop over the number of selected nodes
      for (j=0;j<nsnmas[i];j++){
        av[cnsnmas[i][j]]++;
      }
    }
    
    //  the total number of selected nodes
    tnsn=0;
    for (i=0;i<tnbn;i++){
      if (av[i]>0)
      tnsn++;
    }
    
    //  number of multiplicity of selected boundary/interface nodes
    if (snicmultip!=NULL){
      delete [] snicmultip;
    }
    snicmultip = new long [tnsn];
    
    k=0;
    for (i=0;i<tnbn;i++){
      if (av[i]>0){
        snicmultip[k]=icmultip[i];
        k++;
      }
    }
  }
  
  if (myrank==0){
    if (md==all_nodes){
      
      //  global numbers of boundary/interface nodes appropriate to selected coarse node
      if (sngnbncn!=NULL){
        for (i=0;i<tnsn;i++){
          delete [] sngnbncn[i];
        }
        delete [] sngnbncn;
      }
      sngnbncn = new long* [tnsn];
      for (i=0;i<tnsn;i++){
        sngnbncn[i]=new long [snicmultip[i]];
      }
      
      //  subdomain id of interface/boundary nodes appropriate to coarse node
      if (snsid!=NULL){
        for (i=0;i<tnsn;i++){
          delete [] snsid[i];
        }
        delete [] snsid;
      }
      snsid = new long* [tnsn];
      for (i=0;i<tnsn;i++){
        snsid[i] = new long [snicmultip[i]];
      }
      
      //  loop over the number of all boundary/interface nodes
      k=0;
      for (i=0;i<tnbn;i++){
        if (av[i]>0){
          for (j=0;j<snicmultip[k];j++){
            sngnbncn[k][j]=gnbncn[i][j];
            snsid[k][j]=sid[i][j];
          }
          k++;
        }
      }
      
      for (i=0;i<tnbn;i++){
        delete [] gnbncn[i];
      }
      delete [] gnbncn;
      gnbncn = NULL;
      
      for (i=0;i<tnbn;i++){
        delete [] sid[i];
      }
      delete [] sid;
      sid = NULL;
      
    }
    delete [] icmultip;
    icmultip = NULL;
    delete [] av;
  }
  fprintf (out,"\n\n the total number of selected nodes (tnsn)  %ld",tnsn);
}



/**
   function searches for all possible DOFs in selected nodes
   prescribed DOFs are included

   snndof  is determined
   maxsnndof is determined
   array snndofmas is assembled
   
   @param top - pointer to the general topology
   @param out - output file (for auxiliary output)
   
   JK, 30.7.2009
*/
void selnodes::number_all_dofs (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k;
  MPI_Status stat;
  
  //  the number of DOFs in selected nodes
  snndof=0;
  //  loop over the number of selected nodes
  for (i=0;i<nsn;i++){
    j=lnsn[i];
    snndof+=top->give_ndofn (j);
  }
  
  if (myrank==0){
    //  number of all dofs on subdomains
    if (snndofmas!=NULL)
      delete [] snndofmas;
    snndofmas=new long [nproc];
    
    //  the maximum number of selected DOFs on subdomains
    maxsnndof=snndof;
    
    //  master contributions
    j=domproc[0];
    snndofmas[j]=snndof;

    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      j=domproc[stat.MPI_TAG];
      snndofmas[j]=k;
      if (maxsnndof<k)
        maxsnndof=k;
    }

    for (i=1;i<nproc;i++){
      MPI_Send (&maxsnndof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
    
  }
  else{
    MPI_Send (&snndof,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&maxsnndof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n the maximum number of selected DOFs on subdomains (maxsnndof)   %ld\n",maxsnndof);
    fprintf (out,"\n\n numbers of DOFs on subdomains (snndofmas)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %6ld   %ld",i,snndofmas[i]);
    }
  }
}



/**
   function collects numbers of DOFs on selected nodes
   function allocates and assembles the array snndofnmas
   
   the following arrays are assembled:
   
   snndofnmas - array of numbers of DOFs in selected nodes on the master

   @param top - pointer to the genral topology
   @param domproc - correspondence between processors and subdomains
   @param out - output file (for auxiliary output)

   JK, 30.7.2009
*/
void selnodes::ndofn_on_master (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;

  buffsize=maxsnndof;
  buff = new long [buffsize];
  //  loop over the number of selected nodes
  for (i=0;i<nsn;i++){
    //  node id
    j=lnsn[i];
    buff[i]=top->give_ndofn (j);
  }
  
  if (myrank==0){
    //  the number of DOFs in selected nodes
    snndofnmas = new long* [nproc];
    for (i=0;i<nproc;i++){
      snndofnmas[i] = new long [nsnmas[i]];
    }
    
    //  master contribution
    k=domproc[0];
    for (j=0;j<nsnmas[k];j++){
      snndofnmas[k][j]=buff[j];
    }
    
    //  slave contribution
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      for (j=0;j<nsnmas[k];j++){
        snndofnmas[k][j]=buff[j];
      }
    }
    
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  delete [] buff;

  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n\n numbers of DOFs on selected nodes (snndofnmas)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %6ld",i);
      for (j=0;j<nsnmas[i];j++){
        fprintf (out,"\n selected node %6ld   %ld",j,snndofnmas[i][j]);
      }
    }
  }
}



/**
   function assembles DOFs indicators on master
   
   array sndofmas is assembled
   
   @param top - pointer to the genral topology
   @param domproc - correspondence between processors and subdomains
   @param out - output file (for auxiliary output)
   
   JK, 30.7.2009
*/
void selnodes::dof_indicators (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,l,m,ndofn,buffsize;
  long *buff;
  MPI_Status stat;

  buffsize=maxsnndof;
  buff = new long [buffsize];

  //  index
  l=0;
  //  loop over the number of selected nodes
  for (i=0;i<nsn;i++){
    //  node id
    j=lnsn[i];
    //  the number of DOFs in selected node
    ndofn=top->give_ndofn (j);
    //  loop over the number of DOFs in selected node
    for (k=0;k<ndofn;k++){
      buff[l]=top->give_dof (j,k);
      l++;
    }
  }
  
  if (myrank==0){
    //  array of DOFs or indicators at selected nodes
    if (sndofmas!=NULL){
      for (i=0;i<nproc;i++){
        for (j=0;j<nsnmas[i];j++){
          delete [] sndofmas[i][j];
        }
        delete [] sndofmas[i];
      }
      delete [] sndofmas;
    }
    
    sndofmas = new long** [nproc];
    for (i=0;i<nproc;i++){
      sndofmas[i] = new long* [nsnmas[i]];
      for (j=0;j<nsnmas[i];j++){
        sndofmas[i][j] = new long [snndofnmas[i][j]];
      }
    }
    
    //  master contribution
    k=domproc[0];
    m=0;
    for (j=0;j<nsnmas[k];j++){
      for (l=0;l<snndofnmas[k][j];l++){
        sndofmas[k][j][l]=buff[m];
        m++;
      }
    }
    
    //  slave contribution
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      m=0;
      for (j=0;j<nsnmas[k];j++){
        for (l=0;l<snndofnmas[k][j];l++){
          sndofmas[k][j][l]=buff[m];
          m++;
        }
      }
    }
  }
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  delete [] buff;
  
  printf ("\n selnodes %d %p ",myrank,(void*)out);
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n\n DOFs indicators on master (sndofmas)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %6ld",i);
      for (j=0;j<nsnmas[i];j++){
        fprintf (out,"\n selected node %6ld    ",j);
        for (k=0;k<snndofnmas[i][j];k++){
          fprintf (out,"  %ld",sndofmas[i][j][k]);
        }
      }
    }
  }
}



/**
   function generates ordering of selected nodes with respect
   to the Schur complement method
   
   the function generates code numbers / numbers of unknowns
   on the coarse problem
   
   array ldof is assembled

   @param top - pointer to the general topology
   @param out - output file (for auxiliary output)
   
   JK, 31.7.2007
*/
void selnodes::schur_ordering (gtopology */*top*/,long **dofind,FILE *out)
{
  long i,j,k,l,m,g,ndofn;
  long *av,*aux;
  
  if (myrank==0){
    
    //  array of code numbers of coarse nodes
    if (codensn!=NULL){
      for (i=0;i<tnsn;i++){
        delete [] codensn[i];
      }
      delete [] codensn;
    }
    codensn = new long* [tnsn];
    for (i=0;i<tnsn;i++){
      codensn[i] = NULL;
    }
    
    //  array of numbers of DOFs for selected nodes
    if (ndofnsn!=NULL){
      delete [] ndofnsn;
    }
    ndofnsn = new long [tnsn];

    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      //  loop over the number of selected nodes
      for (j=0;j<nsnmas[i];j++){
        //  global glued number of selected node
        g=gnsnmas[i][j];
        //  coarse number of the selected node
        l=cnsnmas[i][j];
        //  number of DOFs on the selected node
        ndofn=snndofnmas[i][j];
        if (codensn[l]==NULL){
          codensn[l] = new long [ndofn];
          ndofnsn[l]=ndofn;
        }
        //  loop over the number of DOFs on node
        for (k=0;k<ndofn;k++){
          if (dofind[g][k]==1){
            //  this DOF is boundary/interface DOF
            codensn[l][k]=sndofmas[i][j][k];
          }
          else{
            //  this DOF is internal DOF
            codensn[l][k]=0;
          }
        }
      }
    }
    // ***************************************************
    //  generation of code numbers in the coarse problem
    // ***************************************************
    
    //  searching of maximum code number
    tndofsn=0;
    //  loop over the number of selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the number of DOFs on selected node
      for (j=0;j<ndofnsn[i];j++){
        if (tndofsn<codensn[i][j])
          tndofsn=codensn[i][j];
      }
    }
    tndofsn--;
    if (tndofsn<0)  tndofsn=0;
    aux = new long [tndofsn];
    for (i=0;i<tndofsn;i++){
      aux[i]=-1;
    }
    //  total number of DOFs on selected nodes
    tndofsn=1;
    //  loop over the number of selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the number of DOFs on selected node
      for (j=0;j<ndofnsn[i];j++){
        k=codensn[i][j];
        if (k==1){
          codensn[i][j]=tndofsn;
          tndofsn++;
        }
        if (k>1){
          if (aux[k-2]==-1){
            codensn[i][j]=tndofsn;
            aux[k-2]=tndofsn;
            tndofsn++;
          }
          else{
            codensn[i][j]=aux[k-2];
          }
        }
      }
    }
    tndofsn--;
    delete [] aux;
    // **********************************************************
    //  end of generation of code numbers in the coarse problem
    // **********************************************************
    aux = new long [tndofsn];
    //  this array will be recalculated
    //  at this moment, it contains number of all
    for (i=0;i<nproc;i++){
      snndofmas[i]=0;
    }
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      for (j=0;j<tndofsn;j++){
        aux[j]=0;
      }
      //  loop over the number of selected nodes
      for (j=0;j<nsnmas[i];j++){
        //  coarse number of the selected node
        l=cnsnmas[i][j];
        //  number of DOFs on the selected node
        ndofn=snndofnmas[i][j];
        //  loop over the number of DOFs on node
        for (k=0;k<ndofn;k++){
          if (codensn[l][k]>0){
            if (aux[codensn[l][k]-1]==0){
              snndofmas[i]++;
              aux[codensn[l][k]-1]=1;
            }
          }
        }
      }
    }
    //  code numbers on master
    if (cnmas!=NULL){
      for (i=0;i<nproc;i++){
        delete [] cnmas[i];
      }
      delete [] cnmas;
    }
    cnmas = new long* [nproc];
    for (i=0;i<nproc;i++){
      cnmas[i] = new long [snndofmas[i]];
    }
    //  auxiliary array
    av = new long [nproc];
    for (i=0;i<nproc;i++){
      av[i]=0;
    }
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      for (j=0;j<tndofsn;j++){
        aux[j]=0;
      }
      //  loop over the number of selected nodes
      for (j=0;j<nsnmas[i];j++){
        //  coarse number of the selected node
        l=cnsnmas[i][j];
        //  number of DOFs on the selected node
        ndofn=snndofnmas[i][j];
        //  loop over the number of DOFs on node
        for (k=0;k<ndofn;k++){
          m=codensn[l][k];
          if (m>0){
            if (aux[m-1]==0){
              cnmas[i][av[i]]=m;
              av[i]++;
              aux[m-1]=1;
            }
          }
        }
      }
    }
    delete [] av;
    delete [] aux;
  }
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n\n\n tndofsn - total number of DOFs in selected nodes  %ld",tndofsn);
    
    fprintf (out,"\n\n code numbers of Schur ordering (cnmas)\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n subdomain number %6ld",i);
      for (j=0;j<snndofmas[i];j++){
        fprintf (out,"\n dof %6ld   %ld",j,cnmas[i][j]);
      }
    }
  }
}



/**
   function generates ordering of selected nodes with respect
   to the Schur complement method
   
   the function generates code numbers / numbers of unknowns
   on the coarse problem
   
   array ldof is assembled

   @param top - pointer to the general topology
   @param out - output file (for auxiliary output)
   
   JK, 31.7.2007
*/
void selnodes::schur_ordering (gtopology */*top*/,FILE *out)
{
  long i,j,k,l,m,ndofn;
                   
  if (myrank==0){
    //  determination of numbers of DOFs at selected nodes
    if (ndofnsn!=NULL)
      delete [] ndofnsn;
    ndofnsn = new long [tnsn];
    for (i=0;i<tnsn;i++){
      ndofnsn[i]=0;
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nsnmas[i];j++){
        k=snndofnmas[i][j];
        l=cnsnmas[i][j];
        if (ndofnsn[l]>0){
          if (ndofnsn[l]!=k){
            fprintf (stderr,"\n\n incompatible number of DOFs at selected node number %ld (file %s, line %d)\n",l,__FILE__,__LINE__);
          }
        }
        else{
          ndofnsn[l]=k;
        }
      }
    }
    
    //  code number indicators
    codensn = new long* [tnsn];
    for (i=0;i<tnsn;i++){
      codensn[i] = new long [ndofnsn[i]];
    }
    
    for (i=0;i<nproc;i++){
      for (j=0;j<nsnmas[i];j++){
        l=cnsnmas[i][j];
        for (k=0;k<ndofnsn[l];k++){
          codensn[l][k]=sndofmas[i][j][k];
        }
      }
    }
    

    //  code number generation
    tndofsn=1;
    for (i=0;i<tnsn;i++){
      for (j=0;j<ndofnsn[i];j++){
        if (codensn[i][j]>0){
          codensn[i][j]=tndofsn;
          tndofsn++;
        }
      }
    }
    tndofsn--;
    
    
    //  computation of real number of DOFs on subdomains
    //  supports are not included now
    for (i=0;i<nproc;i++){
      snndofmas[i]=0;
      for (j=0;j<nsnmas[i];j++){
        l=cnsnmas[i][j];
        for (k=0;k<ndofnsn[l];k++){
          if (codensn[l][k]>0)
            snndofmas[i]++;
        }
      }
    }
    
    cnmas = new long *[nproc];
    for (i=0;i<nproc;i++){
      cnmas[i] = new long [snndofmas[i]];
    }
    
    for (i=0;i<nproc;i++){
      m=0;
      for (j=0;j<nsnmas[i];j++){
        l=cnsnmas[i][j];
        for (k=0;k<ndofnsn[l];k++){
          if (codensn[l][k]>0){
            cnmas[i][m]=codensn[l][k];
            m++;
          }
        }
      }
    }
  }
  
  //  zde by se dalo smazat pole codensn
  if (myrank==0){
    fprintf (out,"\n\n\n tndofsn  %ld \n",tndofsn);
    fprintf (out,"\n\n\n array cnmas \n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n");
      fprintf (out,"\n snndofmas   %6ld",snndofmas[i]);
      for (j=0;j<snndofmas[i];j++){
        fprintf (out,"\n cnmas %6ld    %6ld",j,cnmas[i][j]);
      }
    }
  }
  /*
  //  number of contributions to the coarse problem
  //  number of DOFs which contribute to the coarse problem
  ndof=0;
  for (i=0;i<nsn;i++){
    l=lnsn[i];
    ndofn=top->give_ndofn (l);
    for (k=0;k<ndofn;k++){
      if (top->give_dof (l,k)>0)
        ndof++;
    }
  }
  */

  //  list of code numbers which are extracted from domain in the FETI method
  //  some code numbers are positive and some of them are negative
  //  it depends on number of subdomain
  //  these code numbers are used on subdomains, master processor contains corresponding array with coarse code numbers
  
/*
  if (ldof!=NULL)
    delete [] ldof;
  ldof = new long [ndof];
  
  m=0;
  for (i=0;i<nsn;i++){
    l=lnsn[i];
    ndofn=top->give_ndofn (l);
    for (j=0;j<ndofn;j++){
      k=top->give_dof (l,j);
      if (k>0){
        ldof[m]=k;
        m++;
      }
    }
  }
*/
}
