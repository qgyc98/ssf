#include <mpi.h>
#include "partop.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include "galias.h"
#include "intp.h"

/**
   constructor
   
   @param np - number of processors
   @param mr - my rank
   @param nd - number of subdomain
   @param meshd - type of mesh description
   @param *nameproc - name of processor
   @param nl - length of processor name  
   
   JK, 3.5.2004
*/
partop::partop(int np,int mr,long nd,meshdescription meshd,char *nameproc, int nl)
{
  //  number of processors
  nproc=np;
  //  my rank
  myrank=mr;
  //  mesh description
  md = meshd;

  //  number of nodes on subdomain
  nn=0;
  //  the number of DOFs on subdomain
  ndof=0;
  //  the number of coupled DOFs
  ncdof=0;
  //  the number of boundary/interface coupled DOFs
  nbcdof=0;

  //  the maximum number of nodes on one subdomain
  maxnn=0;
  
  //  the number of internal nodes
  nin=0;
  //  the number of boundary/interface nodes on subdomain
  nbn=0;
  //  the maximum number of boundary nodes on subdomain
  maxnbn = 0;
  //  the total number of all boundary/interface nodes in the problem
  tnbn=0;
  //  the total number of nodes in the whole problem
  tnnp=0;
  
  //  maximum number of all DOFs on one subdomain
  maxndof=0;

  bltg = NULL;
  
  lnindom = NULL;
  
  lnbndom = NULL;
  
  //  array containing nodal multiplicity
  nodmultip = NULL;
  //  array containing coarse numbers of interface/boundary nodes on one subdomain
  icnbn = NULL;
  //  global numbers of boundary/interface nodes on subdomains
  gnbndom = NULL;
  //  array containing the numbers of internal nodes on subdomains
  nind = NULL;
  
  
  if (myrank==0){
    //  array containing the numbers of nodes on subdomains
    nnsd = NULL;
    //  array containing number of subdomains which each boundary/interface node belongs to
    bmultip = NULL;
    //  array containing number of subdomains which each node belongs to
    amultip = NULL;
    //  array containing the numbers of boundary/interface nodes on subdomains
    nbnd = NULL;
    //  array containing global node numbers
    allnodes = NULL;
    //  array containing coarse node numbers
    bnodes = NULL;
    //  array containing indicators of coupled DOFs
    coupdof = NULL;
    //  array containing suspicious indicators of coupled DOFs
    coupdofmas = NULL;
    //  array containing numbers of all degrees of freedom on subdomains
    nalldof = NULL;
    //  array containing the numbers of DOF on nodes
    ndofnmas = NULL;
    //  array containing DOF indicators
    dofindmas = NULL;
    //  array containing DOF indicators
    dofind = NULL;
    //  array containing coarse numbers of interface/boundary nodes
    icnbnmas = NULL;
    //  array containing node multiplicity of boundary/interface nodes
    icmultip = NULL;
    //  global numbers of boundary/interface nodes on the master
    gnbn = NULL;
    //  global numbers of internal nodes on the master
    gnin = NULL;
    //  global numbers of boundary/interface nodes appropriate to coarse node
    gnbncn = NULL;
    //  subdomain id of interface/boundary nodes appropriate to coarse node
    sid = NULL;
    
  }

  


  // ****************************************************
  //  to co je nize neni zkontrolovano
  // ****************************************************


  //  number of subdomain
  ndom=nd;
  // name of processor
  strcpy (procName,nameproc);
  nameLength=nl;
  
  
  // number of elements on subdomain
  ne=0;

  //  local numbers of boundary nodes
  lnbn = NULL;
  //  local numbers of internal nodes
  lnin = NULL;
  //
  lgnbn = NULL;
  // aray containing DOF multiplicity
  dofmultip = NULL;
    
  if (myrank==0){
    multip = NULL;
    cnbn = NULL;
    nbdofnd = NULL;
    nbdofd = NULL;
    lbcn = NULL;
    bnmultip = NULL;
    llnbn = NULL;
    ldn = NULL;
    
    
    //gnn = NULL;
    //pgcn = NULL;
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

  delete [] bltg;
  delete [] lnindom;
  delete [] lnbndom;

  delete [] nodmultip;
  delete [] icnbn;
  delete [] gnbndom;

  delete [] lnbn;
  delete [] lnin;
  delete [] lgnbn;
  delete [] dofmultip;
  delete [] nind; 

  if (myrank==0){
    delete [] nnsd;
    delete [] bmultip; 
    delete [] amultip; 
    delete [] nbnd;

    for (i=0;i<nproc;i++){
      if (allnodes)
        delete [] allnodes[i];
      if (bnodes)
        delete [] bnodes[i];
      if (coupdof)
        delete [] coupdof[i];
    }
    delete [] allnodes;
    delete [] bnodes;
    delete [] coupdof;

    delete [] coupdofmas;

    delete [] nalldof;

    delete [] ndofnmas;

    for (i=0;i<tnnp;i++){
      if (dofindmas)
        delete [] dofindmas[i];
      if (dofind)
        delete [] dofind[i];
    }
    delete [] dofindmas;
    delete [] dofind;

    if (icnbnmas)
    {
      for (i=0;i<nproc;i++){
        delete [] icnbnmas[i];
      }
    }
    delete [] icnbnmas;

    delete [] icmultip;

    for (i=0;i<nproc;i++){
      if (gnbn)
        delete [] gnbn[i];
      if (gnin)
        delete [] gnin[i];
    }
    delete [] gnbn;
    delete [] gnin;

    for (i=0;i<tnbn;i++){
      if (gnbncn)
        delete [] gnbncn[i];
      if (sid)
        delete [] sid[i];
    }
    delete [] gnbncn;
    delete [] sid;

    delete [] multip;

    for (i=0;i<nproc;i++){
      if (cnbn)
        delete [] cnbn[i];
      if (nbdofnd)
        delete [] nbdofnd[i];
    }
    delete [] cnbn;
    delete [] nbdofnd;

    delete [] nbdofd;

    delete [] bnmultip;
    
    for (i=0;i<tnbn;i++){
      if (llnbn)
        delete [] llnbn[i];
      if (ldn)
        delete [] ldn[i];
    }
    delete [] llnbn;
    delete [] ldn;
    
    //delete [] gnn;
    //delete [] pgcn;

// !!!! maze se primo v jdenom bloku po alokaci, zvazit odstraneni jakozto atributu tridy
//  delete [] gcnbn;

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
  
  //  number of nodes on subdomain
  nn = top->nn;
  //  number of elements on subdomain
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
      top->gnodes[i].ai = ltg[i];
    }
    break;
  }
  default:{
    par_print_err(myrank,"unknown type of mesh description", __FILE__, __LINE__, __func__);
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
   maxnbn - maximum number of boundary nodes
   
   nbn - number of boundary nodes on each subdomain
   
   tnbn - total number of boundary nodes
   
   nbnd - array containing numbers of boundary nodes on subdomains
   
   tnnp - total number of nodes in problem
   
   multip - array of boundary node multiplicity (on the master)
   
   nodmultip - array of node multiplicity (on all processors)
   
   
   @param top - pointer to subdomain topology
   @param domproc - array containing domain-processor correspondence
   
   JK, 30.6.2005
*/
void partop::compute_multiplicity (long *ltg,long *domproc,FILE *out)
{
  long i,j,k,buffsize;
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
    nbn=0;
    maxnbn = 0;
    tnbn=0;
    for (i=0;i<nn;i++){
      if (tnbn<ltg[i])
        tnbn=ltg[i];
      if (ltg[i]>-1)
        nbn++;
    }
    
    buff[0]=tnbn;
    buff[1]=nbn;
    
    if (myrank==0){
      //  list of number of boundary nodes on subdomains
      if (nbnd != NULL)
        delete [] nbnd;
      nbnd = new long [nproc];
      memset(nbnd, 0, sizeof(*nbnd)*nproc);

      //  master contribution
      j=domproc[0];
      if (tnbn<buff[0])  tnbn=buff[0];
      if (maxnbn<buff[1])  maxnbn=buff[1];
      nbnd[j]=buff[1];
      
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        j=domproc[stat.MPI_TAG];
        //  slave contributions
        if (tnbn<buff[0])
          tnbn=buff[0];
        if (maxnbn<buff[1])
          maxnbn=buff[1];
        nbnd[j]=buff[1];
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
    fprintf (stderr,"\n function compute_multiplicity (file %s, line %d)\n",__FILE__,__LINE__);
  }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  fprintf(out,"\n\n\n Multiplicity of nodes on subdomain\n\n");  
  for (i=0;i<nn;i++)
    fprintf(out,"%ld   %ld\n",i,nodmultip[i]);  
  

  fprintf(out,"\n\n\n total number of boundary nodes in whole problem is %ld\n\n",tnbn);
  fprintf (out,"\n\n\n maximum number of boundary nodes on subdomain is %ld\n\n",maxnbn);
  if (myrank==0){
    fprintf (out,"\n\n\n numbers of boundary nodes on each subdomain \n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nbnd[i]);
    }
  }  
  fprintf(out,"\n\n\n total number of boundary nodes in whole problem is %ld\n\n",tnbn);
  fprintf (out,"\n\n\n maximum number of boundary nodes on subdomain is %ld\n\n",maxnbn);
  if (myrank==0){
    fprintf (out,"\n\n\n numbers of boundary nodes on each subdomain \n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nbnd[i]);
    }
  }  
  
}

/**
   function assembles multiplicity of degrees of freedom
   it is different from the function compute_multiplicity, which
   computes multiplicity of nodes
   
   JK, 14.6.2007
*/
void partop::dof_multiplicity (gtopology *top,FILE */*out*/)
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
   
   lgnbn - global numbers of boundary nodes
   lnbn - local numbers of boundary nodes
   nin - number of interanl nodes
   lnin - local numbers of internal nodes
   
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
  /*  
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
  */
  
  
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
    fprintf (out,"\n %ld   %ld",i,lnbn[i]);
  }

  fprintf (out,"\n\n\n numbers of internal nodes on subdomain \n\n");
  fprintf (out,"\n %ld",nin);
  fprintf (out,"\n\n\n local numbers of internal nodes on subdomain");
  for (i=0;i<nin;i++){
    fprintf (out,"\n %ld   %ld",i,lnin[i]);
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
void partop::boundary_nodes_on_master (gtopology *top,long *domproc,FILE */*out*/)
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
  
  /*
  if (myrank==0){
    fprintf (out,"\n\n\n\n coarse numbering of boundary nodes  - cnbn\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain number %ld\n",i+1);
      for (j=0;j<nbnd[i];j++){
        fprintf (out,"%ld %ld\n",j,cnbn[i][j]);
      }
    }
  }
  fprintf(out,"\ntnbn je %ld\n\n",tnbn);
  */
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
void partop::sort_nodes (FILE */*out*/)
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
  
  /*
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
    
    
    
    //fprintf (out,"\n\n  ldn and llnbn array\n");
    // llnbn ukazuje na pozici v poli lnbn
    //for (i=0;i<tnbn;i++){
    //fprintf (out,"\n boundary node %4ld  ",i+1);
    //for (j=0;j<bnmultip[i];j++){
    //fprintf (out," ldn =  %ld  llnbn =  %ld  lnbn = %ld,",ldn[i][j],llnbn[i][j],lnbn[llnbn[i][j]]);
    //}
    //}
  }
  */
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
void partop::number_of_bdofs_on_nodes (gtopology *top,long *domproc,FILE */*out*/)
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
  
  /*
  if (myrank==0){
    fprintf (out,"\n\n\n numbers of DOFs at boundary nodes");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain %ld",i);
      for (j=0;j<nbnd[i];j++){
        fprintf (out,"\n %ld %ld",j+1,nbdofnd[i][j]);
      }
    }
  }
  */
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
void partop::code_numbers_on_master (gtopology *top,long *domproc,FILE */*out*/)
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
  
  /*
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
  */
  
  
  MPI_Barrier (MPI_COMM_WORLD);
  
}


































/*
  ZACATEK ROZSAHLYCH UPRAV
  27.7.2009
*/

/**
   function collects numbers of all nodes on subdomains
   
   the following function has to be called before this function:
   void initiation (gtopology *top,long *ltg);
   
   
   function establishes:

   maxnn - maximum number of nodes on one subdomain
   nnsd - (M) array containing numbers of nodes on subdomains
   
   @param domproc - array containing domain-processor correspondence
   @param out - output file (for auxiliary output)

   JK, 30.6.2005
*/
void partop::numbers_of_all_nodes_on_subdomains (long *domproc,FILE *out)
{
  long i,j,k;
  MPI_Status stat;
  
  //  the maximum number of nodes defined on subdomain (will be sent to the master)
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
    fprintf (out,"\n\n\n partop::numbers_of_all_nodes_on_subdomains");
    fprintf (out,"\n\n the maximum number of nodes on subdomains (maxnn) %ld",maxnn);
    fprintf (out,"\n\n array nnsd");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nnsd[i]);
    }
  }
  
}

/**
   function assembles the array bmultip or amultip
   function computes the variables tnnp or tnbn
   
   this function can assemble also arrays nbnd and nind
   but the coupled DOFs generate additional boundary/interface nodes
   which are not known in this function
   therefore, there is function assemble_nbnd which does it
   
   function determines:
   
   maxnbn - maximum number of boundary/interface nodes on subdomain
   nbn - number of boundary/interface nodes
   tnbn - total number of boundary/interface nodes in the problem
   bnodes - (M)
   bmultip - (M)
   maxnbn - maximum number of boundary/interface nodes on subdomain
   tnnp - (M) total number of nodes in the problem
   allnodes - (M)
   amultip - (M) 

   @param ltg - local to global correspondence
   @param domproc - array containing domain-processor correspondence
   @param out - output file for auxiliary print
   @param proc_name - processor name
   
   27.7.2009, JK
*/
void partop::assemble_multip (long *ltg,long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;
  
  if (maxnn==0){
    par_print_err(myrank,proc_name,"maximum number of nodes on one subdomain is equal to zero",__FILE__, __LINE__, __func__);
    par_print_err(myrank,proc_name,"function numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
  }
  if (myrank==0){
    if (nnsd==NULL){
      par_print_err(myrank,proc_name,"array nnsd (number of nodes on subdomains) is not assembled",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
    }
  }
  
  switch (md){
  case bound_nodes:{
    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes
    
    /*
    buffsize=2;
    buff = new long [buffsize];
    
    //  the number of boundary/interface nodes on subdomain
    nbn = 0;
    //  the maximum number of boundary/interface nodes on subdomain
    maxnbn = 0;
    //  the total number of boundary/interface nodes in the problem
    tnbn = 0;

    //  loop over the number of all nodes on subdomain
    for (i=0;i<nn;i++){
      if (tnbn<ltg[i])
        tnbn=ltg[i];
      if (ltg[i]>-1)
        nbn++;
    }
    
    buff[0]=tnbn;
    buff[1]=nbn;
    
    if (myrank==0){
      //  master contribution
      j=domproc[0];
      if (tnbn<buff[0])  tnbn=buff[0];
      if (maxnbn<buff[1])  maxnbn=buff[1];
      nbnd[j]=buff[1];
      
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        j=domproc[stat.MPI_TAG];
        //  slave contributions
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
    */



    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes
    
    buffsize=2;
    buff = new long [buffsize];
    
    //  the maximum number of boundary/interface nodes on subdomain
    maxnbn = 0;
    //  the number of boundary/interface nodes on subdomain
    nbn = 0;
    //  the total number of boundary/interface nodes in the problem
    tnbn = 0;

    //  loop over the number of all nodes on subdomain
    for (i=0;i<nn;i++){
      if (tnbn<ltg[i])
        tnbn=ltg[i];
      if (ltg[i]>-1)
        nbn++;
    }
    
    buff[0]=tnbn;
    buff[1]=nbn;
    
    if (myrank==0){
      //  master contribution
      j=domproc[0];
      if (tnbn<buff[0])  tnbn=buff[0];
      if (maxnbn<buff[1])  maxnbn=buff[1];
      
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        j=domproc[stat.MPI_TAG];
        //  slave contributions
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
    
    
    
    buffsize = maxnn;
    buff = new long [buffsize];
    for (i=0;i<maxnn;i++){
      buff[i]=-1;
    }
    
    for(i=0;i<nn;i++){
      buff[i] = ltg[i];
    }
    
    
    if(myrank == 0){

      if (bnodes != NULL){
        for (i=0;i<nproc;i++){
          delete [] bnodes[i];
        }
        delete [] bnodes;
      }
      bnodes = new long* [nproc];
      for (i=0;i<nproc;i++){
        bnodes[i] = new long [nnsd[i]];
      }
      
      if (bmultip!=NULL){
        delete [] bmultip;
      }
      bmultip = new long [tnbn];
      for(i=0;i<tnbn;i++){
        bmultip[i]=0;
      }
      
      // master contribution
      k=domproc[0];
      for(j=0;j<nnsd[k];j++){
        bnodes[k][j]=buff[j];
        if (buff[j]>-1){
          bmultip[buff[j]]++;
        }
      }
      
      //  slave contributions
      for (i=1;i<nproc;i++){
        MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

        k=domproc[stat.MPI_TAG];
        for(j=0;j<nnsd[k];j++){
          bnodes[k][j]=buff[j];
          if (buff[j]>-1){
            bmultip[buff[j]]++;
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
        buff[i]=bmultip[i];
      }
      
      for (i=1;i<nproc;i++){
        MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    delete [] buff;
       
    break;
  }



  case all_nodes:
  case neg_bound_nodes:{
    
    buffsize=maxnn+1;
    buff = new long [buffsize];

    //  the total number of nodes in the problem
    tnnp=0;
    //  loop over the number of nodes on subdomain
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
      
      if (amultip != NULL)
        delete [] amultip;
      amultip = new long [tnnp];
      for (i=0;i<tnnp;i++){
        amultip[i]=0;
      }
      
      //  computation of node incidences
      for (i=0;i<nproc;i++){
        for (j=0;j<nnsd[i];j++){
          amultip[allnodes[i][j]]++;
        }
      }
    }
    delete [] buff;
    break;
  }
  default:{
    par_print_err(myrank,proc_name,"unknown type of mesh description",__FILE__, __LINE__, __func__);
  }
  }
  MPI_Barrier (MPI_COMM_WORLD);

  
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    if (allnodes!=NULL){
      fprintf (out,"\n\n array allnodes \n");
      for (i=0;i<nproc;i++){
        fprintf (out,"\n");
        for (j=0;j<nnsd[i];j++){
          fprintf (out,"\n nproc %3ld  node %6ld   allnode %6ld",i,j,allnodes[i][j]);
        }
      }
    }
    if (amultip!=NULL){
      fprintf(out,"\n\n\n tnnp - total number of nodes in whole problem is %ld\n\n",tnnp);
      for (i=0;i<tnnp;i++){
        fprintf (out,"\n node %6ld  amultip %3ld",i,amultip[i]);
      }
    }
    if (bmultip!=NULL){
      for (i=0;i<tnbn;i++){
        fprintf (out,"\n boundary/interface node %6ld  bmultip %3ld",i,bmultip[i]);
      }
    }
    
    if (bnodes!=NULL){
      fprintf(out,"\n\n\n array bnodes \n\n");
      for (i=0;i<nproc;i++){
        fprintf (out,"\n domain %ld\n",i);
        for (j=0;j<nnsd[i];j++){
          fprintf (out,"\n bnodes  %6ld %6ld",j,bnodes[i][j]);
        }
      }
    }

    fprintf(out,"\n\n\n total number of boundary/interface nodes in whole problem (tnbn) is %ld\n\n",tnbn);
    fprintf (out,"\n\n\n maximum number of boundary/interface nodes on subdomain (maxnbn) is %ld\n\n",maxnbn);
  }
  
}



/**
   function searches for coupled DOFs which are shared
   by interfaces

   ncdof - the number of coupled DOFs
   coupdof - (M) array containing number of coupled DOFs on subdomains
   coupdofmas - (M) array containing suspicious indicators of coupled DOFs
   nbcdof - the number of boundary/interface coupled DOFs
   
   @param top - pointer to the general topology
   @param domproc - array containing domain-processor correspondence
   @param out - output file
   
   28.7.2009, JK
   TKo, 20.5.2010
*/
void partop::coupled_dofs (gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,l,ndofn;
  long *buff,buffsize;
  MPI_Status stat;
  
  //  the number of coupled DOFs
  //  
  ncdof=0;

  //  loop over the number of nodes on subdomain
  for (i=0;i<nn;i++){
    //  the number of DOFs on node
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs on node
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>1){
        if (ncdof<k)
          ncdof=k;
      }
    }
  }
  
  
  // *********
  //  master
  // *********
  if (myrank==0){
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //  slave contributions
      if (ncdof<k)  ncdof=k;
      
    }
    
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      MPI_Send (&ncdof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
    
  }
  // ********
  // slaves
  // ********
  else{
    MPI_Send (&ncdof,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    MPI_Recv (&ncdof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  buffsize=ncdof;
  buff=new long [buffsize];
  for (i=0;i<buffsize;i++){
    buff[i]=0;
  }
  
  //  loop over the number of nodes on subdomain
  for (i=0;i<nn;i++){
    //  the number of DOFs on node
    ndofn=top->give_ndofn (i);
    //  loop over the number of DOFs on node
    for (j=0;j<ndofn;j++){
      k=top->give_dof (i,j);
      if (k>1){
        buff[k-1]=k;
      }
    }
  }
  
  //fprintf (out,"\n\n");
  //for (i=0;i<buffsize;i++){
  //fprintf (out,"\n buff %3ld  %ld",i,buff[i]);
  //}


  // *********
  //  master
  // *********
  if (myrank==0){
    
    //  array containing number of coupled DOFs on subdomains
    //  the indicator is greater than one, otherwise the
    //  code number is not coupled
    if (coupdof!=NULL){
      for (i=0;i<nproc;i++){
        delete [] coupdof[i];
      }
      delete [] coupdof;
    }
    coupdof = new long* [nproc];
    for (i=0;i<nproc;i++){
      coupdof[i] = new long [ncdof];
      for (j=0;j<ncdof;j++){
        coupdof[i][j]=0;
      }
    }

    //  master contribution
    k=domproc[0];
    for (j=0;j<ncdof;j++){
      l=buff[j];
      if (l>1)
        coupdof[k][l-1]++;
    }
    
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //  slave contributions
      k=domproc[stat.MPI_TAG];
      for (j=0;j<ncdof;j++){
        l=buff[j];
        if (l>1)
          coupdof[k][l-1]++;
      }
    }
    
    //  array containing suspicious indicators of coupled DOFs
    if (coupdofmas!=NULL){
      delete [] coupdofmas;
    }
    coupdofmas = new long [ncdof];
    for (i=0;i<ncdof;i++){
      coupdofmas[i]=0;
    }
    
    //  loop over the number of coupled DOFs
    for (i=0;i<ncdof;i++){
      k=0;
      //  loop over the number of subdomains
      for (j=0;j<nproc;j++){
        if (coupdof[j][i]>0)
          k++;
      }
      if (k>1){
        //  coupled DOFs described by the i-th indicator occur at least on two subdomains
        //  such coupled DOFs have to be treated as the boundary/interface DOFs
        coupdofmas[i]=1;
      }
    }
    
    //  the number of boundary/interface coupled DOFs
    nbcdof=0;
    for (i=0;i<ncdof;i++){
      if (coupdofmas[i]==1)
        nbcdof++;
    }
    
  }
  // ********
  // slaves
  // ********
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  // *********
  //  master
  // *********
  if (myrank==0){
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      MPI_Send (&nbcdof,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  // ********
  // slaves
  // ********
  else
  {
    MPI_Recv (&nbcdof,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n\n\n partop::coupled_dofs");
    fprintf (out,"\n\n ncdof - the number of coupled DOFs %ld",ncdof);
    fprintf (out,"\n nbcdof - the number of boundary/interface coupled DOFs %ld",nbcdof);
    fprintf (out,"\n\n array coupdof");
    for (i=0;i<nproc;i++){
      for (j=0;j<ncdof;j++){
        fprintf (out,"\n coupdof %3ld %3ld   %ld",i,j,coupdof[i][j]);
      }
    }
    fprintf (out,"\n\n");
    for (i=0;i<ncdof;i++){
      fprintf (out,"\n coupdofmas %6ld   %ld",i,coupdofmas[i]);
    }
  }
  delete [] buff;
}



/**
   function collects the numbers of all degrees of freedom on subdomains
   constrained and prescribed DOFs are taken into account, it means
   they are added to the number of DOFs
   
   function establishes:
   
   maxndof - maximum number of degrees of freedom on one subdomain
             it summarizes number of DOFs at nodes belonging to one subdomain
             it differs from the actual degrees of freedom of subdomains in constraints
             
   nalldof - (M) array of numbers of degrees of freedom on subdomains


   @param top - pointer to the sequential general topology
   @param domproc - array containing domain-processor correspondence
   @param out - output file (for auxiliary output)
   
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
    fprintf (out,"\n\n\n the maximum number of DOFs on subdomain (maxndof) %ld\n",maxndof);
    fprintf (out,"\n\n\n the numbers of all DOFs on each subdomain (alldof)\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nalldof[i]);
    }
  }

}

/**
   function assembles array ndofnmas
   it contains the number of DOFs on nodes
   
   @param top - pointer to the general topology
   @param domproc - array containing domain-processor correspondence
   @param out - output file
   @param proc_name - processor name
   
   JK, 28.7.2009
*/
void partop::ndofn_master (gtopology *top,long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,l;
  long *buff,buffsize;
  MPI_Status stat;
  
  /*
  if (md!=all_nodes){
    par_print_err(myrank,proc_name,"wrong mesh description is used",__FILE__, __LINE__, __func__);
    abort ();
  }
  */

  if (maxnn==0){
    par_print_err(myrank,proc_name,"maximum number of nodes on one subdomain is equal to zero",__FILE__, __LINE__, __func__);
    par_print_err(myrank,proc_name,"function numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
    abort ();
  }
  
  buffsize=maxnn;
  buff = new long [buffsize];
  
  //  loop over the number of nodes on subdomain
  for (i=0;i<nn;i++){
    //  number of DOFs on nodes
    buff[i] = top->give_ndofn (i);
  }
  
  // *********
  //  master
  // *********
  if (myrank==0){
    
    if (tnnp==0){
      par_print_err(myrank,proc_name,"total number of nodes in the problem (tnnp) is zero",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
    }
    if (nnsd==NULL){
      par_print_err(myrank,proc_name,"array nnsd (number of nodes on subdomains) is not assembled",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
    }
    if (allnodes==NULL){
      par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
    }
    
    //  array containing the numbers of DOF on nodes
    if (ndofnmas!=NULL){
      delete [] ndofnmas;
    }
    ndofnmas = new long [tnnp];
    memset(ndofnmas, 0, sizeof(*ndofnmas)*tnnp);
    
    //  master contribution
    
    //  subdomain id
    k=domproc[0];
    //  loop over the number of nodes
    for (j=0;j<nnsd[k];j++){
      //  global node number
      l=allnodes[k][j];
      //  the number of DOFs on the node
      ndofnmas[l]=buff[j];
    }
    
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      
      //  subdomain id
      k=domproc[stat.MPI_TAG];
      //  loop over the number of nodes
      for (j=0;j<nnsd[k];j++){
        //  global node number
        l=allnodes[k][j];
        //  the number of DOFs on the node
        ndofnmas[l]=buff[j];
      }
      
    }
    
  }
  // *********
  //  slaves
  // *********
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf(out,"\n\n\n the numbers of DOFs on nodes (ndofnmas)\n\n");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n node %6ld  ndofnmas %3ld",i,ndofnmas[i]);
    }
  }  
  
  delete [] buff;  
}

/**
   function collects id of DOFs on the master
   
   function establishes:
   
   difind_mas - 

   @param top - pointer to the general topology
   @param domproc - array containing domain-processor correspondence
   @param out - output file
   @param proc_name - processor name
   
   JK, 28.7.2009
*/
void partop::dofind_master (gtopology *top,long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,l,m,ii,ndofn;
  long *buff,buffsize;
  MPI_Status stat;
  
  /*
  if (md!=all_nodes){
    par_print_err(myrank,proc_name,"wrong mesh description is used",__FILE__, __LINE__, __func__);
    abort ();
  }
  */

  if (maxndof==0){
    par_print_err(myrank,proc_name,"maximum number of DOFs is zero",__FILE__, __LINE__, __func__);
    par_print_err(myrank,proc_name,"function partop::numbers_of_all_dofs_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
  }
  
  buffsize=maxndof;
  buff = new long [buffsize];
  
  //  index of array buff
  k=0;
  //  loop over the number of nodes on subdomain
  for (i=0;i<nn;i++){
    //  number of DOFs on nodes
    ndofn = top->give_ndofn (i);
    //  loop over the number of DOFs on node
    for (j=0;j<ndofn;j++){
      buff[k]=top->give_dof (i,j);
      k++;
    }
  }
  
  // *********
  //  master
  // *********
  if (myrank==0){
    
    //  array containing DOF indicators
    if (dofindmas!=NULL){
      for (i=0;i<tnnp;i++){
        delete [] dofindmas[i];
      }
      delete [] dofindmas;
    }
    dofindmas = new long* [tnnp];
    for (i=0;i<tnnp;i++){
      dofindmas[i] = new long [ndofnmas[i]];
    }

    
    //  master contribution
    
    ii=0;
    //  subdomain id
    k=domproc[0];
    //  loop over the number of nodes
    for (j=0;j<nnsd[k];j++){
      //  global node number
      l=allnodes[k][j];
      //  loop over the number of DOFs in node
      for (m=0;m<ndofnmas[l];m++){
        dofindmas[l][m]=buff[ii];
        ii++;
      }
    }
    
    
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      //  slave contributions
      
      ii=0;
      //  subdomain id
      k=domproc[stat.MPI_TAG];
      //  loop over the number of nodes
      for (j=0;j<nnsd[k];j++){
        //  global node number
        l=allnodes[k][j];
        //  loop over the number of DOFs in node
        for (m=0;m<ndofnmas[l];m++){
          dofindmas[l][m]=buff[ii];
          ii++;
        }
      }
      
    }
    
  }
  // *********
  //  slaves
  // *********
  else{
    MPI_Send (buff,buffsize,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf(out,"\n\n\n DOF id (dofindmas)\n\n");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n node %6ld  ",i);
      for (j=0;j<ndofnmas[i];j++){
        fprintf (out,"  dofindmas %3ld",dofindmas[i][j]);
      }
    }
    fprintf (out,"\n\n");
  }  
  delete [] buff;
}




/**
   function updates node multiplicity

   
   function establishes:
   
   tnbn - (M) total number of boundary nodes
   amultip - (M) array of all node multiplicity, it is assembled if md=all_nodes
   dofind - (M) array containing DOF indicators
   
   @param top - pointer to the general topology
   @param out - output file for auxiliary print
   @param proc_name - processor name
   
   JK, 28.7.2009
*/
void partop::update_multip (gtopology */*top*/,FILE *out,char */*proc_name*/)
{
  long i,j,k,l,ndofn,dofid;
  
  /*
  if (md!=all_nodes){
    par_print_err(myrank,proc_name,"wrong mesh description is used",__FILE__, __LINE__, __func__);
    abort ();
  }
  */

  if (myrank==0){
    
    //  array containing DOF indicators
    if (dofind!=NULL){
      for (i=0;i<tnnp;i++){
        delete [] dofind[i];
      }
      delete [] dofind;
    }
    dofind = new long* [tnnp];
    for (i=0;i<tnnp;i++){
      dofind[i] = new long [ndofnmas[i]];
      for (j=0;j<ndofnmas[i];j++){
        dofind[i][j]=0;
      }
    }

    // *********************************************************************************
    //  in the case of coupled DOFs, treatment with nodes is not enough
    //  DOFs have to be split to group of internal DOFs and boundary/interface DOFs
    //
    //  there are three types of boundary/interface DOFs:
    //
    //  all DOFs belonging to the boundary/interface nodes are boundary/interface DOFs
    //
    //  coupled DOFs connected to any boundary/interface node are boundary/interface DOFs
    //
    //  coupled DOFs not connected to any boundary/interface node but shared by at least
    //  two subdomains are also boundary/interface DOFs
    //
    //  all remaining DOFs are internal DOFs and can be eliminated on subdomains
    //
    // *********************************************************************************
    
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      //  loop over the number of nodes on subdomains
      for (j=0;j<nnsd[i];j++){
        //  global node number
        l=allnodes[i][j];
        //  the number of DOFs in node
        ndofn = ndofnmas[l];
        //  loop over the number of DOFs in node
        for (k=0;k<ndofn;k++){
          //  DOF indicator
          dofid=dofindmas[l][k];
          if (amultip[l]>1){
            //  the i-th node is boundary/interface node, therefore this is boundary/interface DOF
            dofind[l][k]=1;
            if (dofid>1){
              //  there is coupled DOF defined in this boundary/interface node
              //  this information has to be stored in the array coupdofmas
              //  it will be used in the following loop
              if (coupdofmas[dofid-1]==0)
                coupdofmas[dofid-1]=1;
            }
          }
          if (dofid>1){
            
            //  there is coupled DOF connected to at least two subdomains
            if (coupdofmas[dofid-1]==1){
              //  this coupled DOF has to be taken into account
              //  this is boundary/interface DOF
              dofind[l][k]=1;
            }

          }
        }
      }
    }
    
    fprintf (out,"\n\n\n");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n dofind %6ld ",i);
      for (j=0;j<ndofnmas[i];j++){
        fprintf (out,"  %3ld",dofind[i][j]);
      }
    }
    
    
    // *********************************************************************************
    //  modification of the array amultip
    //  the array has to be modified with respect to the array dofind
    //  if any DOF in node is denoted as the boundary/interface DOF, the appropriate node
    //  has to be denoted as the boundary/interface node, this fact is
    //  assured with the help of the array amultip, where value 1 is rewritten to the value 2
    // *********************************************************************************
    
    //  total number of boundary/interface nodes
    tnbn=0;
    //  loop over the number of subdomain
    for (i=0;i<nproc;i++){
      //  loop over the number of nodes on subdomain
      for (j=0;j<nnsd[i];j++){
        //  global node number
        l=allnodes[i][j];
        //  the number of DOFs in node
        ndofn = ndofnmas[l];
        //  loop over the number of DOFs in node
        for (k=0;k<ndofn;k++){
          if (dofind[l][k]==1){
            //  the DOF is boundary/interface DOF
            if (amultip[l]==1){
              //  the node is internal node
              //  it has to be changed from internal to boundary/interface
              amultip[l]=2;
            }
          }
        }
      }
    }
    // *********************************************************************************
    //  end of modification of the array amultip
    // *********************************************************************************
    
  }
  
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    if (amultip!=NULL){
      fprintf(out,"\n\n\n modified array amultip with respect to the coupled DOFs\n");
      for (i=0;i<tnnp;i++){
        fprintf (out,"\n node %6ld  amultip %3ld",i,amultip[i]);
      }
    }
    fprintf(out,"\n\n\n DOF id (dofind)\n\n");
    for (i=0;i<tnnp;i++){
      fprintf (out,"\n node %6ld  ",i);
      for (j=0;j<ndofnmas[i];j++){
        fprintf (out,"  dofind %3ld",dofind[i][j]);
      }
    }
    fprintf (out,"\n\n");
  }
  
}

/**
   function assembles the arrays nbnd a nind
   
   these arrays are assembled separately in this function
   they can be assembled in the function assemble_multip
   but the coupled DOFs generate additional boundary/interface
   nodes, these nodes are not known in the function assemble_multip
   
   function assembles:
   nbnd - (M) array of the numbers of boundary/interface nodes
   nind - array of the numbers of internal nodes
   tnbn - the total number of all boundary/interface nodes
          it is computed once again because additional boundary/interface nodes
          may be added due to the coupled DOFs
   maxnbn - maximum number of boundary/interface nodes
   nbn - the number of boundary/interface nodes
   nin - the number of internal nodes
   

   @param ltg - local to global correspondence
   @param domproc - array containing domain-processor correspondence
   @param out - output file (for auxiliary output)
   @param proc_name - processor name

   JK, 29.7.2009
   TKo, 20.5.2010
*/
void partop::assemble_nbnd_nind (long *ltg,long *domproc,FILE *out,char *proc_name)
{
  long i,j,k;
  long *buff,buffsize;
  MPI_Status stat;
  
  switch (md){
  case bound_nodes:{
    
    if (maxnn==0){
      par_print_err(myrank,proc_name,"maximum number of nodes on one subdomain is equal to zero",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
    }
    
    //  determination of total number of boundary nodes
    nbn=0;
    for (i=0;i<nn;i++){
      if (ltg[i]>-1)
        nbn++;
    }
    
    if (myrank==0){
      //  list of number of boundary nodes on subdomains
      if (nbnd != NULL)
        delete [] nbnd;
      nbnd = new long [nproc];
      
      maxnbn=0;
      
      //  master contribution
      j=domproc[0];
      nbnd[j]=nbn;
      if (maxnbn<nbn)
        maxnbn=nbn;
      
      for (i=1;i<nproc;i++){
        MPI_Recv (&k,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        
        //  slave contributions
        j=domproc[stat.MPI_TAG];
        nbnd[j]=k;
        if (maxnbn<k)
          maxnbn=nbnd[j];
        
      }
    }
    else{
      MPI_Send (&nbn,1,MPI_LONG,0,myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (myrank==0){
      for (i=1;i<nproc;i++){
        MPI_Send (&maxnbn,1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
    }
    else
      MPI_Recv (&maxnbn,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

    MPI_Barrier (MPI_COMM_WORLD);

    break;
  }
    
    
    
    
  case all_nodes:
  case neg_bound_nodes:{
    
    buffsize=2;
    buff = new long [buffsize];
    
    // *********
    //  master
    // *********
    if (myrank==0){
      
      if (nnsd==NULL){
        par_print_err(myrank,proc_name,"array nnsd is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      }
      if (amultip==NULL){
        par_print_err(myrank,proc_name,"array amultip is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      if (allnodes==NULL){
        par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }

      //  list of the numbers of boundary nodes on subdomains
      if (nbnd != NULL)
        delete [] nbnd;
      nbnd = new long [nproc];
      
      //  maximum number of boundary/interface nodes on subdomain
      maxnbn=0;
      for (i=0;i<nproc;i++){
        k=0;
        for (j=0;j<nnsd[i];j++){
          if (amultip[allnodes[i][j]]>1)  k++;
        }
        nbnd[i]=k;
        if (maxnbn<nbnd[i])
          maxnbn=nbnd[i];
      }
      
      //  the total number of all boundary nodes
      //  it has to be computed once again because additional
      //  boundary/interface nodes have been defined
      tnbn=0;
      for (i=0;i<tnnp;i++){
        if (amultip[i]!=1)  tnbn++;
      }

      //  loop over the number of slaves
      for (i=1;i<nproc;i++){
        k=domproc[i];
        buff[0]=nbnd[k];
        buff[1]=maxnbn;
        MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
      
      k=domproc[0];
      nbn=nbnd[k];
    }
    else{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      nbn = buff[0];
      maxnbn = buff[1];
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    delete [] buff;
    
    break;
  }
  default:{
    par_print_err(myrank,proc_name,"unknown type of mesh description is required",__FILE__, __LINE__, __func__);
  }
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  //  number of internal nodes
  nin = nn - nbn;
  
  if (nind != NULL)
    delete [] nind;
  nind = new long [nproc];
  tnin = 0;
  buff = new long[nproc+1];
  if (myrank==0){
    //  list of the numbers of internal nodes on subdomains
    for (i=0;i<nproc;i++){
      nind[i]=nnsd[i]-nbnd[i];
      tnin += nind[i];
      buff[i] = nind[i];
    }
    buff[nproc] = tnin;
    for (i=1;i<nproc;i++){
      MPI_Send(buff,nproc+1,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(buff,nproc+1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    for (i=0;i<nproc;i++)
      nind[i] = buff[i];
    tnin = buff[nproc];
  }
  MPI_Barrier (MPI_COMM_WORLD);

  delete [] buff;
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n the number of boundary/interface nodes on subdomain (nbn)             %ld\n",nbn);
  fprintf (out,"\n the number of internal nodes on subdomain (nin)                       %ld\n",nin);
  fprintf (out,"\n the number of all nodes on subdomain (nn)                             %ld\n",nn);
  fprintf (out,"\n the maximum number of all nodes on subdomain (maxnn)                  %ld\n",maxnn);
  fprintf (out,"\n the maximum number of boundary/interface nodes on subdomain (maxnbn)  %ld\n",maxnbn);
  if (myrank==0){
    fprintf (out,"\n\n\n the numbers of boundary/interface nodes on each subdomain (nbnd)\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nbnd[i]);
    }
    fprintf (out,"\n\n\n the numbers of internal nodes on each subdomain (nind)\n\n");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n %ld",nind[i]);
    }
  }  

}




/**
   function assembles the array nodmultip
   
   the array nodmultip can be assembled in the function assemble_multip
   but the coupled DOFs can generate additional boundary/interface nodes
   
   @param ltg - local to global correspondence
   @param domproc - array containing domain-processor correspondence
   @param out - output file for auxiliary print
   @param proc_name - processor name
   
   27.7.2009, JK
*/
void partop::assemble_nodmultip (long *ltg,long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,buffsize;
  long *buff;
  MPI_Status stat;
  
  switch (md){
  case bound_nodes:{
    //  ltg=-1 for the internal nodes
    //  ltg>-1 for the interface/boundary nodes
    
    buffsize=tnbn;
    buff = new long [buffsize];
    
    if (myrank==0){
      for (i=0;i<tnbn;i++){
        buff[i]=bmultip[i];
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
    
    if (maxnn==0){
      par_print_err(myrank,proc_name,"maximum number of nodes on one subdomain is equal to zero",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      abort ();
    }
    if (myrank==0){
      if (nnsd==NULL){
        par_print_err(myrank,proc_name,"array nnsd is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      }
      if (amultip==NULL){
        par_print_err(myrank,proc_name,"array amultip is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      if (allnodes==NULL){
        par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
    }
    
    buffsize=maxnn;
    buff = new long [buffsize];
    
    if (myrank==0){
      for (i=1;i<nproc;i++){
        k=domproc[i];
        for (j=0;j<nnsd[k];j++){
          buff[j]=amultip[allnodes[k][j]];
        }
        MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
      
      k=domproc[0];
      for (j=0;j<nnsd[k];j++){
        buff[j]=amultip[allnodes[k][j]];
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
    par_print_err(myrank,proc_name,"unknown type of mesh description is required",__FILE__, __LINE__, __func__);
  }
  }
  MPI_Barrier (MPI_COMM_WORLD);

  
  // *******************
  //  auxiliary output
  // *******************
  fprintf(out,"\n\n\n Multiplicity of nodes on subdomain (nodmultip)\n\n");  
  for (i=0;i<nn;i++){
    fprintf(out,"%ld   %ld\n",i,nodmultip[i]);  
  }
  
  
}


/**
   function assembles global numbers of nodes
   
   the following arrays are assembled:
   
   gnbn - (M) global numbers of interface/boundary nodes
   gnin - (M) global numbers of internal nodes
   gnbndom -  array containing global numbers of interface/boundary nodes on subdomain

   @param domproc - array containing domain-processor correspondence
   @param out - output file for auxiliary print
   @param proc_name - processor name
   
   JK, 30.7.2009
*/
void partop::node_global_numbers (long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,ii,jj;
  long *buff,buffsize;
  MPI_Status stat;
  
  if (myrank==0){
    //  array containing global numbers of internal nodes
    if (gnin != NULL){
      for (i=0;i<nproc;i++){
        delete [] gnin[i];
      }
      delete [] gnin;
    }
    gnin = new long* [nproc];
    for (i=0;i<nproc;i++){
      gnin[i] = new long [nind[i]];
      memset(gnin[i], 0, sizeof(*gnin[i])*nind[i]);
    }
    //  array containing global numbers of interface/boundary nodes
    if (gnbn != NULL){
      for (i=0;i<nproc;i++){
        delete [] gnbn[i];
      }
      delete [] gnbn;
    }
    gnbn = new long* [nproc];
    for (i=0;i<nproc;i++){
      gnbn[i] = new long [nbnd[i]];
      memset(gnbn[i], 0, sizeof(*gnbn[i])*nbnd[i]);
    }
  }
  
  //  array containing global numbers of interface/boundary nodes on subdomain
  if (gnbndom != NULL){
    delete [] gnbndom;
  }
  gnbndom = new long [nbn];
  memset(gnbndom, 0, sizeof(*gnbndom)*nbn);

  switch (md){
  case bound_nodes:{
    par_print_err(myrank,proc_name,"mesh description bound_nodes cannot be used",__FILE__, __LINE__, __func__);
    break;
  }
  case metis:
  case all_nodes:{
    
    if (myrank==0){
      
      if (nnsd==NULL){
        par_print_err(myrank,proc_name,"array nnsd is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      }
      if (amultip==NULL){
        par_print_err(myrank,proc_name,"array amultip is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      if (allnodes==NULL){
        par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }

      //  loop over the number of subdomains
      for (i=0;i<nproc;i++){
        //  index in the array gnbn
        ii=0;
        //  index in the array gnin
        jj=0;
        //  loop over the number of all nodes on subdomain
        for (j=0;j<nnsd[i];j++){
          if (amultip[allnodes[i][j]]>1){
            //  boundary/interface node
            
            //  global number of boundary/interface node
            gnbn[i][ii]=allnodes[i][j];
            ii++;
          }
          if (amultip[allnodes[i][j]]==1){
            //  internal node
            
            //  global number of internal node
            gnin[i][jj]=allnodes[i][j];
            jj++;
          }
          
        }
        //printf ("\n ii %ld     jj %ld",ii,jj);
      }
      
    }
    break;
  }
  case neg_bound_nodes:{
    
    if (myrank==0){
      
      if (nnsd==NULL){
        par_print_err(myrank,proc_name,"array nnsd is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      }
      if (allnodes==NULL){
        par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }

      //  loop over the number of subdomains
      for (i=0;i<nproc;i++){
        //  index in the array gnbn
        ii=0;
        //  index in the array gnin
        jj=0;
        //  loop over the number of all nodes on subdomain
        for (j=0;j<nnsd[i];j++){
          if (allnodes[i][j]<0){
            //  boundary/interface node
            gnbn[i][ii]=0-allnodes[i][j]-2;
            ii++;
          }
          if (allnodes[i][j]>-1){
            //  internal node
            gnin[i][jj]=allnodes[i][j];
            jj++;
          }
        }
      }
    }
    
    break;
  }
  default:{
    par_print_err(myrank,proc_name,"unknown type of mesh description is required",__FILE__, __LINE__, __func__);
  }
  }
  
  buffsize=maxnbn;
  buff = new long [buffsize];
  
  if (myrank==0){
    
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      //  subdomain id
      k=domproc[i];
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[k];j++){
        buff[j]=gnbn[k][j];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }

    //  subdomain id
    k=domproc[0];
    //  loop over the number of boundary/interface nodes
    for (j=0;j<nbnd[k];j++){
      buff[j]=gnbn[k][j];
    }
    
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  for (i=0;i<nbn;i++){
    gnbndom[i]=buff[i];
  }
  
  delete [] buff;
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n global numbers of interface/boundary nodes on subdomain (gnbndom)");
  for (i=0;i<nbn;i++){
    fprintf (out,"\n %6ld   %6ld",i,gnbndom[i]);
  }
  if (myrank==0){
    fprintf (out,"\n\n\n global numbers of interface/boundary nodes on subdomain (gnbn)");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain %ld",i);
      for (j=0;j<nbnd[i];j++){
        fprintf (out,"\n %6ld   %6ld",j,gnbn[i][j]);
      }
    }
    
    fprintf (out,"\n\n\n global numbers of internal nodes on subdomain (gnin)");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain %ld",i);
      for (j=0;j<nind[i];j++){
        fprintf (out,"\n %6ld   %6ld",j,gnin[i][j]);
      }
    }
  }
}

/**
   function assembles local numbers of nodes
   
   the following arrays are assembled:
   
   lnbndom - local numbers of interface/boundary nodes
   lnindom - local numbers of internal nodes

   @param out - output file for auxiliary print
   
   JK, 16.9.2009
*/
void partop::node_local_numbers (FILE *out)
{
  long i,j;
  
  //  array containing local numbers of internal nodes on subdomain
  if (lnindom != NULL){
    delete [] lnindom;
  }
  lnindom = new long [nin];
  //  array containing local numbers of interface/boundary nodes on subdomain
  if (lnbndom != NULL){
    delete [] lnbndom;
  }
  lnbndom = new long [nbn];
  
  j=0;
  for (i=0;i<nn;i++){
    if (nodmultip[i]>1){
      lnbndom[j]=i;
      j++;
    }
  }
  j=0;
  for (i=0;i<nn;i++){
    if (nodmultip[i]==1){
      lnindom[j]=i;
      j++;
    }
  }
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n local numbers of interface/boundary nodes on subdomain (lnbndom)");
  for (i=0;i<nbn;i++){
    fprintf (out,"\n %6ld   %6ld",i,lnbndom[i]);
  }
  fprintf (out,"\n\n\n local numbers of internal nodes on subdomain (lnindom)");
  for (i=0;i<nin;i++){
    fprintf (out,"\n %6ld   %6ld",i,lnindom[i]);
  }
}


/**
   function assembles coarse numbers of nodes
   
   the following arrays are assembled:
   
   icnbnmas - (M) coarse numbers of interface/boundary nodes
   icmultip - (M) number of multiplicity of boundary/interface nodes
   icnbn - array containing coarse numbers of interface/boundary nodes on one subdomain

   @param domproc - array containing domain-processor correspondence
   @param out - output file for auxiliary print
   @param proc_name - processor name
   
   JK, 29.7.2009
*/
void partop::node_coarse_numbers (long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,l,buffsize;
  long *av,*buff;
  av=NULL;
  MPI_Status stat;
  
  if (myrank==0){
    
    
    //  array containing coarse numbers of interface/boundary nodes
    if (icnbnmas != NULL){
      for (i=0;i<nproc;i++){
        delete [] icnbnmas[i];
      }
      delete [] icnbnmas;
    }
    icnbnmas = new long* [nproc];
    for (i=0;i<nproc;i++){
      icnbnmas[i] = new long [nbnd[i]];
      memset(icnbnmas[i], 0, sizeof(*icnbnmas[i])*nbnd[i]);
    }
    

    switch (md){
    case bound_nodes:{
      
      if (nnsd==NULL){
        par_print_err(myrank,proc_name,"array nnsd is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      }
      if (bnodes==NULL){
        par_print_err(myrank,proc_name,"array bnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      
      //  loop over the number of subdomains
      for (i=0;i<nproc;i++){
        //  index
        l=0;
        //  loop over the number of boundary/interface nodes
        for (j=0;j<nnsd[i];j++){
          k=bnodes[i][j];
          if (k>-1){
            icnbnmas[i][l]=k;
            l++;
          }
        }
      }
      
      break;
    }
    case metis:
    case all_nodes:
    case neg_bound_nodes:{
      
      if (tnnp==0){
        par_print_err(myrank,proc_name,"total number of nodes in the problem (tnnp) is zero",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      if (amultip==NULL){
        par_print_err(myrank,proc_name,"array amultip is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      if (allnodes==NULL){
        par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      
      j=0;
      av = new long [tnnp];
      for (i=0;i<tnnp;i++){
        av[i]=amultip[i];
        if (av[i]==1)
          av[i]=-1;
        if (av[i]>1){
          av[i]=j;
          j++;
        }
      }
      break;
    }
    default:{
      print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
    }
    }
    
    switch (md){
    case bound_nodes:{
      break;
    }
    case metis:
    case all_nodes:{
      //  loop over the number of subdomains
      for (i=0;i<nproc;i++){
        //  index
        l=0;
        //  loop over the number of boundary/interface nodes
        for (j=0;j<nnsd[i];j++){
          k=allnodes[i][j];
          if (av[k]>-1){
            icnbnmas[i][l]=av[k];
            l++;
          }
        }
      }
      break;
    }
    case neg_bound_nodes:{
      //  loop over the number of subdomains
      for (i=0;i<nproc;i++){
        //  index
        l=0;
        //  loop over the number of boundary/interface nodes
        for (j=0;j<nnsd[i];j++){
          k=allnodes[i][j];
          if (k<0){
            k=0-k-2;
            icnbnmas[i][j]=av[k];
            l++;
          }
        }
      }
      break;
    }  
    default:{
      print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
    }
    }
    
    
    //  number of multiplicity of boundary/interface nodes
    if (icmultip!=NULL){
      delete [] icmultip;
    }
    icmultip = new long [tnbn];
    for (i=0;i<tnbn;i++){
      icmultip[i]=0;
    }
    
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[i];j++){
        icmultip[icnbnmas[i][j]]++;
      }
    }

    if (av!=NULL)
      delete [] av;
  }

  buffsize=maxnbn;
  buff = new long [buffsize];
  
  //  array containing coarse numbers of interface/boundary nodes on one subdomain
  if (icnbn!=NULL){
    delete [] icnbn;
  }
  icnbn = new long [nbn];
  
  if (myrank==0){
    
    //  loop over the number of slaves
    for (i=1;i<nproc;i++){
      //  subdomain id
      k=domproc[i];
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[k];j++){
        buff[j]=icnbnmas[k][j];
      }
      MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
    }

    //  subdomain id
    k=domproc[0];
    //  loop over the number of boundary/interface nodes
    for (j=0;j<nbnd[k];j++){
      buff[j]=icnbnmas[k][j];
    }
    
  }
  else{
    MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  for (i=0;i<nbn;i++){
    icnbn[i]=buff[i];
  }
  
  delete [] buff;
  

  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n\n\n coarse numbers of interface/boundary nodes on subdomain");
    for (i=0;i<nproc;i++){
      fprintf (out,"\n\n domain %ld",i);
      for (j=0;j<nbnd[i];j++){
        fprintf (out,"\n %6ld   %6ld",j,icnbnmas[i][j]);
      }
    }
    
    fprintf (out,"\n\n\n coarse-local correspondence \n");
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n coarse node %6ld, number of nodes %2ld",i,icmultip[i]);
    }
  }
  fprintf (out,"\n\n\n coarse numbers of interface/boundary nodes on subdomain");
  for (i=0;i<nbn;i++){
    fprintf (out,"\n %6ld   %6ld",i,icnbn[i]);
  }
}



/**
   function assembles coarse - global numbers map
   
   the following arrays are assembled:
   
   gnbncn - (M) global numbers of boundary/interface nodes of the coarse node
   sid - (M) subdomain id of interface/boundary nodes of the coarse node
   
   @param out - output file for auxiliary print
   @param proc_name - processor name
   
   JK, 30.7.2009
*/
void partop::node_coarse_global_map (FILE *out,char *proc_name)
{
  long i,j,ln,cn;
  
  if (myrank==0){
    
    if (tnbn==0){
      par_print_err(myrank,proc_name,"total number of boundary/interface nodes (tnbn) is zero",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::assemble_nbnd_nind or partop::update_multip or partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
    }
    if (icmultip==NULL){
      par_print_err(myrank,proc_name,"array icmultip is not assembled",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::node_coarse_numbers has to be called first",__FILE__, __LINE__, __func__);
    }
    if (icnbnmas==NULL){
      par_print_err(myrank,proc_name,"array icnbnmas is not assembled",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::node_coarse_numbers has to be called first",__FILE__, __LINE__, __func__);
    }
    if (gnbn==NULL){
      par_print_err(myrank,proc_name,"array gnbn is not assembled",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function partop::node_global_numbers has to be called first",__FILE__, __LINE__, __func__);
    }

    //  global numbers of boundary/interface nodes of the coarse node
    if (gnbncn!=NULL){
      for (i=0;i<tnbn;i++){
        delete [] gnbncn[i];
      }
      delete [] gnbncn;
    }
    gnbncn = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      gnbncn[i]=new long [icmultip[i]];
    }
    //  subdomain id of interface/boundary nodes of the coarse node
    if (sid!=NULL){
      for (i=0;i<tnbn;i++){
        delete [] sid[i];
      }
      delete [] sid;
    }
    sid = new long* [tnbn];
    for (i=0;i<tnbn;i++){
      sid[i] = new long [icmultip[i]];
    }
    
    //  loop over the number of subdomains
    for (i=0;i<tnbn;i++){
      icmultip[i]=0;
    }
    
    //  loop over the number of subdomains
    for (i=0;i<nproc;i++){
      //  loop over the number of boundary/interface nodes
      for (j=0;j<nbnd[i];j++){
        //  coarse number
        cn=icnbnmas[i][j];
        //  global number
        ln=gnbn[i][j];
        
        //  global number of boundary/interface node
        gnbncn[cn][icmultip[cn]]=ln;
        //  subdomain id if boundary/interface node
        sid[cn][icmultip[cn]]=i;
        icmultip[cn]++;
      }
    }
  }
  
  // *******************
  //  auxiliary output
  // *******************
  if (myrank==0){
    fprintf (out,"\n\n\n coarse-global correspondence \n");
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n coarse node %6ld, number of nodes %2ld",i,icmultip[i]);
      for (j=0;j<icmultip[i];j++){
        fprintf (out,"   node %ld, global number %6ld, subdomain id %4ld",j,gnbncn[i][j],sid[i][j]);
      }
    }
    fprintf (out,"\n\n\n coarse-global correspondence \n");
    for (i=0;i<tnbn;i++){
      fprintf (out,"\n coarse node %6ld",i);
      for (j=0;j<icmultip[i];j++){
        fprintf (out,"   ggn %6ld, sub.%4ld",gnbncn[i][j],sid[i][j]);
      }
    }
  }
}




/**
   function generates the Schur complement ordering
   it takes into account the coupled DOFs
   
   ordering of subdomains is performed
   ordering of the coarse problem is in selnodes::schur_ordering
   
   function determines:
   ndof - the number of all DOFs (unknowns)
   nidof - the number of internal DOFs (unknowns)
   nbdof - the number of boundary/interface DOFs (unknowns)
   
   @param top - pointer to general topology
   @param domproc - array containing domain-processor correspondence
   @param out - output file
   @param proc_name - processor name
   
   JK, 29.7.2009
*/
long partop::schur_ordering (gtopology *top,long *domproc,FILE *out,char *proc_name)
{
  long i,j,k,l,m,ndofn,nid,buffsize;
  long *aux,**auxdofind,*buff;
  MPI_Status stat;
  
  switch (md){
  case bound_nodes:{
    
    // ********************************
    //  ordering of internal unknowns
    // ********************************
    ndof=1;
    for (i=0;i<nin;i++){
      //  node id (in local ordering)
      nid=lnindom[i];
      //  number of DOFs at required node
      ndofn = top->give_ndofn (nid);
      
      for (j=0;j<ndofn;j++){
        k=top->give_dof (nid,j);
        if (k>0){
          top->save_dof (nid,j,ndof);
          ndof++;
        }
      }
    }
    //  number of internal DOFs (unknowns)
    nidof = ndof-1;
    
    
    // ********************************
    //  ordering of boundary unknowns
    // ********************************
    for (i=0;i<nbn;i++){
      //  node id (in local ordering)
      nid=lnbndom[i];
      //  number of DOFs at required node
      ndofn = top->give_ndofn (nid);
      
      for (j=0;j<ndofn;j++){
        k=top->give_dof (nid,j);
        if (k>0){
          top->save_dof (nid,j,ndof);
          ndof++;
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
    
    //  number of boundary DOFs (unknowns)
    nbdof = ndof-nidof;
    
    break;
  }
  case metis:
  case all_nodes:{
    
    // *********************************************
    //  distribution of the array dofind to slaves
    // *********************************************
    if (maxndof==0){
      par_print_err(myrank,proc_name,"maximum number of DOFs on one subdomain is equal to zero",__FILE__, __LINE__, __func__);
      par_print_err(myrank,proc_name,"function numbers_of_all_dofs_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
    }
    
    buffsize=maxndof;
    buff = new long [buffsize];
    for (i=0;i<buffsize;i++){
      buff[i]=0;
    }
    
    if (myrank==0){
      
      if (nnsd==NULL){
        par_print_err(myrank,proc_name,"array nnsd is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::numbers_of_all_nodes_on_subdomains has to be called first",__FILE__, __LINE__, __func__);
      }
      if (allnodes==NULL){
        par_print_err(myrank,proc_name,"array allnodes is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::assemble_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      if (ndofnmas==NULL){
        par_print_err(myrank,proc_name,"array ndofnmas is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::ndofn_master has to be called first",__FILE__, __LINE__, __func__);
      }
      if (dofind==NULL){
        par_print_err(myrank,proc_name,"array dofind is not assembled",__FILE__, __LINE__, __func__);
        par_print_err(myrank,proc_name,"function partop::update_multip has to be called first",__FILE__, __LINE__, __func__);
      }
      
      
      //  loop over the number of slaves
      for (i=1;i<nproc;i++){
        //  subdomain id
        k=domproc[i];
        
        //  index in the array buff
        m=0;
        //  loop over the number of all nodes on the k-th subdomain
        for (j=0;j<nnsd[k];j++){
          //  global node number
          nid=allnodes[k][j];
          //printf ("  nid %ld ",nid);
          //  the number of DOFs in node
          ndofn=ndofnmas[nid];
          //  loop over the number of DOFs in node
          for (l=0;l<ndofn;l++){
            buff[m]=dofind[nid][l];
            m++;
          }
        }
        //printf ("\n m %ld",m);
        
        MPI_Send (buff,buffsize,MPI_LONG,i,myrank,MPI_COMM_WORLD);
      }
      //  master part
      
      //  subdomain id
      k=domproc[0];
      
      //  index in the array buff
      m=0;
      //  loop over the number of all nodes on the k-th subdomain
      for (j=0;j<nnsd[k];j++){
        //  global node number
        nid=allnodes[k][j];
        //  the number of DOFs in node
        ndofn=ndofnmas[nid];
        //  loop over the number of DOFs in node
        for (l=0;l<ndofn;l++){
          buff[m]=dofind[nid][l];
          m++;
        }
      }
      
    }
    else{
      MPI_Recv (buff,buffsize,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    auxdofind = new long* [nn];
    for (i=0;i<nn;i++){
      auxdofind[i] = new long [top->give_ndofn(i)];
    }
    
    //  index in the array buff
    k=0;
    //  loop over the number of all nodes on subdomain
    for (i=0;i<nn;i++){
      //  the number of DOFs in node
      ndofn=top->give_ndofn(i);
      //  loop over the number of DOFs in node
      for (j=0;j<ndofn;j++){
        auxdofind[i][j]=buff[k];
        k++;
      }
    }
    
    delete [] buff;
    
    // ****************************************************
    //  end of distribution of the array dofind to slaves
    // ****************************************************
    
    //  searching of maximum code number indicator
    ndof=0;
    //  loop over the number of all nodes in the problem
    for (i=0;i<nn;i++){
      //  the number of DOFs in node
      ndofn=top->give_ndofn (i);
      //  loop over the number of DOFs in node
      for (j=0;j<ndofn;j++){
        k=top->give_dof (i,j);
        if (k>ndof)  ndof=k;
      }
    }
    ndof--;
    if (ndof<0)  ndof=0;
    aux = new long [ndof];
    for (i=0;i<ndof;i++){
      aux[i]=-1;
    }
    
    
    
    //  number of actual DOF
    ndof=1;
    
    // ******************************
    //  generation of internal DOFs
    // ******************************
    
    //  loop over the number of nodes on subdomains
    for (j=0;j<nn;j++){
      //  the number of DOFs
      ndofn=top->give_ndofn (j);
      //  loop over the number of DOFs in node
      for (l=0;l<ndofn;l++){
        if (auxdofind[j][l]==0){
          //  this is internal DOF
          
          k=top->give_dof (j,l);
          if (k<0)  continue;
          if (k==0)  continue;
          if (k==1){
            top->gnodes[j].cn[l]=ndof;  ndof++;
          }
          if (k>1){
            if (aux[k-2]==-1){
              top->gnodes[j].cn[l]=ndof;
              aux[k-2]=ndof;
              ndof++;
            }
            else{
              top->gnodes[j].cn[l]=aux[k-2];
            }
          }
        }
      }
    }
    
    //  the number of internal DOFs (unknowns)
    nidof=ndof-1;
    
    // ****************************************
    //  generation of boundary/interface DOFs
    // ****************************************
    
    //  loop over the number of all nodes on subdomains
    for (j=0;j<nn;j++){
      //  the number of DOFs
      ndofn=top->give_ndofn (j);
      //  loop over the number of DOFs in node
      for (l=0;l<ndofn;l++){
        if (auxdofind[j][l]==1){
          //  this is boundary/interface DOF
          
          k=top->give_dof (j,l);
          if (k<0)  continue;
          if (k==0)  continue;
          if (k==1){
            top->gnodes[j].cn[l]=ndof;  ndof++;
          }
          if (k>1){
            if (aux[k-2]==-1){
              top->gnodes[j].cn[l]=ndof;
              aux[k-2]=ndof;
              ndof++;
            }
            else{
              top->gnodes[j].cn[l]=aux[k-2];
            }
          }
        }
      }
    }
    ndof--;
    
    //  number of boundary DOFs (unknowns)
    nbdof = ndof-nidof;
    
    if (aux != NULL)
      delete [] aux;
    
    for (i=0;i<nn;i++){
      delete [] auxdofind[i];
    }
    delete [] auxdofind;
    
    break;
  }
  default:{
    print_err("unknown type of mesh description is required", __FILE__, __LINE__, __func__);
  }
  }
  
  //  state of code numbers is changed
  //  code numbers are generated
  top->cnstate=1;
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n nidof   %6ld",nidof);
  fprintf (out,"\n\n nbdof   %6ld",nbdof);
  fprintf (out,"\n\n ndof    %6ld",ndof);
  fprintf (out,"\n\n\n kontrola kodovych cisel \n\n");
  //  loop over the number of all nodes on subdomains
  for (j=0;j<nn;j++){
    fprintf (out,"\n node %6ld  ",j);
    //  the number of DOFs
    ndofn=top->give_ndofn (j);
    //  loop over the number of DOFs in node
    for (l=0;l<ndofn;l++){
      fprintf (out,"  %6ld",top->gnodes[j].cn[l]);
    }
  }
  
  
  return ndof;
}

/**
   function prepares all data needed

   @param domproc - array containing domain-processor correspondence
   @param ltg - local to global correspondence
   @param top - pointer to the sequential general topology
   @param out - output file (used for auxiliary output)
   @param proc_name - processor name
   
   JK, 16.9.2009
   TKo, 20.5.2010
*/
void partop::prepare_data (long *domproc,long *ltg,gtopology *top,FILE *out,char *proc_name)
{
  long i, k, delta;

  initiation (top,ltg);
  numbers_of_all_nodes_on_subdomains (domproc,out);
  coupled_dofs (top,domproc,out);
  if ((md==bound_nodes) && nbcdof!=0)
  {
    fprintf(stdout, "\n Coupled dofs detected on boundary/interface nodes\n");
    fprintf(stdout, " Global node numbers will be regenerated\n");
    assemble_nbnd_nind(ltg,domproc,out,proc_name);
    delta = 0;
    for(i=0; i<ndom; i++)
      delta += nind[i];
    k = delta;
    bltg = new long[top->nn];
    for(i=0; i< top->nn; i++)
      bltg[i] = ltg[i];
    for(i=0; i< top->nn; i++)
    {
      if (ltg[i]<0)
      {
        ltg[i] = k;
        k++;
      }
      else
        ltg[i] += tnin;
    }
    md = all_nodes;
    initiation (top,ltg);    
    fprintf(out, "\n\n\n Coupled dofs on boundary/interface nodes were detected");
    fprintf(out, "\n global node numbers (ltg) has been regenerated");
    fprintf(out, "\n mesh description changed to all_nodes, delta=%ld\n", delta);
    for(i=0; i< top->nn; i++)
      fprintf(out, "%ld %ld %ld\n", i+1, bltg[i], ltg[i]);
  }

  assemble_multip (ltg,domproc,out,proc_name);
  numbers_of_all_dofs_on_subdomains (top,domproc,out);
  
  //  pridano kvuli vypoctu castecnych skalarnich soucinu








  long j,buffsize,*buff;
  buffsize=maxnbn;
  buff = new long [buffsize];


  // *******************************
  //  indication of boundary nodes 
  // *******************************
  //  local numbers of boundary nodes
  if (lnbn != NULL)
    delete [] lnbn;
  lnbn = new long [nbn];
  
  switch (md){
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
  
  j=0;
  for (i=0;i<nn;i++){
    if (ltg[i]==-1){
      lnin[j]=i;
      j++;
    }
  }

  MPI_Barrier (MPI_COMM_WORLD);
  


































  if (md==all_nodes || md==neg_bound_nodes){
    ndofn_master (top,domproc,out,proc_name);
    dofind_master (top,domproc,out,proc_name);
    update_multip (top,out,proc_name);
  }

  assemble_nbnd_nind (ltg,domproc,out,proc_name);
  assemble_nodmultip (ltg,domproc,out,proc_name);
  if (md==all_nodes || md==neg_bound_nodes){
    node_global_numbers (domproc,out,proc_name);
  }
  node_local_numbers (out);

  node_coarse_numbers (domproc,out,proc_name);
  
  if (md==all_nodes || md==neg_bound_nodes){
    node_coarse_global_map (out,proc_name);
  }
  
}




/**
   function assembles global code numbers on master processor 
   for using of class parcongrad
   
   function assembles array gcnd
   gcnd[i][j] = k - the j-th DOF on i-th subdomain has global code number k
   
   @param top - pointer to gtopology of subdomain
   @param domproc - domain-processor correspondence
   @param out - output file
   
   10.5.2007, JB
*/
void partop::assemble_gcnd(gtopology *top,long *domproc,FILE *out)
{
  long i,j,k,ii,m;
  long a,b,buffsize;
  long adr,ndofn;
  long *buff;
  long *gnn = NULL;
  long *pgcn = NULL;
  MPI_Status stat;
  
  if (md==bound_nodes){
    // creation of array allnodes
    
    buffsize = maxnbn;
    buff = new long [buffsize];
    
    for(i = 0; i < nbn; i++){
      buff[i] = lnbndom[i];
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
        allnodes[0][buff[k]]=icnbnmas[0][k];
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
          allnodes[j][buff[k]]=icnbnmas[j][k];
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
        
    //fprintf (stdout,"\n\n\n tndof je %ld\n\n\n",tndof);

    
    
    
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
