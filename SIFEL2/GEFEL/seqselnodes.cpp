#include "seqselnodes.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**
   constructor
   
   @param nd - number of subdomains
   @param nnsd - number of all nodes on subdomains
   @param jj - array containing list of selected nodes
               it contains nd,innsd[i] components

   @param nodmultip - array containing number of subdomains which share the nodes
   @param ggnbn - array containing global glued numbers of boundary/interface nodes
   @param icnbnmas - array containing coarse numbers of interface/boundary nodes
                  it contains nd,nbnd[i] components

   @param itnbn - total number of boundary/interface nodes
   @param iicmultip - array containing node multiplicity of boundary/interface nodes
   @param ilnbncn - local numbers of boundary/interface nodes appropriate to coarse node
   @param iggnbncn - global glued numbers of boundary/interface nodes appropriate to coarse node
   @param isid - subdomain id of interface/boundary nodes of the coarse node

   @param d - type of mesh description
   @param out - output stream
   @param mespr - message printing indicator
   
   jj[i][k]>-1 - the k-th node on the i-th subdomain is selected
   
   
   ns - number of subdomains is determined
   tnbn - total number of boundary nodes is determined
   md - mesh description is determined
   nsnmas - array of numbers of selected nodes on subdomains is assembled
   ggnsn - array of global glued numbers of selected nodes is assembled
   cnsnmas - array of coarse numbers of selected boundary/interface nodes is assembled
   icmultip - array of multiplicity of all boundary/interface nodes is assembled
   ggnbncn - local numbers of boundary/interface nodes appropriate to coarse node
   sid - subdomain id of interface/boundary nodes appropriate to coarse node

   JK, 14.9.2007, 5.7.2009
*/
seqselnodes::seqselnodes(long nd,long *nnsd,long **jj,
			 long **nodmultip,long **ggnbn,long **icnbnmas,
			 long itnbn,long *iicmultip,long **ilnbncn,long **iggnbncn,long **isid,
			 meshdescription d,FILE *out,long mespr)
{
  long i,j,k,ii,kk;
  
  //  number of subdomains
  ns=nd;
  
  //  total number of boundary/interafce nodes
  tnbn=itnbn;
  
  //  mesh description
  md = d;
  
  //  total number of selected nodes
  tnsn=0;
  
  //  total number of DOFs on selected nodes
  tndofsn=0;
  
  
  //  number of selected nodes
  nsnmas = new long [ns];
  for (i=0;i<ns;i++){
    nsnmas[i]=0;
    for (j=0;j<nnsd[i];j++){
      if (jj[i][j]>-1)
        nsnmas[i]++;
    }
  }
  
  for (i=0;i<ns;i++){
    if (mespr>0)
      fprintf (stdout,"\n subdomain %4ld, number of selected nodes %ld",i,nsnmas[i]);
    if (nsnmas[i]==0){
      fprintf (stderr,"\n\n wrong number of selected nodes in constructor (file %s, line %d)\n",__FILE__,__LINE__);
      //abort();
    }
  }
  
  //  local numbers of selected nodes
  lnsn = new long* [ns];
  for (i=0;i<ns;i++){
    lnsn[i]=new long [nsnmas[i]];
  }
  //  global glued numbers of selected nodes
  ggnsn = new long* [ns];
  for (i=0;i<ns;i++){
    ggnsn[i]=new long [nsnmas[i]];
  }
  //  coarse numbers of selected nodes
  cnsnmas = new long* [ns];
  for (i=0;i<ns;i++){
    cnsnmas[i]=new long [nsnmas[i]];
  }
  
  //  loop over the number of subdomain
  for (i=0;i<ns;i++){
    //  index of selected nodes
    //  k=nsnmas[i] at the end
    k=0;
    //  index of internal nodes
    //  ii=nind[i] (number of internal nodes on subdomain) at the end
    ii=0;
    //  index of boundary/interface nodes
    //  kk=nbnd[i] (number of boundary/interface nodes on subdomain) at the end
    kk=0;

    //  loop over all nodes on the i-th subdomain
    for (j=0;j<nnsd[i];j++){
      if (nodmultip[i][j]>1){
	//  boundary/interface node
	
	if (jj[i][j]>-1){
	  //  the node is selected
	  
	  //  local number of selected node
	  lnsn[i][k]=j;
	  //  global glued number of selected node
	  ggnsn[i][k]=ggnbn[i][kk];
	  //  coarse numbers of selected nodes
	  cnsnmas[i][k]=icnbnmas[i][kk];
	  k++;
	}
	kk++;
      }
      else{
	//  internal node
	
	if (jj[i][j]>-1){
	  //  the node is selected
	  
	  print_err("internal node has been selected", __FILE__, __LINE__, __func__);
	  
	  ii++;
	}
      }
    }
  }



  //  number of multiplicity of boundary/interface nodes
  icmultip = new long [tnbn];
  for (i=0;i<tnbn;i++){
    icmultip[i]=iicmultip[i];
  }
  //  local numbers of boundary/interface nodes appropriate to coarse node
  lnbncn = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    lnbncn[i]=new long [icmultip[i]];
  }
  //  global glued numbers of boundary/interface nodes appropriate to coarse node
  ggnbncn = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    ggnbncn[i]=new long [icmultip[i]];
  }
  //  subdomain id of interface/boundary nodes of the coarse node
  sid = new long* [tnbn];
  for (i=0;i<tnbn;i++){
    sid[i] = new long [icmultip[i]];
  }
  
  for (i=0;i<tnbn;i++){
    icmultip[i]=iicmultip[i];
  }
  for (i=0;i<tnbn;i++){
    for (j=0;j<icmultip[i];j++){
      lnbncn[i][j]=ilnbncn[i][j];
      ggnbncn[i][j]=iggnbncn[i][j];
      sid[i][j]=isid[i][j];
    }
  }
  
  
  //fprintf (out,"\n\n\n list of local and global/group numbers of selected nodes (lsnl, sncnbn)\n");
  fprintf (out,"\n\n\n list of global glued numbers of selected nodes (ggnsn)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain %4ld,  number of selected nodes is nsnmas %ld",i,nsnmas[i]);
    for (j=0;j<nsnmas[i];j++){
      fprintf (out,"\n global glued node number %6ld",ggnsn[i][j]+1);
      //fprintf (out,"\n local node number %6ld       global/group node number %6ld",lsnl[i][j],lsng[i][j]);
    }
  }
  

  //  number of multiplicity of selected boundary/interface nodes
  snicmultip = NULL;

  //  local numbers of boundary/interface nodes appropriate to selected coarse node
  snlnbncn = NULL;

  //  global glued numbers of selected boundary/interface nodes appropriate to coarse node
  snggnbncn = NULL;
  
  //  subdomain id of selected interface/boundary nodes of the coarse node
  snsid = NULL;
  
  //  array of numbers of DOFs on subdomains at selected nodes
  snndofmas = NULL;
  
  //  array of numbers of DOF at selected nodes
  snndofnmas = NULL;
  
  //  array of DOFs or indicators at selected nodes
  sndofmas = NULL;
  
  //  array of numbers of DOFs for selected nodes
  ndofnsn = NULL;
  
  //  code numbers at selected nodes on master
  codensn = NULL;
  
  //  code numbers on master
  cndofmas=NULL;
  
  //  type of FETI implementation
  //  this means no type is selected
  fetiimpl = no_impl;

  doffeti=NULL;
}

/**
   destructor
   
   JK, 14.9.2007
*/
seqselnodes::~seqselnodes()
{
  long i,j;

  //  global glued numbers of selected nodes
  if (ggnsn!=NULL){
    for (i=0;i<ns;i++){
      delete [] ggnsn[i];
    }
    delete [] ggnsn;
  }
  
  //  coarse numbers of selected nodes
  if (cnsnmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] cnsnmas[i];
    }
    delete [] cnsnmas;
  }
  
  //  number of multiplicity of boundary/interface nodes
  if (icmultip!=NULL)
    delete [] icmultip;

  //  global glued numbers of boundary/interface nodes appropriate to coarse node
  if (ggnbncn!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] ggnbncn[i];
    }
    delete [] ggnbncn;
  }
  
  //  subdomain id of interface/boundary nodes of the coarse node
  if (sid!=NULL){
    for (i=0;i<tnbn;i++){
      delete [] sid[i];
    }
    delete [] sid;
  }
  
  //  number of multiplicity of selected boundary/interface nodes
  if (snicmultip!=NULL)
    delete [] snicmultip;
  
  //  global glued numbers of boundary/interface nodes appropriate to selected coarse node
  if (snggnbncn!=NULL){
    for (i=0;i<tnsn;i++){
      delete [] snggnbncn[i];
    }
    delete [] snggnbncn;
  }
  
  //  subdomain id of interface/boundary nodes appropriate to coarse node
  if (snsid!=NULL){
    for (i=0;i<tnsn;i++){
      delete [] snsid[i];
    }
    delete [] snsid;
  }
  
  //  array of numbers of DOFs on subdomains at selected nodes
  if (snndofmas!=NULL)
    delete [] snndofmas;
  

  //  array of numbers of DOF at selected nodes
  if (snndofnmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] snndofnmas[i];
    }
    delete [] snndofnmas;
  }

  //  array of DOFs or indicators at selected nodes
  if (sndofmas!=NULL){
    for (i=0;i<ns;i++){
      for (j=0;j<nsnmas[i];j++){
	delete [] sndofmas[i][j];
      }
      delete [] sndofmas[i];
    }
    delete [] sndofmas;
  }

  //  array of code numbers of coarse nodes
  if (codensn!=NULL){
    for (i=0;i<tnsn;i++){
      delete [] codensn[i];
    }
    delete [] codensn;
  }

  //  array of numbers of DOFs for selected nodes
  if (ndofnsn!=NULL){
    delete [] ndofnsn;
  }

  //  code numbers on master
  if (cndofmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] cndofmas[i];
    }
    delete [] cndofmas;
  }

  
  //  number of selected nodes on subdomains (array nsnmas)
  if (nsnmas!=NULL)
    delete [] nsnmas;
  
}



/**
   function collects node numbers on master
   
   function assembles arrays:

   snicmultip - number of multiplicity of selected boundary/interface nodes
   snggnbncn - global glued numbers of selected boundary/interface nodes appropriate to coarse node
   snsid - subdomain id of selected interface/boundary nodes of the coarse node

   tnsn - total number of selected nodes
   
   @param out - output file
   
   5.7.2009, JK
*/
void seqselnodes::node_coarse_numbers (FILE *out)
{
  long i,j,k;
  long *av;

  
  av = new long [tnbn];
  for (i=0;i<tnbn;i++){
    av[i]=0;
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      av[cnsnmas[i][j]]++;
    }
  }
  
  //  total number of selected nodes
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
  
  //  local numbers of boundary/interface nodes appropriate to selected coarse node
  if (snlnbncn!=NULL){
    for (i=0;i<tnsn;i++){
      delete [] snlnbncn[i];
    }
    delete [] snlnbncn;
  }
  snlnbncn = new long* [tnsn];
  for (i=0;i<tnsn;i++){
    snlnbncn[i]=new long [snicmultip[i]];
  }

  //  global glued numbers of boundary/interface nodes appropriate to selected coarse node
  if (snggnbncn!=NULL){
    for (i=0;i<tnsn;i++){
      delete [] snggnbncn[i];
    }
    delete [] snggnbncn;
  }
  snggnbncn = new long* [tnsn];
  for (i=0;i<tnsn;i++){
    snggnbncn[i]=new long [snicmultip[i]];
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
      snicmultip[k]=icmultip[i];
      for (j=0;j<snicmultip[k];j++){
	snlnbncn[k][j]=lnbncn[i][j];
	snggnbncn[k][j]=ggnbncn[i][j];
	snsid[k][j]=sid[i][j];
      }
      k++;
    }
  }
  
  /*
 //  node sorting
  for (i=0;i<tnsn;i++){
    for (j=0;j<nodmultip[i];j++){
      min=ns;
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
  */




  
  for (i=0;i<tnbn;i++){
    delete [] lnbncn[i];
  }
  delete [] lnbncn;
  lnbncn = NULL;
  
  for (i=0;i<tnbn;i++){
    delete [] ggnbncn[i];
  }
  delete [] ggnbncn;
  ggnbncn = NULL;
  
  for (i=0;i<tnbn;i++){
    delete [] sid[i];
  }
  delete [] sid;
  sid = NULL;
  
  delete [] icmultip;
  icmultip = NULL;
  
  delete [] av;


  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n total number of selected nodes (tnsn)  %ld",tnsn);
  fprintf (out,"\n\n\n local numbers of boundary/interface nodes appropriate to selected coarse node (snlnbncn)\n");
  fprintf (out,"\n global glued numbers of boundary/interface nodes appropriate to selected coarse node (snggnbncn)\n");
  fprintf (out,"\n subdomain id of interface/boundary nodes appropriate to coarse node (snsid)\n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n sel. n. %6ld   snicmultip  %6ld\n",i,snicmultip[i]);
    for (j=0;j<snicmultip[i];j++){
      fprintf (out,"    %6ld  %6ld  %6ld\n",snlnbncn[i][j],snggnbncn[i][j],snsid[i][j]);
    }
  }
  
}


/**
   function searches for all possible DOFs in selected nodes
   prescribed DOFs are included
   
   function assembles the following array:
   snndofmas - array of numbers of DOFs on subdomains at selected nodes

   @param top - pointer to the general topology
   @param out - output file

   JK, 31.7.2007
*/
void seqselnodes::number_all_dofs (gtopology *top,FILE *out)
{
  long i,j,k;
  
  if (snndofmas!=NULL)
    delete [] snndofmas;
  snndofmas=new long [ns];
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    snndofmas[i]=0;
    
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      k=ggnsn[i][j];
      snndofmas[i]+=top->give_ndofn (k);
    }
  }
  
  fprintf (out,"\n\n numbers of DOFs on subdomains (snndofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld   %ld",i+1,snndofmas[i]);
  }

}


/**
   function collects numbers of DOFs on selected nodes
   numbers are stored on master in the array snndofnmas
   
   array snndofnmas is assembled
   
   @param top - pointer to the general topology
   @param out - output file
   
   JK, 31.7.2007
*/
void seqselnodes::ndofn_on_master (gtopology *top,FILE *out)
{
  long i,j,k;

  //  array of numbers of DOF at selected nodes
  if (snndofnmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] snndofnmas[i];
    }
    delete [] snndofnmas;
  }
  snndofnmas = new long* [ns];
  for (i=0;i<ns;i++){
    snndofnmas[i] = new long [nsnmas[i]];
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      k=ggnsn[i][j];
      snndofnmas[i][j]=top->give_ndofn (k);
    }
  }
  
  fprintf (out,"\n\n numbers of DOFs on selected nodes (snndofnmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<nsnmas[i];j++){
      fprintf (out,"\n selected node %6ld   %ld",j+1,snndofnmas[i][j]);
    }
  }

}


/**
   function assembles DOFs indicators at selected nodes
   
   function assembles the array
   sndofmas - array of DOFs or indicators at selected nodes

   @param top - pointer to the general topology
   @param out - output file
   
   JK, 14.9.2007
*/
void seqselnodes::dof_indicators (gtopology *top,FILE *out)
{
  long i,j,k,l;
  
  //  array of DOFs or indicators at selected nodes
  if (sndofmas!=NULL){
    for (i=0;i<ns;i++){
      for (j=0;j<nsnmas[i];j++){
	delete [] sndofmas[i][j];
      }
      delete [] sndofmas[i];
    }
    delete [] sndofmas;
  }
  sndofmas = new long** [ns];
  for (i=0;i<ns;i++){
    sndofmas[i] = new long* [nsnmas[i]];
    for (j=0;j<nsnmas[i];j++){
      sndofmas[i][j] = new long [snndofnmas[i][j]];
    }
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  loop over the number of DOFs on nodes
      for (k=0;k<snndofnmas[i][j];k++){
	//  global glued number of the selected node
	l=ggnsn[i][j];
	sndofmas[i][j][k]=top->give_dof (l,k);
      }
    }
  }
  
  fprintf (out,"\n\n DOFs indicators on master (sndofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<nsnmas[i];j++){
      fprintf (out,"\n selected node %6ld    ",j+1);
      for (k=0;k<snndofnmas[i][j];k++){
	fprintf (out,"  %ld",sndofmas[i][j][k]);
      }
    }
  }

}


/**
   function generates the Schur complement ordering
   
   array codensn and ndofnsn are assembled
   
   @param out - output file
   
   JK, 10.7.2009
*/
void seqselnodes::schur_ordering (long **dofind,FILE *out)
{
  long i,j,k,l,m,ndofn,g;
  long *av,*aux;
  
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
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  global glued number of selected node
      g=ggnsn[i][j];
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
  
  fprintf (out,"\n\n pole codensn \n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n %ld %ld",codensn[i][0],codensn[i][1]);
  }
  fprintf (out,"\n\n");

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

  fprintf (out,"\n\n pole codensn \n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n %ld %ld",codensn[i][0],codensn[i][1]);
  }
  fprintf (out,"\n\n");
  // **********************************************************
  //  end of generation of code numbers in the coarse problem
  // **********************************************************
  
  aux = new long [tndofsn];

  //  this array will be recalculated
  //  at this moment, it contains number of all
  for (i=0;i<ns;i++){
    snndofmas[i]=0;
  }
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
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
  if (cndofmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] cndofmas[i];
    }
    delete [] cndofmas;
  }
  cndofmas = new long* [ns];
  for (i=0;i<ns;i++){
    cndofmas[i] = new long [snndofmas[i]];
  }
  
  //  auxiliary array
  av = new long [ns];
  for (i=0;i<ns;i++){
    av[i]=0;
  }
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
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
	    cndofmas[i][av[i]]=m;
	    av[i]++;
	    aux[m-1]=1;
	  }
	}
      }
      
    }
  }
  delete [] av;
  delete [] aux;

  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n tndofsn - total number of DOFs in selected nodes  %ld",tndofsn);

  fprintf (out,"\n\n code numbers of Schur ordering (cndofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<snndofmas[i];j++){
      fprintf (out,"\n dof %6ld   %ld",j+1,cndofmas[i][j]);
    }
  }
  
}

/**
   function generates the Schur complement ordering
   
   array codensn and ndofnsn are assembled
   
   @param out - output file
   
   JK, 10.7.2009
*/
/*
void seqselnodes::schur_ordering_old (long **dofind,FILE *out)
{
  long i,j,k,l,m,ndofn,g;
  long *av,*aux;
  
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
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  global glued number of selected node
      g=ggnsn[i][j];
      //  coarse number of the selected node
      l=cnsnmas[i][j];
      //  number of DOFs on the selected node
      ndofn=snndofnmas[i][j];
      
      if (codensn[l]==NULL){
	codensn[l] = new long [ndofn];
	ndofnsn[l]=ndofn;
      }
      //  loo over the number of DOFs on node
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
  
  fprintf (out,"\n\n pole codensn \n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n %ld %ld",codensn[i][0],codensn[i][1]);
  }
  fprintf (out,"\n\n");

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

  fprintf (out,"\n\n pole codensn \n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n %ld %ld",codensn[i][0],codensn[i][1]);
  }
  fprintf (out,"\n\n");
  // **********************************************************
  //  end of generation of code numbers in the coarse problem
  // **********************************************************
  
  
  //  this array will be recalculated
  //  at this moment, it contains number of all
  for (i=0;i<ns;i++){
    snndofmas[i]=0;
  }
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  coarse number of the selected node
      l=cnsnmas[i][j];
      //  number of DOFs on the selected node
      ndofn=snndofnmas[i][j];
      //  loop over the number of DOFs on node
      for (k=0;k<ndofn;k++){
	if (codensn[l][k]>0)
	  snndofmas[i]++;
      }
      
    }
  }


  //  code numbers on master
  if (cnmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] cnmas[i];
    }
    delete [] cnmas;
  }
  cnmas = new long* [ns];
  for (i=0;i<ns;i++){
    cnmas[i] = new long [snndofmas[i]];
  }
  
  //  auxiliary array
  av = new long [ns];
  for (i=0;i<ns;i++){
    av[i]=0;
  }
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
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
	  cnmas[i][av[i]]=m;
	  av[i]++;
	}
      }
      
    }
  }
  delete [] av;

  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n tndofsn - total number of DOFs in selected nodes  %ld",tndofsn);

  fprintf (out,"\n\n code numbers of Schur ordering (cnmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<snndofmas[i];j++){
      fprintf (out,"\n dof %6ld   %ld",j+1,cnmas[i][j]);
    }
  }
  
}
*/



/**
   function generates the Schur complement ordering
   
   array codensn and ndofnsn are assembled
   
   @param out - output file
   
   JK, 10.7.2009
*/
 /*
void seqselnodes::schur_ordering_old_old (FILE *out)
{
  long i,j,k,l,ndofn;
  long *av;

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
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  coarse number of the selected node
      l=cnsnmas[i][j];
      //  number of DOFs on the selected node
      ndofn=snndofnmas[i][j];
      
      if (codensn[l]==NULL){
	codensn[l] = new long [ndofn];
	ndofnsn[l]=ndofn;
      }
      //  loo over the number of DOFs on node
      for (k=0;k<ndofn;k++){
	codensn[l][k]=sndofmas[i][j][k];
      }
    }
  }
  
  //  total number of DOFs on selected nodes
  tndofsn=1;
  //  loop over the number of selected nodes
  for (i=0;i<tnsn;i++){
    //  loop over the number of DOFs on selected node
    for (j=0;j<ndofnsn[i];j++){
      if (codensn[i][j]==1){
	codensn[i][j]=tndofsn;
	tndofsn++;
      }
    }
  }
  tndofsn--;
  
  
  //  code numbers on master
  if (cnmas!=NULL){
    for (i=0;i<ns;i++){
      delete [] cnmas[i];
    }
    delete [] cnmas;
  }
  cnmas = new long* [ns];
  for (i=0;i<ns;i++){
    cnmas[i] = new long [snndofmas[i]];
  }
  
  //  auxiliary array
  av = new long [ns];
  for (i=0;i<ns;i++){
    av[i]=0;
  }
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  coarse number of the selected node
      l=cnsnmas[i][j];
      //  number of DOFs on the selected node
      ndofn=snndofnmas[i][j];
      //  loo over the number of DOFs on node
      for (k=0;k<ndofn;k++){
	cnmas[i][av[i]]=codensn[l][k];
	av[i]++;
      }
      
    }
  }
  delete [] av;

  
  fprintf (out,"\n\n code numbers of Schur ordering (cnmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<snndofmas[i];j++){
      fprintf (out,"\n dof %6ld   %ld",j+1,cnmas[i][j]);
    }
  }
  
  
  //  odsud to bylo zakomentovano

  long i,j,k,l,m,ii,ndofn;
  long *aux;
  
  //  determination of numbers of DOFs at selected nodes
  if (ndofnsn!=NULL)
    delete [] ndofnsn;
  ndofnsn = new long [tnsn];
  for (i=0;i<tnsn;i++){
    ndofnsn[i]=0;
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<nsnmas[i];j++){
      k=ndofnmas[i][j];
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
  
  
  //  code number indicators
  cnm = new long* [tnsn];
  for (i=0;i<tnsn;i++){
    cnm[i] = new long [ndofnsn[i]];
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<nsnmas[i];j++){
      l=gnn[i][j];
      for (k=0;k<ndofnsn[l];k++){
	cnm[l][k]=dofmas[i][j][k];
      }
    }
  }
  
  
  
  //  added 4.5.2009
  //  searching of maximum code number
  tndof=0;
  for (i=0;i<tnsn;i++){
    for (j=0;j<ndofnsn[i];j++){
      if (cnm[i][j]>tndof)
	tndof=cnm[i][j];
    }
  }
  tndof--;
  if (tndof<0)  tndof=0;
  aux = new long [tndof];
  for (i=0;i<tndof;i++){
    aux[i]=-1;
  }
  
  
  //  code number generation
  tndof=1;
  for (i=0;i<tnsn;i++){
    for (j=0;j<ndofnsn[i];j++){
      k=cnm[i][j];
      if (k<0)  continue;
      if (k==0)  continue;
      if (k==1){
	cnm[i][j]=tndof;
	tndof++;
      }
      if (k>1){
	if (aux[k-2]==-1){
	  cnm[i][j]=tndof;
	  aux[k-2]=tndof;
	  tndof++;
	}
	else{
	  cnm[i][j]=aux[k-2];
	}
      }
    }
  }
  tndof--;
  
  delete [] aux;
  
  //  computation of real number of DOFs on subdomains
  //  supports are not included now
  for (i=0;i<ns;i++){
    ndofdom[i]=0;
    for (j=0;j<nsnmas[i];j++){
      l=gnn[i][j];
      for (k=0;k<ndofnsn[l];k++){
	if (cnm[l][k]>0)
	  ndofdom[i]++;
      }
    }
  }
  
  cnmas = new long *[ns];
  for (i=0;i<ns;i++){
    cnmas[i] = new long [ndofdom[i]];
  }
  
  for (i=0;i<ns;i++){
    m=0;
    for (j=0;j<nsnmas[i];j++){
      l=gnn[i][j];
      for (k=0;k<ndofnsn[l];k++){
	if (cnm[l][k]>0){
	  cnmas[i][m]=cnm[l][k];
	  m++;
	}
      }
    }
  }
  
  
  //  zde by se dalo smazat pole cnm
  
  
  
  
  
  //  number of contributions to the coarse problem
  //  number of DOFs which contribute to the coarse problem
  for (ii=0;ii<ns;ii++){
    ndofdom[ii]=0;
    for (i=0;i<nsnmas[ii];i++){
      l=lsngg[ii][i];
      ndofn=top->give_ndofn (l);
      for (k=0;k<ndofn;k++){
	if (top->give_dof (l,k)>0)
	  ndofdom[ii]++;
      }
    }
  }
  
  //  list of code numbers which are extracted from domain in the FETI method
  //  some code numbers are positive and some of them are negative
  //  it depends on number of subdomain
  //  these code numbers are used on subdomains, master processor contains corresponding array with coarse code numbers
  
  if (ldof!=NULL){
    for (i=0;i<ns;i++){
      delete [] ldof[i];
    }
    delete [] ldof;
  }
  ldof = new long* [ns];
  for (i=0;i<ns;i++){
    ldof[i] = new long [ndofdom[i]];
  }
  
  for (ii=0;ii<ns;ii++){
    m=0;
    for (i=0;i<nsnmas[ii];i++){
      l=lsngg[ii][i];
      ndofn=top->give_ndofn (l);
      for (j=0;j<ndofn;j++){
	k=top->give_dof (l,j);
	if (k>0){
	  ldof[ii][m]=k;
	  m++;
	}
      }
    }
  }
  

  fprintf (out,"\n array ldof");
  for (ii=0;ii<ns;ii++){
    for (i=0;i<ndofdom[ii];i++){
      fprintf (out,"\n %6ld %6ld   %ld",ii+1,i+1,ldof[ii][i]);
    }
  }
  
  fprintf (out,"\n\n actual numbers of DOFs on subdomains (ndofdom)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld   %ld",i+1,ndofdom[i]);
  }
  
  fprintf (out,"\n\n code numbers of Schur ordering (cnmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<ndofdom[i];j++){
      fprintf (out,"\n dof %6ld   %ld",j+1,cnmas[i][j]);
    }
  }



}
 */

/**
   function prepares all data for the Schur complement method
   
   @param gt - pointer to the general topology
   @param out - output file
   
   JK, 21.5.2009
*/
void seqselnodes::prepare_schur (gtopology *top,FILE *out)
{
  node_coarse_numbers (out);
  number_all_dofs (top,out);
  ndofn_on_master (top,out);
  dof_indicators (top,out);
  schur_ordering (top->stop->dofind,out);
}































/**
   function computes multiplicity of selected nodes
   
   function assembles array nodmultip
   
   JK, 14.9.2007
*/
/*
void seqselnodes::node_multiplicity (FILE *out)
{
  long i,j;
  
  if (nodmultip!=NULL)
    delete [] nodmultip;
  nodmultip = new long [tnsn];
  for (i=0;i<tnsn;i++){
    nodmultip[i]=0;
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<nsnmas[i];j++){
      nodmultip[gnn[i][j]]++;
    }
  }
  
  fprintf (out,"\n\n node multiplicity (nodmultip)\n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n selected node %6ld    %ld",i+1,nodmultip[i]);
  }

}
*/



/**
   function assembles lists of connected nodes to coarse nodes
   and numbers of subdomains which contain connected nodes
   
   arrays ljn and lsn are assembled
   
   JK, 14.9.2007
*/
     /*
void seqselnodes::group_local_nodes (FILE *out)
{
  long i,j,k,m,n,min;
  
  ljn = new long* [tnsn];
  lsn = new long* [tnsn];
  for (i=0;i<tnsn;i++){
    ljn[i] = new long [nodmultip[i]];
    lsn[i] = new long [nodmultip[i]];
    nodmultip[i]=0;
  }
  
  //  list of local node numbers of connected nodes
  //  list of numbers of subdomains which contains connected nodes
  for (i=0;i<ns;i++){
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
      min=ns;
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
  
  fprintf (out,"\n\n list of local node numbers of connected nodes (ljn)\n");
  fprintf (out,"\n\n list of numbers of subdomains which contains connected nodes (lsn)\n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n coarse node %6ld   multip %6ld   ",i+1,nodmultip[i]);
    for (j=0;j<nodmultip[i];j++){
      fprintf (out,"    %ld %ld",ljn[i][j],lsn[i][j]);
    }
  }

}
     */

/**
   function assembles indicators of code numbers and then generates
   code numbers
   
   JK, 14.9.2007
*/
      /*
void seqselnodes::dof_feti (FILE *out)
{
  long i,j,k,l,m;
  
  //  determination of numbers of DOFs at selected nodes
  if (ndofnsn!=NULL)
    delete [] ndofnsn;
  ndofnsn = new long [tnsn];
  for (i=0;i<tnsn;i++){
    ndofnsn[i]=0;
  }
  
  for (i=0;i<ns;i++){
    for (j=0;j<nsnmas[i];j++){
      k=ndofnmas[i][j];
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
	doffeti[i][j][k]=dofmas[m][l][k];
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
      */

/**
   function determines the multiplicity of unknowns
   the multiplicity of unknown is obtained from the multiplicity of nodes

   array dofmultip is allocated and assembled

   @param out - output file
   
   JB
*/
/*
void seqselnodes::dof_multiplicity (FILE *out)
{
  long i,j,k,n;
  long ndofn;
  long nm;
  
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

  //fprintf (out,"\n\n\n Multiplicity of coarse DOF n = %ld\n\n",n); 
  //for(i = 0; i < tndof; i++){
  //fprintf (out,"coarse DOF %ld has multiplicity %ld\n",i+1,dofmultip[i]); 
  //}
  
  //ldofmultip = new long*[ns];
  //for(i = 0; i < ns; i++){
  //ldofmultip[i] = new long[ndofdom[i]];
  //for(j = 0; j < ndofdom[i]; j++){
  //ldofmultip[i][j] =  dofmultip[cnmas[i][j]-1];
  //}
  //}
  
  //fprintf (out,"\n\n\n local multiplicity of coarse DOF\n\n"); 
  //for(i = 0; i < ns; i++){
  //fprintf (out,"domain %ld\n",i+1);
  //for(j = 0; j < ndofdom[i]; j++){
  //fprintf (out,"%ld %ld %ld\n",j+1,cnmas[i][j],ldofmultip[i][j]); 
  //}
  //}
  

}
*/

/**
   function determines number of contributions to coarse problem
   from particular subdomains
   
   maxndof is rewritten
   ndofdom is rewritten
   
   array ncn is assembled
   
   JK, 14.9.2007
*/
	/*
void seqselnodes::number_contrib (FILE *out)
{
  long i,j,k,l,m,ii;
  
  if (ndofdom!=NULL)
    delete [] ndofdom;
  ndofdom = new long [ns];
  for (i=0;i<ns;i++){
    ndofdom[i]=0;
  }
  
  if (ncnmas!=NULL)
    delete [] ncnmas;
  ncnmas = new long [ns];
  for (i=0;i<ns;i++){
    ncnmas[i]=0;
  }
  
  //  number of contributing DOFs from subdomains to coarse problem
  for (i=0;i<ns;i++){
    for (j=0;j<nsnmas[i];j++){
      ii=gnn[i][j];
      for (k=0;k<nodmultip[ii];k++){
	if (lsn[ii][k]==i){
	  if (k==0){
	    for (l=0;l<nodmultip[ii]-1;l++){
	      ncnmas[i]++;
	      for (m=0;m<ndofnsn[ii];m++){
		if (doffeti[ii][l][m]>0)
		  ndofdom[i]++;
	      }
	    }
	  }
	  else{
	    ncnmas[i]++;
	    for (m=0;m<ndofnsn[ii];m++){
	      if (doffeti[ii][k-1][m]>0)
		ndofdom[i]++;
	    }
	  }
	  break;
	}
      }
    }
  }
    
  fprintf (out,"\n\n number of contributing DOFs to coarse problem (ndofdom)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld   %ld",i+1,ndofdom[i]);
  }
  fprintf (out,"\n\n number of nodes which contribute to coarse FETI problem (ncn)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n domain %ld   %ld",i,ncnmas[i]);
  }
  
}

	*/


/**
   function assmebles nodes contributing to the coarse proble in the FETI method
      
   JK, 14.9.2007
*/
	 /*
void seqselnodes::contrib_dofs (gtopology *top,FILE *out)
{
  long i,j,k,l,m,ii,jj,kk,ndofn,min;
  long *buff;
  
  //ndofdom = new long [ns];
  ldof = new long* [ns];

  //  array of local numbers of nodes contributing to the coarse problem
  for (i=0;i<ns;i++){
    
    buff = new long [ncnmas[i]];
    for (j=0;j<ncnmas[i];j++){
      buff[j]=0;
    }
    
    m=0;
    for (j=0;j<tnsn;j++){
      for (k=0;k<nodmultip[j];k++){
	if (lsn[j][k]==i){
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
    
    
    for (ii=0;ii<ncnmas[i];ii++){
      if (buff[ii]<0)
	min=0-buff[ii]-1;
      else
	min=buff[ii];
      
      l=ii;
      for (j=ii+1;j<ncnmas[i];j++){
	if (buff[j]<0)
	  k=0-buff[j]-1;
	else
	  k=buff[j];
	
	if (min>k){
	  l=j;
	  min=k;
	}
      }
      
      
      k=buff[ii];
      buff[ii]=buff[l];
      buff[l]=k;
    }
    
    fprintf (out,"\n\n kontrola bufferu na vsech uzlech");
    for (j=0;j<ncnmas[i];j++){
      fprintf (out,"\n %ld   %ld",j+1,buff[j]);
    }
    //fprintf (out,"\n konec kontroly\n");
    //fprintf (out,"\n\n kontrola lsnl na vsech uzlech");
    //for (i=0;i<nsn;i++){
    //fprintf (out,"\n %ld   %ld",i+1,lsnl[i]);
    //}
    //fprintf (out,"\n konec kontroly\n");
    
    
    
    //  number of contributions to the coarse problem
    //  number of DOFs which contribute to the coarse problem
    //ndofdom[i]=0;
    //for (ii=0;ii<ncnmas[i];ii++){
    //j=buff[ii];
    //if (j<0){
    //j=0-j-1;
    //}
      //l=lsng[i][j];
      //l=lsnl[i][j];
      //ndofn=top->give_ndofn (l);
      //for (k=0;k<ndofn;k++){
    //if (top->give_dof (l,k)>0)
    //ndofdom[i]++;
    //}
    //}
    
    //  list of code numbers which are extracted from domain in the FETI method
    //  some code numbers are positive and some of them are negative
    //  it depends on number of subdomain
    //  these code numbers are used on subdomains, master processor contains corresponding array with coarse code numbers
    
    ldof[i] = new long [ndofdom[i]];
    
    m=0;
    for (kk=0;kk<ncnmas[i];kk++){
      ii=buff[kk];
      if (ii<0){
	ii=0-ii-1;
	//jj=lsng[i][ii];
	jj=lsnl[i][ii];
	ndofn=top->give_ndofn (jj);
	for (j=0;j<ndofn;j++){
	  l=top->give_dof (jj,j);
	  if (l>0){
	    ldof[i][m]=0-l;
	    m++;
	  }
	}
      }
      else{
	//jj=lsng[i][ii];
	jj=lsnl[i][ii];
	ndofn=top->give_ndofn (jj);
	for (j=0;j<ndofn;j++){
	  l=top->give_dof (jj,j);
	  if (l>0){
	    ldof[i][m]=l;
	    m++;
	  }
	}
      }
    }
    
    delete [] buff;
  }

  
  cnmas = new long* [ns];
  for (i=0;i<ns;i++){
    cnmas[i] = new long [ndofdom[i]];
  }
  
  for (i=0;i<ns;i++){
    l=0;
    for (j=0;j<nsnmas[i];j++){
      m=gnn[i][j];
      for (k=0;k<nodmultip[m];k++){
	if (lsn[m][k]==i){
	  if (k==0){
	    for (ii=0;ii<nodmultip[m]-1;ii++){
	      for (jj=0;jj<ndofnsn[m];jj++){
		if (doffeti[m][ii][jj]>0){
		  cnmas[i][l]=doffeti[m][ii][jj];
		  l++;
		}
	      }
	    }
	  }
	  else{
	    for (jj=0;jj<ndofnsn[m];jj++){
	      if (doffeti[m][k-1][jj]>0){
		cnmas[i][l]=doffeti[m][k-1][jj];
		l++;
	      }
	    }
	  }
	  break;
	}
      }
    }
  }
  
  fprintf (out,"\n array ldof");
  for (i=0;i<ns;i++){
    fprintf (out,"\n\n %ld   number of DOFs which contribute to coarse FETI problem %ld",i,ndofdom[i]);
    for (j=0;j<ndofdom[i];j++){
      fprintf (out,"\n %6ld   %ld",j+1,ldof[i][j]);
    }
  }
  
  fprintf (out,"\n\n code numbers on master (cnmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld (number of DOFs %ld)",i+1,ndofdom[i]);
    for (j=0;j<ndofdom[i];j++){
      fprintf (out,"\n %6ld   %ld",j+1,cnmas[i][j]);
    }
  }

}
	 */






/**
   function defines type of implementation of the FETI method
   possible types can be found in the file galias.h
   
   @param fi - type of the FETI implementation
   
   JK, 29.9.2009
*/
void seqselnodes::define_feti_implementation (fetiimplem fi)
{
  fetiimpl = fi;
}

/**
   function searches for all possible contributing DOFs in selected nodes
   prescribed DOFs are included
   
   if a node is shared by n subdomains, it may contribute n-1 times
   this strategy is used in the FETI method
   in the case of the Schur complement method, each node contributes only once
   
   function assembles the following array:
   snndofmas - array of numbers of contributing DOFs on subdomains in selected nodes

   @param top - pointer to the general topology
   @param out - output file

   JK, 29.9.2009
*/
void seqselnodes::number_contrib_dofs (gtopology *top,FILE *out)
{
  long i,j,k,sdid,ggn,ndofn,mult,cn;
  
  if (snicmultip==NULL)
    print_err("array snicmultip has not been allocated", __FILE__, __LINE__, __func__);
  if (snsid==NULL)
    print_err("array snsid has not been allocated", __FILE__, __LINE__, __func__);
  if (snggnbncn==NULL)
    print_err("array snggnbncn has not been allocated", __FILE__, __LINE__, __func__);
    
    
  if (snndofmas!=NULL)
    delete [] snndofmas;
  snndofmas=new long [ns];
  for (i=0;i<ns;i++){
    snndofmas[i]=0;
  }
  
  // ****************************************************************
  //  number of contributing DOFs from subdomains to coarse problem
  // ****************************************************************
  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	//  global glued number of node shared by the coarse node
	ggn=snggnbncn[i][j];
	//  number of DOFs in the selected node
	ndofn = top->give_ndofn (ggn);
	if (j==0){
	  //  subdomain with the least number
	  //  master subdomain
	  snndofmas[sdid]+=ndofn*(snicmultip[i]-1);
	}
	else{
	  //  subdomain with the non-least number
	  //  slave subdomains
	  snndofmas[sdid]+=ndofn;
	}
      }
    }
    
    break;
  }
  case redundant:{
    
    //  loop over the number of subdomains
    for (i=0;i<ns;i++){
      snndofmas[i]=0;
      
      //  loop over the number of selected nodes
      for (j=0;j<nsnmas[i];j++){
	//  coarse number of selected node
	cn = cnsnmas[i][j];
	//  multiplicity of the selected node
	mult = snicmultip[cn];
	//  global glued number of the selected node
	ggn = ggnsn[i][j];      
	//  number of DOFs in the selected node
	ndofn = top->give_ndofn (ggn);
	
	//  loop over the node multiplicity
	for (k=0;k<mult-1;k++){
	  snndofmas[i]+=ndofn;
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }
  




  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n numbers of contributing DOFs on subdomains (snndofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld   %ld",i+1,snndofmas[i]);
  }
  
}


/**
   function assembles indicators of code numbers and then generates
   code numbers
   
   @param top - pointer to the general topology
   @param out - output file
   
   JK, 14.9.2007, revised 29.9.2009
*/
void seqselnodes::dof_feti (gtopology *top,FILE *out)
{
  long i,j,k,l,ggn,ggnj,ggnk,ndofn,dofj,dofk;
  
  //  determination of numbers of DOFs at selected nodes
  if (ndofnsn!=NULL)
    delete [] ndofnsn;
  ndofnsn = new long [tnsn];
  for (i=0;i<tnsn;i++){
    ndofnsn[i]=0;
  }
  
  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{

    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  global glued number of node shared by the coarse node
      ggn=snggnbncn[i][0];
      //  number of DOFs in the selected node
      ndofn = top->give_ndofn (ggn);
      
      //ndofnsn[i]=ndofn*(snicmultip[i]-1);
      ndofnsn[i]=ndofn;
    }

    break;
  }
  case redundant:{

    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  global glued number of node shared by the coarse node
      ggn=snggnbncn[i][0];
      //  number of DOFs in the selected node
      ndofn = top->give_ndofn (ggn);
      
      //ndofnsn[i]=ndofn*snicmultip[i]*(snicmultip[i]-1)/2;
      ndofnsn[i]=ndofn;
    }
    
    break;
  }
  default:{
    print_err("unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }
  //  array ndofnsn is assembled
  
  
  //  array allocation
  doffeti = new long*** [tnsn];
  for (i=0;i<tnsn;i++){
    doffeti[i] = new long** [snicmultip[i]];
    for (j=0;j<snicmultip[i];j++){
      doffeti[i][j] = new long* [snicmultip[i]];
      for (k=0;k<snicmultip[i];k++){
	doffeti[i][j][k] = new long [ndofnsn[i]];
      }
    }
  }
  
  
  // ****************************************
  //  assembling of code numbers indicators
  // ****************************************
  
  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  global glued number of node shared by the coarse node
	ggnj=snggnbncn[i][j];
	for (k=0;k<snicmultip[i];k++){
	  if (k==j){
	    //  node has to be connected with different node
	    
	    //  loop over the number of DOFs in the selected node
	    for (l=0;l<ndofnsn[i];l++){
	      doffeti[i][j][k][l]=0;
	    }
	    continue;
	  }

	  //  global glued number of node shared by the coarse node
	  ggnk=snggnbncn[i][k];
	  
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    //  DOF in the j-th node connected to the coarse node
	    dofj=top->give_dof (ggnj,l);
	    //  DOF in the k-th node connected to the coarse node
	    dofk=top->give_dof (ggnk,l);
	    
	    if (dofj!=dofk)
	      print_err("incompatible DOFs in the j-th unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
	    
	    if (j==0 || k==0){
	      //  only the first row and column are assembled
	      doffeti[i][j][k][l]=dofj;
	    }else{
	      //  all remaining components are equal to zero
	      doffeti[i][j][k][l]=0;
	    }
	  }
	}
      }
    }
    
    break;
  }
  case redundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  global glued number of node shared by the coarse node
	ggnj=snggnbncn[i][j];
	
	//  loop over the node multiplicity
	for (k=0;k<snicmultip[i];k++){
	  if (k==j){
	    //  the node has to be connected to different node

	    //  loop over the number of DOFs in the selected node
	    for (l=0;l<ndofnsn[i];l++){
	      doffeti[i][j][k][l]=0;
	    }
	    continue;
	  }
	  
	  //  global glued number of node shared by the coarse node
	  ggnk=snggnbncn[i][k];
	  
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    //  DOF in the j-th node connected to the coarse node
	    dofj=top->give_dof (ggnj,l);
	    //  DOF in the k-th node connected to the coarse node
	    dofk=top->give_dof (ggnk,l);
	    
	    if (dofj!=dofk)
	      print_err("incompatible DOFs in the j-th unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
	    
	    doffeti[i][j][k][l]=dofj;
	  }
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }
  //  array doffeti is assembled
  

  // **************************
  //  code numbers generation
  // **************************
  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{
    tndofsn=1;
    //  loop over the number of selected nodes
    for (i=0;i<tnsn;i++){
      j=0;
      //  loop over the node multiplicity
      for (k=1;k<snicmultip[i];k++){
	//  loop over the number of DOFs in the node
	for (l=0;l<ndofnsn[i];l++){
	  if (doffeti[i][j][k][l]>0){
	    doffeti[i][j][k][l]=tndofsn;
	    doffeti[i][k][j][l]=tndofsn;
	    tndofsn++;
	  }
	}
      }
    }
    tndofsn--;
    
    break;
  }
  case redundant:{
    tndofsn=1;
    //  loop over the number of selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  loop over the node multiplicity
	for (k=j+1;k<snicmultip[i];k++){
	  //  loop over the number of DOFs in the node
	  for (l=0;l<ndofnsn[i];l++){
	    if (doffeti[i][j][k][l]>0){
	      doffeti[i][j][k][l]=tndofsn;
	      doffeti[i][k][j][l]=tndofsn;
	      tndofsn++;
	    }
	  }
	}
      }
    }
    tndofsn--;
    
    break;
  }
  default:{
    print_err("unknown type of FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }
  
  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n total number of DOFs in selected nodes (tndofsn)  %ld",tndofsn);
  fprintf (out,"\n\n code numbers for FETI method (doffeti)\n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n\n coarse node %6ld",i+1);
    for (j=0;j<snicmultip[i];j++){
      fprintf (out,"\n ");
      for (k=0;k<snicmultip[i];k++){
	fprintf (out,"   ");
	for (l=0;l<ndofnsn[i];l++){
	  fprintf (out," %5ld",doffeti[i][j][k][l]);
	}
      }
    }
  }
  
}

/**
   function assembles local number of DOFs contributing to the coarse problem in the FETI method
   
   @param top - pointer to the topology
   @param out - output file
   
   JK, 14.9.2007, revised 2.10.2009
*/
void seqselnodes::contrib_dofs_ln (gtopology *top,FILE *out)
{
  long i,j,k,l,ndofn,sdid,ggnj,dofj,sdidk;
  
  for (i=0;i<ns;i++){
    snndofmas[i]=0;
  }

  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	//  global glued number of node shared by the coarse node
	ggnj=snggnbncn[i][j];
	if (j==0){
	  //  node belonging to the subdomain with the least number
	  
	  //  loop over the node multiplicity
	  for (k=1;k<snicmultip[i];k++){
	    //  loop over the number of DOFs in the selected node
	    for (l=0;l<ndofnsn[i];l++){
	      if (doffeti[i][j][k][l]>0){
		snndofmas[sdid]++;
	      }
	    }
	  }
	}
	else{
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    if (doffeti[i][j][0][l]>0){
	      snndofmas[sdid]++;
	    }
	  }
	}
      }
    }
    
    break;
  }
  case redundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	//  global glued number of node shared by the coarse node
	ggnj=snggnbncn[i][j];
	//  number of DOFs in the selected node
	ndofn = top->give_ndofn (ggnj);
	//  loop over the node multiplicity
	for (k=0;k<snicmultip[i];k++){
	  if (k==j){
	    //  node has to be connected with different node
	    continue;
	  }
	  
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    if (doffeti[i][j][k][l]>0){
	      snndofmas[sdid]++;
	    }
	  }
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }
  
  
  
  //  list of local code numbers which contribute to the coarse problem
  lndofmas = new long* [ns];
  for (i=0;i<ns;i++){
    lndofmas[i] = new long [snndofmas[i]];
    snndofmas[i]=0;
  }

  

  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	//  global glued number of node shared by the coarse node
	ggnj=snggnbncn[i][j];
	//  number of DOFs in the selected node
	ndofn = top->give_ndofn (ggnj);
	if (j==0){
	  //  node belonging to the subdomain with the least number
	  
	  //  loop over the node multiplicity
	  for (k=1;k<snicmultip[i];k++){
	    //  loop over the number of DOFs in the selected node
	    for (l=0;l<ndofnsn[i];l++){
	      if (doffeti[i][j][k][l]>0){
		//  DOF in the j-th node connected to the coarse node
		dofj=top->give_dof (ggnj,l);
		lndofmas[sdid][snndofmas[sdid]]=0-dofj;
		snndofmas[sdid]++;
	      }
	    }
	  }
	}
	else{
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    if (doffeti[i][j][0][l]>0){
	      //  DOF in the j-th node connected to the coarse node
	      dofj=top->give_dof (ggnj,l);
	      lndofmas[sdid][snndofmas[sdid]]=dofj;
	      snndofmas[sdid]++;
	    }
	  }
	}
      }
    }
    
    break;
  }
  case redundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	//  global glued number of node shared by the coarse node
	ggnj=snggnbncn[i][j];
	//  number of DOFs in the selected node
	ndofn = top->give_ndofn (ggnj);
	//  loop over the node multiplicity
	for (k=0;k<snicmultip[i];k++){
	  if (k==j){
	    //  node has to be connected with different node
	    continue;
	  }
	  
	  //  id of subdomain which shares the node belonging to the coarse node
	  sdidk=snsid[i][k];

	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    //  DOF in the j-th node connected to the coarse node
	    dofj=top->give_dof (ggnj,l);

	    if (doffeti[i][j][k][l]>0){
	      if (sdid<sdidk){
		lndofmas[sdid][snndofmas[sdid]]=0-dofj;
	      }else{
		lndofmas[sdid][snndofmas[sdid]]=dofj;
	      }
	      snndofmas[sdid]++;
	    }
	  }
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }


  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n\n array of numbers of DOFs on subdomains at selected nodes (snndofmas)");
  fprintf (out,"\n\n\n list of code numbers which contribute to the coarse problem (lndofmas)");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain %6ld   snndofmas %6ld\n",i,snndofmas[i]);
    for (j=0;j<snndofmas[i];j++){
      fprintf (out,"   %ld",lndofmas[i][j]);
    }
  }
  
}

/**
   function assembles coarse code numbers of DOFs contributing to the coarse proble in the FETI method
   
   @param out - output file
   
   JK, 14.9.2007, revised 5.10.2009
*/
void seqselnodes::contrib_dofs_cn (FILE *out)
{
  long i,j,k,l,sdid;
  
  //  list of coarse code numbers which contribute to the coarse problem
  cndofmas = new long* [ns];
  for (i=0;i<ns;i++){
    cndofmas[i] = new long [snndofmas[i]];
    snndofmas[i]=0;
  }
  
  switch (fetiimpl){
  case boolean_matrices:{
    break;
  }
  case nonredundant:{

    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	if (j==0){
	  //  node belonging to the subdomain with the least number
	  
	  //  loop over the node multiplicity
	  for (k=1;k<snicmultip[i];k++){
	    //  loop over the number of DOFs in the selected node
	    for (l=0;l<ndofnsn[i];l++){
	      if (doffeti[i][j][k][l]>0){
		cndofmas[sdid][snndofmas[sdid]]=doffeti[i][j][k][l];
		snndofmas[sdid]++;
	      }
	    }
	  }
	}
	else{
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    if (doffeti[i][j][0][l]>0){
	      cndofmas[sdid][snndofmas[sdid]]=doffeti[i][j][0][l];
	      snndofmas[sdid]++;
	    }
	  }
	}
      }
    }
    
    break;
  }
  case redundant:{
    
    //  loop over the number of all selected nodes
    for (i=0;i<tnsn;i++){
      //  loop over the node multiplicity
      for (j=0;j<snicmultip[i];j++){
	//  id of subdomain which shares the node belonging to the coarse node
	sdid=snsid[i][j];
	//  loop over the node multiplicity
	for (k=0;k<snicmultip[i];k++){
	  if (k==j){
	    //  node has to be connected with different node
	    continue;
	  }
	  
	  //  loop over the number of DOFs in the selected node
	  for (l=0;l<ndofnsn[i];l++){
	    
	    if (doffeti[i][j][k][l]>0){
	      cndofmas[sdid][snndofmas[sdid]]=doffeti[i][j][k][l];
	      snndofmas[sdid]++;
	    }
	  }
	}
      }
    }
    
    break;
  }
  default:{
    print_err("unknown type of the FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }

  
  // *******************
  //  auxiliary output
  // *******************
  fprintf (out,"\n\n coarse code numbers on master (cndofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld (number of DOFs %ld)",i,snndofmas[i]);
    for (j=0;j<snndofmas[i];j++){
      fprintf (out,"\n %6ld   %ld",j,cndofmas[i][j]);
    }
  }
  
}









/**
   function prepares all data for the FETI method
   
   @param fi - type of FETI implementation, possibilities are in the file galias.h
   @param gt - pointer to the general topology
   @param out - output file
   
   JK, 18.11.2007
*/
long seqselnodes::prepare_feti (fetiimplem fi,gtopology *top,FILE *out)
{
  long ndof;
  
  define_feti_implementation (fi);
  node_coarse_numbers (out);
  number_contrib_dofs (top,out);
  dof_feti (top,out);
  top->cngen=1;
  ndof = top->codenum_generation (out);
  contrib_dofs_ln (top,out);
  contrib_dofs_cn (out);
  
  return ndof;

  //contrib_dofs_ln (top,out);
  //contrib_dofs_cn (out);

  /*
  //nodes_on_master (gt->nn,out);
  //node_multiplicity (out);
  number_all_dofs (gt,out);
  ndofn_on_master (gt,out);
  dof_indicators (gt,out);
  group_local_nodes (out);
  dof_feti (out);
  //dof_multiplicity (out);
  number_contrib (out);
  contrib_dofs (gt,out);
  dof_multiplicity (out);
  */
}



































































/**
   function assembles lists of connected nodes to coarse nodes
   and numbers of subdomains which contain connected nodes
   
   arrays ljn and lsn are assembled
   
   JK, 14.9.2007, revised 23.9.2009
*/
/*
void seqselnodes::coarse_local_nodes (FILE *out)
{
  long i,j,k,m,n,min;
  
  ljn = new long* [tnsn];
  lsn = new long* [tnsn];
  for (i=0;i<tnsn;i++){
    ljn[i] = new long [nodmultip[i]];
    lsn[i] = new long [nodmultip[i]];
    nodmultip[i]=0;
  }
  
  //  list of local node numbers of connected nodes
  //  list of numbers of subdomains which contains connected nodes
  for (i=0;i<ns;i++){
    for (j=0;j<nsndom[i];j++){
      k=gnn[i][j];
      ljn[k][nodmultip[k]]=j;
      lsn[k][nodmultip[k]]=i;
      nodmultip[k]++;
    }
  }
  
  
  //  node sorting
  for (i=0;i<tnsn;i++){
    for (j=0;j<nodmultip[i];j++){
      min=ns;
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
  
  fprintf (out,"\n\n list of local node numbers of connected nodes (ljn)\n");
  fprintf (out,"\n\n list of numbers of subdomains which contains connected nodes (lsn)\n");
  for (i=0;i<tnsn;i++){
    fprintf (out,"\n coarse node %6ld   multip %6ld   ",i+1,nodmultip[i]);
    for (j=0;j<nodmultip[i];j++){
      fprintf (out,"    %ld %ld",ljn[i][j],lsn[i][j]);
    }
  }

}
*/

/**
   function assembles indicators of code numbers and then generates
   code numbers
   
   @param out - output file
   
   JK, 14.9.2007, revised 29.9.2009
*/
/*
void seqselnodes::dof_feti (FILE *out)
{
  long i,j,k,l,m;
  
  //  determination of numbers of DOFs at selected nodes
  if (ndofnsn!=NULL)
    delete [] ndofnsn;
  ndofnsn = new long [tnsn];
  for (i=0;i<tnsn;i++){
    ndofnsn[i]=0;
  }
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  number of DOFs in the selected node
      k=ndofnmas[i][j];
      //  global glued number of the selected node
      l=ggnsn[i][j];
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
  //  array ndofnsn is assembled

  //  array allocation
  doffeti = new long** [tnsn];
  for (i=0;i<tnsn;i++){
    doffeti[i] = new long* [snicmultip[i]];
    for (j=0;j<snicmultip[i];j++){
      doffeti[i][j] = new long [ndofnsn[i]];
    }
  }
  
  //  assembling of code numbers indicators
  //  loop over the number of all selected nodes
  for (i=0;i<tnsn;i++){
    //  loop over the node multiplicity
    for (j=0;j<snicmultip[i];j++){
      //  local number of node shared by the coarse node
      //l=ljn[i][j];
      l=snlnbncn[i][j];
      //  id of subdomain which shares the node belonging to the coarse node
      //m=lsn[i][j];
      m=snsid[i][j];
      //  loop over the number of DOFs in the selected node
      for (k=0;k<ndofnsn[i];k++){
	doffeti[i][j][k]=dofmas[m][l][k];
      }
    }
  }
  //  array doffeti is assembled
  
  // **************************
  //  code numbers generation
  // **************************
  switch (fetiimpl){
  case nonredundant:{
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
    
    break;
  }
  case redundant:{
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
    
    break;
  }
  default:{
    print_err("unknown type of FETI implementation is required", __FILE__, __LINE__, __func__);
  }
  }
  
  
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

*/







/**
   function searches for all possible contributing DOFs in selected nodes
   prescribed DOFs are included
   
   if a node is shared by n subdomains, it contributes n-1 times
   this strategy is used in the FETI method
   in the case of the Schur complement method, each node contributes only once
   
   function assembles the following array:
   snndofmas - array of numbers of DOFs on subdomains at selected nodes

   @param top - pointer to the general topology
   @param out - output file

   JK, 22.9.2009
*/
/*
void seqselnodes::number_contrib_dofs (gtopology *top,FILE *out)
{
  long i,j,k;
  
  if (snicmultip==NULL)
    print_err("array snicmultip has not been allocated", __FILE__, __LINE__, __func__);
  if (ndofnsn==NULL)
    print_err("array ndofnsn has not been allocated", __FILE__, __LINE__, __func__);
  if (snsid==NULL)
    print_err("array snsid has not been allocated", __FILE__, __LINE__, __func__);
  if (doffeti==NULL)
    print_err("array doffeti has not been allocated", __FILE__, __LINE__, __func__);
    
    
  if (snndofmas!=NULL)
    delete [] snndofmas;
  snndofmas=new long [ns];
  
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    snndofmas[i]=0;
    
    //  loop over the number of selected nodes
    for (j=0;j<nsnmas[i];j++){
      //  coarse number of selected node
      cn = cnsnmas[i][j];
      //  multiplicity of the selected node
      mult = snicmultip[cn];
      //  global glued number of the selected node
      ggn = ggnsn[i][j];      
      //  number of DOFs in the selected node
      ndofn = top->give_ndofn (ggn);
      
      //  loop over the node multiplicity
      for (k=0;k<mult-1;k++){
	snndofmas[i]+=ndofn;
      }
    }
  }
  

  //  number of contributing DOFs from subdomains to coarse problem
  //  loop over the number of subdomains
  for (i=0;i<ns;i++){
    //  loop over the number of selected nodes on subdomains
    for (j=0;j<nsnmas[i];j++){
      //  global glued number of the selected node
      ii=ggnsn[i][j];
      //  coarse number of selected node
      cn = cnsnmas[i][j];
      //  loop over the multiplicity
      for (k=0;k<snicmultip[cn];k++){
	if (snsid[cn][k]==i){
	  if (k==0){
	    for (l=0;l<snicmultip[cn]-1;l++){
	      ncnmas[i]++;
	      for (m=0;m<ndofnsn[ii];m++){
		if (doffeti[ii][l][m]>0)
		  snndofmas[i]++;
	      }
	    }
	  }
	  else{
	    ncnmas[i]++;
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
 




  fprintf (out,"\n\n numbers of contributing DOFs on subdomains (snndofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld   %ld",i+1,snndofmas[i]);
  }
  
}
*/



/**
   function assembles DOFs indicators on master
   
   array dofmas is assembled
   
   JK, 14.9.2007
*/
/*
void seqselnodes::dof_indicators (gtopology *top,FILE *out)
{
  long i,j,k,l;
  
  dofmas = new long** [ns];
  for (i=0;i<ns;i++){
    dofmas[i] = new long* [nsndom[i]];
    for (j=0;j<nsndom[i];j++){
      dofmas[i][j] = new long [ndofnmas[i][j]];
    }
  }
    
  for (i=0;i<ns;i++){
    for (j=0;j<nsndom[i];j++){
      for (k=0;k<ndofnmas[i][j];k++){
        //l=lsng[i][j];
        l=lsnl[i][j];
        dofmas[i][j][k]=top->give_dof (l,k);
      }
    }
  }

  fprintf (out,"\n\n DOFs indicators on master (dofmas)\n");
  for (i=0;i<ns;i++){
    fprintf (out,"\n subdomain number %6ld",i+1);
    for (j=0;j<nsndom[i];j++){
      fprintf (out,"\n selected node %6ld    ",j+1);
      for (k=0;k<ndofnmas[i][j];k++){
        fprintf (out,"  %ld",dofmas[i][j][k]);
      }
    }
  }

}

*/
