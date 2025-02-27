#include "sselnodes.h"

sselnodes::sselnodes (long nd,long kk,long *j)
{
  long i;
  
  //  number of subdomains/aggregates
  ndom=nd;
  
  //  number of all nodes
  nn = kk;
  
  //  number of selected nodes
  nsn=0;
  for (i=0;i<nn;i++){
    if (j[i]>-1)
      nsn++;
  }
  
  if (nsn==0){
    fprintf (stderr,"\n\n wrong number of selected nodes in constructor (file %s, line %d)\n",__FILE__,__LINE__);
    //abort();
  }
  
  //  list of selected nodes - local numbers
  lsnl=new long [nsn];
  //  list of selected nodes - global numbers
  lsng=new long [nsn];
  
  k=0;
  for (i=0;i<nn;i++){
    if (j[i]>-1){
      lsnl[k]=i;
      lsng[k]=j[i];
      k++;
    }
  }

}

sselnodes::sselnodes (long nd,long *ii,long **jj)
{
  //  number of subdomains/aggregates
  ndom=nd;
  
  nsndom = new long [ndom];
  for (i=0;i<ndom;i++){
    nsndom[i]=ii[i];
  }
  
  gnn = new long* [ndom];
  for (i=0;i<ndom;i++){
    gnn[i] = new long [nsndom[i]];
  }
  
  for (i=0;i<ndom;i++){
    for (j=0;j<nsndom[i];j++){
      gnn[i][j]=jj[i][j];
    }
  }
  
}

sselnodes::~sselnodes ()
{
  delete [] lsnl;
  delete [] lsng;

}


void sselnodes::assemble_list_unknowns (gtopology *gt)
{
  long i,j,k,l,ann,nu,ndofn;
  
  if (ndofdom!=NULL){
    delete [] ndofdom;
  }
  ndofdom = new long [ndom];
  
  maxnu=0;
  for (i=0;i<ndom;i++){
    nu=0;
    for (j=0;j<nsndom[i];j++){
      //  actual node number
      ann=gnn[i][j];
      //  number of unknowns on node
      ndofn = gt->give_ndofn (ann);
      for (k=0;k<ndofn;k++){
	l=gt->give_dof (ann,k);
	if (l>0)
	  nu++;
      }
      
    }
    ndofdom[i]=nu;
    if (nu>maxnu)
      maxnu=nu;
  }
  
  if (cndom!=NULL){
    for (i=0;i<ndom;i++){
      delete [] cndom[i];
    }
    delete [] cndom;
  }
  cndom = new long* [ndom];
  for (i=0;i<ndom;i++){
    cndom[i] = new long [ndofdom[i]];
  }
  
  
  for (i=0;i<ndom;i++){
    ndofdom[i]=0;
    for (j=0;j<nsndom[i];j++){
      //  actual node number
      ann=gnn[i][j];
      //  number of unknowns on node
      ndofn = gt->give_ndofn (ann);
      for (k=0;k<ndofn;k++){
	l=gt->give_dof (ann,k);
	if (l>0){
	  cndom[i][ndofdom[i]]=l;
	  ndofdom[i]++;
	}
      }
      
    }
  }

}
