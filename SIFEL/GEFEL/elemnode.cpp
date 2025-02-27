#include "elemnode.h"

elemnode::elemnode ()
{
  //  number of selected elements
  nse=0;
  //  list of selected elements
  lse=NULL;

  //  number of selected nodes
  nsn=0;
  //  list of selected nodes
  lsn=NULL;
  
  //  number of influenced elements
  nie=0;
  
  
  elnod = NULL; 
}

elemnode::~elemnode ()
{
  long i;
  
  delete [] lse;
  delete [] lsn;
  
  for (i=0;i<nie;i++){
    delete [] elnod[i];
  }
  delete [] elnod;
}

/**
   @param nselelem - number of selected elements
   @param lselem - list of selected elements
   
   JK, 20.10.2007
*/
void elemnode::selelem (long nselelem,long *lselem)
{
  long i;
  
  //  number of selected elements
  nse = nselelem;
  
  //  list of selected elements
  lse = new long [nse];
  
  for (i=0;i<nse;i++){
    lse[i]=lselem[i];
  }
}

/**
   @param nselnod - number of selected nodes
   @param lselnod - list of selected nodes
   
   JK, 20.10.2007
*/
void elemnode::selnode (long nselnod,long *lselnod)
{
  long i;
  
  //  number of selected nodes
  nsn = nselnod;
  
  //  list of selected nodes
  lsn = new long [nsn];
  
  for (i=0;i<nsn;i++){
    lsn[i]=lselnod[i];
  }
}

/**
   function assembles array elnod
   
   @param gt - pointer to general topology
   
   JK, 20.10.2007
*/
void elemnode::elemnodes (gtopology *gt)
{
  long i,j,k,l,ne,an;
  long *nne,**aux;
  ivector nod;
  
  //  number of elements in whole mesh
  ne = gt->ne;
  
  aux = new long* [ne];
  //  number of nodes on elements
  nne = new long [ne];

  //  number of influenced elements
  nie=0;
  
  for (i=0;i<ne;i++){
    //  number of nodes on the i-th element
    nne[i] = gt->give_nne (i);
    //  allocation of auxiliary array
    aux[i] = new long [nne[i]];
    //  allocation of memory
    allocv (nne[i],nod);
    //  nodes on the i-th element
    gt->give_nodes (i,nod);
    l=0;
    for (j=0;j<nne[i];j++){
      aux[i][j]=-1;
      //  number of actual node
      an=nod[j];
      for (k=0;k<nsn;k++){
	if (an==lsn[k]){
	  aux[i][j]=k;
	  l++;
	  break;
	}
      }
    }
    destrv (nod);
    
    if (l==0){
      nne[i]=0;
    }
    else{
      nie++;
    }
  }
  
  lse = new long [nie];
  elnod = new long* [nie];
  nie=0;
  for (i=0;i<ne;i++){
    if (nne[i]>0){
      lse[nie]=i;
      elnod[nie] = new long [nne[i]];
      for (j=0;j<nne[i];j++){
	elnod[nie][j]=aux[i][j];
      }
      nie++;
    }
  }
  
  
  for (i=0;i<ne;i++){
    delete [] aux[i];
  }
  delete [] aux;
  
  delete [] nne;
}
