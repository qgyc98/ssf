#include "floatsub.h"
#include <stdlib.h>

floatsub::floatsub (void)
{
  //  the number of subdomains
  nsd=0;
  
  //  the number of nodes and first node on subdomains
  nnsd=NULL;  fnnsd=NULL;
  //  the number of elements and first element on subdomains
  nesd=NULL;  fensd=NULL;
  //  the number of DOFs and first DOF on subdomain
  ndofd = NULL;  fdofd = NULL;

}

floatsub::~floatsub (void)
{
  delete [] nnsd;   delete [] fnnsd;
  delete [] nesd;   delete [] fensd;
  delete [] ndofd;  delete [] fdofd;
}

/**
   function assembles arrays with numbers of DOFs and first numbers of DOFs on subdomains
   
   JK, 6.11.2005
*/
void floatsub::ndof_subdom (gnode *gnodes)
{
  long i,j,k;
  
  //  allocation of arrays
  if (ndofd==NULL)
    ndofd = new long [nsd];
  if (fdofd==NULL)
    fdofd = new long [nsd+1];
  
  //  loop over subdomains
  fdofd[0]=0;
  for (i=0;i<nsd;i++){
    j=fnnsd[i];  k=fnnsd[i+1];
    ndofd[i] = gnodes[k].cn[0] - gnodes[j].cn[0];
    fdofd[i+1] = gnodes[k].cn[0];
  }
  
}

/**
   function computes number of Lagrange multipliers in the problem
   
   @param nlgn - number of layered general nodes
   @param lgnodes - layered general nodes
   @param gnodes - general nodes
   
   JK, 18.11.2005
*/
long floatsub::number_of_lagr_mult (long nlgn,lgnode *lgnodes,gnode *gnodes)
{
  long i,j,nn,n1,n2,ndofn;
  
  //  number of Lagrange multipliers
  nlm=0;
  
  //  loop over layered nodes
  for (i=0;i<nlgn;i++){
    //  the number of nodes connected via layered node
    nn = lgnodes[i].nl;
    if (nn!=2){
      fprintf (stderr,"\n\n wrong number of connected nodes in layered node %ld (file %s, line %d).\n",i+1,__FILE__,__LINE__);
    }
    
    n1=lgnodes[i].nodes[0];
    n2=lgnodes[i].nodes[1];
    
    ndofn = gnodes[n1].ndofn;
    j = gnodes[n2].ndofn;
    
    if (j!=ndofn){
      fprintf (stderr,"\n\n wrong numbers of DOFs at nodes (file %s, line %d).\n",__FILE__,__LINE__);
    }
    
    for (j=0;j<ndofn;j++){
      if (gnodes[n1].cn[j]!=gnodes[n2].cn[j])
	nlm++;
    }
  }
  
  return nlm;
}


/**
   function assembles extraction %matrix E for problems with floating subdomains
   
   JK, 18.11.2005
*/
void floatsub::displ_extract (long nlgn,lgnode *lgnodes,gnode *gnodes)
{
  long i,j,k,nn,n1,n2,cn1,cn2,ndofn;
  
  //  allocation of memory
  extrtab = new long* [nlm];
  for (i=0;i<nlm;i++){
    extrtab[i] = new long [2];
  }
  
  k=0;
  //  loop over layered nodes
  for (i=0;i<nlgn;i++){
    //  number of nodes connected via layered node
    nn = lgnodes[i].nl;
    if (nn!=2){
      fprintf (stderr,"\n\n wrong number of connected nodes in layered node %ld (file %s, line %d).\n",i+1,__FILE__,__LINE__);
    }
    
    n1=lgnodes[i].nodes[0];
    n2=lgnodes[i].nodes[1];
    
    ndofn = gnodes[n1].ndofn;
    j = gnodes[n2].ndofn;
    
    if (j!=ndofn){
      fprintf (stderr,"\n\n wrong numbers of DOFs at nodes (file %s, line %d).\n",__FILE__,__LINE__);
    }
    
    for (j=0;j<ndofn;j++){
      cn1=gnodes[n1].cn[j];
      cn2=gnodes[n2].cn[j];
      if (cn1!=cn2){
	if (cn1>0 && cn2>0){
	  extrtab[k][0]=0-cn1;
	  extrtab[k][1]=cn2;
	}
	if (cn1<1 && cn2>0){
	  extrtab[k][0]=0;
	  extrtab[k][1]=cn2;
	}
	if (cn1>0 && cn2<1){
	  extrtab[k][0]=0-cn1;
	  extrtab[k][1]=0;
	}
	if (cn1<1 && cn2<1){
	  extrtab[k][0]=0;
	  extrtab[k][1]=0;
	}
	k++;
      }
    }
  }
}

/**
   function extracts appropriate components from local %vector to coarse %vector
   
   @param cv - coarse %vector
   @param lv - local %vector
   
   JK, 18.11.2005
*/
void floatsub::local_coarse (double *cv,double *lv)
{
  long i,cn1,cn2;
  
  for (i=0;i<nlm;i++){
    cn1=extrtab[i][0];
    cn2=extrtab[i][1];
    
    
    if (cn1<0 && cn2>0){
      cv[i]-=lv[0-cn1-1];
      cv[i]+=lv[cn2-1];
    }
    if (cn1==0 && cn2>0){
      cv[i]+=lv[cn2-1];
    }
    if (cn1<0 && cn2==0){
      cv[i]-=lv[0-cn1-1];
    }
  }
  
}

/**
   function localizes appropriate components from coarse %vector to local %vector
   
   @param cv - coarse %vector
   @param lv - local %vector
   
   JK, 18.11.2005
*/
void floatsub::coarse_local (double *cv,double *lv)
{
  long i,cn1,cn2;
  
  for (i=0;i<nlm;i++){
    cn1=extrtab[i][0];
    cn2=extrtab[i][1];
    
    
    if (cn1<0 && cn2>0){
      lv[0-cn1-1]-=cv[i];
      lv[cn2-1]+=cv[i];
    }
    if (cn1==0 && cn2>0){
      lv[cn2-1]+=cv[i];
    }
    if (cn1<0 && cn2==0){
      lv[0-cn1-1]-=cv[i];
    }
  }
  
}
