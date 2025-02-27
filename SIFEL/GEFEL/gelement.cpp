#include "gelement.h"
#include <string.h>
#include "gnode.h"

gelement::gelement (void)
{
  //  the number of nodes on element
  nne = 0;
  //  the number of master nodes on element
  nmne = 0;
  //  the number of degrees of freedom
  ndofe = 0;
  //  indicator of code numbers on element (code numbers are on nodes)
  cne = 0;
  //  type of element (with respect to geometry)
  get=noelem;
  //  number of additional degrees of freedom (especially Lagrange multipliers)
  nmult=0;

  //  array of node numbers
  nodes=NULL;
  //  array containing node numbers with respect to graphic purposes
  //  this array is used only when the hanging nodes are used
  master_nodes=NULL;
  //  auxiliary information about the element
  auxinf=0;
  
  //edgn = NULL;
  //  array for code numbers, if necessary
  cn = NULL;
  //  general time function
  tgf=0;
}

gelement::~gelement (void)
{
  //  array of node numbers
  delete [] nodes;
  //  array containing the master node numbers
  //  this array is used only when the hanging nodes are used
  delete [] master_nodes;
  //delete [] edgn;
  delete [] cn;
}

/**
   function reads informations from input file
   
   @param[in] in - input stream
   @param[in] m - the number nodes on element
   @param[in] n - the number of degrees of freedom on element
   @param[in] et - type of element
   @param[in] maxnn - the maximum number of nodes in the problem (for the checking purpose)
   
   16.7.2002
*/
void gelement::read (XFILE *in,long m,long n,gelemtype et, long maxnn)
{
  long i;

  //  type of finite element with respect to geometry
  get=et;
  //  the number nodes on element
  nne=m;
  //  the number of degrees of freedom on element
  ndofe=n;
  
  nodes = new long [nne];
  for (i=0;i<nne;i++){
    xfscanf (in,"%ld",nodes+i);
    nodes[i]--;
    if ((nodes[i] >= 0) || (nodes[i] < maxnn))
      continue;
    else{
      print_err("invalid node number %ld of the %ld-th element node, it must be in [1;%ld]", __FILE__, __LINE__, __func__, nodes[i]+1, i+1, maxnn);
      abort();
    }
  }
  
  xfscanf (in,"%ld",&cne);
  if (cne!=0 && cne!=1 && cne!=2){
    print_err("identification of code numbers on element (cne) must be equal to 0, 1 or 2", __FILE__, __LINE__, __func__);
  }
  
  // spring element type
  if ((ndofe < -1) && (nne == 1))
    return;

  if (cne==1 || cne==2){
    cn = new long [ndofe];
    for (i=0;i<ndofe;i++){
      xfscanf (in,"%ld",cn+i);
    }
  }
  
}

/**
   function prints informations into output file
   
   @param out - output stream
   @param m - the number nodes on element
   @param n - the number of degrees of freedom on element
   @param et - type of element
   
   19/12/2012, Tkr
*/
void gelement::print (FILE *out,long /*m*/,long /*n*/,gelemtype /*et*/)
{
  long i;

  for (i=0;i<nne;i++){
    fprintf (out," %ld ",nodes[i]+1);
  }
  
  fprintf (out," %ld ",cne);

  if (cne==1 || cne==2){
    for (i=0;i<ndofe;i++){
      fprintf (out," %ld ",cn[i]+1);
    }
  }
}



/**
   function reads informations from input file
   
   @param[in] in - input stream
   @param[in] m - number nodes on element
   @param[in] n - number of degrees of freedom on element
   @param[in] et - type of element
   @param[in] maxnn - the maximum node number }for the checking purpose)
   
   16.7.2002
*/
void gelement::read_gf (XFILE *in,long m,long n,gelemtype et,long maxnn)
{
  long i;

  //  type of finite element with respect to geometry
  get=et;
  //  number nodes on element
  nne=m;
  //  number of degrees of freedom on element
  ndofe=n;
  
  nodes = new long [nne];
  for (i=0;i<nne;i++){
    xfscanf (in,"%ld",nodes+i);
    nodes[i]--;
    if ((nodes[i] >= 0) || (nodes[i] < maxnn))
      continue;
    else{
      print_err("invalid node number %ld of the %ld-th element node, it must be in [1;%ld]", __FILE__, __LINE__, __func__, nodes[i]+1, i+1, maxnn);
      abort();
    }
  }
  
  xfscanf (in,"%ld",&cne);
  if (cne!=0 && cne!=1 && cne!=2){
    print_err("identification of code numbers on element (cne) must be equal to 0, 1 or 2", __FILE__, __LINE__, __func__);
  }
  
  // spring element type
  if ((ndofe < -1) && (nne == 1))
    return;

  if (cne==1 || cne==2){
    cn = new long [ndofe];
    for (i=0;i<ndofe;i++){
      xfscanf (in,"%ld",cn+i);
    }
  }
  
  xfscanf (in,"%ld",&tgf);
  tgf--;
}



/**
   function prints informations into output file
   
   @param out - output stream
   @param m - number nodes on element
   @param n - number of degrees of freedom on element
   
   19/12/2012, Tkr
*/
void gelement::print_gf (FILE *out,long /*m*/,long /*n*/)
{
  long i;

  for (i=0;i<nne;i++){
    fprintf (out," %ld ",nodes[i]+1);
  }
  
  fprintf (out," %ld ",cne);

  if (cne==1 || cne==2){
    for (i=0;i<ndofe;i++){
      fprintf (out," %ld ",cn[i]+1);
    }
  }
  
  fprintf (out," %ld ",tgf+1);
}


/**
   function initiates element values
   function is used in parallel computation
   
   @param icn - array containing element code numbers
   @param m - number of degrees of freedom on element
   
   24.7.2002
*/
void gelement::initiate (long *icn,long m)
{
  long i;
  cne=1;
  ndofe=m;
  
  cn = new long [ndofe];
  for (i=0;i<ndofe;i++){
    cn[i]=icn[i];
  }
  
}

/**
   function returns number of nodes on element
   
   JK
*/
long gelement::give_nne ()
{
  return nne;
}

/**
   function returns the number of master nodes on element
   
   JK
*/
long gelement::give_nmne ()
{
  return nmne;
}

/**
   function returns number of degrees of freedom of element
   
   JK
*/
long gelement::give_ndofe ()
{
  return ndofe;
}

/**
   function returns number of Lagrange multipliers defined on element
   
   JK
*/
long gelement::give_nmult ()
{
  return nmult;
}

/**
   function selects nodes of element
   
   @param nod - %vector containing selected nodes
   
   JK
*/
void gelement::give_nodes (ivector &nod)
{
  long i;
  if (nod.n!=nne){
    print_err("wrong size of array is required", __FILE__, __LINE__, __func__);
  }
  
  for (i=0;i<nne;i++){
    nod[i]=nodes[i];
  }
}

/**
   function selects nodes of element
   
   @param nod - %vector containing selected nodes
   
   JK, 17.6.2012
*/
void gelement::give_master_nodes (ivector &nod)
{
  long i;
  if (nod.n!=nmne){
    print_err("wrong size of array is required", __FILE__, __LINE__, __func__);
  }
  
  for (i=0;i<nmne;i++){
    nod[i]=master_nodes[i];
  }
}

/**
   function returns code number indicator of element
   if cne=1, code numbers are defined on element,
   otherwise code numbers are defined on nodes
   
   JK
*/
long gelement::give_cne ()
{
  return cne;
}

/**
   Function computes centroid of element.
   
   @param dim - dimension of solved problem
   @param gnodes - array of all gnodes
   @param - answer, array of coordinates of centroid
   
   created  11.12.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
void gelement::centroid (long dim,gnode *gnodes,double *coord)
{
  long i;
  
  switch (auxinf){
  case 312:
  case 412:{
    coord[0]=coord[1]=0.0;
    for (i=0;i<nne;i++){
      coord[0] += gnodes[nodes[i]].x;
      coord[1] += gnodes[nodes[i]].y;
    }
    for (i=0;i<dim;i++)
      coord[i] /= nne;
    
    break;
  }
  default:{
    print_err("wrong dimdegnne", __FILE__, __LINE__, __func__);
    break;
  }
  }
}
