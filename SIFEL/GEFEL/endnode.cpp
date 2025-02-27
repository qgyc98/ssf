#include "endnode.h"
#include "gnode.h"
#include "gtopology.h"
#include <math.h>

endnode::endnode (void)
{
  //  number of nodes in end point
  nn=0;
  //  number of first node
  fn=-1;
  //  number of last node
  ln=-1;
  //  node multiplicity
  nm=0;
  //  number of DOFs
  ndofn=0;
  

  //  number of reference element
  re=-1;
  //  list of adjacent elements
  adjel=NULL;

  //  array of code numbers of Lagrange multipliers
  cnm = NULL;

  //  threshold
  threshold=1.0e-6;
}


endnode::~endnode (void)
{
  delete [] adjel;
  delete [] cnm;
}


/**
   function prints
   
   @param out - output stream
   
   JK, 12.10.2008
*/
void endnode::print (FILE *out)
{
  fprintf (out,"\n");
  fprintf (out,"\n number of nodes    %ld",nn);
  fprintf (out,"\n node multiplicity  %ld",nm);
}

/**
   function allocates array for code numbers of Lagrange multipliers
   
   @param nccn - number of components in array cnm
   
   JK, 12.10.2008
*/
void endnode::alloc_cnm (long nccn)
{
  //  number of DOFs
  ndofn=nccn;
  
  if (cnm!=NULL)
    delete [] cnm;
  cnm = new long [ndofn];
}

/**
   function returns node numbers shared in the end node
   
   @param nid - array containing node numbers shared in the end node
   
   JK, 12.10.2008
*/
void endnode::give_node_numbers (long *nid)
{
  nid[0]=fn;
  nid[1]=ln;

}

/**
   function assembles code numbers of Lagrange multipliers
   defined between selected nodes
   
   @param mcn - code numbers of selected multipliers
   
   JK, 12.10.2008
*/
void endnode::give_mult_code_numbers (long *mcn)
{
  long i;
  
  for (i=0;i<ndofn;i++){
    mcn[i]=cnm[i];
  }
}
