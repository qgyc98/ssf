#include "bnodvalt.h"

bnodvalt::bnodvalt ()
{
  //  number of stored components
  nsc=0;
  //  nodal values defined on boundaries
  nodval=NULL;
}

bnodvalt::~bnodvalt ()
{
  delete [] nodval;
}

/**
   function reads input data
   
   @param in - input file
   
   JK, 24.11.2008
*/
void bnodvalt::read (XFILE *in)
{
  long i;
  
  //  number of stored components
  xfscanf (in,"%k%ld", "ncomp", &nsc);
  
  nodval = new gfunct [nsc];
  
  for (i=0;i<nsc;i++){
    nodval[i].read (in);
  }
  
}

/**
   function prints input data
   
   @param out - output file
   
   JK, 24.11.2008
*/
void bnodvalt::print (FILE *out)
{
  long i;
  
  //  number of stored components
  fprintf (out,"%ld\n",nsc);
  
  for (i=0;i<nsc;i++){
    nodval[i].print (out);
  }
  
}

/**
   function computes nodal values
   it assembles nodal values on end nodes, edges
   or surfaces
   
   @param t - input variable, usually time
   @param nv - array of nodal values
   
   JK, 24.11.2008
*/
void bnodvalt::give_val (double t,vector &nv)
{
  long i;
  
  for (i=0;i<nsc;i++){
    nv[i]=nodval[i].getval (t);
  }
}
