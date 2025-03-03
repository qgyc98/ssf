#include <string.h>
#include "lhsrhsc.h"
#include "globalc.h"
#include "gtopology.h"



lhsrhsc::lhsrhsc ()
{
  nlc=Cb->nlc;  lhs=NULL;  lhsi=NULL;  tdlhs=NULL;  rhs=NULL;
}

lhsrhsc::~lhsrhsc ()
{
  delete [] lhs;  delete [] lhsi;  delete [] tdlhs;  delete [] rhs;
}

/**
   function allocates array for left hand side, time derivative of
   left hand side, initial values and right sides of problem
   
   6.4.2003, JK
*/
void lhsrhsc::alloc ()
{
  lhs = new double [nlc*Ndofc];
  memset (lhs,0,(nlc*Ndofc)*sizeof(double));
  lhsi = new double [nlc*Ndofc];
  memset (lhsi,0,(nlc*Ndofc)*sizeof(double));
  tdlhs = new double [nlc*Ndofc];
  memset (tdlhs,0,(nlc*Ndofc)*sizeof(double));
  rhs = new double [nlc*Ndofc];
  memset (rhs,0,(nlc*Ndofc)*sizeof(double));
}

double *lhsrhsc::give_lhs (long lcid)
{
  return (lhs+lcid*Ndofc);
}

double *lhsrhsc::give_lhsi (long lcid)
{
  return (lhsi+lcid*Ndofc);
}

double *lhsrhsc::give_tdlhs (long lcid)
{
  return (tdlhs+lcid*Ndofc);
}

double *lhsrhsc::give_rhs (long lcid)
{
  return (rhs+lcid*Ndofc);
}


/**
   function reads initial conditions
   
   @param in - input data stream 
   17.1.2002
*/
void lhsrhsc::initcond (FILE */*in*/)
{

}

