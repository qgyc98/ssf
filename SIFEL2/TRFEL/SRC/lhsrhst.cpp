#include "globalt.h"
#include "genfile.h"
#include "globmatt.h"
#include "lhsrhst.h"
#include <string.h>

lhsrhst::lhsrhst ()
{
  nlc=Tb->nlc;  lhs=NULL;  rhs=NULL;  lhsi=NULL;  tdlhs=NULL;
  // long vectors deallocation flag - by default the vectors should not be deallocated, flag is set by alloc method
  deallocf = no;
}

lhsrhst::~lhsrhst ()
{
  if (deallocf == yes)
  {
    delete [] lhs;  
    delete [] rhs;  
    delete [] lhsi;  
    delete [] tdlhs;
  }
}

/**
   function allocates array for left and right sides of problem
   29.6.2001
*/
void lhsrhst::alloc ()
{
  switch (Tp->tprob){
  case stationary_problem:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case nonlinear_stationary_problem:{
    lhs = new double [nlc*Ndoft];
    memset (lhs,0,(nlc*Ndoft)*sizeof(double));
    lhsi = new double [nlc*Ndoft];
    memset (lhsi,0,(nlc*Ndoft)*sizeof(double));
    rhs = new double [2*nlc*Ndoft];
    memset (rhs,0,(2*nlc*Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case nonstationary_problem:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    tdlhs = new double [Ndoft];
    memset (tdlhs,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case nonlinear_nonstationary_problem:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    tdlhs = new double [Ndoft];
    memset (tdlhs,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case discont_nonstat_problem:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    tdlhs = new double [Ndoft];
    memset (tdlhs,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case discont_nonlin_nonstat_problem:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    tdlhs = new double [Ndoft];
    memset (tdlhs,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case growing_np_problem:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    tdlhs = new double [Ndoft];
    memset (tdlhs,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  case growing_np_problem_nonlin:{
    lhs = new double [Ndoft];
    memset (lhs,0,(Ndoft)*sizeof(double));
    lhsi = new double [Ndoft];
    memset (lhsi,0,(Ndoft)*sizeof(double));
    tdlhs = new double [Ndoft];
    memset (tdlhs,0,(Ndoft)*sizeof(double));
    rhs = new double [Ndoft];
    memset (rhs,0,(Ndoft)*sizeof(double));
    deallocf = yes;
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

double *lhsrhst::give_lhs (long lcid)
{
  return (lhs+lcid*Ndoft);
}

double *lhsrhst::give_lhsi (long lcid)
{
  return (lhsi+lcid*Ndoft);
}

double *lhsrhst::give_tdlhs (long lcid)
{
  return (tdlhs+lcid*Ndoft);
}

double *lhsrhst::give_rhs (long lcid)
{
  return (rhs+lcid*Ndoft);
}

/**
   function reads initial conditions

   @param in - input data stream 
   
   17.1.2002
*/
void lhsrhst::initcond (XFILE *in)
{
  long i,j,k,l,ndofn;
  double v;
  
  switch (Tp->tprob){
  case stationary_problem:{
    break;
  }
  case nonlinear_stationary_problem:
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:
  case discont_nonstat_problem:
  case discont_nonlin_nonstat_problem:{
    for (i=0;i<Tt->nn;i++){
      xfscanf (in,"%ld",&l);
      ndofn=Gtt->give_ndofn (l-1);
      if (ndofn < 0)
        continue; 
      for (j=0;j<ndofn;j++){
	k=Gtt->give_dof (l-1,j)-1;
	xfscanf (in,"%lf",&v);
	if (k>-1) {  lhsi[k]=v;  }
      }
    }
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

/**
   function prints initial conditions

   @param out - output data stream 

   TKr 3.1.2006
*/
void lhsrhst::initcondprint (FILE *out)
{
  long i,j,k,ndofn;
  
  fprintf (out,"\n\n");
  switch (Tp->tprob){
  case nonlinear_stationary_problem:
  case nonstationary_problem:
  case nonlinear_nonstationary_problem:
  case growing_np_problem:
  case growing_np_problem_nonlin:
  case discont_nonlin_nonstat_problem:{
    for (i=0;i<Tt->nn;i++){
      ndofn=Gtt->give_ndofn (i);
      fprintf (out,"\n %ld  ",i+1);
      for (j=0;j<ndofn;j++){
	k=Gtt->give_dof (i,j)-1;
	if (k>-1){ fprintf (out," %lf",lhsi[k]);}
	else {fprintf (out," 0.0");}
      }
    }
    break;
  }
  default:{
    print_err("unknown problem type is required",__FILE__,__LINE__,__func__);
  }
  }
}

