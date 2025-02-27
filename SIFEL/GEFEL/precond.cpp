#include "precond.h"
#include "galias.h"
#include "gtopology.h"
#include "gmatrix.h"
#include "aggregator.h"

precond::precond ()
{
  //  type of preconditioning
  pt = noprecond;
  //  number of unknowns/equations
  n = 0;
  //  treshold for incomplete factorization
  incompltresh=0.0;

  agg=NULL;
  agm=NULL;
}

precond::~precond ()
{
  
}

/**
   function reads input data
   
   @param in - input stream
   @param mespr - type of message printing
   
   JK, 15.3.2007
*/
void precond::read (gtopology *gt,XFILE *in,long mespr)
{
  //  type of preconditioner
  xfscanf (in, "%k%m", "precond_type", &precondtype_kwdset, (int*)&pt);
  
  switch (pt){
  case noprecond:{
    if (mespr==1)  fprintf (stdout,"\n no preconditioner will be used");
    break;
  }
  case diagprec:{
    if (mespr==1)  fprintf (stdout,"\n the diagonal (Jacobi) preconditioner will be used");
    break;
  }
  case ssorprec:{
    if (mespr==1)  fprintf (stdout,"\n the SSOR preconditioner will be used");
    xfscanf (in,"%lf",&ssoromega);
    break;
  }
  case incomdec:{
    if (mespr==1)  fprintf (stdout,"\n preconditioner based on incomplete factorization of matrix will be used");
    xfscanf (in,"%lf",&incompltresh);
    break;
  }
  case boss:{
    agg = new aggregator [1];
    if (mespr==1)  fprintf (stdout,"\n preconditioner based on the BOSS (smoothed aggregations) will be used");
    agg->read (gt,in,mespr);
    break;
  }
  default:{
    print_err("unknown type of preconditioner of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
  
}

/**
   function prints input data
   
   @param out - output stream
   
   JK, 15.3.2007
*/
void precond::print (FILE *out)
{
  //  type of preconditioner
  fprintf (out,"%d\n",pt);
  
  switch (pt){
  case noprecond:{
    break;
  }
  case diagprec:{
    break;
  }
  case ssorprec:{
    fprintf (out,"%f\n",ssoromega);
    break;
  }
  case incomdec:{
    fprintf (out,"%f\n",incompltresh);
    break;
  }
  case boss:{
    agg->print (out);
    break;
  }
  default:{
    print_err("unknown type of preconditioner of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }

}



/**
   function initiates necessary vectors, matrices, etc.
   e.g., it computes incomplete factorization or prepare aggregates
   in the case of BOSS preconditioning
   
   @param gt - pointer to general topology
   @param gm - pointer to general %matrix
   @param out - output stream
   
   JK, 14.3.2007
*/
void precond::initiation (gtopology *gt,gmatrix *gm,FILE *out)
{
  //  number of unknowns
  n=gm->n;
  
  switch (pt){
  case noprecond:{
    break;
  }
  case incomdec:{
    agm = new gmatrix [1];
    gm->copygm (*agm);
    agm->incomplete_fact (incompltresh);
    break;
  }
  case boss:{
    //agg = new aggregator [1];
    //agg->prepare_boss (gt,naggr,gm,tlinsol,out);
    agg->prepare_boss (gt,gm,out);
    agm = new gmatrix;
    gm->copygm (*agm);
    break;
  }
  default:{
    print_err("unknown type of preconditioner of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }

}

/**
   function computes preconditioning
   function modifies the residual %vector
   
   @param r - residual %vector
   @param h - preconditioned residual %vector
   @param rhs - right hand side %vector
   
   JK, 15.3.2007
*/
void precond::preconditioning (double *r,double *h,double */*rhs*/)
{
  long i;
  
  switch (pt){
  case noprecond:{
    for (i=0;i<n;i++){
      h[i]=r[i];
    }
    break;
  }
  case boss:{
    agg->boss (agm,r,h);
    break;
  }
  case incomdec:{
    agm->back_incomplete_fact (h,r,incompltresh);
    break;
  }
  default:{
    print_err("unknown type of preconditioner of solver of linear equations is required",__FILE__,__LINE__,__func__);
  }
  }
}
