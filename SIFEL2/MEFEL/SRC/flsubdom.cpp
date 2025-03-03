#include "flsubdom.h"
#include "gtopology.h"
#include "global.h"
#include "probdesc.h"
#include "gmatrix.h"

flsubdom::flsubdom ()
{
  //  the number of changes
  nch=0;
  //  array of Lagrange multipliers
  lambda = NULL;
  //  array of total Lagrange multipliers
  totlambda = NULL;
  //  array of compliances
  compli=NULL;
}

flsubdom::~flsubdom ()
{
  delete [] lambda;
  delete [] totlambda;
  delete [] compli;
}

/**
   function initializes variables and arrays
   
   JK, 3. 8. 2012
*/
void flsubdom::initiation (gmatrix *gm,FILE *out)
{
  Mp->ssle->feti.matrices_assembl (gm,out);
}


/**
   function solves system of linear algebraic equations
   by the modified conjugate gradient method
   Lagrange multipliers are computed first and then displacements are obtained
   for more details see Kruis: Domain Decomposition Methods for
   Distributed Computing. Saxe-Coburg, 2006.
   
   @param lhs - left hand side %vector (%vector of nodal unknowns)
   @param rhs - right hand side
   
   JK, 20.6.2006, modified 2.8.2012
*/
void flsubdom::solve_lin_alg_system (double *lhs,double *rhs,FILE *out)
{
  //  assembling of the vectors of the right hand side for each subdomain
  Mp->ssle->feti.assemble_ff (rhs);

  //  assembling of the e vector
  Mp->ssle->feti.evector (out);

  //  preconditioned modified conjugate gradient method
  Mp->ssle->feti.mpcg (out);
  
  //  computation of nodal unknowns from Lagrange multipliers
  Mp->ssle->feti.nodalunknowns (lhs,out);
}

/**
   function adds increments of Lagrange multipliers to
   the %vector of total multipliers in nonlinear problems
   
   @param dlambda - increment from arclength method
   
   JK, 22.6.2006
*/
/*
void flsubdom::add_mult (double dlambda)
{
  long i;
  
  //for (i=0;i<nlm;i++){
    //tw[i]+=dlambda*w[i];
  //}

  for (i=0;i<nlm;i++){
    w[i]*=dlambda;
    tw[i]+=w[i];
  }
}
*/
/**
   function corrects Lagrange multipliers with respect to
   defined law
   
   JK, 24.6.2006
*/
 /*
void flsubdom::mult_correction ()
{
  long i;
  
  //cv = new double [nlm];
  //Gtm->flsub.local_coarse (cv,Lsrs->lhs);
  
  nch=0;
  for (i=0;i<nlm;i++){
    
    if (tw[i]>1.0 && tw[i]<-1.0){
      compli[i]=5.0;
    }

 */
    /*
    if (mu[i]>100.0){
      //tw[i]=10.0;
      //compli[i]=20000.0;
      //compli[i]=2.0;
      compli[i]=(mu[i]-100.0)/10000.0;
      nch++;
    }
    */
    /*
    if (tw[i]>10.0 && cv[i]<10.0){
      //tw[i]=10.0;
      //compli[i]=20000.0;
      compli[i]=2.0;
      //compli[i]=tw[i]-10.0;
      nch++;
    }
    if (tw[i]>10.0 && cv[i]>10.0){
      //tw[i]=10.0;
      compli[i]=20000.0;
      //compli[i]=2.0;
      //compli[i]=tw[i]-10.0;
      nch++;
    }
    */
    /*
    if (tw[i]>100.0)
      tw[i]=100.0;
    */

    /*
  }
  
  //delete [] cv;
}

    */
