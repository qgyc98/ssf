#include "lhsrhs.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechbclc.h"
#include "mechprint.h"
#include "gtopology.h"
#include "dloadcase.h"
#include "inicd.h"
#include <string.h>


/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
lhsrhs::lhsrhs ()
{
  //  number of DOFs
  ndof=Ndofm;
  //  number of load case
  nlc=Mb->nlc;
  
  //  array containing unknowns/DOFs
  lhs=NULL;
  //  array containing first time derivatives of unknowns/DOFs
  tdlhs=NULL;
  //  array containing second derivatives of unknowns/DOFs
  stdlhs=NULL;
  //  array containing the right hand side
  rhs=NULL;
  //  array containing eigenvalues
  w=NULL;
  // long vectors deallocation flag - by default the vectors should not be deallocated, flag is set by alloc method
  deallocf = no;
}




/**
  Destructor releases allocated memory of the lhsrhs object.

  Created by JK,
*/
lhsrhs::~lhsrhs ()
{
  if (deallocf == yes)
  {
    delete [] lhs;
    delete [] tdlhs;
    delete [] stdlhs;
    delete [] rhs;
    delete [] w;
  }
}



/**
  Function allocates arrays for left and right hand sides of problem.
   
  @return The function does not return anything.

  Created by JK, 29.6.2001
*/
void lhsrhs::alloc ()
{
  switch (Mp->tprob){
  case linear_statics:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case linear_stability:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    w = new double [nlc*ndof];
    memset (w,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case eigen_dynamics:{
    lhs = new double [Mp->eigsol.nev*ndof];
    memset (lhs,0,(Mp->eigsol.nev*ndof)*sizeof(double));
    
    switch (Mp->eigsol.teigsol){
    case inv_iteration:{
      w = new double [1];
      break;
    }
    case subspace_it_jacobi:{
      w = new double [Mp->eigsol.nev];
      break;
    }
    case subspace_it_gsortho:{
      w = new double [Mp->eigsol.nev];
      break;
    }
    case shifted_subspace_it_gsortho:{
      w = new double [Mp->eigsol.nev];
      break;
    }
    default:{
      print_err("unknown eigenvalue solver is required", __FILE__, __LINE__, __func__);
    }
    }
    deallocf = yes;
    break;
  }
  case forced_dynamics:{
    if (Mp->tforvib==modal_analysis){
      lhs = new double [ndof+Mp->eigsol.nev*ndof];
      memset (lhs,0,(ndof+Mp->eigsol.nev*ndof)*sizeof(double));

      switch (Mp->eigsol.teigsol){
      case inv_iteration:{
	w = new double [1];
	break;
      }
      case subspace_it_jacobi:{
	w = new double [Mp->eigsol.nev];
	break;
      }
      case subspace_it_gsortho:{
	w = new double [Mp->eigsol.nev];
	break;
      }
      case shifted_subspace_it_gsortho:{
	w = new double [Mp->eigsol.nev];
	break;
      }
      default:{
	print_err("unknown eigenvalue solver is required", __FILE__, __LINE__, __func__);
      }
      }
      
    }
    else{
      lhs = new double [nlc*ndof];
      memset (lhs,0,(nlc*ndof)*sizeof(double));
      tdlhs = new double [nlc*ndof];
      memset (tdlhs,0,(nlc*ndof)*sizeof(double));
      stdlhs = new double [nlc*ndof];
      memset (stdlhs,0,(nlc*ndof)*sizeof(double));
      rhs = new double [nlc*2*ndof];
      memset (rhs,0,(nlc*2*ndof)*sizeof(double));
    }
    deallocf = yes;
    break;
  }
  case mat_nonlinear_statics:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*2*ndof];
    memset (rhs,0,(nlc*2*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case geom_nonlinear_statics:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*2*ndof];
    memset (rhs,0,(nlc*2*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case mech_timedependent_prob:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case growing_mech_structure:{
    lhs = new double [ndof];
    memset (lhs,0,(ndof)*sizeof(double));
    rhs = new double [ndof];
    memset (rhs,0,(ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case earth_pressure:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*2*ndof];
    memset (rhs,0,(nlc*2*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case layered_linear_statics:{    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case lin_floating_subdomain:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case nonlin_floating_subdomain:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [2*nlc*ndof];
    memset (rhs,0,(2*nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case hemivar_inequal:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
  case load_balancing:{
    lhs = new double [nlc*ndof];
    memset (lhs,0,(nlc*ndof)*sizeof(double));
    rhs = new double [nlc*ndof];
    memset (rhs,0,(nlc*ndof)*sizeof(double));
    deallocf = yes;
    break;
  }
    
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }

}



/**
  Function returns pointer to %vector of unknowns (lhs %vector).
  This %vector contains nodal values.
   
  @param i - load case id
   
  @return The function returns pointer to the array containing %vetcor of unknowns.

  Created by JK,
*/
double *lhsrhs::give_lhs (long i)
{
  return (lhs+i*ndof);
}



/**
  Function returns pointer to %vector of first time derivatives of unknowns (tdlhs %vector).
  This %vector contains first time derivatives of nodal values.
   
  @param i - load case id
   
  @return The function returns pointer to the array containing %vetcor of first time drivatives of unknowns.

  Created by JK,
*/
double *lhsrhs::give_tdlhs (long i)
{
  return (tdlhs+i*ndof);
}



/**
  Function returns pointer to %vector of second time derivatives of unknowns (stdlhs %vector)
  this %vector contains second time derivatives of nodal values
   
  @param i - load case id
   
  @return The function returns pointer to the array containing %vetcor of second time drivatives of unknowns.

  Created by JK,
*/
double *lhsrhs::give_stdlhs (long i)
{
  return (stdlhs+i*ndof);
}



/**
   Function returns pointer to %vector of right hand side (rhs %vector)
   This %vector contains prescribed nodal values.
   
   @param i - load case id
   
   @return The function returns pointer to the array containing %vetcor of precribed nodal values.
   
   Created by JK,
*/
double *lhsrhs::give_rhs (long i)
{
  return (rhs+i*ndof);
}



/**
  Function initializes the solution vectors by initial conditions.
   
  @return The function does not return anything.

  Created by JK, 17.1.2002
  Modified by TKo & JK, 14.6.2017
*/
void lhsrhs::initcond ()
{
  long i,j,k,l,ndofn;

  if (Mp->tprob == forced_dynamics && Mb->nico && Mb->dlc[0].tdl!=responsespectrum){
    memset(lhs, 0, sizeof(*lhs)*Ndofm);
    memset(tdlhs, 0, sizeof(*tdlhs)*Ndofm);
    for (i=0;i<Mt->nn;i++) 
    {
      ndofn=Mt->give_ndofn (i);
      //  nodal values
      for (j=0;j<ndofn;j++)
      {
	k=Gtm->give_dof (i,j) - 1;
	if ((k>-1) && (Mb->ico[i].type & inidisp))  
          lhs[k]=Mb->ico[i].val[j];
      }
      if (Mb->ico[i].type & inidisp)
        l = ndofn;            
      else
        l = 0L;
      for (j=0;j<ndofn;j++)
      {
	k=Gtm->give_dof (i,j) - 1;
	if ((k>-1) && (Mb->ico[i].type & iniderdisp))  
          tdlhs[k]=Mb->ico[i].val[l];
        l++;
      }            
      /* old version
      xfscanf (in,"%ld",&l);      
      //  nodal values
      for (j=0;j<ndofn;j++){
	k=Gtm->give_dof (l-1,j) - 1;
	xfscanf (in,"%lf",&v);
	if (k>-1)  lhs[k]=v;
      }    
      //  time derivatives of nodal values
      for (j=0;j<ndofn;j++){
	k=Gtm->give_dof (l-1,j) - 1;
	xfscanf (in,"%lf",&v);
	if (k>-1)  tdlhs[k]=v;
      }
      */
    }
  }
  if ((Mp->tprob == mat_nonlinear_statics) && Mb->nico)
  {
    memset(lhs, 0, sizeof(*lhs)*Ndofm);
    for (i=0;i<Mt->nn;i++) 
    {
      ndofn=Mt->give_ndofn (i);
      //  nodal values
      for (j=0;j<ndofn;j++)
      {
	k=Gtm->give_dof (i,j) - 1;
	if ((k>-1) && (Mb->ico[i].type & inidisp))  
          lhs[k]=Mb->ico[i].val[j];
      }            
    }
  }
}



/**
  Function cleans %vector of unknowns.
   
  @return The function does not return anything.

  Created by JK,
*/
void lhsrhs::clean_lhs ()
{
  long i;
  for (i=0;i<ndof;i++){
    lhs[i]=0.0;
  }
}

