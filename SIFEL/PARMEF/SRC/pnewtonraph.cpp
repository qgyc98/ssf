#include "mpi.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "pglobal.h"
#include "pnewtonraph.h"
#include "nssolver.h"
#include "backupsol.h"
#include "global.h"
#include "globmat.h"
#include "mechprint.h"
#include "gmatrix.h"
#include "gtopology.h"
#include "mechprint.h"
#include "mathem.h"
#include "vector.h"
#include "matrix.h"




/**
  Function performs calculation of one parallel load/time step of the Newton-Raphson 
  method for the given load case. Solved equation system does not contain 
  time variable.
  
  @param lcid[in]   - load case id
  @param nlman[in]  - pointer to structure conatining setup of the solver
  @param ra[in/out] - %vector of attained displacements
  @param fa[in]     - attained load %vector
  @param fb[ut]     - residual %vector or load %vector increment - right-hand side %vector
  @param dr[out]    - %vector of displacement increment - left-hand side %vector
  @param fi[out]    - %vector of internal forces
  @param dtr[out]   - reduction factor of step length obtained from material models
  @param istep[in]  - time/load step id
  @param j[in/out]  - inner loop step id (output parameter, set to -1 if no inner loop was performed) 
  @param li[in]     - initial value of time/load step id
  @param ierr[in]    - required normed error of residual %vector

  @return The function returns reached lambda parameter.

  TKr, 19/02/2013 according to Tomas Koudelka
*/
double par_gnewton_raphson_one_step(long lcid, nonlinman *nlman, double *fa, double *ra, double *fb, double *dr, double *fi, 
                                    double &dtr, long istep, long &j, long li, long ini, double ierr)
{
  long n = Ndofm;
  double norf, norfa;
  matrix lsm_a(3,3); // least square matrix for the divergence detection
  vector lsm_r(3), lsm_l(3); // right hand and left hand side for the divergence detection
  resnormt normt;    // type of residual vector norm
  
  // type of residual vector norm
  normt = nlman->rnormtnr;

  //  assembling of tangent stiffness matrix
  assemble_stiffness_matrix(lcid,istep,-1,li,no);
    
  //  solution of K(r).v=F
  Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);//Mp->ssle->solve_system (Gtm,Smat,dr,fb,Out);

  //  update total displacement vector
  //  ra[j]+=dr[j];
  addv(ra, dr, n);
    
  //  computation of internal forces
  internal_forces (lcid,fi);

  // check time step size requirements which originate in material models, 
  // after check, the dtr will contain the minimum ratio of required time step size to the actual one 
  // ( dtr \in (0.0; 1.0) ) collected from all integration points or dtr = 1.0 for no change 
  // in time step size
  dtr = par_check_dtr();  
  if (dtr < 1.0) // reduction of time step size was required by some material model
  {
    j=-1;
    return (2.0*ierr);
  }
  
    
  //  vector of unbalanced forces
  //  fb[j] = fa[j] - fi[j];
  subv(fa, fi, fb, n);

  //  norm of vector of unbalanced forces
  norf = par_compute_res_norm_nr(lcid, normt, fb, fa, n, norfa);

  if ((Myrank==0)&&(Mespr==1)) 
    fprintf (stdout,"\n increment %ld:  init. norf %le, req. err %le",istep, norf, ierr);

  if (Psolm->compare_on_master(norf, ierr) <= 0){ // norf <= ierr
    j = -1;
    return norf;
  }
    
  //  internal iteration loop
  fillm(0.0, lsm_a);
  fillv(0.0, lsm_r);
  for (j=0; j<ini; j++)
  {
    Mp->jstep = j;

    // re-assemble stifness matrix if it is required
    assemble_stiffness_matrix(lcid,istep,j,li,no);
    
    // solution of K(r).v=F
    Psolm->par_linear_solver (Gtm,Smat,dr,fb,Out,Mespr);//Mp->ssle->solve_system(Gtm,Smat,dr,fb,Out);
    
    // ra[k]+=dr[k];
    addv(ra, dr, n);
      
    // computation of internal forces
    internal_forces(lcid,fi);
      
    // check time step size requirements which originate in material models, 
    // after check, the dtr will contain the minimum ratio of required time step size to the actual one 
    // ( dtr \in (0.0; 1.0) ) collected from all integration points or dtr = 1.0 for no change 
    // in time step size
    dtr = par_check_dtr();    
    if (dtr < 1.0){ // reduction of time step size was required by some material model
      j=-1;
      return (2.0*ierr);
    }

    // vector of unbalanced forces
    // fb[k]=fa[k]-fi[k]
    subv(fa, fi, fb, n);
      
    // norm of vector of unbalanced forces
    norf = par_compute_res_norm_nr(lcid, normt, fb, fa, n, norfa);

    if ((Myrank==0) && (Mespr==1))
      fprintf (stdout,"\n increment %ld   inner loop j %ld     norf %e", istep, j, norf);
      
    if (Psolm->compare_on_master(norf, ierr) <= 0) // norf <= ierr
      // equilibrium was attained
      return norf;

    // divergence detection with help of the least square method
    if (check_divergency(nlman, lsm_a, lsm_r, lsm_l, j, norf))
      return norf;
  }

  // equilibrium was not attained in this step
  return norf;  
}


/**
  The function check time step size requirements which originate in material models, 
  after check, the dtr will contain the minimum ratio of required time step size to the actual one 
  ( dtr \in (0.0; 1.0) ) collected from all integration points and all domains or dtr = 1.0 for no change 
  in time step size

  Created by Tomas Koudelka, 08.2021
*/
double par_check_dtr(void)
{
  long k;
  double dtr, aux;
  MPI_Status stat;
  
  dtr = Mm->dstep_red_mat();
  if (Myrank == 0){
    for(k=1; k<Nproc; k++){
      MPI_Recv(&aux, 1, MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      if (aux < dtr)
        dtr = aux;
    }    
  }
  else{   
    MPI_Send(&dtr, 1, MPI_DOUBLE, 0, Myrank, MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  if (Myrank==0){
    for (k=1; k<Nproc; k++){
      MPI_Send(&dtr, 1, MPI_DOUBLE, k, Myrank, MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv(&dtr, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

  return dtr;
}



/**
  The function computes norm of the residual vector with respect to the given norm type.
  The norm is computed according to the argument normt which specifies whether the norm is either
  relative to the attained load vector fb or reactions or it may be absolute.

  @param[in] lcid - load case id
  @param[in] normt - type of the residual norm
  @param[in] fb - residual vector
  @param[in] fa - load vector
  @param[in] n  - dimension of vectors fb and fa
  @param[out] norfa - attained norm of the load vector

  @return The function returns norm of the residual vector fb with respect to settings given by normt.
*/
double par_compute_res_norm_nr(long lcid, resnormt normt, double *fb, double *fa, long n, double &norfa)
{
  double norf=0.0;
  norfa=0.0;

  //  norm of vector of unbalanced forces (residual norm)
  norf = Psolm->pss(fb, fb, Out);
  norf = sqrt(norf);
  switch(normt){
    case rel_load_norm: // norm is relative to attained load vector
      // compute norm of attained load vector
      norfa = Psolm->pss(fa, fa, Out);
      norfa = sqrt(norfa);
      if (norfa != 0.0)
        norf /= norfa;
      break;
    case rel_react_norm: // norm is relative to attained reactions
      compute_reactions(lcid);
      norfa = par_compute_react_norm();
      norfa = sqrt(norfa);
      if (norfa != 0.0)
        norf /= norfa;
      break;
    case rel_loadreact_norm: // norm is relative to attained reactions and load
      // compute norm of attained load vector
      compute_reactions(lcid);
      norfa = par_compute_react_norm();
      norfa += Psolm->pss(fa, fa, Out);
      norfa  = sqrt(norfa);
      if (norfa != 0.0)
        norf /= norfa;
      break;
    case absol_norm: //absolute norm of the residual vector
      break;
    default:
      print_err("unknown type of residual norm %d is required", __FILE__, __LINE__, __func__, normt);
      abort();
  }

  return norf;
}



/**
  The function computes norm of all reactions in the problem. It is used for the scaling
  of the residual %vector norm.

  @retval The function returns the resulting norm of reaction component %vector.

  Created by Tomas Koudelka, 08.2021
*/
double par_compute_react_norm()
{
  long i, j, k, nid, ndofn;
  double normi;
  double totnorm;
  long nin = Psolm->ptopjk->nin; // number of internal nodes
  long nbn = Psolm->ptopjk->nbn; // number of boundary/interface nodes
  long n = Psolm->ptopjk->maxnbnwcd*Psolm->ptopjk->maxndofbn + 1;
  vector lv(n);

  // compute contributions from internal nodes with constraints on the given subdomain
  normi = 0.0;
  for (i=0L; i<nin; i++){
    nid = Psolm->ptopjk->inid[i];
    if (Mt->nodes[nid].react == 1){
      ndofn = Mt->give_ndofn(nid);
      for (j=0L; j<ndofn; j++){
        normi += sqr(Mt->nodes[nid].r[j]);
      }
    }
  }
  // store norm contribution from internal nodes
  k = 0;
  lv(k) = normi;
  k++;
  // reaction components at boundary nodes with constraints on the given subdomain
  for(i=0; i<nbn; i++){
    nid = Psolm->ptopjk->bnid[i];
    if (Mt->nodes[nid].react == 1){
      for (j=0L; j<ndofn; j++){
        lv(k) = Mt->nodes[nid].r[j];
        k++;
      }
    }    
  }

  totnorm = Psolm->compute_quantfluxres_norm(lv, Out);
  
  return totnorm;  
}
