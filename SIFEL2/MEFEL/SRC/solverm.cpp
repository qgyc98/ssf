#include "solverm.h"
#include "lssolver.h"
#include "edsolver.h"
#include "fdsolver.h"
#include "nssolver.h"
#include "mtsolver.h"
#include "llssolver.h"
#include "epsolver.h"
#include "slsolver.h"
#include "lfssolver.h"
#include "nfssolver.h"
#include "cpsolver.h"
#include "hvisolver.h"
#include "lbsolver.h"
#include "homog.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "lhsrhs.h"
#include "homogmech.h"
#include "globmat.h"
#include "elemswitch.h"



void solve_mefel_problem ()
{
  if (Mp->stochasticcalc==0){
    solve_mefel_deterministic_problem ();
  }
  if (Mp->stochasticcalc==1 || Mp->stochasticcalc==2 || Mp->stochasticcalc==3){
    solve_mefel_stochastic_problem ();
  }
}



void solve_mefel_deterministic_problem ()
{
  if(Mm->elisopd){
    double *rhs = Lsrs->give_rhs(0);
    double *lhs = Lsrs->give_lhs(0);
    double *fl;
    //  vector of prescribed force loads, it does not contain forces caused by temperature, etc.
    fl   = new double [Ndofm];
    if (Smat == NULL){ // matrix has not yet been assembled (changed number of DOFs) or forced update is required
      stiffness_matrix (0);
    }
    double old_time = Mp->time;
    Mp->time = 1.0;
    mefel_right_hand_side (0, rhs, fl);
    //  solution of K(r).v=F
    Mp->ssle->solve_system (Gtm,Smat,lhs,rhs,Out);
    Mp->time = old_time;
    compute_ipstrains(0);
    //stress_initdispl(0);   
    if (Mp->tprob == linear_statics)
      Mm->initmaterialmodels(0, false);

    nullv(lhs, Ndofm);
    nullv(rhs, Ndofm);
    delete [] fl;
  }
  switch (Mp->tprob){
  case linear_statics:{
    solve_linear_statics ();
    break;
  }
  case eigen_dynamics:{
    solve_eigen_dynamics (Lsrs->lhs,Lsrs->w);
    break;
  }
  case linear_stability:{
    solve_linear_stability (Lsrs->lhs,Lsrs->w);
    break;
  }
  case forced_dynamics:{
    solve_forced_dynamics ();
    break;
  }
  case mat_nonlinear_statics:{
    solve_nonlinear_statics ();
    break;
  }
  case mech_timedependent_prob:{
    solve_time_dep_prob ();
    break;
  }
  case growing_mech_structure:{
    solve_prob_constr_phases ();
    //solve_prob_constr_phases_old ();
    break;
  }
  case earth_pressure:{
    solve_epressure ();
    break;
  }
  case layered_linear_statics:{
    solve_layered_linear_statics ();
    break;
  }
  case lin_floating_subdomain:{
    //solve_linear_statics ();
    solve_linear_floating_subdomains ();
    //solve_nonlinear_floating_subdomains ();
    break;
  }
  case nonlin_floating_subdomain:{
    solve_incremental_floating_subdomains ();
    break;
  }
  case hemivar_inequal:{
    solve_hemivariational_inequalities ();
    break;
  }
  case load_balancing:{
    solve_load_balancing ();
    break;
  }
  default:{
    print_err("unknown problem type is required", __FILE__, __LINE__, __func__);
  }
  }
  
  //  homogenization has different strategy than in TRFEL
  if (Mp->homog!=0){
    switch (Mp->homog){
    case 1:{//new version 07/01/2015
      homogenization (Out,0);//mechanical_homogenization();
      break;
    }
    case 2:{
      homogenization (Out,0);
      break;
    }
    case 3:
    case 4:
    case 9:{
      break;
    }
    case 5:
    case 6:
    case 7:
    case 8:{
      mechanical_homogenization();
      
      // only for testing:
      /* double *matrix_k;
	 matrix_k = new double [6];
	 
	 paral_mechanical_homogenization (NULL, matrix_k);
      */
      break;
    }
    default:{
      print_err("unknown homogenization problem type is required",__FILE__,__LINE__,__func__);
    }   
    }
  } 
}



void solve_mefel_stochastic_problem ()
{
  long i,j;
  stochdriver *stochd=St;
    
  for (i=0;i<stochd->nsampl;i++){
    //fprintf (Out,"\n\n %ld\n",i);
    fprintf (stdout,"\n sample %ld",i);
    //Mp->ns=i;
    stochd->changevalues (i);
    Lsrs->clean_lhs ();
    Mm->clean_ip ();
    Mt->clean_nodes ();
    solve_mefel_deterministic_problem ();
    stochd->extractor ();
    stochd->save_results (i);
    //  presunuto do extractoru
    //stochd->update_auxparam ();
  }
  


  if (Mp->stochasticcalc==1){
    stochd->writetable ();
  }
  if (Mp->stochasticcalc==3){
    //for (i=0;i<stochd->nprunknowns;i++){
    //fprintf (Out,"\n\n fuzzy number %ld",i);
    //stochd->fn[i].print (Out);
    //}
    
    long k=0;
    for (i=0;i<stochd->npnd;i++){
      for (j=0;j<stochd->pnd[i].num;j++){
	//fprintf (Out,"\n\n fuzzy posun %ld",i);
	stochd->fn[k].print (stochd->datout);
	k++;
      }
    }
    
    /*
    long k;
    for (i=0;i<stochd->npev;i++){
      fprintf (Out,"\n\n");
      
      for (k=0;k<6;k++){
	fprintf (Out,"\n");
	stochd->fn[j].print (Out);
	j++;
      }
    }
    */
    
    
    /*
    fprintf (Out,"\n\n\n\n");
    for (i=0;i<stochd->npev*6;i+=3){
      fprintf (Out,"\n");
      stochd->fn[j+i].print (Out);
    }
    fprintf (Out,"\n\n\n\n");
    for (i=1;i<stochd->npev*6;i+=3){
      fprintf (Out,"\n");
      stochd->fn[j+i].print (Out);
    }
    fprintf (Out,"\n\n\n\n");
    for (i=2;i<stochd->npev*6;i+=3){
      fprintf (Out,"\n");
      stochd->fn[j+i].print (Out);
    }

    */
    
    
  }

  //stochd->diagpostproc ();

}
