#include "cpnnpsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "elemswitcht.h"
#include "transprint.h"
#include <string.h>

/**
   Function solves system of non-linear TRFEL algebraic equations by Newton-Raphson method
   for time-dependent problems.
   
   TKr 6.4.2007
*/
void solve_nonstationary_growing_problem_nonlin ()
{
  long i,j,k,nt,tncd,tnce,ini;
  double zero,dt,prev_timet,end_time,alpha;
  double *d,*p,*lhst,*tdlhst,*rhst;
  double *fbt,*fit,*lhst_last;
  //double s;
  double norf_last,err,norf;

  //  maximum number of iterations in inner loop
  ini = Tp->nii;
  //  required norm of vector of unbalanced forces
  err = Tp->err;
  //  number of transport degrees of freedom
  nt=Ndoft;

  //  left hand side
  lhst = Lsrst->give_lhs (0);
  tdlhst = Lsrst->give_tdlhs (0);
  //  right hand side
  rhst = Lsrst->give_rhs (0);
  
  //  array containing predictors
  d = new double [nt];
  //  auxiliary vector 
  p = new double [nt];
  fbt = new double [nt];
  fit = new double [nt];
  lhst_last = new double [nt];

  //  initial values
  nullv (lhst,nt);
  nullv (tdlhst,nt);
  nullv (d,nt);
  nullv (p,nt);
  nullv (fbt,nt);
  nullv (fit,nt);

  alpha=Tp->alpha;
  zero=Tp->zero;
  
  //  starting time
  Tp->time=Tp->timecont.starttime ();
  //  time increment
  dt=Tp->timecont.initialtimeincr ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initiation of the variable prev_timet
  prev_timet=Tp->time;

  //  initiation of transport material models
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation ();  

  // **************************
  //  main iteration loop  ***
  // **************************
  i=0;
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);
  print_flusht();  

  //saving values for the first iteration
  Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);

  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);
    
    //  searching for new elements
    //  only new elements are switched on
    //  former elements are temporarily switched off
    tnce = Gtt->search_newelem (Tp->time,prev_timet);
    
    //  marking of elements switched on
    Gtt->update_elem (Tp->time);
    
    //  searching for new DOFs
    //  only new DOFs are switched on
    //  former DOFs are temporarily switched off
    tncd = Gtt->search_newdofs (Tp->time,prev_timet);
 
    //  marking of DOFs switched on
    Gtt->update_dofs (Tp->time);
    
    //  generation of new code numbers
    Ndoft = Gtt->codenum_generation (Outt);
    nt=Ndoft;
            
    Tm->updateipval ();
    i++;

    if (tnce!=0 || tncd!=0){
      //  assembling of displacement vector for new code numbers
      Tt->lhs_restore (lhst,Lsrst->lhsi,tdlhst);
    }

    //  cleaning of the conductivity matrix
    if (Kmat != NULL){
      delete Kmat;
      Kmat=NULL;
    }
    //  cleaning of the capacity matrix
    if (Cmat != NULL){
      delete Cmat;
      Cmat=NULL;
    }

 
    //  assignment of initial nodal values
    //  it works for problem with growing number of elements
    //  it does not work generally for problems with increasing and decreasing number of elements
    if (tnce>0){
      Tt->initial_nodval ();
    }

    //  conductivity matrix
    conductivity_matrix (0);
    
    //  capacity matrix
    capacity_matrix (0);
    
    for (j=0;j<nt;j++){
      p[j] = lhst[j] + (1.0-alpha)*dt*tdlhst[j];
    }
    
    //  right hand side of transport part
    trfel_right_hand_side (0,rhst,nt);
    
    //  auxiliary vector
    Kmat->gmxv (p,d);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    for (j=0;j<nt;j++){
      rhst[j] = rhst[j] - d[j];
    }
    
    //  solution of the system of algebraic equations
    Tp->ssle->solve_system (Gtt,Kmat,tdlhst,rhst,Outt);

    for (j=0;j<nt;j++){
      //verze1:
      //s=(1.0-alpha)*d[j]+alpha*tdlhst[j];
      //lhst[j]+=dt*s;
      //verze2:
      lhst[j] = p[j] + alpha*dt*tdlhst[j];
    }
    
    //  physically corrected solution
    solution_correction ();    

    Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);

    //  approximation of nodal values into ontegration points
    approximation ();

    norf_last = 1.0e20;
    //inner iteration loop//////////////////
    for (j=0;j<ini;j++){
      
      //  physically corrected solution
      solution_correction ();
      
      Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
      
      //  approximation of nodal values into ontegration points
      approximation ();
      
      // full newton-raphson
      if (Tp->trsolv == fullnewtont){
	// matrices are computing in each inner iteration
	//  capacity matrix
	capacity_matrix (0);
	
	//  conductivity matrix
	conductivity_matrix (0);
	
	//  auxiliary vector  K (d+(1-alpha)*dt*v)
	Kmat->gmxv (p,d);
	
	//  matrix of the system of equations
	//  C + alpha.dt.K
	Kmat->scalgm (dt*alpha);
	Kmat->addgm (1.0,*Cmat);
      }
      
      if (Tp->trestype==lrhst){
	//  right hand side of transport part
	trfel_right_hand_side (0,rhst,nt);
	// Solver computes residuum from system of equations
	Kmat->gmxv (tdlhst,fit);
	// Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
	for (k=0;k<nt;k++){
	  fbt[k] = rhst[k] - d[k] - fit[k];
	}
      }
      
      if (Tp->trestype==fluxest){
	// Solver computes unbalanced fluxes
	internal_fluxes (fit,nt);
	//  vector of unbalanced fluxes
	for (k=0;k<nt;k++){
	  fbt[k]=fit[k];
	}
      }
      
      norf = ss (fbt,fbt,nt);
      
      if (Mesprt==1)  fprintf (stdout,"\n inner iteration %ld   error %e",j,norf);
      
      if (norf<err) break;
      
      Tp->ssle->solve_system (Gtt,Kmat,d,fbt,Outt);
      
      for (k=0;k<nt;k++){
	tdlhst[k]+=d[k];
	lhst[k]+=alpha*dt*d[k];
      }
      
      //  physically corrected solution
      solution_correction ();
      
      Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
      
      //  approximation of nodal values into ontegration points
      approximation ();
      
      //  convergence control for newton-raphson
      //  condition of decreesing of error
      if (Tp->convergcontrolt==yes){
	if (norf > norf_last){
	  for (k=0;k<nt;k++){
	    lhst[k] = lhst_last[k]; //storing results from previous inner step
	  }
	  
	  //  physically corrected solution
	  solution_correction ();
	  
	  Tt->lhs_save (lhst,Lsrst->lhsi,tdlhst);
	  
	  //  approximation of nodal values into ontegration points
	  approximation ();	  

	  if (Mesprt==1)  fprintf (stdout,"\n\n convergence control: inner iteration skiped %ld error %e\n\n",j,norf);	  

	  break;
	}
	for (k=0;k<nt;k++){
	  lhst_last[k] = lhst[k]; //storing results from previous inner step
	}
	
	norf_last = norf; //storing norf from previous inner step
      }
    }
    ////////////////////////////////////////
    //  new time and new time increment
    prev_timet=Tp->time;
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();  
    
  }while(Tp->time<=end_time);
  
  print_closet ();
  
  delete [] p;
  delete [] d;

  delete [] fit;
  delete [] fbt;
  delete [] lhst_last;
}
