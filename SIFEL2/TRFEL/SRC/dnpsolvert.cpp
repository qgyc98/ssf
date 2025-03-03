#include "dnpsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "transprint.h"
#include "saddpoint.h"
#include <string.h>

/**
   function solves nonstationary transport problem with
   discontinuities across material interfaces
   
   function selects a suitable numerical scheme for solution
   
   JK
*/
void solve_discont_nonstationary_problem ()
{
  lin_nonstat_dform ();

  //lin_nonstat_dform_resistance ();

}

/**
   function solves nonstationary transport problem with discontinuity
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   JK, 6.8.2008
*/
void lin_nonstat_dform ()
{
  long i,j,n,nbdof,nsad,lcid;
  double dt,end_time,alpha,zero,totsalt;
  double *f,*d,*p,*lhs,*tdlhs,*rhs;
  
  //  load case id must be equal to zero in this type of problem
  lcid=0;
  
  //  the number of primal unknowns / primal degrees of freedom
  //  it does not contain the number of multipliers
  n=Ndoft;
  //  the number of primal boundary/interface unknowns
  nbdof=Gtt->nbdof;
  //  the number of unknowns in reduced problem
  //  the sum of the number of primal boundary/interface unknowns and the number of Lagrange multipliers
  nsad=Gtt->nsad;
  
  //  computer zero
  zero=Tp->zero;
  
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  //  step number
  i=0;
  
  //  nodes - integration points interpolation
  approximation ();

  // vlozil JM 15.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  
  //  initiation of material models
  Tm->initmaterialmodels();

  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime (dt);
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    //  predictor d
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j] = lhs[j]+dt*(1.0-alpha)*tdlhs[j];
    }
    
    //  C.dd
    Cmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed fluxes f_{n+1}
    trfel_right_hand_side (lcid,f,n);
    
    //  C.d + alpha dt f_{n+1}
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
    }

    Tt->compute_jumps (rhs);

    
    Gtt->mult_localization (Kmat);
    
    
    Kmat->diag_check (zero,rhs);

    //  solution of the system of algebraic equations
    //  (C + alpha dt K) d_{n+1} = alpha dt f_n+1 + C dd
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);


    //  first time derivatives of nodal values
    // v_{n+1} = 1/alpha/dt (d_{n+1}-dd)
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();
    
    totsalt = total_integral_ip (5);
    fprintf (stdout,"\n total salt content  %lf",totsalt);
    totsalt = total_integral_ip (6);
    fprintf (stdout,"\n integral of cf      %lf",totsalt);
    totsalt = total_integral_ip (7);
    fprintf (stdout,"\n integral of cc      %lf",totsalt);
    totsalt = total_integral_ip (8);
    fprintf (stdout,"\n integral of cb      %lf",totsalt);
    
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
  }while(Tp->time<end_time);
  
  delete [] p;
  delete [] d;
  delete [] f;
  
  print_closet();
}














/**
   trial subroutine devoted to solution of problems with resistance factor
   
   function solves nonstationary transport problem with discontinuity
   the d-form version of the generalized trapezoidal method
   see T.J.R. Hughes: The Finite Element Method; Linear Static
   and Dynamic Finite Element Analysis
   
   JK, 2.2.2011
*/
void lin_nonstat_dform_resistance ()
{
  long i,j,n,ii,ipp,nbdof,nsad,lcid;
  double dt,end_time,alpha,zero;
  double *f,*d,*p,*lhs,*tdlhs,*rhs;
  
  //  load case id must be equal to zero in this type of problem
  lcid=0;
  
  //  the number of primal unknowns / primal degrees of freedom
  //  it does not contain the number of multipliers
  n=Ndoft;
  //  the number of primal boundary/interface unknowns
  nbdof=Gtt->nbdof;
  //  the number of unknowns in reduced problem
  //  the sum of the number of primal boundary/interface unknowns and the number of Lagrange multipliers
  nsad=Gtt->nsad;
  
  //  computer zero
  zero=Tp->zero;
  
  
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  nullv (lhs,n);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  nullv (tdlhs,n);
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  nullv (rhs,n);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  
  //  coefficient in generalized trapezoidal integration rule
  alpha=Tp->alpha;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  //  step number
  i=0;
  
  //  nodes - integration points interpolation
  approximation ();

  // vlozil JM 15.1.2008
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  
  
  //  initiation of material models
  Tm->initmaterialmodels();

  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    //  determination of new time instance
    Tp->time=Tp->timecont.newtime (dt);
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    //  update of step number
    i++;
    
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    //  predictor d
    //  dd = d_n + (1-alpha) dt v_n
    for (j=0;j<n;j++){
      d[j] = lhs[j]+dt*(1.0-alpha)*tdlhs[j];
    }
    
    //  C.dd
    Cmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    
    //  prescribed fluxes f_{n+1}
    trfel_right_hand_side (lcid,f,n);
    
    //  C.d + alpha dt f_{n+1}
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
    }

    Tt->compute_jumps (rhs);

    
    Gtt->mult_localization (Kmat);
    
    
    //  solution of the system of algebraic equations
    //  (C + alpha dt K) d_{n+1} = alpha dt f_n+1 + C dd
    //  nodal values d_{n+1} are solved
    Tp->ssle->solve_system (Gtt,Kmat,lhs,rhs,Outt);


    //  first time derivatives of nodal values
    // v_{n+1} = 1/alpha/dt (d_{n+1}-dd)
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    //  physically corrected solution
    solution_correction ();    
    //  approximation of nodal values into ontegration points
    approximation ();
    
    
    for (ii=0;i<Tt->ne;ii++){
      //  integration point id
      ipp=Tt->elements[ii].ipp[0][0];
      
      Tm->computenlfluxes (0,ipp);
      Tm->computenlfluxes (1,ipp);
    }
    Tt->compute_resistance_factor (rhs);
    




    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();
    
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
    
  }while(Tp->time<end_time);
  
  delete [] p;
  delete [] d;
  delete [] f;
  
  print_closet();
}

