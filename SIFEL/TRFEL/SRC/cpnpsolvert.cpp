#include "cpnpsolvert.h"
#include "globalt.h"
#include "globmatt.h"
#include "transprint.h"
#include <string.h>

/**
   function solves nonstationary problems with growing number of elements
   
   JK, 5.3.2006
   corrected by TKr 11.10.2006
*/
void solve_nonstationary_growing_problem ()
{
  long i,j,n,nce,ncd,lcid;
  double s,dt,prev_time,end_time,alpha,*d,*p,*lhs,*tdlhs,*rhs;
  
  n=Ndoft;
  lcid=0;
  
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  p = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  
  alpha=Tp->alpha;

  //  initial time
  Tp->time=Tp->timecont.starttime ();
  prev_time=Tp->time;
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;

  //  nodes - integration points interpolation
  approximation ();
  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  //saving values for the first iteration
  Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

  do{

    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);

    //  searching for new elements
    //  only new elements are switched on
    //  former elements are temporarily switched off
    nce = Gtt->search_newelem (Tp->time,prev_time);
    
    //  marking of elements switched on
    Gtt->update_elem (Tp->time);

    //  searching for new DOFs
    //  only new DOFs are switched on
    //  former DOFs are temporarily switched off
    ncd = Gtt->search_newdofs (Tp->time,prev_time);
     
    //  marking of DOFs switched on
    Gtt->update_dofs (Tp->time);
    
    //  generation of new code numbers
    Ndoft = Gtt->codenum_generation (Outt);
    n=Ndoft;
        
    Tm->updateipval ();
    i++;

    if (nce!=0 || ncd!=0){
      //  assembling of displacement vector for new code numbers
      Tt->lhs_restore (lhs,Lsrst->lhsi,tdlhs);
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
    if (nce>0){
      Tt->initial_nodval ();
    }

    //  capacity matrix    
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    

    //fprintf (Outt,"\n\n\n CAPACITY MATRIX \n\n\n");
    //Cmat->printmat (Outt);
    //Cmat->printdiag (Outt);
    //fprintf (Outt,"\n\n\n CONDUCTIVITY MATRIX \n\n\n");
    //Kmat->printmat (Outt);
    //Kmat->printdiag (Outt);

    
    for (j=0;j<n;j++){
      p[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    Kmat->gmxv (p,d);
    
    Kmat->scalgm (dt*alpha);
    
    Kmat->addgm (1.0,*Cmat);
    
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    for (j=0;j<n;j++){
      rhs[j] = rhs[j] - d[j];
      d[j]=tdlhs[j];
    }
    
    //Kmat->solve_system (Gtt,tdlhs,rhs);
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,rhs,Outt);
    
    for (j=0;j<n;j++){
      s=(1.0-alpha)*d[j]+alpha*tdlhs[j];
      lhs[j]+=dt*s;
    }
    

    //  nodes - integration points interpolation
    solution_correction ();

    Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

    approximation ();
    
    //  time increment
    prev_time=Tp->time;
    Tp->time=Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] p;  delete [] d;
  print_closet();
}



/**
   function solves nonstationary problems with growing number of elements
   
   TKr, 26.3.2006
*/
void solve_nonstationary_growing_vform ()
{
  long i,j,n,nce,ncd,lcid;
  double dt,prev_time,end_time,alpha,*d,*p,*lhs,*tdlhs,*rhs;
  
  n=Ndoft;
  lcid=0;
  
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  p = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (d,n);
  nullv (p,n);
  
  alpha=Tp->alpha;

  //  initial time
  Tp->time=Tp->timecont.starttime ();
  prev_time=Tp->time;
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  i=0;

  //  nodes - integration points interpolation
  approximation ();
  Tm->initmaterialmodels();
  print_initt(-1, "wt");
  print_stept(0,i,Tp->time,NULL);

  //saving values for the first iteration
  Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

  do{

    if (Mesprt != 0)  fprintf (stdout,"\n\n increment number %ld,   time %f,   time step = %f\n",i,Tp->time,dt);

    //  searching for new elements
    //  only new elements are switched on
    //  former elements are temporarily switched off
    nce = Gtt->search_newelem (Tp->time,prev_time);
    
    //  marking of elements switched on
    Gtt->update_elem (Tp->time);

    //  searching for new DOFs
    //  only new DOFs are switched on
    //  former DOFs are temporarily switched off
    ncd = Gtt->search_newdofs (Tp->time,prev_time);
     
    //  marking of DOFs switched on
    Gtt->update_dofs (Tp->time);
    
    //  generation of new code numbers
    Ndoft = Gtt->codenum_generation (Outt);
    n=Ndoft;    

    Tm->updateipval ();
    i++;
    
    if (nce!=0 || ncd!=0){
      //  assembling of displacement vector for new code numbers
      Tt->lhs_restore (lhs,Lsrst->lhsi,tdlhs);
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
    if (nce>0){
      Tt->initial_nodval ();
    }

    //  capacity matrix    
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    

    //fprintf (Outt,"\n\n\n CAPACITY MATRIX \n\n\n");
    //Cmat->printmat (Outt);
    //Cmat->printdiag (Outt);
    //fprintf (Outt,"\n\n\n CONDUCTIVITY MATRIX \n\n\n");
    //Kmat->printmat (Outt);
    //Kmat->printdiag (Outt);

    
    for (j=0;j<n;j++){
      d[j]=lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  auxiliary vector
    Kmat->gmxv (d,p);
    
    Kmat->scalgm (dt*alpha);
    
    Kmat->addgm (1.0,*Cmat);
    
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    for (j=0;j<n;j++){
      rhs[j] = rhs[j] - p[j];
    }
    
    //Kmat->solve_system (Gtt,tdlhs,rhs);
    Tp->ssle->solve_system (Gtt,Kmat,tdlhs,rhs,Outt);
    
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }

    //  nodes - integration points interpolation
    solution_correction ();

    Tt->lhs_save (lhs,Lsrst->lhsi,tdlhs);

    approximation ();
    
    //  time increment
    prev_time=Tp->time;
    Tp->time=Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

  }while(Tp->time<end_time);
  
  delete [] p;  delete [] d;
  print_closet();
}


