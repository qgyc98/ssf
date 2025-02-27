#include "ndsolver.h"
#include "edsolver.h"
#include "global.h"
#include "probdesc.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "gmatrix.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "mechprint.h"
#include "elemswitch.h"

void solve_nonlinear_dynamics ()
{
  switch (Mp->tforvib){
  case newmark:{
    nonlin_newmark_method ();
    break;
  }
  default:{
    fprintf (stderr,"\n\n unknown solver of nonlinear dynamics is required");
    fprintf (stderr,"\n in function solve_nonlinear_dynamics () (file %s, line %d).\n",__FILE__,__LINE__);
  }
  }
}

/**
   Nonlinear Newmark method for solution of nonlinear dynamic problems
   the method is implicit, equilibrium is solved by the Newton-Raphson method
   
   tam co byla alpha je nyni beta
   tam co byla delta je nyni gamma
   
   5.3.2007, JK
*/
void nonlin_newmark_method ()
{
  long i,j,k,n,ini,lcid;
  double beta,gamma,time,dt,end_time,err,norres;
  double *a,*v,*p,*dd,*vv,*lhs,*rhs,*fs,*brhs;
  
  //  number of degrees of freedom
  n = Ndofm;
  //  coefficient of the method
  beta=Mp->alphafvn;
  //  coefficient of the method
  gamma=Mp->deltafvn;
  
  
  //  maximum number of iterations in one increment
  ini = Mp->nlman->niilnr;
  //  required error in inner loop
  err = Mp->nlman->errnr;
  
  
  //  initial time
  time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  first time increment
  dt = Mp->timecon.initialtimeincr ();
  
  //  load case id
  //  this solver enables only one load case with respect to nonlinearity
  lcid=0;
  
  //  nodal values - displacements
  //  (assigment of initial values)
  lhs = Lsrs->give_lhs (lcid);
  //  time derivatives of nodal values - velocities
  //  (assigment of initial values)
  v = Lsrs->give_tdlhs (lcid);
  //  vector of dynamic load
  rhs = Lsrs->give_rhs (2*lcid);
  //  vector of static load
  fs = Lsrs->give_rhs (2*lcid+1);
  //  backup of the right hand side
  brhs = new double [n];

  //  vector of acceleration
  a = Lsrs->give_stdlhs (lcid);
  nullv (a,n);
  //  predictor of displacements
  dd = new double [n];
  //  predictor of velocities
  vv = new double [n];
  //  auxiliary vector
  p = new double [n];
  
  //  assembling of static loads  
  nullv (fs+lcid*n, n);
  Mb->lc[lcid].assemble (lcid,fs,NULL,1.0);
  
  //  number of step
  i=0;
  print_init(-1, "wt");  
  print_step(lcid, i, Mp->time, rhs);
  print_flush();
  
  
  //  time initialization
  Mp->time = Mp->timecon.newtime ();
  //  time increment
  dt = Mp->timecon.actualbacktimeincr ();
  
  
  //  loop over time interval
  do{
    
    //  stiffness matrix
    stiffness_matrix (lcid);
    //  mass matrix
    mass_matrix (lcid);
    //  damping matrix
    damping_matrix (lcid);
    
    if (Mespr==1)  fprintf (stdout,"\n                                time  %f",Mp->time);
    
    //  assembling of dynamic loads
    nullv (rhs+lcid*n, n);
    Mb->dlc[lcid].assemble (lcid,rhs,NULL,n,Mp->time);
    
    //  total load
    addv (rhs,fs,n);
    
    //  modification of vector of prescribed nodal forces
    cmulv (dt*dt*beta,rhs,n);
    //  copy of the actual prescribed forces
    copyv (rhs,brhs,n);

    
    //  predictor vectors
    for (j=0;j<n;j++){
      dd[j] = lhs[j] + dt*v[j] + dt*dt*(0.5-beta)*a[j];
      vv[j] = v[j] + dt*(1.0-gamma)*a[j];
    }
    
    //  M.dd
    Mmat->gmxv (dd,p);
    
    //  contribution to the right hand side
    addv (rhs,p,n);
    
    //  C.dd
    Dmat->gmxv (dd,p);
    cmulv (dt*gamma,p,n);
    addv (rhs,p,n);
    
    //  C.vv
    Dmat->gmxv (vv,p);
    cmulv (-1.0*dt*dt*beta,p,n);
    addv (rhs,p,n);
    
    //  dynamic stiffness matrix
    //  (stored in array for stiffness matrix)
    //  M + gamma.dt.C + beta.dt.dt.K
    Smat->scalgm (dt*dt*beta);
    Smat->addgm (1.0,*Mmat);
    Smat->addgm (dt*gamma,*Dmat);
    
    
    //  solution of system of algebraic equations
    Mp->ssle->solve_system (Gtm,Smat,lhs,rhs,Out);
    
    for (j=0;j<n;j++){
      a[j] = 1.0/dt/dt/beta*(lhs[j]-dd[j]);
      v[j] = vv[j] + gamma*dt*a[j];
    }
    
    
    //  restoration of the prescribed forces
    copyv (brhs,rhs,n);
    
    
    //  inertial forces
    //  M.a
    Mmat->gmxv (a,p);
    cmulv (-1.0*dt*dt*beta,p,n);
    addv (rhs,p,n);
    
    //  damping forces
    //  C.v
    Dmat->gmxv (v,p);
    cmulv (-1.0*dt*dt*beta,p,n);
    addv (rhs,p,n);
    
    //  internal forces
    internal_forces (lcid,p);
    cmulv (-1.0*dt*dt*beta,p,n);
    addv (rhs,p,n);
    
    //  norm of the residual vector
    norres = ss (rhs,rhs,n);
    fprintf (stdout,"\n norm of residual vector  %le",norres);
    if (norres>err){
      //  norm of residual vector is greater than required tolerance
      //  the Newton-Raphson method has to be used for reaching equilibrium
      for (j=0;j<ini;j++){
	
	//  solution of system of algebraic equations
	Mp->ssle->solve_system (Gtm,Smat,lhs,rhs,Out);
	
	for (k=0;k<n;k++){
	  lhs[k]+=p[k];
	  a[k] = 1.0/dt/dt/beta*(lhs[k]-dd[k]);
	  v[k] = vv[j] + gamma*dt*a[k];
	}
	
	//  restoration of the prescribed forces
	copyv (brhs,rhs,n);
	
	//  inertial forces
	//  M.a
	Mmat->gmxv (a,p);
	cmulv (-1.0*dt*dt*beta,p,n);
	addv (rhs,p,n);
	
	//  damping forces
	//  C.v
	Dmat->gmxv (v,p);
	cmulv (-1.0*dt*dt*beta,p,n);
	addv (rhs,p,n);
	
	//  internal forces
	internal_forces (lcid,p);
	cmulv (-1.0*dt*dt*beta,p,n);
	addv (rhs,p,n);
	
	//  norm of the residual vector
	norres = ss (rhs,rhs,n);
	fprintf (stdout,"\n norm of residual vector  %le",norres);
	if (norres<err){
	  break;
	}
      }
      
    }
    
    if (j==ini){
      //  time correction, time step is reduced two times
      //Mp->timecon.reduce_step (2.0);
      //  new time
      Mp->time = Mp->timecon.actualtime ();
      //  new time increment
      dt = Mp->timecon.actualforwtimeincr ();
    }
    else{
      compute_req_val (lcid);
      print_step(lcid, i, Mp->time, rhs);
      print_flush();
      
      //  time increment
      Mp->time = Mp->timecon.newtime ();
      dt = Mp->timecon.actualbacktimeincr ();
      
      //  step number increment
      i++;
    }
    
  }while(Mp->time<end_time);
  
  print_close();  
  delete [] p;
  delete [] dd;
  delete [] vv;
}


