#include "fdsolver.h"
#include "edsolver.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "mechtop.h"
#include "mechbclc.h"
#include "lhsrhs.h"
#include "globmat.h"
#include "gmatrix.h"
#include "loadcase.h"
#include "dloadcase.h"
#include "mechprint.h"
#include "elemswitch.h"



#include "intpoints.h"
#include "element.h"


void solve_forced_dynamics ()
{
  long i;
  
  for (i=0;i<Lsrs->nlc;i++){
    switch (Mp->tforvib){
    case newmark:{
      //newmark_method (i);
      optim_newmark_method (i);
      break;
    }
    case findiff:{
      difference_method (i);
      break;
    }
    case explfindiff:{
      //explicit_difference_method (i);
      difference_method_2 (i);
      break;
    }
    case modal_analysis:{
      response_spectrum_method (i);
      break;
    }
    default:{
      print_err("unknown solver of forced vibration is required", __FILE__, __LINE__, __func__);
    }
    }
  }
  
}

/**
   Newmark method for solution of forced dynamic problems
   the method is implicit
   
   @param lcid - load case id
   
   JK, revised 26.7.2005
*/
void newmark_method (long lcid)
{
  long i,j,n;
  double alpha,delta,time,dt,end_time;
  double *a,*v,*p,*dd,*vv,*lhs,*rhs,*fs;
  
  //  number of degrees of freedom
  n = Ndofm;
  //  coefficient of the method
  alpha=Mp->alphafvn;
  //  coefficient of the method
  delta=Mp->deltafvn;
  
  //  initial time
  time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  first time increment
  dt = Mp->timecon.initialtimeincr ();
  
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
    Mb->dlc[lcid].assemble (lcid,rhs,NULL,Ndofm,Mp->time);
    
    //  total load
    addv (rhs,fs,n);
    
    //  modification of vector of prescribed nodal forces
    cmulv (dt*dt*alpha,rhs,n);

    //  predictor vectors
    for (j=0;j<n;j++){
      dd[j] = lhs[j] + dt*v[j] + dt*dt*(0.5-alpha)*a[j];
      vv[j] = v[j] + dt*(1.0-delta)*a[j];
    }
    
    //  M.dd
    Mmat->gmxv (dd,p);
    
    //  contribution to the right hand side
    addv (rhs,p,n);
    
    //  C.dd
    Dmat->gmxv (dd,p);
    cmulv (dt*delta,p,n);
    addv (rhs,p,n);
    
    //  C.vv
    Dmat->gmxv (vv,p);
    cmulv (-1.0*dt*dt*alpha,p,n);
    addv (rhs,p,n);
    
    //  dynamic stiffness matrix
    //  (stored in array for stiffness matrix)
    //  M + delta.dt.C + alpha.dt.dt.K
    Smat->scalgm (dt*dt*alpha);
    Smat->addgm (1.0,*Mmat);
    Smat->addgm (dt*delta,*Dmat);
    
    
    //  solution of system of algebraic equations
    //Smat->solve_system (Gtm,lhs,rhs);
    Mp->ssle->solve_system (Gtm,Smat,lhs,rhs,Out);

    for (j=0;j<n;j++){
      v[j] = vv[j] + delta/dt/alpha*(lhs[j]-dd[j]);
      a[j] = 1.0/dt/dt/alpha*(lhs[j]-dd[j]);
    }
    
    compute_req_val (lcid);
    print_step(lcid, i, Mp->time, rhs);
    print_flush();
    
    //  time increment
    Mp->time = Mp->timecon.newtime ();
    dt = Mp->timecon.actualbacktimeincr ();
    i++;

  }while(Mp->time<end_time);

  print_close();  
  delete [] p;
  delete [] dd;
  delete [] vv;
}


/**
   Newmark method for solution of forced dynamic problems
   the method is implicit
   
   @param lcid - load case id
   
   JK, 11. 4. 2024
*/
void optim_newmark_method (long lcid)
{
  long i,j,n;
  double alpha,delta,time,dt,pdt,end_time,zero=1.0e-10;
  double *a,*v,*p,*dd,*vv,*lhs,*rhs,*fs;
  
  //  number of degrees of freedom
  n = Ndofm;
  //  coefficient of the method
  alpha=Mp->alphafvn;
  //  coefficient of the method
  delta=Mp->deltafvn;
  
  //  initial time
  time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  first time increment
  dt = Mp->timecon.initialtimeincr ();
  
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
  

  i=0;
  print_init(-1, "wt");  
  print_step(lcid, i, Mp->time, rhs);
  print_flush();
  
  
  //  time initialization
  Mp->time = Mp->timecon.newtime ();
  //  time increment
  dt = Mp->timecon.actualbacktimeincr ();
  pdt = dt;
  
  //  stiffness matrix
  stiffness_matrix (lcid);
  //  mass matrix
  mass_matrix (lcid);
  //  damping matrix
  damping_matrix (lcid);
  
  
  //  loop over time interval
  do{
    
    if (Mespr==1)  fprintf (stdout,"\n                                time  %f",Mp->time);
    
    if (fabs(pdt-dt)>zero){
      //  the time step was changed
      //  therefore the matrices have to be recalculated
      
      pdt = dt;
      
      //  stiffness matrix
      stiffness_matrix (lcid);
      //  mass matrix
      mass_matrix (lcid);
      //  damping matrix
      damping_matrix (lcid);

      //  dynamic stiffness matrix
      //  (stored in array for stiffness matrix)
      //  M + delta.dt.C + alpha.dt.dt.K
      Smat->scalgm (dt*dt*alpha);
      Smat->addgm (1.0,*Mmat);
      Smat->addgm (dt*delta,*Dmat);
      
    }
    
    //  assembling of dynamic loads
    nullv (rhs+lcid*n, n);
    Mb->dlc[lcid].assemble (lcid,rhs,NULL,Ndofm,Mp->time);
    
    //  total load
    addv (rhs,fs,n);
    
    //  modification of vector of prescribed nodal forces
    cmulv (dt*dt*alpha,rhs,n);
    
    //  predictor vectors
    for (j=0;j<n;j++){
      dd[j] = lhs[j] + dt*v[j] + dt*dt*(0.5-alpha)*a[j];
      vv[j] = v[j] + dt*(1.0-delta)*a[j];
    }
    
    //  M.dd
    Mmat->gmxv (dd,p);
    
    //  contribution to the right hand side
    addv (rhs,p,n);
    
    //  C.dd
    Dmat->gmxv (dd,p);
    cmulv (dt*delta,p,n);
    addv (rhs,p,n);
    
    //  C.vv
    Dmat->gmxv (vv,p);
    cmulv (-1.0*dt*dt*alpha,p,n);
    addv (rhs,p,n);
    
    
    //  solution of system of algebraic equations
    //Smat->solve_system (Gtm,lhs,rhs);
    Mp->ssle->solve_system (Gtm,Smat,lhs,rhs,Out);

    for (j=0;j<n;j++){
      v[j] = vv[j] + delta/dt/alpha*(lhs[j]-dd[j]);
      a[j] = 1.0/dt/dt/alpha*(lhs[j]-dd[j]);
    }
    
    compute_req_val (lcid);
    print_step(lcid, i, Mp->time, rhs);
    print_flush();
    
    //  time increment
    Mp->time = Mp->timecon.newtime ();
    dt = Mp->timecon.actualbacktimeincr ();
    i++;

  }while(Mp->time<end_time);

  print_close();  
  delete [] p;
  delete [] dd;
  delete [] vv;
}


/**
   difference method for solution of forced dynamic problems
   the method is explicit
   
   @param lcid - load case id
   
   JK, 26.7.2005
*/
void difference_method (long lcid)
{
  long i,j,n;
  double time,dt,end_time;
  double *a,*v,*p,*da,*dp,*lhs,*rhs,*fs;
  
  //  number of degrees of freedom
  n = Ndofm; 
  //  initial time
  time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  first time increment
  dt = Mp->timecon.initialtimeincr ();
  
  //  nodal values - displacements at next step
  //  (assigment of initial values)
  lhs = Lsrs->give_lhs (lcid);
  //  time derivatives of nodal values - velocities
  //  (assigment of initial values)
  v = Lsrs->give_tdlhs (lcid);
  //  second time derivatives of nodal values - accelerations
  a = Lsrs->give_stdlhs (lcid);
  //  vector of dynamic load
  rhs = Lsrs->give_rhs (2*lcid);
  //  vector of static load
  fs = Lsrs->give_rhs (2*lcid+1);
  

  //  vector of displacements at previous step
  dp = new double [n];
  //  vector of displacements at actual step
  da = new double [n];
  //  auxiliary vector
  p = new double [n];
  
  //  assembling of static loads  
  nullv (fs+lcid*n, n);
  Mb->lc[lcid].assemble (lcid,fs,NULL,1.0);
  

  i=0;
  print_init(-1, "wt");  
  print_step(lcid, i, Mp->time, rhs);
  print_flush();
  
  
  //  initial values
  for (j=0;j<n;j++){
    da[j] = lhs[j];
    dp[j] = da[j] - dt*v[j];
  }

  
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
    Mb->dlc[lcid].assemble (lcid,rhs,NULL,Ndofm,Mp->time);
    
    //  total load
    addv (rhs,fs,n);
    
    //  K.da
    Smat->gmxv (da,p);
    cmulv (-1.0,p,n);
    addv (rhs,p,n);
    
    //  M.da.2/dt/dt
    Mmat->gmxv (da,p);
    cmulv (2.0/dt/dt,p,n);
    addv (rhs,p,n);
    
    //  M.dp/dt/dt
    Mmat->gmxv (dp,p);
    cmulv (-1.0/dt/dt,p,n);
    addv (rhs,p,n);
    
    //  C.vv
    Dmat->gmxv (dp,p);
    cmulv (1.0/2.0/dt,p,n);
    addv (rhs,p,n);
    
    //  dynamic stiffness matrix
    //  (stored in array for mass matrix)
    //  M/dt/dt + C/2/dt
    Mmat->scalgm (1.0/dt/dt);
    Mmat->addgm (1.0/2.0/dt,*Dmat);
    
    //  solution of system of algebraic equations
    //Mmat->solve_system (Gtm,lhs,rhs);
    Mp->ssle->solve_system (Gtm,Mmat,lhs,rhs,Out);

    for (j=0;j<n;j++){
      v[j]=(lhs[j]-dp[j])/2.0/dt;
      a[j]=(lhs[j]-2.0*da[j]+dp[j])/dt/dt;
      
      dp[j]=da[j];
      da[j]=lhs[j];
    }
    
    //  time increment
    Mp->time = Mp->timecon.newtime ();
    dt = Mp->timecon.actualbacktimeincr ();
    i++;

    compute_req_val (lcid);
    print_step(lcid, i, Mp->time, rhs);
    print_flush();
    

  }while(Mp->time<end_time);

  print_close();  
  delete [] da;
  delete [] dp;
  delete [] p;
}




/**
   function solves system of ordinary differential equations of second order
   such equations are used e.g. in dynamics
   function is intended for molecular dynamics
   
   @param lcid - load case id
   
   JK, 24.7.2005
*/
void verlet_method (long lcid)
{
  long i,n,ni;
  double mm,zero,dt,end_time;
  double *m,*dp,*da,*dn,*fp,*fi;
  
  //  number of rows of the matrix
  n = Ndofm;
  //  number of iterations
  ni = 0;
  //  computer zero
  zero = Mp->zero;
  
  //  vector of displacements at previous step
  dp = new double [n];
  //  vector of displacements at new step
  dn = new double [n];
  //  vector of displacements at actual step
  da = Lsrs->give_lhs (lcid);
  //  vector of prescribed forces
  fp = Lsrs->give_rhs (lcid);
  //  vector of internal forces
  fi = new double [n];
  //  vector of reciprocal values of weights
  m = new double [n];
  
  //  initial time
  Mp->time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  initial time increment
  dt = Mp->timecon.initialtimeincr ();
  
  //  extraction of weights
  for (i=0;i<n;i++){
    mm=Mmat->give_entry (i,i);
    if (mm<zero){

    }
    else{
      m[i]=1.0/mm;
    }
  }
  
  
  //  first step
  //  prescribed forces
  mefel_right_hand_side (i,fp);
  
  //  internal forces
  internal_forces (lcid,fi);
  
  for (i=0;i<n;i++){
    //  initial displacements
    da[i]=Lsrs->lhsi[i];
    //  initial velocities (temporarily stored in vector dn)
    dn[i]=Lsrs->tdlhsi[i];
    //  computation of vector of displacements at time t_init - dt
    dp[i]=da[i]-dt*dn[i]+dt*dt/2.0*m[i]*(fp[i]+fi[i]);
  }
  
  
  print_init(-1, "wt");
  //print_step(lcid, i, Mp->time, f);
  print_flush();
  
  
  do{
    fprintf (stdout,"\n time %e",Mp->time);
    
    //  prescribed forces
    mefel_right_hand_side (i,fp);
    
    //  internal forces
    internal_forces (lcid,fi);
    
    for (i=0;i<n;i++){
      dn[i]=2.0*da[i]-dp[i]+dt*dt*m[i]*(fi[i]+fp[i]);
      dp[i]=da[i];
      da[i]=dn[i];
    }
    
    
    //print_step(lcid, i, Mp->time, f);
    print_flush();
    
    //  time increment
    Mp->time = Mp->timecon.newtime ();
    ni++;
    
    
    
  }while(Mp->time<end_time);
  
  print_close ();

  delete [] dp;
  delete [] dn;
  delete [] da;
  delete [] fp;
  delete [] fi;
  delete [] m;
}

/**
   function solves seismic analysis using response spectrum
   
   @param lcid - load case id
   
   JK, 20.8.2005
*/
void response_spectrum_method (long /*lcid*/)
{
  long i,j,n,m;
  double coeff,per;
  double *v,*auxv;
  
  //  number of DOFs
  n=Ndofm;
  //  number of eigenvectors used in analysis
  m=Mp->eigsol.neigv;
  
  //  computation of eigenvalues and eigenmodes
  solve_eigen_dynamics (Lsrs->lhs+n,Lsrs->w);
  
  v = new double [n];
  auxv = new double [n];
  
  nullv (Lsrs->lhs,n);
  
  for (i=0;i<m;i++){
    //  period (computed from circular frequency)
    per=2.0*3.14159265358979/sqrt(Lsrs->w[i]);
    
    nullv (v,n);
    
    //  vector of directions multiplied by response spectrum values
    Mb->dlc[0].stool.assemble (v,per);
    
    //  mass matrix multiplication
    Mmat->gmxv (v,auxv);
    
    //  dot product with eigenmode
    coeff = ss (auxv,Lsrs->lhs+(i+1)*n,n)/(Lsrs->w[i]);
    
    //fprintf (Out,"\n response spectrum - circ. frequency  %le   period %le     contribution to the %ld-th mode    %le",Lsrs->w[i],per,i,coeff);

    for (j=0;j<n;j++){
      Lsrs->lhs[j]-=coeff*Lsrs->lhs[(i+1)*n+j];
    }
    
  }
  
  
  for (j=0;j<n;j++){
    Lsrs->lhs[j]=Lsrs->lhs[n+j];
  }
  
  /*
  fprintf (Out,"\n\n\n displacements");
  for (i=0;i<n;i++){
    fprintf (Out,"\n %5ld   %le",i,Lsrs->lhs[i]);
  }
  */
  
  compute_ipstresses (0);
  
  /*
  long ipp;
  fprintf (Out,"\n\n");
  for (i=0;i<Mt->ne;i++){
    ipp=Mt->elements[i].ipp[0][0];
    fprintf (Out,"\n %le %le %le",Mm->ip[ipp].stress[0],Mm->ip[ipp].stress[1],Mm->ip[ipp].stress[2]);
    ipp++;
    fprintf (Out,"     %le %le %le",Mm->ip[ipp].stress[0],Mm->ip[ipp].stress[1],Mm->ip[ipp].stress[2]);
  }
  */

  delete [] auxv;
  delete [] v;
}


/**
   difference method for explicit solution of particle dynamic problems
   
   @param lcid - load case id
   
   JK, 6. 7. 2019
*/
void explicit_difference_method (long lcid)
{
  long i,j,n;
  double time,dt,end_time;
  double *a,*v,*p,*da,*dp,*lhs,*rhs,*fs;
  
  //  number of degrees of freedom
  n = Ndofm;  
  //  initial time
  time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  first time increment
  dt = Mp->timecon.initialtimeincr ();
  
  //  nodal values - displacements at next step
  //  (assigment of initial values)
  lhs = Lsrs->give_lhs (lcid);
  //  time derivatives of nodal values - velocities
  //  (assigment of initial values)
  v = Lsrs->give_tdlhs (lcid);
  //  second time derivatives of nodal values - accelerations
  a = Lsrs->give_stdlhs (lcid);
  //  vector of dynamic load
  rhs = Lsrs->give_rhs (2*lcid);
  //  vector of static load
  fs = Lsrs->give_rhs (2*lcid+1);
  

  //  vector of displacements at previous step
  dp = new double [n];
  //  vector of displacements at actual step
  da = new double [n];
  //  auxiliary vector
  p = new double [n];
  
  //  assembling of static loads  
  nullv (fs+lcid*n, n);
  Mb->lc[lcid].assemble (lcid,fs,NULL,1.0);
  

  i=0;
  print_init(-1, "wt");  
  print_step(lcid, i, Mp->time, rhs);
  print_flush();
  
  
  //  initial values
  for (j=0;j<n;j++){
    da[j] = lhs[j];
    dp[j] = da[j] - dt*v[j];
  }

  //  mass matrix
  if (Mp->tstormm != diag_mat){
    print_err("the mass matrix is not diagonal in explicit algorithm", __FILE__, __LINE__, __func__);
    abort ();
  }
  mass_matrix (lcid);
  

  //  damping matrix
  //damping_matrix (lcid);
  //  stiffness matrix
  stiffness_matrix (lcid);
 
  
  //  loop over time interval
  do{
    
    if (Mespr==1)  fprintf (stdout,"\n                                time  %f",Mp->time);
    
    //  assembling of dynamic loads
    nullv (rhs+lcid*n, n);
    Mb->dlc[lcid].assemble (lcid,rhs,NULL,Ndofm,Mp->time);
    
    //  total load = static + dynamic
    addv (rhs,fs,n);

    //  internal forces
    //nullv (p,n);
    //internal_forces (lcid, p);
    //cmulv (-1.0,p,n);
    //addv (rhs,p,n);

    //  K.da
    Smat->gmxv (da,p);
    cmulv (-1.0,p,n);
    addv (rhs,p,n);


    cmulv (dt*dt,rhs,n);

    //  M.da.2
    Mmat->gmxv (da,p);
    cmulv (2.0,p,n);
    addv (rhs,p,n);
    
    //  M.dp
    Mmat->gmxv (dp,p);
    cmulv (-1.0,p,n);
    addv (rhs,p,n);
    
    //  C.dp.dt/2
    //Dmat->gmxv (dp,p);
    //cmulv (0.5*dt,p,n);
    //addv (rhs,p,n);
    
    //  dynamic stiffness matrix
    //  (stored in array for mass matrix)
    //  M + C dt/2
    //Mmat->addgm (0.5*dt,*Dmat);
    
    
    
    //  solution of system of algebraic equations
    //  Gtm
    //  Mp->ssle->prec - preconditioning
    //  lhs - left hand side
    //  rhs - right hand side
    Mp->ssle->solve_system (Gtm,Mmat,lhs,rhs,Out);
    //Mmat->solve_system (Gtm,Mp->ssle->prec,lhs,rhs,Out);
    
    for (j=0;j<n;j++){
      //  velocity
      v[j]=(lhs[j]-dp[j])/2.0/dt;
      //  acceleration
      a[j]=(lhs[j]-2.0*da[j]+dp[j])/dt/dt;
      
      dp[j]=da[j];
      da[j]=lhs[j];
    }
    
    //  time increment
    Mp->time = Mp->timecon.newtime ();
    dt = Mp->timecon.actualbacktimeincr ();
    i++;

    compute_req_val (lcid);
    print_step(lcid, i, Mp->time, rhs);
    print_flush();
    

  }while(Mp->time<end_time);

  print_close();  
  delete [] da;
  delete [] dp;
  delete [] p;
}


/**
   difference method for solution of forced dynamic problems
   the method is explicit
   it serves for linear problems and benchmarks
   
   @param lcid - load case id
   
   JK, 5. 7. 2019
   modification 3. 3. 2021
*/
void difference_method_2 (long lcid)
{
  long i,j,n;
  double time,dt,end_time;
  double *a,*v,*p,*da,*dp,*lhs,*rhs,*fs;
  
  //  number of degrees of freedom
  n = Ndofm;
  //  initial time
  time=Mp->timecon.starttime ();
  //  end time
  end_time = Mp->timecon.endtime ();
  //  first time increment
  dt = Mp->timecon.initialtimeincr ();
  
  //  nodal values - displacements at next step
  //  (assignment of initial values)
  lhs = Lsrs->give_lhs (lcid);
  //  time derivatives of nodal values - velocities
  //  (assignment of initial values)
  v = Lsrs->give_tdlhs (lcid);
  //  second time derivatives of nodal values - accelerations
  a = Lsrs->give_stdlhs (lcid);
  //  vector of dynamic load
  rhs = Lsrs->give_rhs (2*lcid);
  //  vector of static load
  fs = Lsrs->give_rhs (2*lcid+1);
  

  //  vector of displacements at previous step
  dp = new double [n];
  //  vector of displacements at actual step
  da = new double [n];
  //  auxiliary vector
  p = new double [n];
  
  //  assembling of static loads  
  nullv (fs+lcid*n, n);
  Mb->lc[lcid].assemble (lcid,fs,NULL,1.0);
  

  i=0;
  print_init(-1, "wt");  
  print_step(lcid, i, Mp->time, rhs);
  print_flush();
  
  
  //  initial values
  for (j=0;j<n;j++){
    da[j] = lhs[j];
    dp[j] = da[j] - dt*v[j];
  }

  //  mass matrix
  mass_matrix (lcid);
  //  damping matrix
  //damping_matrix (lcid);
  
  
  
  //  loop over time interval
  do{
    
    if (Mespr==1)  fprintf (stdout,"\n                                time  %f",Mp->time);
    
    //  assembling of dynamic loads
    nullv (rhs+lcid*n, n);
    Mb->dlc[lcid].assemble (lcid,rhs,NULL,Ndofm,Mp->time);
    
    //  total load
    addv (rhs,fs,n);
    
    //  internal forces
    internal_forces (lcid,lhs);
    
    for (j=0;j<n;j++){
      p[j] = 2.0*da[j] - dp[j];
      rhs[j]=dt*dt*(rhs[j]-lhs[j]);
    }
    
    Mmat->gmxv (p,lhs);
    
    for (j=0;j<n;j++){
      rhs[j]+=lhs[j];
    }
    
    Mp->ssle->solve_system (Gtm,Mmat,lhs,rhs,Out);
    
    for (j=0;j<n;j++){
      v[j]=(lhs[j]-dp[j])/2.0/dt;
      a[j]=(lhs[j]-2.0*da[j]+dp[j])/dt/dt;
      
      dp[j]=da[j];
      da[j]=lhs[j];
    }
    
    //  time increment
    Mp->time = Mp->timecon.newtime ();
    dt = Mp->timecon.actualbacktimeincr ();
    i++;

    compute_req_val (lcid);
    print_step(lcid, i, Mp->time, rhs);
    print_flush();
    

  }while(Mp->time<end_time);

  print_close();  
  delete [] da;
  delete [] dp;
  delete [] p;
}

