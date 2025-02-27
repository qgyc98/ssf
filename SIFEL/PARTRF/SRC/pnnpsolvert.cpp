#include "mpi.h"
#include "backupsolt.h"
#include "pnnpsolvert.h"
#include "pglobalt.h"
#include "seqfilest.h"
#include <string.h>
#include <math.h>

/**
   JK, 12.7.2011
*/
void par_solve_nonlinear_nonstationary_problem_dform ()
{
  long i,j,k,n,ani,ini,lcid,stop,nsts;
  double dt,end_time,alpha,zero,dtmin,dtmax;
  double *d,*p,*f,*v,*z,*lhs,*lhsi,*lhsb,*tdlhs,*tdlhsb,*rhs,*gz,*grhs,*err,*thresh;
  double norfb;

  //  load case id must be equal to zero in this type of problem
  lcid=0;
  zero=1.0e-20;
  
  
  n=Ndoft;
    
  //  nodal values
  lhs = Lsrst->give_lhs (0);
  //  time derivatives of nodal values
  tdlhs = Lsrst->give_tdlhs (0);
  //  initial nodal values
  lhsi = Lsrst->lhsi;
  //  prescribed fluxes / right hand side
  rhs = Lsrst->give_rhs (0);
  
  //  vector of prescribed fluxes (right hand side)
  f = new double [n];
  nullv (f,n);
  //  predictor
  d = new double [n];
  nullv (d,n);
  //  auxiliary vector
  p = new double [n];
  nullv (p,n);
  //  auxiliary vector
  v = new double [n];
  nullv (v,n);
  //  auxiliary vector
  z = new double [n];
  nullv (z,n);
  //  backup of nodal values
  lhsb = new double [n];
  nullv (lhsb,n);
  //  backup of time derivatives of nodal values
  tdlhsb = new double [n];
  nullv (tdlhs,n);
  
  gz = new double [Psolt->schcom->ndofcp];
  grhs = new double [Psolt->schcom->ndofcp];

  //  nodes - integration points interpolation
  //  the intial values are interpolated from nodes to integration points
  approximation ();
  
  if (Tp->nvs==1 && Tp->pnvs==1){
    //  if the actual nodal values and nodal values from the previous time
    //  steps are needed, the following function has to be called
    actual_previous_nodval ();
  }
  
  //  parameter of the trapezoidal method
  alpha=Tp->alpha;
  //  prescribed norm of residual
  err=Tp->errarr;
  //  maximum number of iterations in one time step
  ini=Tp->nii;
  //  threshold for the size of the right hand side
  thresh=Tp->threshrhs;

  //  minimum time increment
  dtmin=Tp->timecont.dtmin;
  //  maximum time increment
  dtmax=Tp->timecont.dtmax;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  //  function initiates some materials models, e.g.
  //  cemhydmat, saltmat3, saltmat4
  Tm->initmaterialmodels();
  //  nodes - integration points interpolation
  approximation (); 
  
  //  toto vyhodit
  // nulleqother ();
  if (Tp->nvs==1 && Tp->pnvs==1)
    actual_previous_nodval ();
  // -----*******
  //  az sem
  
  double *uf,*buff;
  buff = new double [10];
  uf = new double [10];
  
  
  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  //  number of time step
  i=0;
  //  number of successful time steps (number of steps without inner iterations)
  nsts=0;

  stop=0;

  
  //  inital printing of output and graphical informations
  if (Tp->hdbcont.restore_stat()){
    solvert_restore (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    print_initt(-1, "at",Ptp->fni,Ptp->fei);
    print_stept(0,i,Tp->time,NULL);
  }
  else{
    compute_req_valt (0);
    print_initt(-1, "wt",Ptp->fni,Ptp->fei);
    print_stept(0,i,Tp->time,rhs);
  }
  print_flusht();
  
  //  toto vyhodit
  Tm->initmaterialmodels();

  do{
    
    //  time update
    Tp->time=Tp->timecont.newtime (dt);
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
    
    Tm->updateipval ();
    i++;

    //  backup of attained nodal values and their time derivatives
    for (j=0;j<n;j++){
      lhsb[j]=lhs[j];
      tdlhsb[j]=tdlhs[j];
    }

    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    //  predictor d
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }

    fprintf (Outt,"\n\n\n\n kontrola vektoru lhs a tdlhs a d\n");
    for (long ijk=0;ijk<Gtt->nn;ijk++){
      long cn1,cn2;
      fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
      cn1=Gtt->gnodes[ijk].cn[0];
      cn2=Gtt->gnodes[ijk].cn[1];
      if (cn1>0 && cn2>0){
	fprintf (Outt,"  lhs % 12.9le    tdlhs % 12.9le    d % 12.9le",lhs[cn1-1],tdlhs[cn1-1],d[cn1-1]);
	fprintf (Outt,"  lhs % 12.9le    tdlhs % 12.9le    d % 12.9le",lhs[cn2-1],tdlhs[cn2-1],d[cn2-1]);
      }
    }

    
    //  C.d
    Cmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);

    //  right hand side - prescribed fluxes f_{n+1}
    trfel_right_hand_side (lcid,f,n);
    
    //  C.d + alpha dt f_{n+1}
    for (j=0;j<n;j++){
      rhs[j] = f[j]*alpha*dt + p[j];
    }


    fprintf (Outt,"\n\n\n\n kontrola vektoru p a f a rhs\n");
    for (long ijk=0;ijk<Gtt->nn;ijk++){
      long cn1,cn2;
      fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
      cn1=Gtt->gnodes[ijk].cn[0];
      cn2=Gtt->gnodes[ijk].cn[1];
      if (cn1>0 && cn2>0){
	fprintf (Outt,"  p % 12.9le    f % 12.9le    rhs % 12.9le",p[cn1-1],f[cn1-1],rhs[cn1-1]);
	fprintf (Outt,"  p % 12.9le    f % 12.9le    rhs % 12.9le",p[cn2-1],f[cn2-1],rhs[cn2-1]);
      }
    }
     
    
    nullv (z,n);
    
    Kmat->diag_check (zero,rhs);
    
    /*
    //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
    //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
    //fprintf (stdout,"\n\n\n\n\n kontrola reseni");
    fprintf (Outt,"\n\n kontrola vektoru pred resenim\n");
    for (long ijk=0;ijk<Gtt->nn;ijk++){
      long cn1,cn2;
      cn1=Gtt->gnodes[ijk].cn[0];
      cn2=Gtt->gnodes[ijk].cn[1];
      //fprintf (stdout,"\n node %6ld   %10.6lf %10.6lf    %6ld %6ld",ijk,Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y,cn1,cn2);
      //fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
      if (cn1>0 && cn2>0){
	fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
	fprintf (Outt,"\n %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",f[cn1-1],p[cn1-1],v[cn1-1],z[cn1-1],rhs[cn1-1]);
	fprintf (Outt," %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",f[cn2-1],p[cn2-1],v[cn2-1],z[cn2-1],rhs[cn2-1]);
	//fprintf (stdout,"\n    %12.9le %12.9le",lhs[cn1-1],lhs[cn2-1]);
	//}
	//if (cn1>0){
	//fprintf (Outt,"\n %6ld  f,p,v,z,rhs  %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",ijk,f[cn1-1],p[cn1-1],v[cn1-1],z[cn1-1],rhs[cn1-1]);
	//}
      }
      else{
	//fprintf (Outt,"    %12.9le %12.9le",0.0,0.0);
	//fprintf (stdout,"\n    %12.9le %12.9le",0.0,0.0);
      }
    }
    */
    
    /*
    if (Kmat->dsky){
      fprintf (Outt,"\n\n\n\n Kontrola pole dsky");
      for (i=0;i<Kmat->dsky->n;i++){
	for (j=Kmat->dsky->adr[i];j<Kmat->dsky->adr[i+1];j++){
	  fprintf (Outt,"\n %13.10le",Kmat->dsky->a[j]);
	}
      }
    }
    if (Kmat->sky){
      fprintf (Outt,"\n\n\n\n Kontrola pole sky");
      for (i=0;i<Kmat->sky->n;i++){
	for (j=Kmat->sky->adr[i];j<Kmat->sky->adr[i+1];j++){
	  fprintf (Outt,"\n %13.10le",Kmat->sky->a[j]);
	}
      }
    }
    */
    
    
    fprintf (Outt,"\n\n\n\n druha kontrola vektoru rhs\n");
    for (long ijk=0;ijk<Gtt->nn;ijk++){
      long cn1,cn2;
      fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
      cn1=Gtt->gnodes[ijk].cn[0];
      cn2=Gtt->gnodes[ijk].cn[1];
      if (cn1>0 && cn2>0){
	fprintf (Outt,"  rhs % 12.9le   % 12.9le",rhs[cn1-1],rhs[cn2-1]);
      }
    }
  
   
    //Kmat->solve_system (tdlhs,fb);
    Psolt->par_linear_solver (Gtt,Kmat,lhs,rhs,Outt,Mesprt);
    
    //fprintf (stdout,"\n\n\n\n\n kontrola reseni");
    fprintf (Outt,"\n\n\n\n\n kontrola reseni");
    //fprintf (Outt,"\n\n\n\n\n ");
    for (long ijk=0;ijk<Gtt->nn;ijk++){
      long cn1,cn2;
      cn1=Gtt->gnodes[ijk].cn[0];
      cn2=Gtt->gnodes[ijk].cn[1];
      //fprintf (stdout,"\n node %6ld   %10.6lf %10.6lf    %6ld %6ld",ijk,Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y,cn1,cn2);
      fprintf (Outt,"\n %10.6lf %10.6lf  ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
      //if (cn1>0 && cn2>0){
      //fprintf (Outt,"    %12.9le %12.9le",lhs[cn1-1],lhs[cn2-1]);
	//fprintf (stdout,"\n    %12.9le %12.9le",lhs[cn1-1],lhs[cn2-1]);
      //}
      if (cn1>0){
	fprintf (Outt,"    % 12.9le",lhs[cn1-1]);
      }
      else{
	fprintf (Outt,"    % 12.9le",0.0);
      }
      if (cn2>0){
	fprintf (Outt,"    % 12.9le",lhs[cn2-1]);
      }
      else{
	fprintf (Outt,"    % 12.9le",0.0);
      }
     }
    
    //fprintf (Outt,"\n\n\n\n\n kontrola leve strany");
    //fprintf (stdout,"\n\n\n\n\n time %le",Tp->time);
    //fprintf (stdout,"\n\n\n\n\n kontrola leve strany - reseni time %le",Tp->time);
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (stdout,"\n lhs %6ld  %le",ijk,lhs[ijk]);
      //fprintf (Outt,"\n %ld %le",ijk,lhs[ijk]);
    //}

    //  computation of time derivatives of the nodal values
    for (j=0;j<n;j++){
      tdlhs[j]=(lhs[j]-d[j])/dt/alpha;
    }
    
    solution_correction ();
    //  approximation of nodal values into ontegration points
    approximation ();
    
    if (Tp->nvs==1 && Tp->pnvs==1)
      actual_previous_nodval ();

    
    //  iteration for equilibrium
    for (j=0;j<ini;j++){
      
      //  capacity matrix
      capacity_matrix (lcid);
      
      //  conductivity matrix
      conductivity_matrix (lcid);
      
      
      //  matrix of the system of equations
      //  C + alpha.dt.K
      Kmat->scalgm (dt*alpha);
      Kmat->addgm (1.0,*Cmat);
      
      //  (C + alpha.dt.K).d
      nullv (v,n);
      Kmat->gmxv (lhs,v);
      
      //  C.d
      nullv (p,n);
      Cmat->gmxv (d,p);
      
      //  computation of residuals
      //norrhs=0.0;
      //norres=0.0;
      
      //fprintf (stdout,"\n\n alpha  %le   dt %le\n\n",alpha,dt);
      fprintf (Outt,"\n\n kontrola vektoru \n");
      for (k=0;k<n;k++){
	rhs[k]=f[k]*alpha*dt+p[k]-v[k];
	

	//fprintf (stdout,"\n rhs %6ld   %le  f  %15.12le   v  %15.12le  p  %15.12le  z %15.12le",k,rhs[k],f[k],v[k],p[k],z[k]);
	//	fprintf (stdout,"\n rhs   %15.10le   k je %d",rhs[k],k);
	//norres+=rhs[k]*rhs[k];
	z[k]+=f[k]*alpha*dt+p[k];

	fprintf (Outt,"\n %6ld f,p,v,z,rhs   %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",k,f[k],p[k],v[k],z[k],rhs[k]);

	//norrhs+=z[k]*z[k];
	//fprintf (stdout,"\n rhs %6ld   %le  f  %15.12le   v  %15.12le  p  %15.12le  z %15.12le",k,rhs[k],f[k],v[k],p[k],z[k]);
      }
      //norres=sqrt(norres);
      //norrhs=sqrt(norrhs);
      
      /*
      //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
      //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
      //fprintf (stdout,"\n\n\n\n\n kontrola reseni");
      fprintf (Outt,"\n\n kontrola vektoru po reseni\n");
      for (long ijk=0;ijk<Gtt->nn;ijk++){
	long cn1,cn2;
	cn1=Gtt->gnodes[ijk].cn[0];
	cn2=Gtt->gnodes[ijk].cn[1];
	//fprintf (stdout,"\n node %6ld   %10.6lf %10.6lf    %6ld %6ld",ijk,Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y,cn1,cn2);
	//fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
	if (cn1>0 && cn2>0){
	  fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
	  fprintf (Outt,"\n %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",f[cn1-1],p[cn1-1],v[cn1-1],z[cn1-1],rhs[cn1-1]);
	  fprintf (Outt," %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",f[cn2-1],p[cn2-1],v[cn2-1],z[cn2-1],rhs[cn2-1]);
	  //fprintf (stdout,"\n    %12.9le %12.9le",lhs[cn1-1],lhs[cn2-1]);
	  //}
	  //if (cn1>0){
	  //fprintf (Outt,"\n %6ld  f,p,v,z,rhs  %12.9le  %12.9le  %12.9le  %12.9le   %12.9le",ijk,f[cn1-1],p[cn1-1],v[cn1-1],z[cn1-1],rhs[cn1-1]);
	  //}
	}
	else{
	  //fprintf (Outt,"    %12.9le %12.9le",0.0,0.0);
	  //fprintf (stdout,"\n    %12.9le %12.9le",0.0,0.0);
	}
      }
      */


      //fprintf (stdout,"\n KONTROLA norres   %15.10le",norres);
      //fprintf (stdout,"\n\n\n\n\n kontrola reseni");
      fprintf (Outt,"\n\n kontrola vektoru rhs a z\n");
      for (long ijk=0;ijk<Gtt->nn;ijk++){
	long cn1,cn2;
	cn1=Gtt->gnodes[ijk].cn[0];
	cn2=Gtt->gnodes[ijk].cn[1];
	fprintf (Outt,"\n %10.6lf %10.6lf   ",Gtt->gnodes[ijk].x,Gtt->gnodes[ijk].y);
	if (cn1>0 && cn2>0){
	  fprintf (Outt,"  rhs % 12.9le   z % 12.9le   rhs % 12.9le  z % 12.9le",rhs[cn1-1],z[cn1-1],rhs[cn2-1],z[cn2-1]);
	}
      }



      
      stop=0;
      for (k=0;k<Tp->ntm;k++){
	if (Myrank==0){
	  //fprintf (stdout,"\n Psolt->schcom->ndofcp %ld",Psolt->schcom->ndofcp);
	  nullv (grhs,Psolt->schcom->ndofcp);
	  nullv (gz,Psolt->schcom->ndofcp);
	}
	//  if the j-th medium is in equilibrium state, the following function returns 1
	stop+= Psolt->selected_norm_calculation (k,err[k],thresh[k],1,Gtt,rhs,grhs,z,gz,n,Outt,norfb);
	fprintf (stdout,"\n i,j,k   %ld %ld %ld    %ld",i,j,k,stop);
      }
      
      if (stop==Tp->ntm){
	break;
      }
      
      if (Mesprt != 0)  fprintf (stdout,"\n iteration number %ld",j);
      
      Kmat->diag_check (zero,rhs);
      
      
      Psolt->par_linear_solver (Gtt,Kmat,z,rhs,Outt,Mesprt);
      
      for (k=0;k<n;k++){
	lhs[k]+=z[k];
      }
      
      nullv (z,n);
      
      //  computation of time derivatives of the nodal values
      for (k=0;k<n;k++){
	tdlhs[k]=(lhs[k]-d[k])/dt/alpha;
      }
      
      //  physically corrected solution
      solution_correction ();    
      //  approximation of nodal values into ontegration points
      approximation ();
      
      
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	actual_previous_nodval ();
      // ***********************************
      
    }//  end of inner iteration loop


    //  actual number of performed iterations
    ani=j;
    
    if (ani==0)
      nsts++;
    
    if (ani==ini){
      //  backup is used for actual values
      for (j=0;j<n;j++){
	lhs[j]=lhsb[j];
	tdlhs[j]=tdlhsb[j];
      }
      
      // *************************
      // vlozil JM 25.4.2008
      
      //  physically corrected solution
      solution_correction ();    
      //  approximation of nodal values into ontegration points
      approximation ();
      
      
      
      // nulleqother ();
      if (Tp->nvs==1 && Tp->pnvs==1)
	actual_previous_nodval ();
      // ***********************************
      
      
      //  reduction of the time increment because
      //  inner loop was not able to enforce equilibrium
      dt/=2.0;
      i--;
      Tp->timecont.oldtime ();
      
      //fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      // fprintf (stdout,"\n time increment is reduced because the inner loop was not able to enforce equilibrium dt = %le",dt);
      //  fprintf (Outt2,"\n time increment is reduced because the inner loop was not able to enforce equilibrium");
      if (dt<dtmin){
	fprintf (stderr,"\n\n time increment is less than minimum time increment");
	fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
	abort ();
      }
    }else{
      //  time increment
      //Tp->time = Tp->timecont.newtime ();
      //dtdef = Tp->timecont.actualforwtimeincr ();
      
      if (nsts==2){
	double times, acttimes;
	times =  Tp->timecont.starttime();
	acttimes = Tp->timecont.actualtime();
	
	if ((times+18000)>acttimes) {
	  dt = dt;
	}
	else{
	  dt*=2;
	  if (dt>dtmax) {
	    dt = dtmax;
	  }
	}
	//dt*=1.2;
	nsts=0;
	//	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	//fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le, dtdef = %le",dt,dtdef);
	//	fprintf (stdout,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt = %le",dt);
	//fprintf (Outt2,"\n time increment is enlarged because no inner loop was neccessary in previous 3 steps dt= %lf",dt);
      }
      
      //if (dt>dtdef)
      //dt=dtdef;
      
      //Tp->time = Tp->timecont.newtime (dt);
      
      
      //time+=dt;
      //Tp->time=time;
      //Tp->timecont.time=time;
      //Tp->timecont.fdt=dt;
      
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      
    }
    
    if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
      solvert_save (lhs,tdlhs,f,i,Tp->time,dt,Tp->timecont,n);
    
  }while(Tp->time<end_time);

  
  /*
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] d;
  */
  //print_closet();
}



void par_solve_nonlinear_nonstationary_problem ()
{
  long i,j,k,n,ini;
  double dt,end_time,alpha,err,totnorfb;
  double *d,*fb,*fi,*p,*lhs,*lhsi,*tdlhs,*rhs,*arhs;

  n=Ndoft;
    
  lhs = Lsrst->give_lhs (0);
  tdlhs = Lsrst->give_tdlhs (0);
  lhsi = Lsrst->lhsi;
  rhs = Lsrst->give_rhs (0);
  
  d = new double [n];
  arhs = new double [n];
  p = new double [n];
  fb = new double [n];
  fi = new double [n];
  
  //  initial values
  nullv (lhs,n);
  nullv (tdlhs,n);
  nullv (rhs,n);
  nullv (arhs,n);
  nullv (d,n);
  nullv (p,n);
  nullv (fb,n);
  nullv (fi,n);
  
  //  nodes - integration points interpolation
  approximation ();
  
  alpha=Tp->alpha;
  err=Tp->err;
  ini=Tp->nii;
  
  //  initial time
  Tp->time=Tp->timecont.starttime ();
  //  end time
  end_time = Tp->timecont.endtime ();
  //  initial time increment
  dt = Tp->timecont.initialtimeincr ();
  
  i=0;
  print_initt(-1, "wt",Ptp->fni,Ptp->fei);
  print_stept(0,i,Tp->time,NULL);
  print_flusht();

  
  // ***************************
  //  main iteration loop  ****
  // ***************************
  do{
    
    if (Mesprt != 0)  fprintf (stdout,"\n\n iteration number %ld, time %f, time step = %f\n",i,Tp->time,dt);
          
    //kontrolni tisk ???!!!
    /* if ((Tp->time >= 304.0))
       {
       fprintf (Outt,"\n\nkontrolni tisk time = %e\n",Tp->time);
       fprintf (Outt,"------------------------\n\n",Tp->time);  
       }
    */
    //  capacity matrix
    capacity_matrix (0);
    
    //  conductivity matrix
    conductivity_matrix (0);
    
    for (j=0;j<n;j++){
      d[j] = lhs[j] + (1.0-alpha)*dt*tdlhs[j];
    }
    
    //  right hand side
    trfel_right_hand_side (0,rhs,n);
    
    //fprintf (Outt,"\n\n\n\n\n kontrola prave strany");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n %ld %le",ijk,rhs[ijk]);
    //}
    
    
    
    //  auxiliary vector  K (d+(1-alpha)*dt*v)
    Kmat->gmxv (d,p);
    
    //  matrix of the system of equations
    //  C + alpha.dt.K
    Kmat->scalgm (dt*alpha);
    Kmat->addgm (1.0,*Cmat);
    
    //  F - K (d+(1-alpha)*dt*v)
    for (j=0;j<n;j++){
      fb[j] = rhs[j] - p[j];
    }
    
    //Kmat->solve_system (tdlhs,fb);
    Psolt->par_linear_solver (Gtt,Kmat,tdlhs,fb,Outt,Mesprt);
    
    //  d_{i+1}=\tilda{d_{i+1}} + alpha.dt.v_{i+1}
    for (j=0;j<n;j++){
      lhs[j] = d[j] + alpha*dt*tdlhs[j];
    }
    
    //fprintf (Outt,"\n\n\n\n\n kontrola leve strany");
    //for (long ijk=0;ijk<n;ijk++){
    //fprintf (Outt,"\n %ld %le",ijk,lhs[ijk]);
    //}

    approximation ();
    internal_fluxes (fi,n);
    
    //  vector of unbalanced fluxes
    for (j=0;j<n;j++){
      fb[j]=fi[j];
      //fprintf (Outt,"\nfb %ld  %e",j,fb[j]);
    }

    
    totnorfb = Psolt->pss (fb, fb, Outt);
    
    
    

    if (Mesprt==1)  fprintf (stdout,"\n%e %e",totnorfb,err);
    //if (Mesprt==1)  fprintf (Outt,"\n\ntotnorfb %e   err %e\n",totnorfb,err);
    
    if (totnorfb<err){
      //  time increment
      Tp->time = Tp->timecont.newtime ();
      dt = Tp->timecont.actualbacktimeincr ();
      i++;

      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      continue;
    }
    

    for (j=0;j<ini;j++){
      

      //fprintf (Outt,"\n\n\n\n\n vnitrni smycka kontrola prave strany");
      //for (long ijk=0;ijk<n;ijk++){
      //fprintf (Outt,"\n %ld %le",ijk,fb[ijk]);
      //}

      Psolt->par_linear_solver (Gtt,Kmat,p,fb,Outt,Mesprt);
      
      for (k=0;k<n;k++){
        tdlhs[k]+=p[k];
        lhs[k]+=alpha*dt*p[k];
      }
      
      //fprintf (Outt,"\n\n\n\n\n vnitrni smycka kontrola leve strany");
      //for (long ijk=0;ijk<n;ijk++){
      //fprintf (Outt,"\n %ld %le",ijk,lhs[ijk]);
      //}

      approximation ();
      internal_fluxes (fi,n);
      
      //  vector of unbalanced fluxes
      for (k=0;k<n;k++){
        fb[k]= fi[k];
      }      
      
      totnorfb = Psolt->pss (fb, fb, Outt);
      
      if (Mesprt==1)  fprintf (stdout,"\n iteration %ld   error %e",j,totnorfb);
      
      if (totnorfb<err){
        break;
      }
    }


    //  time increment
    Tp->time = Tp->timecont.newtime ();
    dt = Tp->timecont.actualbacktimeincr ();
    i++;
    
    print_stept(0,i,Tp->time,NULL);
    print_flusht();

    printf ("\n Tp->time   %e",Tp->time);
    printf ("\n i          %ld",i);
    
  }while(Tp->time<end_time);
  
  delete [] fi;  delete [] fb;  delete [] p;
  delete [] arhs;  delete [] d;
  
  //print_closet();
}

