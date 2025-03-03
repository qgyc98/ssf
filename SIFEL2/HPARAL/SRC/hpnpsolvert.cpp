#include "mpi.h"
#include "hpnpsolvert.h"
#include "hpglobal.h"
#include "seqfilest.h"
#include <string.h>
#include <math.h>
//#include "MPI_Send_mine.h"
//#include "MPI_Recv_mine.h"

/**
   function computes linear nonstationary transport problem
   material parameters are obtained by homogenization
   on microscale, stationary problems are solved
   these microproblems can be defined on every finite element or
   on aggregates of finite elements, the aggregates can send
   data from all aggregated elements or averaged data which
   has the same form as data from a single element
   
   the function parallel_homogenization works and was tested against
   MATLAB code, the function parallel_homogenization_new27_5_2011 contains
   compact subroutines,
   
   if the function parallel_homogenization_new works, the function
   parallel_homogenization will be removed

   JK, 27.5.2011
   TKr 20/07/2015
*/
void parallel_homogenization_lin_nonstat ()
{
  long i,j,k,ii,jj,kk,n,nelidom,buffsize1,buffsize2,ntm,dim;
  long *neldom,*buff;
  double zero,alpha,*d,*p,*f,*lhs,*tdlhs,*rhs;
  double dt,time,end_time,*buff1,*buff2;
  MPI_Status stat;
  vector gr(2);
  
  //  allocation of buffer
  //  buff[0] - the number of connections between the macroproblem and microproblem
  //  buff[1] - the maximum number of connections between the macroproblem and microproblem
  buff = new long [2];
  

  /**************************************************************************/

  if (Myrank==0){  
    switch (Tp->mami){
    case elem_domain:{
      //  each element is connected with its own microproblem
      buff[0]=1;
      buff[1]=1;
            
      break;
    }
    case aggregate_domain:{
      //  generation of auxiliary data needed in parallel computation
      //  elements of macroproblem are aggregated, each finite element
      //  sends its data to microproblem and obtains its matrices
      
      //  map between finite elements of macroproblem and microproblems
      //  this map is read in the subroutine transtop::read_seq_top
      //  eldom[i]=j - the i-th finite element of macroproblem is connected to the j-th microproblem
      
      //  the number of elements on subdomains
      //  it is the number of macroelements connected to one microproblem
      //  the number of subdomains--microproblems must be Nproc-1
      //  neldom[i]=j - j finite elements of the macroproblem are connected to the i-th microproblem
      //  neldom[0]=0 - the master processor (myrank=0) contains no microproblem
      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }
      for (i=0;i<Tt->ne;i++){
	neldom[Gtt->stop->eldom[i]]++;
      }
      buffsize1=0;
      for (i=0;i<Nproc;i++){
	if (buffsize1<neldom[i])
	  buffsize1=neldom[i];
      }
      
      buff[1]=buffsize1;
      
      break;
    }
    case material_aggregate_domain:{
      //one domain represents one material = one processor
      //the aggregate covers elements of one materyal type
      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }
      
      for (i=0;i<Tt->ne;i++)
	neldom[Gtt->stop->eldom[i]] = 1;
      
      buff[0]=1;
      buff[1]=1;
      
      break;
    }
    case eff_aggregate_domain:{

      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }

      for (i=0;i<Tt->ne;i++){
	if(Gtt->stop->eldom[i] != Tt->elements[i].idm[0])
	  neldom[Gtt->stop->eldom[i]]++;
      }
      
      //  the aggregate sends averaged data to microproblem
      //  only aggregate sends and receives data, the same
      //  matrices are used on all elements in the aggregate
      buff[0]=1;
      buff[1]=1;
      
      break;
    }
    default:{
      par_print_err(Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
    }
    }    
  
    for (i=1;i<Nproc;i++){
      if (elem_domain == aggregate_domain)
	buff[0]=neldom[i];

      MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      //MPI_Send_mine (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
    if (elem_domain == aggregate_domain)
      buff[0]=neldom[0];
  }
  else{
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    //MPI_Recv_mine (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  

  //  the number of connections between the macroproblem and microproblem
  //  the number of entities (elements or subdomains) of macroproblem connected with microproblem
  nelidom=buff[0];
  //  the maximum number of connections between the macroproblem and microproblem
  //  the maximum number of entities (elements or subdomains) of macroproblem connected with microproblem
  buffsize1=buff[1];
  
  //  dimension of the problem, one, two or three dimensional case
  //  the dimension is defined with respect to the first finite element
  dim = Tt->give_dimension (0);
  //  the number of transported materials
  ntm = Tp->ntm;
  
  buffsize2=buffsize1;
  //  buffer for values and gradient components
  buffsize1*=ntm*(dim+1);
  //  buffer for the conductivity and capacity matrices
  //  diagonal capacity matrices are assumed
  buffsize2*=((ntm*dim)*(ntm*dim)+ntm*ntm);
  
  fprintf (stdout,"\n procesor %d  nelidom   %ld",Myrank,nelidom);
  fprintf (stdout,"\n procesor %d  buffsize1 %ld",Myrank,buffsize1);
  fprintf (stdout,"\n procesor %d  buffsize2 %ld",Myrank,buffsize2);
  
  
  delete [] buff;
  buff1 = new double [buffsize1];
  nullv (buff1,buffsize1);
  buff2 = new double [buffsize2];
  nullv (buff2,buffsize2);
  
  
  /**************************************************************************/


  if (Myrank==0){
    //  the number of unknonws
    //  it is the number of unknowns in macroproblem in the case of the master processor
    //  otherwise, it is the number of unknowns in microproblems
    n=Ndoft;
    
    //  array for nodal unknowns
    lhs = Lsrst->give_lhs (0);
    nullv (lhs,n);
    //  array for time derivatives of the nodal unknowns
    tdlhs = Lsrst->give_tdlhs (0);
    nullv (tdlhs,n);
    //  array for the right hand side
    rhs = Lsrst->give_rhs (0);
    nullv (rhs,n);
    
    //  auxiliary arrays
    d = new double [n];
    nullv (d,n);
    p = new double [n];
    nullv (p,n);
    //  vector of prescribed fluxes (right hand side)
    f = new double [n];
    nullv (f,n);
        
    //  initial values
    nullv (lhs,n);
    nullv (tdlhs,n);
    nullv (d,n);
    nullv (p,n);
    
    //  coefficient of the trapezoidal rule  
    alpha=Tp->alpha;
    //  computer zero
    zero=Tp->zero;
    
    
    // **************************
    //  main iteration loop  ***
    // **************************
    
    //  loop counter
    i=0;
    
    //  starting time
    Tp->time = Tp->timecont.starttime ();
    //  time increment
    dt = Tp->timecont.initialtimeincr ();
    //  end time
    end_time = Tp->timecont.endtime ();
    
    //  nodes - integration points interpolation
    approximation ();
    actual_previous_change ();
    
    //  initialization of material models (if needed)
    Tm->initmaterialmodels();
    
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }
  /**************************************************************************/

  if (Myrank==0){
    for (ii=1;ii<Nproc;ii++){
      MPI_Send (&end_time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      //MPI_Send_mine (&end_time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    //MPI_Recv_mine (&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  fprintf (stdout,"\n procesor %d end_time %lf\n",Myrank,end_time);
  
  
  do{
    
    /**************************************************************************/

    if (Myrank==0){
      if (Mesprt==1){
	fprintf (stdout,"\n\n --------------------------------------------------------------------------");
	fprintf (stdout,"\n HPTRFEL Time step = %ld,  Time %e,  Time increment = %e",i,Tp->time,dt);
	fprintf (stdout,"\n --------------------------------------------------------------------------\n");
      }
      
      //  determination of new time instance
      Tp->time=Tp->timecont.newtime ();
      time=Tp->time;
      //  computation of backward time step
      dt = Tp->timecont.actualbacktimeincr ();
      
      for (ii=1;ii<Nproc;ii++){
	MPI_Send (&time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
	//MPI_Send_mine (&time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (&time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //MPI_Recv_mine (&time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    /**************************************************************************/
    
    if (Myrank==0){
      //  update of step number
      i++;
      
      // **********************************
      //  solution of microscale problems
      // **********************************
      for (ii=1;ii<Nproc;ii++){
	
	//  function assembles data from the macroproblem for the ii-th microproblem to the array buff1
	// tady pole buff1 musi mit rozmer podle nelidom??!! nebo podle poctu elementu??!!
	higher_to_lower_level (ii,buffsize1,buff1);
	
	MPI_Send (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
	//MPI_Send_mine (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      //MPI_Recv_mine(buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      
      for (ii=0;ii<nelidom;ii++){
	
	//fprintf (Outt,"\n\n *********************************************");
	//fprintf (Outt,"\n Myrank  = %ld",Myrank);
	//fprintf (Outt,"\n nelidom = %ld",ii);
	//fprintf (Outt,"\n\n pred vypoctem paral_transport_homogenization:");
	//for (long ijk=0;ijk<buffsize1;ijk++){
	////for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	//}
	//fprintf (Outt,"\n");
	//for (long ijk=0;ijk<buffsize2;ijk++){
	////for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	//}
	//fflush(Outt);
	
	paral_transport_homogenization (buff1+ii*ntm*(dim+1),buff2+ii*((ntm*dim)*(ntm*dim)+ntm*ntm),
					buff2+ii*((ntm*dim)*(ntm*dim)+ntm*ntm)+(ntm*dim)*(ntm*dim));
	
	//fprintf (Outt,"\n\n po vypoctu paral_transport_homogenization:");
	//for (long ijk=0;ijk<buffsize1;ijk++){
	////for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	//}
	//fprintf (Outt,"\n");
	//for (long ijk=0;ijk<buffsize2;ijk++){
	////for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	//}
	//fflush(Outt);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    /**************************************************************************/

    if (Myrank==0){
      for (ii=1;ii<Nproc;ii++){
	MPI_Recv (buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	//MPI_Recv_mine(buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
      	k=stat.MPI_TAG;
	jj=0;
	
	for (j=0;j<Tt->ne;j++){
	  
	  //fprintf (Outt,"\n Myrank  = %ld",Myrank);
	  //fprintf (Outt,"\n element no. = %ld",j);
	  //fprintf (Outt,"\n Gtt->stop->eldom[j] = %ld",Gtt->stop->eldom[j]);
	  //fprintf (Outt,"\n ii = %ld",ii);
	  //fprintf (Outt,"\n k  =stat.MPI_TAG = %ld",k);
	  
	  if (Gtt->stop->eldom[j]==k){
	    //  only elements connected with the k-th microscale are assumed

	    switch (Tp->mami){
	    case elem_domain:{
	      //  each element is connected with its own microproblem
	      Tm->hommat[j].assemble_matrices (buff2,ntm,dim);
	      break;
	    }
	      
	    case aggregate_domain:{
	      //  generation of auxiliary data needed in parallel computation
	      Tm->hommat[j].assemble_matrices (buff2+jj*((ntm*dim)*(ntm*dim)+ntm*ntm),ntm,dim);
	      jj++;
	      break;
	    }
	    case material_aggregate_domain:{
	      //  aggregate sends averaged data to microproblem
	      //  index jj is not incremented in this case because all elements in aggregate??!!
	      //  obtain the same conductivity and capacity matrices
	      //
	      //  index of subdomain-aggregate = index of material (processor)

	      //fprintf (Outt,"\n Myrank  = %ld",Myrank);
	      //fprintf (Outt,"\n from Myrank  = %ld",ii);
	      //fprintf (Outt,"\n element  = %ld",j);
	      //fprintf (Outt,"\n\n pred ulozenim matic na materialu c. %ld:",j);
	      //for (long ijk=0;ijk<buffsize2;ijk++){
	      //fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	      //}
	      
	      Tm->hommat[j].assemble_matrices (buff2,ntm,dim);

	      //fflush(Outt);
	      break;
	    }
	    case eff_aggregate_domain:{
	      //  aggregate sends averaged data to microproblem
	      //  index jj is not incremented in this case because all elements in aggregate
	      //  obtain the same conductivity and capacity matrices
	      //
	      //  index of subdomain-aggregate = index of element material
	      kk = Tt->elements[i].idm[0];
	      Tm->hommat[kk].assemble_matrices (buff2+kk*((ntm*dim)*(ntm*dim)+ntm*ntm),ntm,dim);
	      break;
	    }
	      
	    default:{
	      par_print_err(Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
	    }
	    }
	    
	  }
	}//  end loop over the number of elements
      }//  end loop over the number of processors
    }//  end of the master processor part
    else{
      MPI_Send (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
      //MPI_Send_mine (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    fflush (stdout);


    /**************************************************************************/

    if (Myrank==0){
      
      // *********************************
      //  solution of macroscale problem
      // *********************************

      //  assembling of conductivity matrix
      conductivity_matrix (0);

      //Kmat->printmat (Outt);

      //  assembling of capacity matrix
      capacity_matrix (0);

      //Cmat->printmat (Outt);

      //  predictor
      //  dd = d_n + (1-alpha) dt v_n
      for (j=0;j<n;j++){
	d[j]=lhs[j]+(1.0-alpha)*dt*tdlhs[j];
      }
      
      //  auxiliary vector
      //  K dd
      Kmat->gmxv (d,p);
      
      //  matrix of the system of equations
      //  C + K alpha dt
      Kmat->scalgm (alpha*dt);
      Kmat->addgm (1.0,*Cmat);

      //  right hand side
      //  prescribed nodal fluxes f_{n+1}
      trfel_right_hand_side (0,rhs,n);
      
      //  computation of the right hand side vector
      //  f_{n+1} - K dd
      for (j=0;j<n;j++){
	f[j] = rhs[j] - p[j];
      }
      
      //  solution of the system of algebraic equations
      //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
      //  time derivatives v_{n+1} are solved
      Tp->ssle->solve_system (Gtt,Kmat,tdlhs,f,Outt);

      //  nodal values computed from nodal derivatives
      //  d_{n+1} = dd + alpha dt v_{n+1}
      for (j=0;j<n;j++){
	lhs[j] = d[j] + alpha*dt*tdlhs[j];
      }
      
      solution_correction ();
      //  nodes - integration points interpolation
      approximation ();
      compute_req_valt (0);
      
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      
      fprintf (stdout,"\n Tp->time   %e",Tp->time);
      fprintf (stdout,"\n i          %ld",i);
    }
    fprintf (stdout,"\n procesor %d time = %lf end_time = %lf",Myrank,time,end_time);
    fflush (stdout);
     
  }while(time<end_time);

  /**************************************************************************/

  if (Myrank==0){
    delete [] p;
    delete [] d;
    delete [] f;
    
    print_closet();
  }
  
  delete [] buff1;
  delete [] buff2;

}


/**
   function computes nonstationary transport problem
   material parameters are obtained by homogenization
   on microscale, stationary problems are solved
   these microproblems can be defined in every integration point,
   on aggregates of integration points or on aggregates of finite elements
   
   the function parallel_homogenization works and was tested against
   MATLAB code, the function parallel_homogenization_new contains
   compact subroutines,
   
   if the function parallel_homogenization_new works, the function
   parallel_homogenization will be removed

   JK, 22.9.2010
*/
void parallel_homogenization_new ()
{
  long i,j,k,ii,jj,n,nelidom,buffsize1,buffsize2,ntm,dim;
  long *neldom,*buff;
  double zero,alpha,*d,*p,*f,*lhs,*tdlhs,*rhs;
  double dt,time,end_time,*buff1,*buff2;
  MPI_Status stat;
  vector gr(2);
  
  //  allocation of buffer
  //  buff[0] - the number of connections between the macroproblem and microproblem
  //  buff[1] - the maximum number of connections between the macroproblem and microproblem
  buff = new long [2];

  if (Myrank==0){  
    
    //  this part is used only for Charles bridge
    //Gtt->stop->eldom = new long [Tt->ne];
    //for (i=0;i<Tt->ne;i++){
    //Gtt->stop->eldom[i] = Tt->elements[i].idm[0];
    //}
    //  end of the part for Charles bridge

    switch (Tp->mami){
    case elem_domain:{
      //  each element is connected with its own microproblem
      buff[0]=1;
      buff[1]=1;
      
      for (i=1;i<Nproc;i++){
	MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      }
      
      break;
    }
      
    case aggregate_domain:{
      
      //  generation of auxiliary data needed in parallel computation
      
      //  map between finite elements of macroproblem and microproblems
      //  this map is read in the subroutine transtop::read_seq_top
      //  eldom[i]=j - the i-th finite element of macroproblem is connected to the j-th microproblem
      
      //  the number of elements on subdomains
      //  it is the number of macroelements connected to one microproblem
      //  the number of subdomains--microproblems must be Nproc-1
      //  neldom[i]=j - j finite elements of the macroproblem are connected to the i-th microproblem
      //  neldom[0]=0 - the master processor (myrank=0) contains no microproblem
      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }
      for (i=0;i<Tt->ne;i++){
	neldom[Gtt->stop->eldom[i]]++;
      }
      buffsize1=0;
      for (i=0;i<Nproc;i++){
	if (buffsize1<neldom[i])
	  buffsize1=neldom[i];
      }
      
      buff[1]=buffsize1;
      
      for (i=1;i<Nproc;i++){
	buff[0]=neldom[i];
	MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      }

      buff[0]=neldom[0];
      
      break;
    }
    default:{
      par_print_err(Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
    }
    }
    
  }
  else{
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  //  the number of connections between the macroproblem and microproblem
  nelidom=buff[0];
  //  the maximum number of connections between the macroproblem and microproblem
  buffsize1=buff[1];
  
  //  dimension of the problem, one, two or three dimensional case
  //  the dimension is defined with respect to the first finite element
  dim = Tt->give_dimension (0);
  //  the number of transported materials
  ntm = Tp->ntm;
  
  buffsize2=buffsize1;
  //  buffer for values and gradient components
  buffsize1*=ntm*(dim+1);
  //  buffer for the conductivity and capacity matrices
  buffsize2*=((ntm*dim)*(ntm*dim)+ntm*dim);
  
  fprintf (stdout,"\n procesor %d  nelidom   %ld",Myrank,nelidom);
  fprintf (stdout,"\n procesor %d  buffsize1 %ld",Myrank,buffsize1);
  fprintf (stdout,"\n procesor %d  buffsize2 %ld",Myrank,buffsize2);
  
  
  delete [] buff;
  buff1 = new double [buffsize1];
  buff2 = new double [buffsize2];
  

  if (Myrank==0){
    
    
    //  the number of unknonws
    //  it is the number of unknowns in macroproblem in the case of the master processor
    //  otherwise, it is the number of unknowns in microproblems
    n=Ndoft;
    
    //  array for nodal unknowns
    lhs = Lsrst->give_lhs (0);
    nullv (lhs,n);
    //  array for time derivatives of the nodal unknowns
    tdlhs = Lsrst->give_tdlhs (0);
    nullv (tdlhs,n);
    //  array for the right hand side
    rhs = Lsrst->give_rhs (0);
    nullv (rhs,n);
    
    //  auxiliary arrays
    d = new double [n];
    nullv (d,n);
    p = new double [n];
    nullv (p,n);
    //  vector of prescribed fluxes (right hand side)
    f = new double [n];
    nullv (f,n);
    
    
    //  initial values
    nullv (lhs,n);
    nullv (tdlhs,n);
    nullv (d,n);
    nullv (p,n);
    
    //  coefficient of the trapezoidal rule  
    alpha=Tp->alpha;
    //  computer zero
    zero=Tp->zero;
    
    
    // **************************
    //  main iteration loop  ***
    // **************************
    
    //  loop counter
    i=0;
    
    //  starting time
    Tp->time = Tp->timecont.starttime ();
    //  time increment
    dt = Tp->timecont.initialtimeincr ();
    //  end time
    end_time = Tp->timecont.endtime ();
    
    
    //  nodes - integration points interpolation
    approximation ();
    actual_previous_change ();
    
    //  initialization of material models (if needed)
    Tm->initmaterialmodels();
    
    
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }
  
  
  if (Myrank==0){
    for (ii=1;ii<Nproc;ii++){
      MPI_Send (&end_time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  fprintf (stdout,"\n procesor %d end_time %lf",Myrank,end_time);
  

  do{
    
    if (Myrank==0){
      
      fprintf (stdout,"\n iteration (step) number %ld",i);
      
      //  determination of new time instance
      Tp->time=Tp->timecont.newtime ();
      time=Tp->time;
      //  computation of backward time step
      dt = Tp->timecont.actualbacktimeincr ();
      
      for (ii=1;ii<Nproc;ii++){
	MPI_Send (&time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (&time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    
    if (Myrank==0){
      
      //  update of step number
      i++;
      

      // **********************************
      //  solution of microscale problems
      // **********************************
      for (ii=1;ii<Nproc;ii++){
	
	//  function assembles data from the macroproblem for the ii-th microproblem to the array buff1
	higher_to_lower_level (ii,buffsize1,buff1);
	
	
	//fprintf (stdout,"\n cas %lf",Tp->time);
	//for (long ijk=0;ijk<buffsize1;ijk++){
	  //fprintf (stdout,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	//}
	
	MPI_Send (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
      for (ii=1;ii<Nproc;ii++){
	MPI_Recv(buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	//
	//fprintf (stdout,"\n cas %lf",Tp->time);
	//for (long ijk=0;ijk<buffsize2;ijk++){
	  //fprintf (stdout,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	//}

      	k=stat.MPI_TAG;
	//fprintf (stdout,"\n jouda  %ld  eldom %ld %ld",k,eldom[0],eldom[1]);
	jj=0;
	for (j=0;j<Tt->ne;j++){
	  if (Gtt->stop->eldom[j]==k){
	    //  only elements connected with the k-th microscale are assumed
	    Tm->hommat[j].assemble_matrices (buff2+jj*20,ntm,dim);
	    jj++;
	  }
	}

      }
      
    }
    else{
      
      MPI_Recv(buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      fprintf (stdout,"\n nelidom  %ld",nelidom);
      for (ii=0;ii<nelidom;ii++){
	//  odkomentovat a zavolat spravnou funkci
	//fprintf (stdout,"\n ii   %ld",ii);
	paral_transport_homogenization (buff1+ii*6,buff2+ii*20,buff2+ii*20+16);
      }
      MPI_Send (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    

    if (Myrank==0){
      
      
      // *********************************
      //  solution of macroscale problem
      // *********************************
      
      //  assembling of conductivity matrix
      conductivity_matrix (0);

      //  assembling of capacity matrix
      capacity_matrix (0);
      //  predictor
      //  dd = d_n + (1-alpha) dt v_n
      for (j=0;j<n;j++){
	d[j]=lhs[j]+(1.0-alpha)*dt*tdlhs[j];
      }
      
      //  auxiliary vector
      //  K dd
      Kmat->gmxv (d,p);
      
      //  matrix of the system of equations
      //  C + K alpha dt
      Kmat->scalgm (alpha*dt);
      Kmat->addgm (1.0,*Cmat);
      
      //  right hand side
      //  prescribed nodal fluxes f_{n+1}
      trfel_right_hand_side (0,rhs,n);
      
      //  computation of the right hand side vector
      //  f_{n+1} - K dd
      for (j=0;j<n;j++){
	f[j] = rhs[j] - p[j];
      }
      
      //  solution of the system of algebraic equations
      //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
      //  time derivatives v_{n+1} are solved
      Tp->ssle->solve_system (Gtt,Kmat,tdlhs,f,Outt);

      //  nodal values computed from nodal derivatives
      //  d_{n+1} = dd + alpha dt v_{n+1}
      for (j=0;j<n;j++){
	lhs[j] = d[j] + alpha*dt*tdlhs[j];
      }
      
      solution_correction ();
      //  nodes - integration points interpolation
      approximation ();
      compute_req_valt (0);
      

      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      
      printf ("\n Tp->time   %e",Tp->time);
      printf ("\n i          %ld",i);
    }
    
  }while(time<end_time);

  
  if (Myrank==0){
    delete [] p;
    delete [] d;
    delete [] f;
    
    print_closet();
  }
  
  delete [] buff1;
  delete [] buff2;

}


/**
   function computes nonstationary transport problem
   material parameters are obtained by homogenization
   on microscale, stationary problems are solved
   these microproblems can be defined in every integration point,
   on aggregates of integration points or on aggregates of finite elements
   
   JK, 22.9.2010
*/
void parallel_homogenization ()
{
  long i,j,k,ii,jj,n,ipp,nelidom,buffsize1,buffsize2;
  long *eldom,*neldom,*buff;
  double zero,alpha,*d,*p,*f,*lhs,*tdlhs,*rhs;
  double dt,time,end_time,*buff1,*buff2;
  MPI_Status stat;
  vector gr(2);
  long dim,ntm;

  buff = new long [2];
  //  dimension of the problem, one, two or three dimensional case
  //  the dimension is defined with respect to the first finite element
  dim = Tt->give_dimension (0);
  //  the number of transported materials
  ntm = Tp->ntm;

  //  generation of auxiliary data needed in parallel computation
  if (Myrank==0){
    //  map between finite elements of macroproblem and microproblems
    //  this map is defined by id of material types, it is read from the input file on the master
    //  eldom[i]=j - the i-th finite element of macroproblem is connected to the j-th microproblem
    eldom = new long [Tt->ne];
    for (i=0;i<Tt->ne;i++){
      eldom[i]=Tt->elements[i].idm[0];
    }
    
    //  the number of elements on subdomains
    //  it is the number of macroelements connected to one microproblem
    //  the number of subdomains--microproblems must be Nproc-1
    //  neldom[i]=j - j finite elements of the macroproblem are connected to the i-th microproblem
    //  neldom[0]=0 - the master processor (myrank=0) contains no microproblem
    neldom = new long [Nproc];
    for (i=0;i<Nproc;i++){
      neldom[i]=0;
    }
    for (i=0;i<Tt->ne;i++){
      neldom[eldom[i]+1]++;
    }
    buffsize1=0;
    for (i=0;i<Nproc;i++){
      if (buffsize1<neldom[i])
	buffsize1=neldom[i];
    }
    
    buff[0]=buffsize1;
    
    for (i=1;i<Nproc;i++){
      buff[1]=neldom[i];
      MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
    }
    
  }
  else{
    MPI_Recv (buff,2,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    buffsize1=buff[0];
    //  number of elements connected with this microscale problem
    nelidom=buff[1];
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  buffsize2=buffsize1;
  buffsize1*=6;
  buffsize2*=20;
  
  fprintf (stdout,"\n procesor %d  nelidom   %ld",Myrank,nelidom);
  fprintf (stdout,"\n procesor %d  buffsize1 %ld",Myrank,buffsize1);
  fprintf (stdout,"\n procesor %d  buffsize2 %ld",Myrank,buffsize2);
  

  delete [] buff;
  buff1 = new double [buffsize1];
  buff2 = new double [buffsize2];
  
  if (Myrank==0){
    
    
    //  the number of unknonws
    //  it is the number of unknowns in macroproblem in the case of the master processor
    //  otherwise, it is the number of unknowns in microproblems
    n=Ndoft;
    
    //  array for nodal unknowns
    lhs = Lsrst->give_lhs (0);
    nullv (lhs,n);
    //  array for time derivatives of the nodal unknowns
    tdlhs = Lsrst->give_tdlhs (0);
    nullv (tdlhs,n);
    //  array for the right hand side
    rhs = Lsrst->give_rhs (0);
    nullv (rhs,n);
    
    //  auxiliary arrays
    d = new double [n];
    nullv (d,n);
    p = new double [n];
    nullv (p,n);
    //  vector of prescribed fluxes (right hand side)
    f = new double [n];
    nullv (f,n);
    
    
    //  initial values
    nullv (lhs,n);
    nullv (tdlhs,n);
    nullv (d,n);
    nullv (p,n);
    
    //  coefficient of the trapezoidal rule  
    alpha=Tp->alpha;
    //  computer zero
    zero=Tp->zero;
    
    
    // **************************
    //  main iteration loop  ***
    // **************************
    
    //  loop counter
    i=0;
    
    //  starting time
    Tp->time = Tp->timecont.starttime ();
    //  time increment
    dt = Tp->timecont.initialtimeincr ();
    //  end time
    end_time = Tp->timecont.endtime ();
    
    
    //  nodes - integration points interpolation
    approximation ();
    actual_previous_change ();
    
    //  initialization of material models (if needed)
    Tm->initmaterialmodels();
    
    
    print_initt(-1, "wt");
    print_stept(0,i,Tp->time,NULL);
    print_flusht();
  }
  
  
  if (Myrank==0){
    for (ii=1;ii<Nproc;ii++){
      MPI_Send (&end_time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&end_time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  fprintf (stdout,"\n procesor %d end_time %lf",Myrank,end_time);
  

  do{
    
    if (Myrank==0){
      
      fprintf (stdout,"\n iteration (step) number %ld",i);
      
      //  determination of new time instance
      Tp->time=Tp->timecont.newtime ();
      time=Tp->time;
      //  computation of backward time step
      dt = Tp->timecont.actualbacktimeincr ();
      
      for (ii=1;ii<Nproc;ii++){
	MPI_Send (&time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (&time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    

    if (Myrank==0){
      
      //  update of step number
      i++;
      
      /*
      for (long ijk=0;ijk<4;ijk++){
	fprintf (stdout,"\n lhsi %ld   %le",ijk,Lsrst->lhsi[ijk]);
      }
      */

      // **********************************
      //  solution of microscale problems
      // **********************************
      for (ii=1;ii<Nproc;ii++){

	jj=0;
	for (j=0;j<Tt->ne;j++){

	  if (eldom[j]==ii-1){
	    //  only elements connected with the ii-th microscale are assumed
	    
	    //  id of integration point
	    ipp=Tt->elements[j].ipp[0][0];
	    fprintf (stdout,"\n cas %lf  prvek %ld   cislo ip bodu  %ld",Tp->time,j,ipp);

	    //
	    buff1[jj]=Tm->ip[ipp].av[0];
	    jj++;

	    //
	    Tm->givegrad (0,ipp,gr);
	    buff1[jj]=gr[0];
	    jj++;
	    buff1[jj]=gr[1];
	    jj++;
	    
	    //
	    buff1[jj]=Tm->ip[ipp].av[1];
	    jj++;
	    
	    //
	    Tm->givegrad (1,ipp,gr);
	    buff1[jj]=gr[0];
	    jj++;
	    buff1[jj]=gr[1];
	    jj++;
	    
	  }

	}
	fprintf (stdout,"\n cas %lf  procesor %d  index j %ld",Tp->time,Myrank,j);
	

	for (long ijk=0;ijk<buffsize1;ijk++){
	  fprintf (stdout,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	}

	
	MPI_Send (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
      for (ii=1;ii<Nproc;ii++){
	MPI_Recv(buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	

	fprintf (stdout,"\n cas %lf  ",Tp->time);
	for (long ijk=0;ijk<buffsize2;ijk++){
	  fprintf (stdout,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	}


      	k=stat.MPI_TAG;
	fprintf (stdout,"\n jouda  %ld  eldom %ld %ld",k,eldom[0],eldom[1]);
	jj=0;
	for (j=0;j<Tt->ne;j++){
	  if (eldom[j]==k-1){
	    //  only elements connected with the k-th microscale are assumed
	    Tm->hommat[j].assemble_matrices (buff2+jj*20,ntm,dim);
	    jj++;
	  }
	}
      }
      
    }
    else{

      MPI_Recv(buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      fprintf (stdout,"\n cas %lf nelidom  %ld",Tp->time,nelidom);
      for (ii=0;ii<nelidom;ii++){
	//  odkomentovat a zavolat spravnou funkci
	fprintf (stdout,"\n ii   %ld",ii);
	paral_transport_homogenization (buff1+ii*6,buff2+ii*20,buff2+ii*20+16);
      }
      MPI_Send (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (Myrank==0){
      
      // *********************************
      //  solution of macroscale problem
      // *********************************
      
      //  assembling of conductivity matrix
      conductivity_matrix (0);

      //  assembling of capacity matrix
      capacity_matrix (0);
      //  predictor
      //  dd = d_n + (1-alpha) dt v_n
      for (j=0;j<n;j++){
	d[j]=lhs[j]+(1.0-alpha)*dt*tdlhs[j];
      }
      
      //  auxiliary vector
      //  K dd
      Kmat->gmxv (d,p);
      
      //  matrix of the system of equations
      //  C + K alpha dt
      Kmat->scalgm (alpha*dt);
      Kmat->addgm (1.0,*Cmat);
      
      //  right hand side
      //  prescribed nodal fluxes f_{n+1}
      trfel_right_hand_side (0,rhs,n);
      
      //  computation of the right hand side vector
      //  f_{n+1} - K dd
      for (j=0;j<n;j++){
	f[j] = rhs[j] - p[j];
      }
      
      //  solution of the system of algebraic equations
      //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
      //  time derivatives v_{n+1} are solved
      Tp->ssle->solve_system (Gtt,Kmat,tdlhs,f,Outt);

      //  nodal values computed from nodal derivatives
      //  d_{n+1} = dd + alpha dt v_{n+1}
      for (j=0;j<n;j++){
	lhs[j] = d[j] + alpha*dt*tdlhs[j];
      }
      
      solution_correction ();
      //  nodes - integration points interpolation
      approximation ();
      compute_req_valt (0);
      
      print_stept(0,i,Tp->time,NULL);
      print_flusht();
      
      printf ("\n Tp->time   %e",Tp->time);
      printf ("\n i          %ld",i);
    }
    
  }while(time<end_time);

  
  if (Myrank==0){
    delete [] p;
    delete [] d;
    delete [] f;
    
    print_closet();
  }
  
  delete [] buff1;
  delete [] buff2;
  
}



/**
   function selects appropriate components of higher level of a multiscale
   problem which are sent to lower level
   
   @param llid - lower level id
   @param ncbuff - the number of components in the array buff
   @param buff - array for selected components
   
   JK, 24.3.2011
*/
void higher_to_lower_level (long llid,long ncbuff,double *buff)
{
  long i,j,nc,*counter,ntm,dim;
  double area,ar,*aux;

  counter = new long [1];
  counter[0]=0;
  
  switch (Tp->mami){
  case elem_domain:{
    //  each element is connected with a microproblem
    
    //  loop over the number of elements
    for (i=0;i<Tt->ne;i++){
      
      if (Gtt->stop->eldom[i]==llid){
	//  only elements connected with the microid-th microscale are assumed
	//  the following function is located in the file TRFEL/SRC/elemswitcht
	higher_to_lower_level_elem (i,counter,buff);
      }
    }//  end of loop over the number of elements
    break;
  }
  case aggregate_domain:{
    //  each element is connected with a microproblem
    //  elements in one aggregate are connected with the same microproblem
    
    //  loop over the number of elements
    for (i=0;i<Tt->ne;i++){
      
      if (Gtt->stop->eldom[i]==llid){
	//  the following function is located in the file TRFEL/SRC/elemswitcht
	higher_to_lower_level_elem (i,counter,buff);
      }
    }//  end of loop over the number of elements
    break;
  }
  case material_aggregate_domain:{
    //  each aggregate (elements of one material) is connected with a microproblem
    //  averaged data is sent to the microproblem
      
    //  dimension of the problem, one, two or three dimensional case
    //  the dimension is defined with respect to the first finite element
    dim = Tt->give_dimension (0);
    //  the number of transported materials
    ntm = Tp->ntm;
    aux = new double [ncbuff];

    for (i=0;i<ncbuff;i++){
      buff[i]=0.0;//erase actual position of array buff
    }

    // total area of the aggregate
    ar = 0.0;

    //  loop over the number of elements
    for (i=0;i<Tt->ne;i++){
      
      //  generally, the averaging should take into account the size
      //  of particular finite elements, not only their number
      //  if needed, the function higher_to_lower_level_elem
      //  will return also area or volume of finite elements
      
      
      if (Gtt->stop->eldom[i]==llid){
	//  the following function is located in the file TRFEL/SRC/elemswitcht
	counter[0]=0;
	for (long ijk=0;ijk<ncbuff;ijk++){
	  aux[ijk]=0.0;//erase array
	}
	higher_to_lower_level_elem (i,counter,aux);//array aux is filled with unknowns and their gradients 

	//averaging according to element area/volume:
	area = give_elemarea(i);
	area = fabs(area);
	//fprintf (Outt,"\n elem %ld area = %lf",i,area);

	for (j=0;j<ncbuff;j++){
	  buff[j]+=aux[j]*area;
	}
	ar += area;
      }
    }//  end of loop over the number of elements
    
    //fprintf (Outt,"\n\n pred prumerovanim buff1:");
    //fprintf (Outt,"\n Proc No.= llid = %ld po: \n",llid);
    //for (long ijk=0;ijk<ncbuff;ijk++){
    //  fprintf (Outt,"\n buff1 %ld   %le",ijk,buff[ijk]);
    //}
  
    //  averaging of the values
    fprintf (stdout,"\n\n Averaging of values of material aggregate on processor No. %ld",llid);
    for (i=0;i<ncbuff;i++){
      buff[i]=buff[i]/ar;
    }
    
    //fprintf (Outt,"\n\n zprumerovane buff1:");
    //fprintf (Outt,"\n Proc No.= llid = %ld po: \n",llid);
    //for (long ijk=0;ijk<ncbuff;ijk++){
    //fprintf (Outt,"\n buff1 %ld   %le",ijk,buff[ijk]);
    //}
    
    delete [] aux;
    break;
  }
  case eff_aggregate_domain:{
    //  each aggregate is connected with a microproblem
    //  averaged data is sent to the microproblem
    
    aux = new double [ncbuff];
    for (i=0;i<ncbuff;i++){
      aux[i]=0.0;
      buff[i]=0.0;
    }
    
    //  the number of contributions = the number of elements in the aggregate
    nc=0;
    
    //  loop over the number of elements
    for (i=0;i<Tt->ne;i++){
      
      //  generally, the averaging should take into account the size
      //  of particular finite elements, not only their number
      //  if needed, the function higher_to_lower_level_elem
      //  will return also area or volume of finite elements
      
      
      if (Gtt->stop->eldom[i]==llid){
	//  the following function is located in the file TRFEL/SRC/elemswitcht
	higher_to_lower_level_elem (i,counter,buff);
	counter[0]=0;
	//nc++;

	//averaging according to element area/volume:
	area = give_elemarea(i);
	for (j=0;j<ncbuff;j++){
	  aux[j]+=buff[j]/area;
	}
      }
    }//  end of loop over the number of elements
    
    //  averaging of the values
    //for (i=0;i<nc;i++){
    //buff[i]=aux[i]/nc;
    //}
    
    delete [] aux;
    break;
  }
  default:{
    par_print_err (Myrank,"unknown type of aggregation is required",__FILE__,__LINE__,__func__);
  }
  }
  
  
  delete [] counter;
  
}
