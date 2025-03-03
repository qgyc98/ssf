#include "hpnnpsolvert.h"
#include "hpnpsolvert.h"
#include "hpglobal.h"
#include "seqfilest.h"
#include <string.h>
#include <math.h>
#include "mpi.h"

/**
   function computes non-linear non-stationary transport problem
   material parameters are obtained by homogenization
   on microscale, stationary problems are solved
   these microproblems can be defined on every finite element or
   on aggregates of finite elements, the aggregates can send
   data from all aggregated elements or averaged data which
   has the same form as data from a single element
   
   JK, 1.10.2011
   TKr 20/07/2015
*/
void parallel_homogenization_nonlin_nonstat ()
{
  long i,j,k,kk,ii,jj,ij,n,nelidom,buffsize1,buffsize2,ntm,dim,ini,nsts,stop,lcid,ani,newton=0;
  long *neldom,*buff;
  double zero,alpha,*fi,*fb,*d,*p,*f,*lhst,*tdlhst,*rhst,*lhstb,*tdlhstb;
  double dt,dtdef,dtmin,dtmax,time,end_time,*buff1,*buff2,*err,*thresh,norf;
  MPI_Status stat;
  vector gr(2);
  
  lcid=0;
  
  //  allocation of buffer
  //  buff[0] - the number of connections between the macroproblem and microproblem
  //  buff[1] - the maximum number of connections between the macroproblem and microproblem
  buff = new long [2];

  if (Myrank==0){  
    
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
      
      for (i=1;i<Nproc;i++){
	buff[0]=neldom[i];
	MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      }

      buff[0]=neldom[0];
      
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
      
      buff[0]=1;//tady jeste popremyslet??!!
      buff[1]=1;
      
      for (i=1;i<Nproc;i++){
	MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      }
      break;
    }
    case eff_aggregate_domain:{
      //tady komentar??!!
      neldom = new long [Nproc];
      for (i=0;i<Nproc;i++){
	neldom[i]=0;
      }
      //tady nahrat do pole neldom pocet subdomen neboli agregatu ze seznamu materialu nebo ze seznamu domem???
      for (i=0;i<Tt->ne;i++){
	if(Gtt->stop->eldom[i] != Tt->elements[i].idm[0])//asi tak nejak to bude??!!
	  neldom[Gtt->stop->eldom[i]]++;
      }
      
      //tady upravit????!!!
      //  the aggregate sends averaged data to microproblem
      //  only aggregate sends and receives data, the same
      //  matrices are used on all elements in the aggregate
      buff[0]=1;
      buff[1]=1;
      
      for (i=1;i<Nproc;i++){
	MPI_Send (buff,2,MPI_LONG,i,Myrank,MPI_COMM_WORLD);
      }
      
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
  

  if (Myrank==0){
    
    //  the number of unknonws
    //  it is the number of unknowns in macroproblem in the case of the master processor
    //  otherwise, it is the number of unknowns in microproblems
    n=Ndoft;
    //  maximum number of iterations in inner loop
    ini = Tp->nii;
    //  required norm of vector of unbalanced forces
    err = Tp->errarr;
    //  threshold for the size of the right hand side
    thresh=Tp->threshrhs;
    
    //  array for nodal unknowns
    lhst = Lsrst->give_lhs (0);
    nullv (lhst,n);
    //  array for time derivatives of the nodal unknowns
    tdlhst = Lsrst->give_tdlhs (0);
    nullv (tdlhst,n);
    //  array for the right hand side
    rhst = Lsrst->give_rhs (0);
    nullv (rhst,n);
    
    //  auxiliary arrays
    d = new double [n];
    nullv (d,n);
    p = new double [n];
    nullv (p,n);
    fb = new double [n];
    nullv (fb,n);
    fi = new double [n];
    nullv (fi,n);
    //  vector of prescribed fluxes (right hand side)
    f = new double [n];
    nullv (f,n);
    //  backup of nodal values
    lhstb = new double [n];
    nullv (lhstb,n);
    //  backup of time derivatives of nodal values
    tdlhstb = new double [n];
    nullv (tdlhst,n);
    
    
    //  coefficient of the trapezoidal rule  
    alpha=Tp->alpha;
    //  computer zero
    zero=Tp->zero;
    //  loop counter
    i=0;
    //  number of successful time steps (number of steps without inner iterations)
    nsts=0;
    
    stop=0;
    
    //  starting time
    Tp->time = Tp->timecont.starttime ();
    //  time increment
    dt = Tp->timecont.initialtimeincr ();
    //  end time
    end_time = Tp->timecont.endtime ();
    //  minimum time increment
    dtmin=Tp->timecont.dtmin;
    //  maximum time increment
    dtmax=Tp->timecont.dtmax;
    
    
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

  if (Myrank==0){
    for (ii=1;ii<Nproc;ii++){
      MPI_Send (&ini,1,MPI_LONG,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&ini,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  fprintf (stdout,"\n proc %ld   end_time  %lf",Myrank,end_time);
  

  //synchronization of full-newton method
  //aded by TKr 16/04/2015
  if (Myrank==0){
    if (Tp->trsolv == fullnewtont)
      newton = 1;
    
    for (ii=1;ii<Nproc;ii++){
      MPI_Send (&newton,1,MPI_LONG,ii,Myrank,MPI_COMM_WORLD);
    }
  }
  else{
    MPI_Recv (&newton,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);

 fprintf (stdout,"\n proc %ld   full-newton %ld",Myrank,newton);
 
  // *********************************
  //  main iteration (time) loop  ****
  // *********************************
  do{
    
    if (Myrank==0){
      
      //  determination of new time instance
      Tp->time=Tp->timecont.newtime (dt);
      time=Tp->time;
      
      for (ii=1;ii<Nproc;ii++){
	MPI_Send (&time,1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv (&time,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);

    //  time synchronization
    fprintf (stdout,"\n proc %ld   time  %lf",Myrank,time);

    if (Myrank==0)
      if (Mesprt==1){
	fprintf (stdout,"\n\n --------------------------------------------------------------------------");
	fprintf (stdout,"\n HPTRFEL Time step = %ld,  Time %e,  Time increment = %e",i,Tp->time,dt);
	fprintf (stdout,"\n --------------------------------------------------------------------------\n");

	fprintf (Outt,"\n\n --------------------------------------------------------------------------");
	fprintf (Outt,"\n HPTRFEL Time step = %ld,  Time %e,  Time increment = %e",i,Tp->time,dt);
	fprintf (Outt,"\n --------------------------------------------------------------------------\n");	
      }
    
    if (Myrank==0){
      //  update of step number
      i++;
      
      // update material properties and auxiliary values
      Tm->updateipval ();

      //  backup of attained nodal values and their time derivatives
      for (j=0;j<n;j++){
	lhstb[j]=lhst[j];
	tdlhstb[j]=tdlhst[j];
      }
      
      // **********************************
      //  solution of microscale problems
      // **********************************
      for (ii=1;ii<Nproc;ii++){
	
	//  function assembles data from the macroproblem for the ii-th microproblem to the array buff1
	higher_to_lower_level (ii,buffsize1,buff1);

	MPI_Send (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv(buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);

      for (ii=0;ii<nelidom;ii++){
	
	//fprintf (Outt,"\n\n *********************************************");
	//fprintf (Outt,"\n Myrank  = %ld",Myrank);
	//fprintf (Outt,"\n elidom = %ld",ii);
	//fprintf (Outt,"\n\n pred vypoctem paral_transport_homogenization:");
	//for (long ijk=0;ijk<buffsize1;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	//}
	//fprintf (Outt,"\n");
	//for (long ijk=0;ijk<buffsize2;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	//}
	//fflush(Outt);
	
        //  this function is in the file TRFEL/SRC/homogtrans.cpp
	paral_transport_homogenization (buff1+ii*ntm*(dim+1),buff2+ii*((ntm*dim)*(ntm*dim)+ntm*ntm),
					buff2+ii*((ntm*dim)*(ntm*dim)+ntm*ntm)+(ntm*dim)*(ntm*dim));

	//fprintf (Outt,"\n\n po vypoctu paral_transport_homogenization:");
	//for (long ijk=0;ijk<buffsize1;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	//}
	//fprintf (Outt,"\n");
	//for (long ijk=0;ijk<buffsize2;ijk++){
	//for (long ijk=0;ijk<1;ijk++){
	//fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	//}
	//fflush(Outt);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    if (Myrank==0){
      for (ii=1;ii<Nproc;ii++){
	MPI_Recv(buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
      	k=stat.MPI_TAG;

	jj=0;
	for (j=0;j<Tt->ne;j++){
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
	      
	      //Tm->hommat[j].assemble_matrices (buff2+(ii-1)*((ntm*dim)*(ntm*dim)+ntm*ntm),ntm,dim);//puvodni
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
    }
    MPI_Barrier (MPI_COMM_WORLD);
    

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
	d[j]=lhst[j]+(1.0-alpha)*dt*tdlhst[j];
      }
      
      //  auxiliary vector
      //  K dd
      Kmat->gmxv (d,p);
      
      //  matrix of the system of equations
      //  C + K alpha dt
      Kmat->scalgm (alpha*dt);
      Kmat->addgm (1.0,*Cmat);
      Kmat->copygm (*Cmat);
      
      //  right hand side
      //  prescribed nodal fluxes f_{n+1}
      trfel_right_hand_side (0,rhst,n);
      
      //  computation of the right hand side vector
      //  f_{n+1} - K dd
      for (j=0;j<n;j++){
	f[j] = rhst[j] - p[j];
      }
      
      //  solution of the system of algebraic equations
      //  (C + alpha dt K) v_{n+1} = f_{n+1} - K dd
      //  time derivatives v_{n+1} are solved
      Tp->ssle->solve_system (Gtt,Cmat,tdlhst,f,Outt);
      
      //  nodal values computed from nodal derivatives
      //  d_{n+1} = dd + alpha dt v_{n+1}
      for (j=0;j<n;j++){
	lhst[j] = d[j] + alpha*dt*tdlhst[j];
      }
      
      solution_correction ();
      //  nodes - integration points interpolation
      approximation ();
      compute_req_valt (0);
      
    }//  end of myrank=0
    
    
    // ****************************
    //  iteration for equilibrium
    // ****************************
    for (j=0;j<ini;j++){
      
      if (Myrank==0){
	// **********************************
	//  solution of microscale problems
	// **********************************
	for (ii=1;ii<Nproc;ii++){
	  
	  if(newton == 1){
	    // full newton-raphson
	    fprintf (stdout,"\n Full newton method -> homogenization is performed in each inner loopnelidom\n");
	    //  function assembles data from the macroproblem for the ii-th microproblem to the array buff1
	    higher_to_lower_level (ii,buffsize1,buff1);
	  }
	  
	  MPI_Send (buff1,buffsize1,MPI_DOUBLE,ii,Myrank,MPI_COMM_WORLD);
	}
      }
      else{
	MPI_Recv(buff1,buffsize1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	
	if(newton == 1){
	  fprintf (stdout,"\n nelidom  %ld",nelidom);
	  for (ii=0;ii<nelidom;ii++){
	    //fprintf (Outt,"\n\n *********************************************");
	    //fprintf (Outt,"\n Myrank  = %ld",Myrank);
	    //fprintf (Outt,"\n elidom = %ld",ii);
	    //fprintf (Outt,"\n\n pred vypoctem paral_transport_homogenization:");
	    //for (long ijk=0;ijk<buffsize1;ijk++){
	      //for (long ijk=0;ijk<1;ijk++){
	      //fprintf (Outt,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	    //}
	    //fprintf (Outt,"\n");
	    //for (long ijk=0;ijk<buffsize2;ijk++){
	      //for (long ijk=0;ijk<1;ijk++){
	      //fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	    //}
	    //fflush(Outt);
	    
	    paral_transport_homogenization (buff1+ii*ntm*(dim+1),buff2+ii*((ntm*dim)*(ntm*dim)+ntm*ntm),
					    buff2+ii*((ntm*dim)*(ntm*dim)+ntm*ntm)+(ntm*dim)*(ntm*dim));
	    
	    //fprintf (Outt,"\n\n po vypoctu paral_transport_homogenization:");
	    //for (long ijk=0;ijk<buffsize1;ijk++){
	    //for (long ijk=0;ijk<1;ijk++){
	    //fprintf (Outt,"\n buff1 %ld   %le",ijk,buff1[ijk]);
	    //}
	    //fprintf (Outt,"\n");
	    //for (long ijk=0;ijk<buffsize2;ijk++){
	    //for (long ijk=0;ijk<1;ijk++){
	    //fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
	    //}
	    //fflush(Outt);
	    
	  }
	}// end of full-newton
	
      }
      MPI_Barrier (MPI_COMM_WORLD);
      
      if (Myrank==0){
	for (ii=1;ii<Nproc;ii++){
	  MPI_Recv(buff2,buffsize2,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	  
	  k=stat.MPI_TAG;
	  
	  // full newton-raphson
	  if(newton == 1){
	    jj=0;
	    for (ij=0;ij<Tt->ne;ij++){
	      if (Gtt->stop->eldom[ij]==k){
		//  only elements connected with the k-th microscale are assumed
		
		switch (Tp->mami){
		case elem_domain:{
		  //  each element is connected with its own microproblem
		  Tm->hommat[ij].assemble_matrices (buff2,ntm,dim);
		  break;
		}
		case aggregate_domain:{
		  //  generation of auxiliary data needed in parallel computation
		  Tm->hommat[ij].assemble_matrices (buff2+jj*((ntm*dim)*(ntm*dim)+ntm*ntm),ntm,dim);
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
		  //fprintf (Outt,"\n element  = %ld",j);
		  //fprintf (Outt,"\n\n pred ulozenim matic na materialu c. %ld:",j);
		  //for (long ijk=0;ijk<buffsize2;ijk++){
		  //fprintf (Outt,"\n buff2 %ld   %le",ijk,buff2[ijk]);
		  //}
		  
		  //Tm->hommat[ij].assemble_matrices (buff2+(ii-1)*((ntm*dim)*(ntm*dim)+ntm*ntm),ntm,dim);
		  Tm->hommat[ij].assemble_matrices (buff2,ntm,dim);
		  
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
	  }// end of full-newton
	}//  end loop over the number of processors
      }//  end of the master processor part
      else{
	MPI_Send (buff2,buffsize2,MPI_DOUBLE,0,Myrank,MPI_COMM_WORLD);
      }
      MPI_Barrier (MPI_COMM_WORLD);
      
      if (Myrank==0){
	
	// full newton-raphson
	if (Tp->trsolv == fullnewtont){
	  // matrices are computing in each inner iteration
	  // capacity matrix
	  capacity_matrix (0);
	  
	  //  conductivity matrix
	  conductivity_matrix (0);
	  
	  //  auxiliary vector  K (d+(1-alpha)*dt*v)
	  Kmat->gmxv (d,p);
	  
	  //  matrix of the system of equations
	  //  C + alpha.dt.K
	  Kmat->scalgm (dt*alpha);
	  Kmat->addgm (1.0,*Cmat);
	  Kmat->copygm (*Cmat);
	} //end of full newton-raphson
	
	if (Tp->trestype==lrhst){
	  //  correct right hand side of transport part
	  trfel_right_hand_side (lcid,rhst,n);
	  
	  // Solver computes residuum from system of equations
	  Kmat->gmxv (tdlhst,fi);//new fi vector
	  // Res. = F - K (d+(1-alpha)*dt*v) - (C + alpha*dt*K)*v 
	  for (k=0;k<n;k++){
	    fb[k] = rhst[k] - p[k] - fi[k];
	  }
	}
	
	if (Tp->trestype==fluxest){
	  // Solver computes unbalanced fluxes
	  internal_fluxes (fi,n);//new fi vector
	  //  vector of unbalanced fluxes
	  for (k=0;k<n;k++){
	    fb[k]=fi[k];
	  }
	}
	
	//stop = norm_computation_vec (fb,f,err,thresh,2,1);//tady rozmyslet??!!
	stop = 0;
	norf = ss (fb,fb,n);//tady pokus??!!
	fprintf (stdout,"\n norf = %e   error = %e",norf,err[0]);

	if (norf <= err[0])
	  stop = 1;

	for (ii=1;ii<Nproc;ii++){
	  MPI_Send (&stop,1,MPI_LONG,ii,Myrank,MPI_COMM_WORLD);
	}
      }//  end of the master part
      else{
	MPI_Recv(&stop,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      }//  end of myrank is not zero
      MPI_Barrier (MPI_COMM_WORLD);
      
      if (stop==1){
	break;
      }

      if (Myrank==0){
	
	if (Mesprt != 0)  fprintf (stdout,"\n inner iteration %ld",j);

	//Kmat->diag_check (zero,fb);//temporarily commented
	
	Tp->ssle->solve_system (Gtt,Cmat,fi,fb,Outt);//fi is output now
	
	for (k=0;k<n;k++){
	  tdlhst[k]+=fi[k];
	  lhst[k]+=alpha*dt*fi[k];
	}
	
	//  physically corrected solution
	solution_correction ();    
	//  approximation of nodal values into ontegration points
	approximation ();
	
	
	// vlozil JM 25.4.2008
	// nulleqother ();
	if (Tp->nvs==1 && Tp->pnvs==1)
	  actual_previous_nodval ();
	// ***********************************
      }//  end of myrank=0
    }//  end of the inner iterations
    

    if (Myrank==0){
      //  actual number of performed iterations
      ani=j;
      
      if (ani==0)
	nsts++;
      
      if (ani==ini){
	//  backup is used for actual values
	for (j=0;j<n;j++){
	  lhst[j]=lhstb[j];//restoring results from previous time step
	  tdlhst[j]=tdlhstb[j];//restoring results from previous time step
	}
	
	//  physically corrected solution
	solution_correction ();    
	//  approximation of nodal values into ontegration points
	approximation ();
	
	if (Tp->nvs==1 && Tp->pnvs==1)
	  actual_previous_nodval ();
	// ***********************************
	
	
	//  reduction of the time increment because
	//  inner loop was not able to enforce equilibrium
	dt/=2.0;
	i--;
	Tp->timecont.oldtime ();
	
	if (dt<dtmin){
	  fprintf (stderr,"\n\n time increment is less than minimum time increment");
	  fprintf (stderr,"\n computation fails (file %s, line %d)\n",__FILE__,__LINE__);
	  compute_req_valt (lcid);
          print_stept_forced(lcid, i, Tp->time, rhst);
          print_flusht();
	  break;
	}
      }else{
	dtdef = Tp->timecont.actualforwtimeincr();
	if (nsts==2)
	  {
	    dt*=2.0;
	    nsts=0;
	    
	    if (Mesprt==1)  
	      fprintf (stdout,"\n\n time increment is enlarged because no inner loop was neccessary in previous 3 steps");
	    
	    if (dt<=dtdef)//increase of time increment according to prescribed one
	      dt = dtdef;
	    
	    if (dt>dtmax)//maximum time increment
	      dt = dtmax;
	  }
	
	//////////////////////////////////////////////////
	//  printing of output and graphical informations
	compute_req_valt (lcid);
	print_stept(0,i,Tp->time,NULL);
	print_flusht();
      }
      
      if ((Tp->timecont.isitimptime ()==1) && Tp->hdbcont.save_stat())
	{
	  if (Mesprt==1)
	    fprintf (stdout,"\n Creating HPTRFEL backup file\n");
	  //npsolvert_save (lhst,tdlhst,f,i,Tp->time,dt,Tp->timecont,n);
	  solvert_save (lhst,tdlhst,f,i,Tp->time,dt,Tp->timecont,n);
      }
    }//  end of myrank=0

  }while(time<end_time);
  

  if (Myrank==0){
    delete [] neldom;
    delete [] p;
    delete [] d;
    delete [] f;
    delete [] fi;
    delete [] fb;
    delete [] lhstb;
    delete [] tdlhstb;
    
    print_closet();
  }
  
  delete [] buff1;
  delete [] buff2;

}


