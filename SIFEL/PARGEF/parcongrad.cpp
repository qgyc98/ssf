#include <mpi.h>
#include "parcongrad.h"
#include <string.h>
#include <math.h>
#include <time.h>

/**
   Constructor of class parcongrad
 **/
parcongrad::parcongrad(int np,int mr,long nd, long mes)
{
  mespr = mes;
  nproc=np;  
  myrank=mr;
  ndom=nd;
  ndof=0;  
  nbdof=0;
  nbn=0;  
  maxnbdof=0;
}

parcongrad::~parcongrad()
{
  
}


void parcongrad::initiate(partop *ptop,gtopology *top,FILE *out)
{
  ndof = top->codenum_generation (out);
  // maximal number of DOF on one subdomain
  maxndof = ptop->maxndof;
  // total number of DOFs in the whole problem
  tndof=ptop->tndof;
  // 
  bcn = ptop->gcnd;
  //
  ndofdom = ptop->nud;
  
  
}


/**
   function scaters boundary DOFs from buffer to coarse %vector
   
   @param lv - local %vector
   @param buff - buffer
   @param nd - number of subdomain
   
   JK, 4.1.2006
*/
void parcongrad::buff_coarse (double *cv,double *buff,long nd)
{
  long i,j;
  
  //for (i=0;i<nbdofdom[nd];i++){
  for (i=0;i<ndofdom[nd];i++){
    j=bcn[nd][i]-1;
    cv[j]+=buff[i];
  }
}

/**
   function scaters boundary DOFs from buffer to coarse %vector
   
   @param lv - local %vector
   @param buff - buffer
   @param nd - number of subdomain
   
   JK, 4.1.2006
*/
void parcongrad::coarse_buff (double *cv,double *buff,long nd)
{
  long i,j;
  
  //for (i=0;i<nbdofdom[nd];i++){
  for (i=0;i<ndofdom[nd];i++){
    j=bcn[nd][i]-1;
    buff[i]=cv[j];
  }
}


/**
   function solves system of algebraic equations by parallel
   conjugate gradient method
   
   @param top - topology
   @param gm - %matrix of the system
   @param domproc - domain-processor correspondence
   @param lhs - array containing solution of the system
   @param rhs - array containing right hand side
   @param out - output stream
   @param iv - initial values indicator
   
   iv=0 - initial vector is zero %vector
   iv=1 - initial vector is taken from x array
   
   JK, 4.1.2006
*/
void parcongrad::solve_system (partop */*ptop*/,gtopology */*top*/,gmatrix *gm,
			       long *domproc,double *lhs,double *rhs,FILE *out,long iv)
{
  long i,j,k,buffsize;
  double nom,denom,alpha,beta;
  double *d,*p,*r,*dd,*dp,*dr,*dx,*buff,*master_rhs;
  MPI_Status stat;

  buffsize=maxndof+1;
  
  //  auxiliary vector
  p = new double [ndof];
  //  residuum vector
  r = new double [ndof];
  //  direction vector
  d = new double [ndof];
  //  buffer
  buff = new double [buffsize];
  

  
  if (myrank==0){
    for (i=1;i<nproc;i++){
      MPI_Send (&ni,1,MPI_LONG,i,myrank,MPI_COMM_WORLD); 
    }
  }
  else{
    MPI_Recv (&ni,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);	 

  
  
  if(myrank == 0){
    master_rhs = new double[tndof];
    nullv (master_rhs,tndof);
    //  master contribution
    k=domproc[0];
    buff_coarse (master_rhs,buff,k);
    //  slave contributions
    for (j=1;j<nproc;j++){
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      buff_coarse (master_rhs,buff,k);
    }
    
    if(mespr == 1) fprintf (out,"\n\n\n\n Check of master RHS\n\n");
    double norm=0.0;
    for(i = 0; i < tndof; i++){
      norm+=master_rhs[i];
    }
    if(mespr == 1){
      for(i = 0; i < tndof; i++){
	fprintf (out,"%ld    %le\n",i,master_rhs[i]);
      }
    }
    if(mespr == 1) fprintf (stdout,"\n\n RHS size is %le",norm);

    delete []master_rhs;
  }
  else{
    MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); 
  }
  MPI_Barrier (MPI_COMM_WORLD);	 
  
  /* konec vektoru zatizeni  */




  //  initial vector
  if (iv==0){
    for (i=0;i<ndof;i++){
      lhs[i]=0.0;
    }
  }
  
  //  matrix multiplied by initial vector
  gm->gmxv (lhs,p);
  
  //  initial residuum
  subv (rhs,p,r,ndof);
  
  nullv (buff,buffsize);
  copyv (r,buff,ndof);
  
  if (myrank==0){
    

    //  vector of coarse residuum
    dr = new double [tndof];
    nullv (dr,tndof);
    //  coarse direction vector
    dd = new double [tndof];
    nullv (dd,tndof);
    //  vector of unknowns
    dx = new double [tndof];
    nullv (dx,tndof);
    //  auxiliary vector
    dp = new double [tndof];
    nullv (dp,tndof);
    
    //  master contribution
    k=domproc[0];
    buff_coarse (dr,buff,k);
    
    //  slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      buff_coarse (dr,buff,k);
    }
    
    nom = ss (dr,dr,tndof);
    
    copyv (dr,dd,tndof);
    
    
    for (i=1;i<nproc;i++){
      k=domproc[i];
      nullv (buff,buffsize);
      coarse_buff (dd,buff,k);
      MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD); 
    }
    
    k=domproc[0];
    nullv (buff,buffsize);
    coarse_buff (dd,buff,k);
  }
  else{
    MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); 
    MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
  }
  MPI_Barrier (MPI_COMM_WORLD);	 
  

  copyv (buff,d,ndof);
  
  //  iteration loop
  for (i=0;i<ni;i++){
    
    //  A s = p
    gm->gmxv (d,p);
    
    copyv (p,buff,ndof);
    

    if (myrank==0){
      nullv (dp,tndof);
      
      //  master contribution
      k=domproc[0];
      buff_coarse (dp,buff,k);
      
      //  slave contributions
      for (j=1;j<nproc;j++){
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	k=domproc[stat.MPI_TAG];
	buff_coarse (dp,buff,k);
      }
      
      denom = ss (dd,dp,tndof);
      
      alpha = nom/denom;
      

      //  new approximation of x and r
      for (j=0;j<tndof;j++){
	dx[j]+=alpha*dd[j];
	dr[j]-=alpha*dp[j];
      }
      
      denom=nom;
      
      nom = ss (dr,dr,tndof);
      
      if(mespr == 1) fprintf (stdout,"\n iteration number   %ld    nom %le",i,nom);

      if (nom<err)
	buff[buffsize-1]=100.0;
      
      beta = nom/denom;

      //  new vector of direction
      for (j=0;j<tndof;j++){
	dd[j]=beta*dd[j]+dr[j];
      }
      

      for (j=1;j<nproc;j++){
	k=domproc[j];
	nullv (buff,buffsize-1);
	coarse_buff (dd,buff,k);
	MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD); 
      }
      
      k=domproc[0];
      nullv (buff,buffsize-1);
      coarse_buff (dd,buff,k);
      
    }
    else{
      MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); 
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    MPI_Barrier (MPI_COMM_WORLD);	 
    
    
    copyv (buff,d,ndof);
    
    if (buff[buffsize-1]>1.0){
      
      if (myrank==0){
	for (j=1;j<nproc;j++){
	  k=domproc[j];
	  nullv (buff,buffsize-1);
	  coarse_buff (dx,buff,k);
	  MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD); 
	}
	k=domproc[0];
	nullv (buff,buffsize-1);
	coarse_buff (dx,buff,k);
      }
      else{
	MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      }
      
      copyv (buff,lhs,ndof);
      
      break;
    }
  }
  MPI_Barrier (MPI_COMM_WORLD);	 
  
  
  
//  time_t t1,t2,t3,t4;
  
/*  
  
  
    



  if (myrank==0){
    
    fprintf (out,"\n\n\n\n\n");
    j=domproc[myrank];
    fprintf (out,"\n\n\n Domain %ld",j);
    fprintf (out,"\n time of condensation             %g",buff[0]);
    fprintf (out,"\n time of solution of red. sys.    %g",buff[1]);
    fprintf (out,"\n time of back substitution        %g",buff[2]);
    
    j=domproc[myrank];
    fprintf (stdout,"\n\n\n Domain %ld",j);
    fprintf (stdout,"\n time of condensation             %g",buff[0]);
    fprintf (stdout,"\n time of solution of red. sys.    %g",buff[1]);
    fprintf (stdout,"\n time of back substitution        %g",buff[2]);

    for (i=1;i<nproc;i++){
      MPI_Recv (buff,3,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      j=domproc[stat.MPI_TAG];
      fprintf (out,"\n\n\n Domain %ld",j);
      fprintf (out,"\n time of condensation             %g",buff[0]);
      fprintf (out,"\n time of solution of red. sys.    %g",buff[1]);
      fprintf (out,"\n time of back substitution        %g",buff[2]);

      fprintf (stdout,"\n\n\n Domain %ld",j);
      fprintf (stdout,"\n time of condensation             %g",buff[0]);
      fprintf (stdout,"\n time of solution of red. sys.    %g",buff[1]);
      fprintf (stdout,"\n time of back substitution        %g",buff[2]);
    }
  }
  else{
    MPI_Send(buff,3,MPI_LONG,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
*/
  


}
