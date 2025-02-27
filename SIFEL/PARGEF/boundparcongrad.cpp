#include <mpi.h>
#include "boundparcongrad.h"
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector.h>
#include "petsc_boundcg.h"

/**
   Constructor of class boundparcongrad
   JB
 **/
boundparcongrad::boundparcongrad(int np,int mr,long nd,long mes)
{
  mespr = mes;
  // number of processors
  nproc=np;
  // rank of processors
  myrank=mr;
  // number of subdomain
  ndom=nd;
  // number of degrees of freedom on subdomain
  ndof=0;
  // number of boundary degrees of freedom
  nbdof=0;
  // number of internal degrees of freedom
  nidof = 0;
  // maximum numbre of boundary degrees of freedom
  maxnbdof=0;  
  
  // number of iterations
  ni = 0;
  // error of computation
  err = 0.0;
  // type of preconditioning
  prec = -1;
  // local numbers of internal degrees of freedomw
  lidof= NULL;
  // local numbers of boundary degrees of freedomw
  lbdof = NULL;
  // array containing boundary global code numbers
  bcn = NULL;
  //  array containing numbers of DOFs on boundary nodes
  nbdofdom = NULL;
  // pointer to PETSc
  petscpoint = new PETSC_CONT;
    
}

/**
   Desctructor of class boundparcongrad
   JB
 **/

boundparcongrad::~boundparcongrad()
{
  if(prec == petscilu){
    MatDestroy(petscpoint->fact);
    VecDestroy(petscpoint->r_rhs);
    VecDestroy(petscpoint->z_lhs);
    // finalize of Petsc
    PetscFinalize();
    
    delete []petscpoint->val;
    delete []petscpoint->ix;
  }
  delete petscpoint;
 
}


/**
   Function initiate boundary version of paralell conjugate gradient solver
   @param selnodschur - pointer to selnodes class (selnodes.cpp)
   @param ptop - pointer to patrop class (patrop.cpp)
   @param top - pointer to gtopology class (../GEFEL/gtopology.cpp)
   @param out - output file
   
   JB 14.08.2007
  
*/
void boundparcongrad::initiate (selnodes *selnodschur,partop *ptop,gtopology *top,FILE *out)
{
  long i,j;
  //long *aux;
  
 
  if(prec == petscilu){
    static char help[] = " ";
    PetscInitialize(&Argc,&Argv,(char *)0,help);
  }
  
  if (myrank==0){
    //  array containing numbers of DOFs on boundary nodes
    nbdofdom =  selnodschur -> snndofmas;
    // array containing boundary global code numbers
    bcn =  selnodschur -> cnmas;
    //  number of global DOFs = total number of boundary DOFs
    ndofcp = selnodschur -> tndofsn;
    
//     fprintf (out,"\n\n boundary global code numbers on master\n\n\n");   
//     for (i=0;i<nproc;i++){
//       fprintf (out,"Domain %ld   snndofmas  %ld\n",i,nbdofdom[i]);
//       for(j = 0; j < nbdofdom[i]; j++){
// 	fprintf (out,"cnmas %ld %ld\n",j,bcn[i][j]);
//       }
//     }

  }
  maxnbdof =  selnodschur->maxsnndof;
  
  ndof=ptop->ndof;
  nidof=ptop->nidof;
  nbdof=ptop->nbdof;
  
}



/*
  Function computes scalar product of internal parts of vectors u anb v
  @param u - the first vector
  @param v - the second vector
  @out sp - scalar product
  JB
*/

double boundparcongrad::v_ixu_i(double *u,double *v)
{
  long i;
  double sp;
  sp = 0.0;
  for(i = 0; i < nidof; i++){
    // kvuli cislovani od 0 - natace C
    //j= lidof[i]-1;
    sp+=u[i]*v[i];
  }
  
  return sp;
}


/*
  Function selects boundary DOFs  from vector v to buff on subdomain
  @param buff - buffer (target vector)
  @param v  - vector with  boundary DOFs (source vector)
  JB
*/
void boundparcongrad::select_bound(double *buff,double *v)
{
  long i;
  
  for(i = nidof; i < ndof; i++){
    buff[i-nidof] = v[i];
  }
}

/*
  Function adds boundary DOFs from buff to v on slaves
  @param buff - buffer (source vector)
  @param v - target vector
*/

void boundparcongrad::add_bound(double *buff,double *v)
{
  long i;
  for(i = 0; i < nbdof; i++){
    v[i+nidof] = buff[i];
  }
}

/*
  Function adds buff into coarse vector on master
  @param buff - buffer (source vector)
  @param cv - coarse vector (target vector)
  @param nd - subdomain number
  
  JB
*/

void boundparcongrad::add_master_buff (double *buff,double *cv,long nd)
{
  long i,j;
  for (i=0;i<nbdofdom[nd];i++){
    // kvuli cislovani od 0 - natace C
    j=bcn[nd][i]-1;
    cv[j]+=buff[i];
  }
  
}

/*
  Function selects appropriate boundary DOFs form cv to buff 
  @param buff - buffer (target vector)
  @param cv - coarse vector (source vector)
  @param nd - subdomain number
  JB
*/

void boundparcongrad::select_master_buff(double *cv,double *buff,long nd)
{
  long i,j;
  for (i=0;i<nbdofdom[nd];i++){
    // kvuli cislovani od 0 - natace C
    j=bcn[nd][i]-1;
    buff[i] =  cv[j];
  }
  
}

/**
   Function computes incomplete factorization of subdomain matrix with PETSC-ILU
   with the help of PETSc library
   @param out - output file for logs
   @param gm - pointer to gmatrix class (../GEFEL/gmatrix.cpp)
   
JB
*/

void boundparcongrad::ilu_mat_petsc(gmatrix *gm,FILE *out)
{
  
  long adr1,adr2;
  
  // PETSC Int
  PetscInt i,j,jj,nz,m;
  PetscInt *nnz;
  IS row,col;
  // PETSC error code 
  //PetscErrorCode ierr;
  // PETSC scalar
  PetscScalar v;
  //
  //PetscTruth flg;

  //fprintf(out,"faktorizace\n");
  
  switch(gm->ts){
  case compressed_rows:{
    //fprintf(out,"faktorizace %ld\n",gm->ts);
    // nnz assembling
    nz = 0;
    nnz = new PetscInt[gm->cr->n];
    //fprintf(out,"n = %ld\n",gm->cr->n);
    //fprintf(out,"negm = %ld\n",gm->cr->negm);
    //fprintf(out,"adr[%ld] = %ld\n",0,gm->cr->adr[0]);
    for(i = 0; i < gm->cr->n; i++){
      //fprintf(out,"adr[%ld] = %ld\n",i,gm->cr->adr[i+1]);
      adr1=gm->cr->adr[i];
      adr2=gm->cr->adr[i+1];
      nnz[i]=adr2-adr1;
    }
    
    // test print
    //fprintf(out,"Test print of nnz array\n");
    //for(i = 0; i < gm->cr->n; i++){
    // fprintf(out,"nnz[%ld]=%ld\n",i,nnz[i]);
    //}
    
    // create of PETSC matrix
    // compressed rows storage
    MatCreateSeqAIJ(PETSC_COMM_SELF,gm->cr->n,gm->cr->n,nz,nnz,&petscpoint->fact);
    
    // set value of matrix
    for(i = 0; i < gm->cr->n; i++){
      adr1=gm->cr->adr[i];
      adr2=gm->cr->adr[i+1];
      for(j = adr1; j < adr2; j++){
	v = gm->cr->a[j];
	jj=gm->cr->ci[j];
	//fprintf(out,"ci[%ld] = %ld\n",j,gm->cr->ci[j]);
	MatSetValues(petscpoint->fact,1,&i,1,&jj,&v,INSERT_VALUES);
      }
    }
    ISCreateStride(PETSC_COMM_WORLD,gm->cr->n,0,1,&row);
    ISCreateStride(PETSC_COMM_WORLD,ndof,0,1,&col);
    break;
  }
  case symm_comp_rows:{
    //fprintf(out,"faktorizace %ld\n",gm->ts);
    nz = 0;
    nnz = new PetscInt[gm->scr->n];
    for(i = 0; i < gm->scr->n; i++){
      nnz[i] = 0;
    }
    //fprintf(out,"n = %ld\n",gm->scr->n);
    //fprintf(out,"negm = %ld\n",gm->scr->negm);
    //fprintf(out,"adr[%ld] = %ld\n",0,gm->scr->adr[0]);
    
    for(i = 0; i < gm->scr->n; i++){
      //fprintf(out,"adr[%ld] = %ld\n",i,gm->scr->adr[i+1]);
      adr1=gm->scr->adr[i];
      adr2=gm->scr->adr[i+1];
      nnz[i] = adr2-adr1;
      for(j = i+1; j < gm->scr->n; j++){
	for(m = adr1=gm->scr->adr[j]; m < gm->scr->adr[j+1]; m++){
	  if(gm->scr->ci[m] == i){
	    nnz[i]++;
	  }
	}
      }
    }
    
    
    // test print
    //fprintf(out,"Test print of nnz array\n");
    //for(i = 0; i < gm->scr->n; i++){
    //fprintf(out,"nnz[%ld]=%ld\n",i,nnz[i]);
    //}
    

    // create of PETSC matrix
    // compressed rows storage
    MatCreateSeqAIJ(PETSC_COMM_SELF,gm->scr->n,gm->scr->n,nz,nnz,&petscpoint->fact);
    
    // set value of matrix
    for(i = 0; i < gm->scr->n; i++){
      adr1=gm->scr->adr[i];
      adr2=gm->scr->adr[i+1];
      for(j = adr1; j < adr2; j++){
	v = gm->scr->a[j];
	jj=gm->scr->ci[j];
	//fprintf(out,"ci[%ld] = %ld\n",j,gm->scr->ci[j]);
	MatSetValues(petscpoint->fact,1,&i,1,&jj,&v,INSERT_VALUES);
	MatSetValues(petscpoint->fact,1,&jj,1,&i,&v,INSERT_VALUES);
      }
    }
    ISCreateStride(PETSC_COMM_WORLD,gm->scr->n,0,1,&row);
    ISCreateStride(PETSC_COMM_WORLD,ndof,0,1,&col);
    break;
  }
  default:{
    //parrerror
    par_print_err(myrank,"Unknown storage type\n", __FILE__, __LINE__, __func__);
    break;
  }
  }
  
  
  MatAssemblyBegin(petscpoint->fact,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(petscpoint->fact,MAT_FINAL_ASSEMBLY);
  //MatView(petscpoint->fact,0);
  
  // set properties of ILU factorization
  MatFactorInfo info;
  info.fill = 1.0;
  //info.dtcol = 1.0;
  //info.zeropivot = 1.e-14;
  //info.pivotinblocks = 1.0;
  //info.shiftpd = PETSC_FALSE;
  info.levels=0;
  
   
  MatILUFactor(petscpoint->fact,row,col,&info);
  
  // creation of PETSC vectors
  i = ndof;
  VecCreateSeq(PETSC_COMM_SELF,i,&petscpoint->r_rhs);
  VecAssemblyBegin(petscpoint->r_rhs);
  VecAssemblyEnd(petscpoint->r_rhs);
  i = ndof;
  VecCreateSeq(PETSC_COMM_SELF,i,&petscpoint->z_lhs);
  VecAssemblyBegin(petscpoint->r_rhs);
  VecAssemblyEnd(petscpoint->r_rhs);
       
  //creation of auxiliary PETSC array
  petscpoint->ix = new PetscInt[ndof];
  petscpoint->val = new PetscScalar[ndof];

  delete []nnz;
}

/**
   Function performs back substitution with the help of PETSc library
   @param r - residual vector
   @param z - preconditoned vector z = PETscILUMAT^-1*r
   @param ndof - number of degrees of freedom on subdomain
   @param out - output file for logs

   JB
 **/

void boundparcongrad::solve_fact_mat_petsc(double *r,double *z,long ndof,FILE *out)
{
  
  PetscInt i,ngv;
  PetscScalar v;
  //time_t tsetmats,tsetmate;
  //double sec1,sec2;
  
  //for(i = 0; i < ndof; i++){
  //fprintf(out,"r %f\n",r[i]);
  //}
  
  
  for(i = 0; i < ndof; i++){
    v = r[i];
    VecSetValue(petscpoint->r_rhs,i,v,INSERT_VALUES);
  }

  //VecView(petscpoint->r_rhs,0);
  MatSolve(petscpoint->fact,petscpoint->r_rhs,petscpoint->z_lhs);
  //VecView(petscpoint->z_lhs,0);
  
  nullv(z,ndof);
  i = ndof;
  //ix = new PetscInt[ndof];
  //val = new PetscScalar[ndof];
  for(i = 0; i < ndof; i++){
    petscpoint->ix[i] = i;
  }
  ngv = ndof;
  VecGetValues(petscpoint->z_lhs,ngv,petscpoint->ix,petscpoint->val);
  
  for(i = 0; i < ndof; i++){
    //fprintf(out,"val %f\n",PetscRealPart(val[ix[i]]));
    z[petscpoint->ix[i]] =petscpoint->val[petscpoint->ix[i]];// PetscRealPart(val[i]);
  }
  

}

/**
   Function initialises diagonal preconditioning
   @param gm - pointer to gmatrix class (../GEFEL/gmatrix.cpp)
   @param precvec - vector for preconditionig
   @param out - output file for logs

   JB
 **/

void boundparcongrad::initialize_diag_prec(gmatrix *gm,double *precvec,FILE *out)
{
  long i,j,jj,adr1,adr2;

  switch(gm->ts){
  case compressed_rows:{
    for(i = 0; i < gm->cr->n; i++){
      adr1=gm->cr->adr[i];
      adr2=gm->cr->adr[i+1];
      for(j = adr1; j < adr2; j++){
	jj=gm->cr->ci[j];
	if(jj == i){
	  precvec[i] = gm->cr->a[j];
	  precvec[i]=1.00/precvec[i];
	  break;
	}
	
      }
      //precvec[i] = 1.0;
    }
    break;
  }
  case symm_comp_rows:{
    for(i = 0; i < gm->scr->n; i++){
      adr1=gm->scr->adr[i];
      adr2=gm->scr->adr[i+1];
      for(j = adr1; j < adr2; j++){
	jj=gm->scr->ci[j];
	if(jj == i){
	  precvec[i] = gm->scr->a[j];
	  precvec[i] = 1.00/precvec[i];
	}
      }
      
      if(mespr == 1) fprintf(out,"prec[%ld]= %le\n",i,precvec[i]);
      //precvec[i] = 1.0;
    }
    break;
  }
  default:{
    printf("Unknown storage type of matrix");
    break;
  }
  }
  
}

/**
   Function performs diagonal preconditioning
   @param precvec - vector for preconditionig
   @param invec - input vector 
   @param outvec - output vector outvec=precvec*inputvec
   @param out - output file for logs
   JB
 **/


void boundparcongrad::jacobi_precondition(double *precvec,double *invec,double *outvec,long ndof,FILE *out)
{
  long i;
  for(i = 0; i < ndof; i++){
    outvec[i] = precvec[i]*invec[i];
    if(mespr == 1) fprintf(out,"precvec[%ld]= %le incvec[%ld] = %le outvec[%ld] = %le\n",i,precvec[i],i,invec[i],i,outvec[i]);
  }
}



/**
   Function solves system of algebraic equations by boundary version of parallel
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
   
  
*/

void boundparcongrad::solve_system (partop *ptop,gtopology *top,gmatrix *gm,long *domproc,double *lhs,double *rhs,FILE *out,long iv)
{
  long i,j,k,iter;
  long buffsize;
  double nom_i,denom_i,alpha,beta,norrhs;
  double *p_dom, *r_dom, *d_dom,*x_dom,*buff,*z_dom,*aux;
  double *precvec;
  double tbacksub;
  double  nom_master_i,nom_master_b,nom,denom;
  double *r_bound, *d_bound, *x_bound, *p_bound,*z_bound;
  double denom_master_i,denom_master_b;
  double tfacts,tfacte,tsols,tsole,titers,titere,ts,te;
  double nory;
  
  MPI_Status stat;

  tbacksub=0.0;
  tsols = clock();
  switch(prec){
  case noprec:{
    break;
  }
  case petscilu:{
    tfacts = clock();
    ilu_mat_petsc(gm,out);
    tfacte = clock();
    if(mespr == 1) fprintf (out,"Time of ILU factorization - ilu_mat_petsc function %le s\n",(tfacte-tfacts)/(double)CLOCKS_PER_SEC);
    break;
  }
  case pardiagprec:{
    precvec = new double[ndof];
    initialize_diag_prec(gm,precvec,out);
    break;
  }
  }
    
  
  // arrays defined on each subdomain
  // ndof - number of degrees of freedom
  //  auxiliary vector
  p_dom = new double [ndof];
  nullv (p_dom,ndof);
  //  residuum vector
  r_dom = new double [ndof];
  nullv (r_dom,ndof);
  //
  z_dom = new double [ndof];
  nullv (z_dom,ndof);
  //  direction vector
  d_dom = new double [ndof];
  nullv (d_dom,ndof);
  //
  x_dom = new double [ndof];
  nullv (x_dom,ndof);
  //
  aux = new double[maxnbdof];
  //  buffer
  buffsize=2*maxnbdof+2;
  // 0 - maxnbdof - r
  // maxnbdof - 2*maxnbdof-1 - z
  //maxnbdof+1 - nom_i or denom_i
  //maxnbdof+2 - error of computation
  buff = new double [buffsize];
  
  //  initial vector
  if (iv==0){
    for (i=0;i<ndof;i++){
      lhs[i]=0.0;
    }
  }
  
  nory = 0.0;
  for(i = 0; i < ndof; i++){
    nory +=rhs[i];
  }
  if(mespr == 1) fprintf (stdout,"\n nory %le",nory);
  if(myrank == 0){
    norrhs = nory;
    for (i=1;i<nproc;i++){
      MPI_Recv (&nory,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      norrhs +=nory;
    }
    if(mespr == 1) fprintf (stdout,"\n norrhs %le",nory);
  }
  else{
    MPI_Send(&nory,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); 
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // vazeni prave strany
  //for(i = 0; i < ndof; i++){
  //rhs[i]/=double(ptop->dofmultip[i]);
  //}
  
  

  //  matrix multiplied by initial vector
  // p_dom,j = K,j*x,j
  switch(prec){
  case noprec:{
    gm->gmxv (x_dom,p_dom); 
    break;
  }
  case petscilu:{
    gm->gmxv (x_dom,p_dom); 
    //petsc_mxv(x_dom,p_dom,ndof,out);
    break;
  }
  case pardiagprec:{
    gm->gmxv (x_dom,p_dom); 
    break;
  }
  }
  
  
  //  initial residuum
  // r_dom,j = rhs,j - p_dom,j
  subv (rhs,p_dom,r_dom,ndof);
  //for(i = 0; i < ndof; i++){
  //r_dom[i] = rhs[i];
  //}
  
  // precondition
  switch(prec){
  case noprec:{
    for(i = 0; i < ndof; i++){
      z_dom[i] = r_dom[i];
    }
    break;
  }
  case petscilu:{
    solve_fact_mat_petsc(r_dom,z_dom,ndof,out);
    //for(i = 0; i < ndof; i++){
    //z_dom[i] = r_dom[i];
    //}
    break;
  }
  case pardiagprec:{
    jacobi_precondition(precvec,r_dom,z_dom,ndof,out);
    break;
  }
  }
  
  // initial direction vector
  // d_dom,j = z_dom,j
  copyv (z_dom,d_dom,ndof);
  
  // test print
  //for(i = 0; i < ndof; i++){
  //fprintf (out,"%ld rhs %le  p %le  r %le z %le  d %le\n",i+1,rhs[i],p_dom[i],r_dom[i],z_dom[i],d_dom[i]);
  //}
  
  // nom_i,j = r_i,j^T*r_i,j
  nom_i = v_ixu_i(r_dom,z_dom);

  // test print
  //fprintf (out,"nom_i %le\n",nom_i);
  
  // r_b,j = Lb,j*r,j
  nullv(aux,maxnbdof);
  select_bound(aux,r_dom);
  for(i = 0; i < nbdof; i++){
    buff[i] = aux[i];
  }
  // z_b,j = Lb,j*z,j
  nullv(aux,maxnbdof);
  select_bound(aux,z_dom);
  k = 0;
  for(i = maxnbdof; i < maxnbdof+nbdof; i++){
    buff[i] = aux[k];
    k++;
  }
  
  buff[2*maxnbdof] = nom_i;
  
  
  if(myrank == 0){
    
    // arrays defined only on master
    //  vector of coarse residuum
    r_bound = new double [ndofcp];
    nullv (r_bound,ndofcp);
    //
    z_bound = new double [ndofcp];
    nullv (z_bound,ndofcp);
    //  coarse direction vector
    d_bound = new double [ndofcp];
    nullv (d_bound,ndofcp);
    //  vector of unknowns
    x_bound = new double [ndofcp];
    nullv (x_bound,ndofcp);
    //  auxiliary vector
    p_bound = new double [ndofcp];
    nullv (p_bound,ndofcp);
    
    //  master contribution
    nom_master_i = 0.0;
    nom_master_i += buff[2*maxnbdof];
    k=domproc[0];
    // r_b=Sum_j=1^m L_b,j^T*r_b,j
    nullv (aux,maxnbdof);
    for(j = 0; j < nbdofdom[k]; j++){
      aux[j] = buff[j];
    }
    add_master_buff (aux,r_bound,k);
    // z_b=Sum_j=1^m L_b,j^T*z_b,j
    nullv (aux,maxnbdof);
    k = 0;
    for(j = maxnbdof; j < maxnbdof+nbdofdom[k]; j++){
      aux[j-maxnbdof] = buff[j];
    }
    add_master_buff (aux,z_bound,k);
    
    // slave contributions
    for (i=1;i<nproc;i++){
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      k=domproc[stat.MPI_TAG];
      //nom_master_i = Sum_j=1^m nom_i,j
      nom_master_i += buff[2*maxnbdof];
      nullv (aux,maxnbdof);
      for(j = 0; j < nbdofdom[k]; j++){
	aux[j] = buff[j];
      }
      add_master_buff (aux,r_bound,k);
      // z_b=Sum_j=1^m L_b,j^T*z_b,j
      nullv (aux,maxnbdof);
      for(j = maxnbdof; j < maxnbdof+nbdofdom[k]; j++){
	aux[j-maxnbdof] = buff[j];
      }
      add_master_buff (aux,z_bound,k);
    }
    
    //test print
    //fprintf (out,"nom_master_i %le\n",nom_master_i);
    
    // test print r_bound
    //for(i = 0; i < ndof; i++){
    //fprintf (out,"%ld r_bound %le\n",i+1,r_bound[i]);
    //}

    //nom_master_b = r_b^T*r_b
    nom_master_b = ss(r_bound,z_bound,ndofcp);
    
    // test print
    //fprintf (out,"nom_master_b %le\n",nom_master_b);
    
    // nom = nom_master_i+nom_master_b
    nom =  nom_master_i+nom_master_b;
    
    //test print nom
    //fprintf (out,"nom %le\n",nom);
    
    // d_b = z_b
    copyv (z_bound,d_bound,ndofcp);
    
    //test print d_bound
    //for(i = 0; i < ndofcp; i++){
    //fprintf (out,"%ld r_bound %lf d_bound %lf\n",i+1,r_bound[i],d_bound[i]); 
    //}
    
    //d_b,j=L_b,s*d_b
    // slaves
    for (i=1;i<nproc;i++){
      k=domproc[i];
      nullv (buff,maxnbdof);
      select_master_buff(d_bound,buff,k);
      MPI_Send (buff,buffsize,MPI_DOUBLE,i,myrank,MPI_COMM_WORLD); 
    }
    
    // master
    k=domproc[0];
    nullv (buff,buffsize);
    select_master_buff(d_bound,buff,k);
    
  }
  else{
    
    MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);  
    MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  // d_dom,j = d_i,j, d_b,j
  add_bound(buff,d_dom);
  
  //test print
  // for(i = 0; i < ndof; i++){
  // fprintf (out,"%ld d_dom %lf\n",i+1,d_dom[i]); 
  // }
  delete []buff;
  delete []aux;
  //  buffer
  buffsize=maxnbdof+2;
  //maxnbdof+1 - nom_i or denom_i
  //maxnbdof+2 - error of computation
  buff = new double [buffsize];
  
  
  // iteration loop
  titers = clock();
  for(i = 0; i < ni; i++){
    //fprintf(out,"%ld iterace\n",i);
    //  p_dom,j = Kj*d_dom,j
    gm->gmxv (d_dom,p_dom);
    
    //  denom_i,j = s_i,j^T*d_i,j
    denom_i =  v_ixu_i(d_dom,p_dom);
    
    // test print
    //fprintf (out,"denom_i %lf\n",denom_i);
    
    nullv (buff,maxnbdof);
    // p_b,j = Lb,j*p,j
    select_bound(buff,p_dom);
    
    buff[maxnbdof] = denom_i;
    
    
    if(myrank == 0){
       
      nullv (p_bound,ndofcp);
       
      //  master contribution
      k=domproc[0];
      denom_master_i = 0.0;
      denom_master_i += buff[maxnbdof];
      // p_b=Sum_j=1^m L_b,j^T*p_b,j
      add_master_buff (buff,p_bound,k);
      
      // slave contributions
      for (j=1;j<nproc;j++){
        MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        k=domproc[stat.MPI_TAG];
        //denom_master_i = Sum_j=1^m nom_i,j
        denom_master_i += buff[maxnbdof];
        // p_b=Sum_j=1^m L_b,j^T*p_b,j
        add_master_buff (buff,p_bound,k);
      }
      
      //test print p_bound
      //for(j = 0; j < ndofcp; j++){
      //fprintf (out,"%ld p_bound %lf\n",j+1,p_bound[j]); 
      //}
      
      //nom_master_b = r_b^T*r_b
      denom_master_b = ss(d_bound,p_bound,ndofcp);
     
      // nom = nom_master_i+nom_master_b
      denom =  denom_master_i+denom_master_b;
      
      //test print
      //fprintf (out,"denom %le\n",denom);
      
      // alpha = nom/denom;
      alpha = nom/denom;
      denom = nom;
      
      
      buff[maxnbdof] = alpha;
      //test print
      //fprintf (out,"alpha %le\n",alpha);
      
      // new approximation of x_b and r_b
      // r_b^k+1 = r_b^k-alpha*p_b
      // x_b^k+1 = x_b^k+alpha*d_b
      for (j=0;j<ndofcp;j++){
        x_bound[j]+=alpha*d_bound[j];
        r_bound[j]-=alpha*p_bound[j];
      }
    
      //test print
      //for(j = 0; j < ndofcp; j++){
      //fprintf (out,"%ld r_bound %lf\n",j+1,r_bound[j]); 
      //}
      //for(j = 0; j < ndofcp; j++){
      //  fprintf (out,"%ld x_bound %lf\n",j+1,x_bound[j]); 
      // }
      
      // slaves
      for (j=1;j<nproc;j++){
	// precondition
	switch(prec){
	case noprec:{
	  break;
	}
	case petscilu:{
	  k=domproc[j];
	  nullv (buff,maxnbdof);
	  select_master_buff(r_bound,buff,k);
	  //  for(k = 0; k < maxnbdof; k++){
	  //fprintf (out,"%ld %ld r_bound %lf\n",j,k+1,buff[k]); 
	  //}
	  break;
	}
	case pardiagprec:{
	  k=domproc[j];
	  nullv (buff,maxnbdof);
	  select_master_buff(r_bound,buff,k);
	  break;
	}
	}
	MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD); 
      }
      // master
      switch(prec){
      case noprec:{
	break;
      }
      case petscilu:{
	k=domproc[0];
	select_master_buff(r_bound,buff,k);
	//for(k = 0; k < maxnbdof; k++){
	//fprintf (out,"%ld %ld r_bound %lf\n",0,k+1,buff[k]); 
	//}
	break;
      }
      case pardiagprec:{
	k=domproc[0];
	nullv (buff,maxnbdof);
	select_master_buff(r_bound,buff,k);
	break;
      }
      }
      
    }
    else{
      MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);   
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat); 
    }
    MPI_Barrier (MPI_COMM_WORLD);
    
    alpha = buff[maxnbdof];
    
    // new approximation of x_dom and r_dom - internal DOFs
    for (j=0;j<ndof;j++){
      x_dom[j]+=alpha*d_dom[j];
      r_dom[j]-=alpha*p_dom[j];
    }
    
    
    
    // precondition
    switch(prec){
    case noprec:{
      for(j = 0; j < ndof; j++){
	z_dom[j] = r_dom[j];
      }
      break;
    }
    case petscilu:{
      // r_dom,j = d_i,j, r_b,j
      add_bound(buff,r_dom);
      // preconditining z_k+1=C^-1*r_k+1
      ts = clock();
      solve_fact_mat_petsc(r_dom,z_dom,ndof,out);
      te = clock();
      if(mespr == 1) fprintf (out,"Time of back substitution %le s\n",(te-ts)/(double)CLOCKS_PER_SEC); 
      tbacksub += ((te-ts)/(double)CLOCKS_PER_SEC);
      //for(j = 0; j < ndof; j++){
      //z_dom[j] = r_dom[j];
      //}
      
      nullv (buff,maxnbdof);    
      // z_b,j = Lb,j*z,j
      select_bound(buff,z_dom);
      break;
    }
    case pardiagprec:{
      // r_dom,j = d_i,j, d_b,j
      add_bound(buff,r_dom);
      // preconditining z_k+1=C^-1*r_k+1
      ts = clock();
      jacobi_precondition(precvec,r_dom,z_dom,ndof,out);
      te = clock();
      if(mespr == 1) fprintf (out,"Time of precondition %le s\n",(te-ts)/(double)CLOCKS_PER_SEC); 
      tbacksub += ((te-ts)/(double)CLOCKS_PER_SEC);
      //for(j = 0; j < ndof; j++){
      //z_dom[j] = r_dom[j];
      //}
      
      nullv (buff,maxnbdof);    
      // z_b,j = Lb,j*z,j
      select_bound(buff,z_dom);
      break;
    }
    }
    
    
    // nom_i,j^k+1 = (r_i,j^k+1)^T*r_i,j^k+1
    nom_i = v_ixu_i(r_dom,z_dom);
    buff[maxnbdof] = nom_i;
   
    
    if(myrank == 0){
      // master contributions
      nom_master_i = buff[maxnbdof];
      k=domproc[0];
      nullv(z_bound,ndofcp);
      switch(prec){
      case noprec:{
	for(j = 0; j < ndofcp; j++){
	  z_bound[j] = r_bound[j];
	}
	break;
      }
      case petscilu:{
	// z_b=Sum_j=1^m L_b,j^T*z_b,j
	add_master_buff (buff,z_bound,k);
	break;
      }
      case pardiagprec:{
	// z_b=Sum_j=1^m L_b,j^T*z_b,j
	add_master_buff (buff,z_bound,k);
	break;
      }
	
      }
      
      // slave contributions
      for (j=1;j<nproc;j++){
        MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
        k=domproc[stat.MPI_TAG];
        nom_master_i += buff[maxnbdof];
	switch(prec){
	case noprec:{
	  break;
	}
	case petscilu:{
	  // z_b=Sum_j=1^m L_b,j^T*z_b,j
	  add_master_buff (buff,z_bound,k);
	  break;
	}
	case pardiagprec:{
	// z_b=Sum_j=1^m L_b,j^T*z_b,j
	add_master_buff (buff,z_bound,k);
	break;
      }
	}
	
      }
      
      
      //for(j = 0; j < ndofcp; j++){
      //fprintf (out,"%ld z_bound %lf r_bound %lf\n",j+1,z_bound[j],r_bound[j]); 
      //}
      
      //nom_master_b = r_b^T*z_b
      nom_master_b = ss(r_bound,z_bound,ndofcp);
      
      // nom = nom_master_i+nom_master_b
      nom =  nom_master_i+nom_master_b; 
      //norrhs = nom;
      //test print
      //fprintf (out,"nom k+1 %le\n",nom);
      if(sqrt(nom/norrhs) < err || i == ni-1){
	//if (nom<err || i == ni-1 ){
        buff[maxnbdof+1]=10000.0;
	iter = i;
      }
      
      //fprintf (stdout,"\n iteration number   %ld    nom %le",i,nom);
      if(mespr == 1) fprintf (stdout,"\n iteration   %ld   norres  %le  norrhs %le   norres/norrhs %e",i,sqrt(nom),sqrt(norrhs),sqrt(nom/norrhs));
      beta = nom/denom;
      
      //test print
      //fprintf (out,"beta %le\n",beta);
      
      //  new aproximation of vector d_b
      for (j=0;j<ndofcp;j++){
        d_bound[j]=beta*d_bound[j]+z_bound[j];
      }
      
      //test print
      //for(j = 0; j < ndofcp; j++){
      //fprintf (out,"%ld d_bound %lf\n",j+1,d_bound[j]); 
      //}
      
      for (j=1;j<nproc;j++){
        k=domproc[j];
        nullv (buff,maxnbdof);
        //d_b,j=L_b,s*d_b
        select_master_buff(d_bound,buff,k);
        
        buff[maxnbdof]=beta;
        
        MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD); 
      }
      
      k=domproc[0];
      nullv (buff,maxnbdof);
      //d_b,j=L_b,s*d_b 
      select_master_buff(d_bound,buff,k);
      
      buff[maxnbdof] = beta;
      
    }
    else{
      MPI_Send (buff,buffsize,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD); 
      MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
    }
    
    // d_dom,j = d_i,j, d_b,j
    add_bound(buff,d_dom);
    
    beta = buff[maxnbdof];
    
    //fprintf (out,"beta_dom %le\n",beta);
    
    //  new aproximation of vector d_dom,j
    for (j=0;j<nidof;j++){
      d_dom[j]=beta*d_dom[j]+z_dom[j];
    }
    
    if (buff[maxnbdof+1]==10000.0){
      if (myrank==0){
        // slaves
        for (j=1;j<nproc;j++){
          k=domproc[j];
          nullv (buff,maxnbdof);
          select_master_buff(x_bound,buff,k);
          MPI_Send (buff,buffsize,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD); 
        }
        // master
        k=domproc[0];
        nullv (buff,buffsize-1);
        select_master_buff(x_bound,buff,k);
        
        // deleting
        delete []r_bound;
        delete []d_bound;
        delete []x_bound;
        delete []p_bound;
	delete []z_bound;
      
      }
      else{
        MPI_Recv (buff,buffsize,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
      }
      // x 
      add_bound(buff,x_dom);
      //
      //for (j=0;j<ndof;j++){
      //fprintf (out,"x_dom[%ld] %le\n",j+1,x_dom[j]); 
      //}
      
      //gm->gmxv (x_dom,p_dom);
      
      //fprintf (out,"kontrola\n"); 
      //for (j=0;j<ndof;j++){
      //fprintf (out,"%ld %le %le\n",j+1,p_dom[j],rhs[j]); 
      //}
      
      // results of computation
      copyv (x_dom,lhs,ndof);
      
      // deleting
      delete []r_dom;
      delete []d_dom;
      delete []x_dom;
      delete []p_dom;
      delete []z_dom;
      
      break;
    }
    
  }
  titere = clock();
  MPI_Barrier (MPI_COMM_WORLD);
 
  
  tsole = clock();
  if(myrank == 0){
    printf ("\n\nTime of computation\n\n\n"); 
    switch(prec){
    case noprec:{
      if(mespr == 1) fprintf (out,"Time of iteration %le s\n",(titere-titers)/(double)CLOCKS_PER_SEC); 
      if(mespr == 1) fprintf (out,"Time of solution of system %le s\n",(tsole-tsols)/(double)CLOCKS_PER_SEC); 
      printf ("Time of iteration %le s\n",(titere-titers)/(double)CLOCKS_PER_SEC); 
      printf ("Time of solution of system %le s\n",(tsole-tsols)/(double)CLOCKS_PER_SEC); 
      break;
    }
    case petscilu:{
      if(mespr == 1) fprintf (out,"Time of ILU factorisation by PETSC %le s\n",(tfacte-tfacts)/(double)CLOCKS_PER_SEC); 
      if(mespr == 1) fprintf (out,"Time of iteration %le s\n",(titere-titers)/(double)CLOCKS_PER_SEC); 
      if(mespr == 1) fprintf (out,"Time of solution of system %le s\n",(tsole-tsols)/(double)CLOCKS_PER_SEC); 
      if(mespr == 1) fprintf (out,"Domain 0:\n"); 
      if(mespr == 1) fprintf (out,"Time of all back substitution %le s\n",tbacksub); 
      if(mespr == 1) fprintf (out,"Time of back substitution per iteration %lf s\n",tbacksub/double(iter)); 

      printf ("Time of ILU factorisation by PETSC %le s\n",(tfacte-tfacts)/(double)CLOCKS_PER_SEC); 
      printf ("Time of iteration %le s\n",(titere-titers)/(double)CLOCKS_PER_SEC); 
      printf ("Time of solution of system %le s\n",(tsole-tsols)/(double)CLOCKS_PER_SEC); 
      printf ("Domain 0:\n"); 
      printf ("Time of all back substitution %le s\n",tbacksub); 
      printf ("Time of back substitution per iteration %le s\n",tbacksub/double(iter)); 

      
      
      for (i=1;i<nproc;i++){
	MPI_Recv (&tbacksub,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	j=domproc[stat.MPI_TAG];
	if(mespr == 1) fprintf (out,"Domain %ld:\n",j); 
	if(mespr == 1) fprintf (out,"Time of all back substitution %le s\n",tbacksub); 
	if(mespr == 1) fprintf (out,"Time of back substitution per iteration %lf s\n",tbacksub/double(iter)); 
	printf ("Domain %ld:\n",j); 
	printf ("Time of all back substitution %le s\n",tbacksub); 
	printf ("Time of back substitution per iteration %le s\n",tbacksub/double(iter)); 
      }
      break;
    }
    }
  }
  else{
    MPI_Send(&tbacksub,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  
  delete []buff;
  switch(prec){
  case pardiagprec:{
    delete []precvec;
    break;
  }
  }
}

