#include "spasol.h"
#include <stdio.h>
#include <stdlib.h>

//extern int omp_get_max_threads();
/* PARDISO prototype. */
//#define PARDISO spasol


spasol::spasol ()
{
  char *var;
  
  nRows=0;  
  nCols=0;
  nNonZeros=0;
  nrhs=1;
  
  mtype=11;
  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  /*for (i = 0; i < 64; i++) {
    pt[i] = 0;
  }*/
  
  pardisoinit_(pt,  &mtype, iparm);

    /* Numbers of processors, value of OMP_NUM_THREADS */
  num_procs=1;
  var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
        sscanf( var, "%d", &num_procs );
  else {
        printf("Assume that environment variable OMP_NUM_THREADS is set to 1");
        //exit(1);
  }
  iparm[2]  = num_procs;
  
  //mtype=-2;

}

spasol::~spasol ()
{
  
}

/**
   function computes symbolic factorization
*/
void spasol::symbfact (double *a,long *ci,long *adr,long ndof)
{
  long i;

  
  nRows=ndof;
  nCols=ndof;
  nNonZeros=adr[ndof];
  
  
  
  for(i=0;i<=ndof;i++) adr[i]++;
  for(i=0;i<nNonZeros;i++) ci[i]++;
  //for(i=0;i<ndof;i++) rhs[i]=1.0;
  //for(i=0;i<ndof;i++) solValues[i]=1.0;  
  
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 1; /* Don't print statistical information in file */
  error = 0; /* Initialize error flag */
  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  
  phase = 11;
  //iparm[3] = 0;
//  printf("JSEM TU 1b\n");
//  fflush(stdout);
  

  // *************************
  //  symbolic factorization
  // *************************
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   (int *)&ndof, a, (int *)adr, (int *)ci, &idum, &nrhs,
	   iparm, &msglvl, &ddum, &ddum, &error);
  

  if (error != 0) {
    fprintf(stderr,"\nERROR during symbolic factorization: %d\n.", error);
    //exit(1);
  }
  
  fprintf(stdout,"\nReordering completed ... ");
  fprintf(stdout,"\nPeak memory in KB = %d", iparm[14]);
  fprintf(stdout,"\nPermament memory in KB = %d", iparm[15]);
  fprintf(stdout,"\nMemory for solve in KB = %d", iparm[16]);        
  fprintf(stdout,"\nNumber of nonzeros in factors = %d", iparm[17]);
  fprintf(stdout,"\nNumber of factorization MFLOPS = %d", iparm[18]);
  
  fflush(stdout);
}

/**
   function computes numerical factorization
*/
void spasol::numfact (double *a,long *ci,long*adr,long ndof)
{
  phase = 22;
  
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   (int *)&ndof, a, (int *)adr, (int *)ci, &idum, &nrhs,
	   iparm, &msglvl, &ddum, &ddum, &error);
  
  if (error != 0) {
    printf("\nERROR during numerical factorization: %d", error);
    //exit(2);
    printf("\nFactorization completed ... ");
    fflush(stdout);    
  }
  
}

/**
   function computes back substitution
*/
void spasol::backsubst (double *a,long *ci,long*adr,long ndof,double *x,double *y)
{
  int i;
  phase = 33;
  iparm[7] = 1; /* Max numbers of iterative refinement steps. */
  
  /*
  if (x==NULL)
  {
    printf("\n X je NULL\n");
    fflush(stdout);
  }
  if (y==NULL)
  {
    printf("\n Y je NULL\n");
    fflush(stdout);
  }*/
  /*
  for(i=0;i<100;i++)
  {
    printf(" X=%g    Y=%g\n",x[i],y[i]);
  }
  for(i=ndof-100;i<ndof;i++)
  {
    printf(" X=%g    Y=%g\n",x[i],y[i]);
  }*/

  /*
  for(i=0;i<ndof;i++)
  {
    if (x[i]!=0.0) printf(" X=%g    Y=%g\n",x[i],y[i]);
  }*/
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   (int *)&ndof, a, (int *)adr, (int *)ci, &idum, &nrhs,
	   iparm, &msglvl, y, x, &error);


  if (error != 0) {
    printf("\nERROR during backward substitution: %d", error);
    //exit(3);
  }
  printf("\nSolve completed ... ");
  printf("\nThe solution of the system is: ");
  /*
  for(i=0;i<100;i++)
  {
    printf(" X=%g    Y=%g\n",x[i],y[i]);
  }
  for(i=ndof-100;i<ndof;i++)
  {
    printf(" X=%g    Y=%g\n",x[i],y[i]);
  }
  
  fflush(stdout);    */
  
}

