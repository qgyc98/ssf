#include "pardiso.h"
#include <stdio.h>

//extern int omp_get_max_threads();
/* PARDISO prototype. */
//#define PARDISO pardiso
extern "C" { int PARDISO (void *, int *, int *, int *, int *, int *,
		    double *, int *, int *, int *, int *, int *,
		    int *, double *, double *, int *);
}


pardiso::pardiso ()
{
  nRows=0;  nCols=0;
  nNonZeros=0;
  nrhs=1;
  
  mtype=11;

}

pardiso::~pardiso ()
{
  
}

/**
   function computes symbolic factorization
*/
void pardiso::symbfact (double *a,long *ci,long*adr,long ndof)
{
  long i;
  
  nRows=ndof;
  nCols=ndof;
  nNonZeros=adr[ndof];
  
  for (i = 0; i < 64; i++) {
    iparm[i] = 0;
  }
  
  iparm[0] = 1; /* No solver default */
  iparm[1] = 2; /* Fill-in reordering from METIS */
  /* Numbers of processors, value of OMP_NUM_THREADS */
  //iparm[2] = omp_get_max_threads();
  iparm[2] = 1;
  iparm[3] = 0; /* No iterative-direct algorithm */
  iparm[4] = 0; /* No user fill-in reducing permutation */
  iparm[5] = 0; /* Write solution into x */
  iparm[6] = 16; /* Default logical fortran unit number for output */
  iparm[7] = 2; /* Max numbers of iterative refinement steps */
  iparm[8] = 0; /* Not in use */
  iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0; /* Not in use */
  iparm[12] = 0; /* Not in use */
  iparm[13] = 0; /* Output: Number of perturbed pivots */
  iparm[14] = 0; /* Not in use */
  iparm[15] = 0; /* Not in use */
  iparm[16] = 0; /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0; /* Output: Numbers of CG Iterations */
  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1; /* Which factorization to use. */
  msglvl = 0; /* Don't print statistical information in file */
  error = 0; /* Initialize error flag */
  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for (i = 0; i < 64; i++) {
    pt[i] = 0;
  }
  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  
  phase = 11;
  iparm[3] = 31;
  

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
void pardiso::numfact (double *a,long *ci,long*adr,long ndof)
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
void pardiso::backsubst (double *a,long *ci,long*adr,long ndof,double *x,double *y)
{
  phase = 33;
  iparm[7] = 2; /* Max numbers of iterative refinement steps. */
  
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
	   (int *)&ndof, a, (int *)adr, (int *)ci, &idum, &nrhs,
	   iparm, &msglvl, x, y, &error);


  if (error != 0) {
    printf("\nERROR during solution: %d", error);
    //exit(3);
  }
  printf("\nSolve completed ... ");
  printf("\nThe solution of the system is: ");
  fflush(stdout);    
  
}

