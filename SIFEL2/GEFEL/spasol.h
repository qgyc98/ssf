/**
   class defining Intel sparse direct solver
   
   IS, 7.12.2005
*/
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif

extern "C" { int PARDISO (void *, int *, int *, int *, int *, int *,
		    double *, int *, int *, int *, int *, int *,
		    int *, double *, double *, int *);
  int pardisoinit_(void *, int *, int *);
		    
}

/*
#ifdef AIX
#define F77_FUNC(func)  func     
#else
#define F77_FUNC(func)  func ## _
#endif

// PARDISO prototype. 
extern  int F77_FUNC(pardisoinit)
    (void *, int *, int *);
*/

class spasol  
{
 public:
  
  spasol ();
  ~spasol ();

  void symbfact (double *a,long *ci,long *adr,long ndof);
  void numfact (double *a,long *ci,long *adr,long ndof);
  void backsubst (double *a,long *ci,long *adr,long ndof,double *x,double *y);
  
  //spasol(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*,int*, int*, int*, double*, double*, int*);

  //void spasol(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*,int*, int*, int*, double*, double*, int*);
  
  ///  number of rows
  int nRows;
  ///  number of columns
  int nCols;
  ///  number of nonzeros
  int nNonZeros;
  ///  number of right hand sides
  int nrhs;
  /* Number of processors. */
  int      num_procs;

  /* Auxiliary variables. */
  // char    *var;

  //double *rhs,*solValues;

  ///  type of solved matrix
  ///  mtype=11 - real unsymmetric matrix
  int mtype;
  
  ///  internal solver memory pointer pt
  ///  32-bit: int pt[64]; 64-bit: long int pt[64]
  ///  or void *pt[64] should be OK on both architectures
  void *pt[64];

  /// spasol control parameters
  int iparm[64];
  int maxfct, mnum, phase, error, msglvl;
  // /*  auxiliary variables. */
  double ddum; /* Double dummy */
  int idum; /* Integer dummy. */
  

};

