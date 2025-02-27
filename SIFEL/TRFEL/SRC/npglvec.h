#ifndef NPGLVEC_H
#define NPGLVEC_H

struct np_glob_vec
{
  /// step id
  long istep;
  /// number of changed DOFs
  long ncd;
  /// number of changed elements
  long nce;

  // Following pointers are not allocated, 
  // they are used as references to arrays allocated somewhere else

  ///  nodal values
  double *lhs;
  ///  time derivatives of nodal values
  double *tdlhs;
  ///  right hand side
  double *rhs;

  // Following pointers to arrays are allocated 
  // allocation of vectors  

  ///  vector of prescribed fluxes (right hand side)
  double *f;
  ///  predictor
  double *d;
  ///  auxiliary vector
  double *p;
  ///  auxiliary vectors for nonlinear problems 
  double *fb,*fi;

  ///  auxiliary vector for nonlinear problems in dform
  double *v;
  ///  auxiliary vector for nonlinear problems in dform
  double *z;
  ///  backup of nodal values
  double *lhsb;
  ///  backup of time derivatives of nodal values
  double *tdlhsb;
  /// %vector of prescribed forces from the previous step due to prescribed force load - upravit??!!
  //double  *flp;
  /// %vector of indicators for nodes on the interface between old and new parts of the structure (used in the growing mechanical problem) - upravit??!!
  long *ifn;


  /// Default constructor 
  np_glob_vec();
  /// function allocates vectors lhs, tdlhs, rhs, f, d, p
  void alloc(long n);
  /// function allocates vectors fb, fi
  void alloc_aux(long n);
  /// function allocates vectors v, z, lhsb and tdlhsb for dform solver
  void alloc_daux(long n);
  /// function deallocates vectors lhs, tdlhs, rhs, f, d, p
  void dealloc();
  /// Destructor deallocates vectors lhs, tdlhs, rhs, f, d, p
  ~np_glob_vec();
};

#endif
