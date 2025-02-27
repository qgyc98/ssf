#ifndef MTGLVEC_H
#define MTGLVEC_H

struct mt_glob_vec
{
  /// step id
  long istep;
  /// number of changed DOFs
  long ncd;
  /// number of changed elements
  long nce;

  // Following pointers are not allocated, 
  // they are used as references to arrays allocated somewhere else

  ///  vector of nodal displacements
  double *r;
  ///  vector of nodal prescribed forces
  double *f;
  
  // Following pointers to arrays are allocated 
  // allocation of vectors  

  ///  %vector of increments of nodal displacements
  double *dr;
  ///  %vector of prescribed forces from the previous step
  double *fp;
  ///  %vector of prescribed force loads, it does not contain forces caused by temperature, etc.
  double *fl;
  ///  vector of internal forces
  double *fi;
  ///  auxiliary force %vector
  double *fb;
  ///  backup of the nodal displacements
  double *lhsb;
  /// %vector of prescribed forces from the previous step due to prescribed force load
  double  *flp;
  /// %vector of indicators for nodes on the interface between old and new parts of the structure (used in the growing mechanical problem)
  long *ifn;

  /// Default constructor 
  mt_glob_vec();
  /// function allocates vectors dr, fp, fl, flp, fi, fb, lhsb and ifn
  void alloc(long n);
  /// function deallocates vectors dr, fp, fl, flp, fi, fb, lhsb and ifn
  void dealloc();
  /// Destructor deallocates vectors dr, fp, fl, flp, fi, fb, lhsb and ifn
  ~mt_glob_vec();
};

#endif
