#ifndef CSV_H
#define CSV_H

/**
   class deals with compressed sparse %vector format
   
   only nonzero %vector components are stored in the array a
   corresponding indices are stored in the array ind
   
   JK, 25.2.2007
*/

class compvect
{
 public:
  compvect ();
  compvect (long *ii,long nonz);
  ~compvect ();
  void copy (compvect *bc);
  void setup_vector (double *b,long nc);
  void test_ordering ();
  void axpy (compvect *x,double aa);
  void axpy_known (compvect *x,double aa);

  void reorder_components ();
  double dotprod (compvect *cv);
  

  
  
  ///  number of all components (including zero components) of the stored %vector
  long n;
  ///  number of nonzero components of the stored %vector
  long nz;
  ///  array containing nonzero components
  double *a;
  ///  array containing indices of nonzero components
  long *ind;
  
  ///  indicator of reordering
  long ordering;
  ///  limit for components extraction
  double limit;
};



#endif
