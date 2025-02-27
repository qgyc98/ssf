#ifndef ADAPTIVITY_H
#define ADAPTIVITY_H

#include "xfile.h"
#include "vector.h"

class adaptivity
{
 private:
  /// characteristics of the problem
  long dim;    ///  dimension
  long ord;    ///  base functions order
  long ncomp;  ///  number of refined components [stress or strain in mefel]
  long nn;     ///  number of nodes
  long ne;     ///  number of elements
  
  ///  type of error smoothing ; 1 = spr_smoothing , 2 = z2_smoothing
  long tad;
  
  ///  sum of printing flags:
  ///  1  - print stdout     - printing of solving information on standard output
  ///  2  - print bgm        - printing of file.bgm necessary for new mesh generation by T3d generator, file contains required element sizes on points(of old mesh) 
  ///  4  - print
  ///  8  - print
  ///  16 - print test       - printing of file.test for testing
  long printflags;
  
  ///  sum of other flags:
  ///  1 - linear (0) X nonlinear(1) compute_refsizel
  ///  2 - spr compute with strain(0) X stress(1)
  ///  4 - (strain X stress) is own(0) X another(2)
  ///  8 - zienkiewiczuv zpusob nelin pocitani(8) X ostatni(0)
  ///    2 - strain is own(0) X another(1)
  ///    4 - stress is own(0) X another(2)
  long otherflags;
  
  ///  correction
  double corr;
  
  /// accuracy
  double adapt_accuracy;
  
  /// decomposed input file name
  char *path;
  char *filename;
  char *suffix;
  
  ///  second suffix of enter file, it contains actual number of iteration
  char *ni;
  
  /// *** COMPUTED VALUES IN ARRAYS ***
  /// array of refined derivatives (stress or strain) in nodes
  vector *refined_ders_in_nodes;
  /// percentual error at elements
  vector elem_error_pct;
  /// absolute new/refined size of elements
  vector refsizel;
  /// relative new/refined size of elements [only for better visualization]
  vector refsizelrel;
  
  
  ///  answer of this class ; 1(0) == remeshing is (is not) recommended
  long answer;
  
 public:
  long give_adaptflag (void) const { return printflags; }
  
 public:
  /// CONSTRUCTOR
  adaptivity (void);
  /// DESTRUCTOR
  ~adaptivity (void);
  
  /// read values from input file
  void readinit (XFILE *in);
  /// print values to input file
  void printinit (FILE *out);
  
  /// main function
  long run (long of,int width,long nincr);
  
 private:
  /// 
  void prepare_ni (int w, long i);
  /// check consistency
  void check_consistency (void) const;
  
  /// average derivatives to nodes
  void spr (void);
  /// compute error...
  void compute_error (void);
  
  ///
  void compute_refsizel_lin (const double *sizel, const double *ei2, double e2, double u2);
  
  ///
  //void print_addat_vtk () const;
  ///
  void print_test (void) const;
};


#endif
