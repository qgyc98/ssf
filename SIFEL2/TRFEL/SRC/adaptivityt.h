#ifndef ADAPTIVITYT_H
#define ADAPTIVITYT_H

#include "xfile.h"
#include "vector.h"
#include "timecontr.h"
#include "gtopology.h"

class adaptivityt
{
 private:
  /// characteristics of the problem
  long dim;    ///  dimension
  long ord;    ///  base functions order
  long ncomp;  ///  number of refined components [stress or strain in mefel]
  long nn;     ///  number of nodes
  long ne;     ///  number of elements
  long ntm;    ///  number of transported matters
  
  
  ///  type of error smoothing ; 1 = spr_smoothing , 2 = z2_smoothing
  long tad;
  
  /// accuracy
  double adapt_accuracy;
  
  /// max enlargement, reduction
  double enlarg;  // > 1.0
  double reduct;  // < 1.0
  
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
  //long otherflags;
  
  ///  correction
  double corr;
  
  /// decomposed input file name
  char *path;
  char *filename;
  char *suffix;
  
  /// internal adaptive loop
  bool interadaloop;
  /// step in adaptive loop
  long step;
  
  ///  second suffix of enter file, it contains actual number of iteration
  int niwidth;
  char *ni;
  
  
  /// *** COMPUTED VALUES IN ARRAYS ***
  /// array of refined derivatives (gradients of fluxes) in nodes, dimension [nn, ncomp]
  vector *refined_ders_in_nodes;
  
  /// percentual error at elements, dimension [ntm, ne]
  vector *elem_error_pct;
  /// relative new/refined size of elements [only for better visualization], dimension [ntm, ne]
  vector *refsizelrel;
  /// absolute new/refined size of elements, dimension [ntm, ne]
  vector *refsizelabs;
  
  
  /// *** BACK UP FOR INTERNAL LOOP ***
  double *r;         /// unknowns in ALL nodes of global problem;  values in array are stored (r[node_1][x] ... r[node_nn][x] , r[node_1][y] ... r[node_nn][y])
  double *rdr;       /// derivatives of unknowns in   --""--    ;   --""--
 public:
  timecontr *tctrl;  /// time controler
 public:
  long istep;        /// istep
  
 public:
  ///  answer of this class ; 1(0) == remeshing is (is not) recommended
  long answer;
  
 private:
  void set_step (long s);
  
  
 public:
  //long give_adaptflag (void) const { return printflags; }
  int give_dim     (void) const { return dim; }
  int give_ntm     (void) const { return ntm; }
  int give_niwidth (void) const { return niwidth; }
  int give_step    (void) const { return step; }
  const char*   give_filename (void) const { return filename; }
  const char*   give_ni       (void) const { return ni; }
  const double* give_r        (void) const { return r; }
  const double* give_rdr      (void) const { return rdr; }
  
 public:
  /// CONSTRUCTOR
  adaptivityt (void);
  /// DESTRUCTOR
  ~adaptivityt (void);
  
  /// read values from input file
  void readinit (XFILE *in);
  /// print values to input file
  void printinit (FILE *out);
  /// intialize atributes
  void initialize (long s);
  
  /// main function
  long run (long of, bool ial);
  
 private:
  /// 
  void prepare_ni (void);
  /// check consistency
  void check_consistency (void) const;
  
  /// average derivatives to nodes
  void spr (int mattid);
  /// compute error...
  void compute_error (int mattid);
  
  ///
  void compute_refsizel_lin (int mattid, const double *sizel, const double *ei2, double e2, double u2);
  
  ///
  void print_addat_vtk () const;
  
  
  // ***   STATE DATA   ***
 public:
  void statedata_backup (void);
  void statedata_transfer (adaptivityt *Adat_old, gtopology *Gtt_old);
  void statedata_restore (void);
};

#endif
