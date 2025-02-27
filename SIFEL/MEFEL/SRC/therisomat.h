#ifndef THERISOMAT_H
#define THERISOMAT_H

#include "iotools.h"
#include "alias.h"
struct matrix;


/**
  Class defines thermal isotropic dilatancy material
*/
class therisomat
{
 public:
  therisomat (void);
  ~therisomat (void);

  /// reads material parameter from file
  void read (XFILE *in);
  void print (FILE *out);
  
  /// assembles thermal dilatancy %matrix 
  void matdilat (matrix &d,strastrestate ssst);
  /// assembles thermal dilatancy %matrix for 1D stress/strain state
  void matdilat_bar (matrix &d);
  /// assembles thermal dilatancy %matrix for 2D beam stress/strain state
  void matdilat_plbeam (matrix &d);
  /// assembles thermal dilatancy %matrix 2D plainstress state
  void matdilat_plstress (matrix &d);
  /// assembles thermal dilatancy %matrix 2D plainstrain state
  void matdilat_plstrain (matrix &d);
  /// assembles thermal dilatancy %matrix axisymmetrical problems
  void matdilat_axi (matrix &d);
  /// assembles thermal dilatancy %matrix for plates
  void matdilat_plate (matrix &d);
  /// assembles thermal dilatancy %matrix for 3D stress/strain state
  void matdilat_spacestr (matrix &d);

  /// computes thermal strains
  void temprstrains (long ipp);

  /// computes thermal strains and accumulates them in Mm->tempstrains
  void cumultemprstrains (long ipp);

  /// marks required non-mechanical quantities
  void give_reqnmq(long *anmq);
  ///  coefficient of thermal dilatancy
  double alpha;
};

#endif
