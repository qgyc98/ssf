#ifndef WINPAST_H
#define WINPAST_H

#include "iotools.h"
#include "alias.h"
#include "matrix.h"
#include "vector.h"

class winpast
{
 public:
  winpast (void);
  ~winpast (void);
  void read (XFILE *in);
  void print (FILE *out);

  void matstiff (matrix &d,long ipp);
  void elmatstiff (matrix &d,long ipp);
  void matstiff_soilplbeam (matrix &d);
  void matstiff_soilbeam (matrix &d);
  void matstiff_soilplate (matrix &d);
  void matstiff_soilspacestr (matrix &d);
  
  void nlstresses (long ipp);
  
  //  material parameters
  vector c1,c2;
  
};

#endif
