#ifndef GEOELASTMAT_H
#define GEOELASTMAT_H

#include "iotools.h"
#include "alias.h"
struct matrix;
struct vector;



/**
   Class geoelastmat defines elastic material model for geotechnical
   purposes. This material has different elastic modulus for loading 
   and unloading. Loading is defined as decreasing of mean stress 
   (i.e. increasing pressure) and unloading is defined as increasing
   of mean stress (i.e. decreasing pressure).
   
   Created by Tomas Koudelka,
*/
class geoelastmat
{
 public:
  geoelastmat (void);
  ~geoelastmat (void);
  void read (XFILE *in);

  void matstiff (matrix &d, long ipp, long ido);
  void nlstresses (long ipp, long im, long ido);
  void updateval (long ipp, long im, long ido);
  
  /// unloading modulus coefficient 
  double keu;
};

#endif
