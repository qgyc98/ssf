#ifndef DLOADEL_H
#define DLOADEL_H

#include "xfile.h"
#include "alias.h"
#include <stdio.h>
#include "vector.h"

class gfunct;

/**
  Class deals with time dependent (dynamic) loads defined on elements.
  It contains edge load, surface load and volume load whose particular 
  components are given by time functions independetly.

  Created by Tomas Koudelka, 26.4.2016
*/
class dloadel
{
 public:
  dloadel();
  ~dloadel();
  /// reads element load from file
  long read(XFILE *in);

  /// reads element volume load from file
  void readvolumeload (XFILE *in);

  /// reads element edge load from file
  void readedgeload (XFILE *in);

  /// reads element surface load from file
  void readsurfaceload (XFILE *in);

  /// reads element load from the preprocessor file
  long read_prep(XFILE *in, long lc);

  /// prints element load to the file
  long print(FILE *out, int ident);

  /// computes nodal values of element load for the actual time
  void compute_load(double t);
  /// computes nodal values from volume load
  void volumeload (double t);

  /// computes nodal values from edge load
  void edgeload (double t);

  /// computes nodal values from surface load
  void surfaceload (double t);

  /// merges element load defined by lel
  long merge(dloadel &lel);
  
  /// sets load type depending on allocated arrays nodval{e|s|v}
  void set_load_type();
  
  ///  the number of loaded element (element id)
  long eid;
  ///  the number of DOFs on element
  long ndofe;
  ///  the type of element load
  elloadtype tel;
  ///  the number of edges
  long ned;
  ///  the number of nodes on one edge
  long nned;
  ///  the number of surfaces
  long nsurf;
  ///  the number of nodes on one surface
  long nnsurf;
  ///  the number of approximated functions on element
  long napfun;
  ///  the number of elements of nodvale array
  long nnve;
  ///  the number of elements of nodvals array
  long nnvs;
  ///  the number of elements of nodvalv array
  long nnvv;
  ///  components of edge load
  gfunct *nodvale;
  ///  components of surface load
  gfunct *nodvals;
  ///  components of volume load
  gfunct *nodvalv;
  ///  components of load %vector - nodal forces
  vector nf;
  ///  indicators of loaded edges
  long *le;
  ///  indicators of loaded surfaces
  long *ls;
  ///  load case id
  long nlc;
  ///  subload case id
  long nslc;
};

#endif
