#ifndef LOADEL_H
#define LOADEL_H

#include <stdio.h>
#include "iotools.h"
#include "alias.h"
#include "vector.h"

/**
  Class deals with loads defined on elements.
  It contains edge load, surface load and volume load.

  Created by JK,
  Modified by Tomas Koudelka, 06.2009
*/
class loadel
{
 public:
  loadel();
  ~loadel();
  /// reads element load from file
  long read(XFILE *in);

  /// reads element volume load from file
  void readvolumeload (XFILE *in);

  /// reads element edge load from file
  void readedgeload (XFILE *in);

  /// reads element surface load from file
  void readsurfaceload (XFILE *in);

  /// reads load types specific to beam elements
  void readbeamload(XFILE *in);

  /// reads element load from the preprocessor file
  long read_prep(XFILE *in, long lc, long *slc);

  /// prints element load to the file
  long print(FILE *out, int ident);

  /// computes nodal values due to volume load
  void volumeload ();

  /// computes nodal values due to edge load
  void edgeload ();

  /// computes nodal values due to surface load
  void surfaceload ();

  /// computes nodal values due to beam load
  void beamload();

  /// merges element load defined by lel
  long merge(loadel &lel);
  
  /// sets load type depending on allocated arrays nodval{e|s|v}
  void set_load_type();
  
  ///  the number of loaded element (element id)
  long eid;
  ///  the number of DOFs on element
  long ndofe;
  ///  the type of element load
  elloadtype tel;
  /// meaning of element load values
  elloadmeaning elm;
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
  double *nodvale;
  ///  components of surface load
  double *nodvals;
  ///  components of volume load
  double *nodvalv;
  ///  offset of beam continuous load from the first node
  double *la;
  ///  length of applied continuous load on beam element
  double *lf;
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
  /// number of point load on beam element
  long npnt;
};

#endif
