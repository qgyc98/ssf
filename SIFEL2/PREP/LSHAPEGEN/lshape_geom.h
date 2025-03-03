#ifndef LSHAPE_GEOM_H
#define LSHAPE_GEOM_H

#include <stdio.h>

enum tmat {nomat=0, sdk=1, sheet=2, insulation=3, steel=4};

class cprof_geom;
class gtopology;
class nodmap;
struct XFILE;

class lshape_geom
{
 public:
  lshape_geom();
  ~lshape_geom();
  int read(XFILE *in);
  int gen_nodes(FILE *out);
  int gen_elements(FILE *out);
  int cprof_hang_nodes_natcoords(gtopology *gt, long wne, long nnmap, nodmap *cpnmap, long ni, double err);
  int print_hang_nodes(FILE *out, gtopology *gt, long nnmap, nodmap *cpnmap);

  double lh; // length of horizontal part of wall (from inner corner)
  double lv; // length of vertical part of wall (from inner corner)

  double th; // thickness of horizontal part of wall (sum of thicknesses of all horizontal part layers)
  double tv; // thickness of vertical part of wall (sum of thicknesses of all vertical part layers)

  int numlh;   // the number of layers in horizontal part of wall
  int numlv;   // the number of layers in vertical part of wall

  int mshdh; // mesh density of the horizontal wall part in horizontal direction (common for all layers)
  int mshdv; // mesh density of the vertical wall part in vertical direction (common for all layers)

  int *mshdhl; // mesh density of layers in the horizontal wall part
  int *mshdvl; // mesh density of layers in the vertical wall part

  double *tlh;    // array of thicknesses of horizontal parts of wall (from outer surface)
  double *tlv;    // array of thicknesses of vertical parts of wall (from outer surface)

  tmat   *tmath;  // array of material types for particular layers of horizontal wall part
  int    *matidh; // array of material id for particular layers of horizontal wall part

  tmat   *tmatv;  // array of material types for particular layers of vertical wall part
  int    *matidv; // array of material id for particular layers of vertical wall part

  int numprof;    // number of steel C-profiles
  double *px;     // array of x-coordinates of left bottom corner for particular C-profiles
  double *py;     // array of y-coordinates of left bottom corner for particular C-profiles
  double *alpha;  // array of angle of rotation for particular C-profiles
  cprof_geom *cprof; // array with the description for particular C-profiles geometry

  int tnn;      // the total number of nodes defined in the problem
  int nn_wall;  // the number of nodes used for the layered L-shape wall domain
  int nn_cprof; // the number of nodes used for all C-profiles in the problem
  
  int tne;      // the total number of elements defined in the problem
  int ne_wall;  // the number of elements used for the layered L-shape wall domain
  int ne_cprof; // the number of elements used for all C-profiles in the problem
};

#endif
