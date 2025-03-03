#ifndef CPROF_GEOM_H
#define CPROF_GEOM_H

#include <stdio.h>

#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif 

struct XFILE;

class cprof_geom
{
 public:
  cprof_geom();
  ~cprof_geom();
  int  read(XFILE *in);
  int  get_nn();
  int  get_ne();
  void gen_nodes(FILE *out, double tx, double ty, double alpha, int inid);
  void gen_elements(FILE *out, int inid, int ieid);

  double wp; // width of profile (width of flange + tw)
  double hp; // height of profile (height of web + 2xtf)

  double tf; // thickness of flange
  double tw; // thickness of web

  int mshdf; // mesh density of the flange
  int mshdw; // mesh density of the web
  int mshdtf; // mesh density of flange across thickness
  int mshdtw; // mesh density of web across thickness
  int matid;  // index of material
};

#endif
