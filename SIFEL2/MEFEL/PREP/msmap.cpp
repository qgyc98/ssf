#include <float.h>
#include "msmap.h"


msmap::msmap()
{
  x = y = z = 0.0;
  xi = eta = zeta = 0.0;
  cerr = DBL_MAX;
  eid = -1;
  nid = -1;
  slmas = 0;
  ent = gentity(-1);
}
