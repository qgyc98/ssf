#include <float.h>
#include "hngen.h"
#include "vector.h"


hngen::hngen()
{
  x = y = z = 0.0;
  xi = eta = zeta = 0.0;
  cerr = DBL_MAX;
  eid = -1;
  et = gtypel(-1);
}
