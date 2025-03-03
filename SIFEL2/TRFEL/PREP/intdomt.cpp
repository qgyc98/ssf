#include <stdlib.h>
#include "intdomt.h"

intdomt::intdomt()
{
  idid = NULL;
  n = 0;
}


intdomt::~intdomt()
{
  delete [] idid;
}
