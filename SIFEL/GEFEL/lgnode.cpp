#include "lgnode.h"

lgnode::lgnode (void)
{
  nl=0;
  nodes=NULL;  cn=NULL;
}

lgnode::~lgnode (void)
{
  delete [] nodes;  delete [] cn;
}

