#include "pointset.h"

/**
  This constructor initializes class attributes to zero values
*/
pointset::pointset()
{
  n      = -1L;
  ncoord =  0L;
  coord  = NULL;
}

/**
  This destructor deallocates used memory
*/
pointset::~pointset()
{
  delete [] coord;
}


long pointset::read(FILE */*in*/)
{
  return(0);
}

long pointset::print(FILE */*out*/)
{
  return(0);
}
