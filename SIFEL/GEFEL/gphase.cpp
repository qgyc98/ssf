#include "gphase.h"

gphase::gphase (void)
{
  f=0.0;
}

gphase::~gphase (void)
{
}

/**
   functions reads basic informations
   
   @param in - input stream
   
   TKr
*/
void gphase::read (FILE *in)
{
  //  fraction
  fscanf (in,"%lf",&f);
}

