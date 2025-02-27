#ifndef DLOADPD_H
#define DLOADPD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iotools.h"

#include "gfunct.h"



/**
  This class defines time dependent prescribed displacements.
*/
class dloadpd
{
  public:
   gfunct gf;

   dloadpd();
   ~dloadpd();
   long read(XFILE *in);
   long print(FILE *out);
   double getval(double t);
};

#endif
