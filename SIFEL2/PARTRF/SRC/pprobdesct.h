#ifndef PPROBDESCT_H
#define PPROBDESCT_H

#include "psolver.h"
#include "pgalias.h"
#include "probdesct.h"
#include "aliast.h"
#include "genfile.h"
#include <iotools.h>
#include <stdio.h>

class pprobdesct
{
 public:
  pprobdesct (void);
  ~pprobdesct (void);
  void shift_indices ();

  //  first node index on subdomain (in global ordering; used for GiD)
  long fni;
  //  first element index on subdomain (in global ordering; used for GiD)
  long fei;
  
};

#endif
