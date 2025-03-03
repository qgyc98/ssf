#ifndef PPROBDESC_H
#define PPROBDESC_H

#include "../../PARGEF/psolver.h"
#include "alias.h"
#include "probdesc.h"
#include "pgalias.h"
#include "../../GEFEL/xfile.h"
#include <stdio.h>

class pprobdesc
{
 public:
  pprobdesc (void);
  ~pprobdesc (void);
  void shift_indices ();
  
  //  first node index on subdomain (in global ordering; used for GiD)
  long fni;
  //  first element index on subdomain (in global ordering; used for GiD)
  long fei;

};

#endif
