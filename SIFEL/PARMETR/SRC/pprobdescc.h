#ifndef PPROBDESCC_H
#define PPROBDESCC_H

#include <stdio.h>
#include "alias.h"
#include "aliast.h"
#include "aliasc.h"
#include "probdesc.h"
#include "probdesct.h"
#include "probdescc.h"
//#include "psolver.h"
#include "xfile.h"

class pprobdescc
{
 public:
  pprobdescc (void);
  ~pprobdescc (void);
  void shift_indices ();

  ///  first node index on subdomain (in global ordering; used for GiD)
  long fnim;
  ///  first element index on subdomain (in global ordering; used for GiD)
  long feim;

  ///  first node index on subdomain (in global ordering; used for GiD)
  long fnit;
  ///  first element index on subdomain (in global ordering; used for GiD)
  long feit;

};

#endif
