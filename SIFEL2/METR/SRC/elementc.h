#ifndef ELEMENTC_H
#define ELEMENTC_H

#include <stdio.h>
#include "aliasc.h"
#include "iotools.h"

/**
   class elementc defines general element for problems of mechanical-transport coupling
*/
class elementc
{
 public:
  elementc (void);
  ~elementc (void);
  void read (XFILE *in,long eid);
  void readmat (long m,XFILE *in);
  void readmatfc (XFILE *in);  
  
  ///  type of element
  elemtypec te;
  ///  number of blocks
  long nb;
  ///  array of integration point pointers (upper part)
  long **ippu;
  ///  array of integration point pointers (lower part)
  long **ippl;
  ///  array of numbers of integration points (upper part)
  long **nipu;
  ///  array of numbers of integration points (lower part)
  long **nipl;
  ///  array of numbers of integration points for fully coupled models
  long **nipc;
  ///  array of integration orders for coupling matrix
  long **intordvum;
  ///  array of integration orders for coupling matrix
  long **intordvlm;

  ///  type of material
  mattypec *tmu;
  ///  number of appropriate material type
  long *idmu;
  ///  type of material
  mattypec *tml;
  ///  number of appropriate material type
  long *idml;
  ///  type of fully coupled material
  mattypec tm;
  ///  number of appropriate fully coupled material type
  long idm;
  
  
  long ipp;
};

#endif
