#ifndef BOUNDFLUXES_H
#define BOUNDFLUXES_H

#include <stdio.h>
#include "iotools.h"
#include "loadelt.h"
#include "elemnode.h"

class boundfluxes
{
 public:
  boundfluxes (void);
  ~boundfluxes (void);
  void read (XFILE *in,long lcid);
  void print (FILE *out,long lcid);
  
  // *************************************************
  //  DESCRIPTION OF BOUNDARY CONDITIONS ON ELEMENTS
  // *************************************************
  ///  number of elements with boundary conditions
  long neb;
  ///  load case id; it is the number of variable which has to be integrated
  long lcid;

  loadelt *elemload;
  
};

#endif
