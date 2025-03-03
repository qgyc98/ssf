#ifndef HANGNODET_H
#define HANGNODET_H

#include <stdio.h>
#include "galias.h"
#include "iotools.h"

/**
   
 This class contains additional information about hanging node
 such as master nodes indeces and natural coordinates.
   
 Created by TKo, 07.2018
*/

class hangnodet
{
 public :
  hangnodet();
  ~hangnodet();
  long read(XFILE *in);
  long print(FILE *out);

  /// number of master nodes
  long nmn;
  ///  array of the master nodes (if the node is hanging node)
  long *mnodes;
  ///  natural coordinates of the hanging node
  double *natcoord;
  ///  type of master entity
  ///  it describes to which quantity (edge, surface, etc.) is the hanging node attached
  gtypel masentity;
};


#endif
