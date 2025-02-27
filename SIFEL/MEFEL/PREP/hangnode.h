#ifndef HANGNODE_H
#define HANGNODE_H

#include <stdio.h>
#include "galias.h"
#include "iotools.h"

struct ivector;
struct vector;
class hngen;

/**
   
 This class contains additional information about hanging node
 such as master nodes indeces and natural coordinates.
   
 Created by TKo, 10.2012
*/

class hangnode
{
 public :
  hangnode();
  hangnode(ivector &mn, vector &natc, gtypel ment);
  ~hangnode();
  long read(XFILE *in);
  long print(FILE *out);
  long compare(const hngen &hng);

  /// number of master nodes
  long nmn;
  ///  array of the master nodes (if the node is hanging node)
  long *mnodes;
  ///  natural coordinates of the hanging node
  double *natcoord;
  ///  type of master element type
  ///  it describes to which general element (linbar, lintriangle, etc.) is the hanging node attached
  gtypel maset;
};


#endif
