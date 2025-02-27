#ifndef HNGEN_H
#define HNGEN_H

#include "galias.h"
#include "vector.h"

class hngen
{
 public:
  double x, y, z;       // global coordinates of hanging node
  double xi, eta, zeta; // natural coordinates of hanging node on element eid
  double cerr; // error of computed natural coordinates if they are not in range [-1.0;1.0] exactly
  long eid;    // master element id on which the hanging node is located eid >=0,
  gtypel et; // general element type used for the hanging node assignment
  ivector mnodes; /// array of master nodes
  
  hngen();
};

#endif
