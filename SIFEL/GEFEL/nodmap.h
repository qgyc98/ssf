#ifndef NODMAP_H
#define NODMAP_H

class nodmap
{
 public:
  long nid; // node id
  double x, y, z;  // global coordinates of node
  double xi, eta, zeta;  // natural coordinates of integration point on element eid
  double cerr;  // error of computed natural coordinates if they are not in range [-1.0;1.0] exactly
  long eid;     // element id on which the node is located

  nodmap();
  ~nodmap(){};
};

#endif
