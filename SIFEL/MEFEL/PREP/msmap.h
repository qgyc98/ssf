#ifndef MSMAP_H
#define MSMAP_H

#include "galias.h"

class msmap
{
 public:
  double x, y, z;       // global coordinates of slave point
  double xi, eta, zeta; // natural coordinates of slave point on element eid
  double cerr; // error of computed natural coordinates if they are not in range [-1.0;1.0] exactly
  long eid;    // master element id on which the slave point is located eid >=0,
               // eid < 0 if the mapping is given by direct master node mapping or not yet specified
  long nid;    // master node id 
               // nid < 0, the the mapping is given by xi, eta, zeta or not yet specified
               // nid >= 0, the direct mapping between points
  long slmas;  // slave/master indicator
               // slmas == -1, master point is set
               // slmas == 1, slave point is set
               // slmas == 0, no mapping is set
               // slmas == 2 if slave point with several master assignments was detected or
               //            if one point is commanded to be both a slave and a master point concurrently
  gentity ent; // entity used for the master/slave point assignment

  msmap();
};

#endif
