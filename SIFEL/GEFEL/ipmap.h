#ifndef IPMAP_H
#define IPMAP_H

class ipmap
{
 public:
  double x, y, z;  // global coordinates of integration point
  double xi, eta, zeta;  // natural coordinates of integration point on element eid
  double cerr;  // error of computed natural coordinates if they are not in range [-1.0;1.0] exactly
  long eid;     // element id on which the point is located
  long ipp;     // integration point pointer 
                // ipp<0, the the mapping is given by xi, eta, zeta
                // ipp >= 0, the direct mapping between integration points
  long app;    // auxiliary point pointer
               // app < 0, the mapping is not set or ipp is given
               // app >= 0, index of auxiliary integration point in the 
               //           array of auxiliary integration points

  ipmap();
  ~ipmap(){};
};

#endif
