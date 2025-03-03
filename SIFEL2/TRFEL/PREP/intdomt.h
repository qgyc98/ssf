#ifndef INTDOMT_H
#define INTDOMT_H

/**
  The class is intended for the description of integration domains (surfaces/edges) for integration of flux 
  resultant on the given element. The length of array idid n is given by the number of boundary objects (surfaces/edges)
  on the given element.
  sid[i] = index of integration domain to which the given element contribution will be added.

  Created by Tomas Koudelka, 02.2018
*/
class intdomt
{
 public:
  intdomt();
  ~intdomt();

  /// array of integration domain indeces
  long *idid;
  /// length of array idid = the number of boundary objects on the given element
  long n;
};

#endif
