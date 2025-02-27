#ifndef DOFRANGE_H
#define DOFRANGE_H

#include <utility>

class dofrange{
 public:
  dofrange() {mindof = deltadof = 0L;}
  static bool compare(const std::pair<dofrange,long> &lhs, const std::pair<dofrange,long> &rhs){
    if ((lhs.first.mindof < rhs.first.mindof) ||
        ((lhs.first.mindof == rhs.first.mindof) && (lhs.first.deltadof < rhs.first.deltadof)))
      return true;
    return false;
  }
  long mindof;   // the minimum dof number
  long deltadof; // maximum dof numer - minimum dof number  
};

#endif
