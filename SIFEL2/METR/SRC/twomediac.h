#ifndef TWOMEDIAC_H
#define TWOMEDIAC_H

#include "genfile.h"

class medc2
{
 public:
  medc2();    //constructor
  ~medc2();   //destructor

  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs_u2 (matrix &d,long ri,long ci,long ipp);
  
 private: 
  
};  

#endif
