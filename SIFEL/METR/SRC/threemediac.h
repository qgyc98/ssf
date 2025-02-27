#ifndef THREEMEDIAC_H
#define THREEMEDIAC_H

#include "genfile.h"

class medc3
{
 public:
  medc3();    //constructor
  ~medc3();   //destructor

  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs_u2 (matrix &d,long ri,long ci,long ipp);
  
 private: 
  
};  

#endif
