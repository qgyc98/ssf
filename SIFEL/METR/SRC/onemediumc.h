#ifndef ONEMEDIUMC_H
#define ONEMEDIUMC_H

#include "genfile.h"

class medc1
{
 public:
  medc1();    //constructor
  ~medc1();   //destructor

  void matcond_u (matrix &d,long ri,long ci,long ipp);
  void matcap_u (matrix &d,long ri,long ci,long ipp);
  void matcond_l (matrix &d,long ri,long ci,long ipp);
  void matcap_l (matrix &d,long ri,long ci,long ipp);
  void rhs_u1 (matrix &d,long ri,long ci,long ipp);
  void rhs_u2 (matrix &d,long ri,long ci,long ipp);

 private: 
  
  double scale_t,scale_u;
  
};  

#endif
