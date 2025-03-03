#ifndef HYDRATIONHEAT_H
#define HYDRATIONHEAT_H

#include <stdio.h>
#include "aliast.h"
#include "genfile.h"

/**
   class hydrationheat describes source of heat due to hydration process in concrete
   
   amount of heat generated by concrete is described by function
   h(t)=r*a((b(t/c)^d)/(e+b(t/c)^d))^f
   where t is in hours
   
   heat source is time derivative of the function h(t)
   time is converted to seconds, because seconds are usually used
   it has the form
   z(t)=aa * (bb*t^d/(a+bb*t^d))^(f-1) * (cc*t^(d-1))/(e+bb*t^d)^2
   where
   aa=a*f*r
   bb=b/(c*3600)^d
   cc=bb*d*e

   JK, 29.12.2009
*/
class hydrationheat
{
 public:
  hydrationheat (void);
  ~hydrationheat (void);
  void read (XFILE *in);
  void print (FILE *out);
  double give_value (double t);
  long compare(hydrationheat &hh);
  
  double a,b,c,d,e,f,r;
  double aa,bb,cc;
};

#endif
