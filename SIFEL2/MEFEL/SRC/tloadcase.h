#ifndef TLOADCASE_H
#define TLOADCASE_H

#include <stdio.h>
#include "alias.h"

class gfunct;
class dloadn;
class dloadpd;
class loadcase;

/**
   class tloadcase defines load cases for time dependent problems
   (dynamic problems, creep problems, etc.)
   
   JK, TKo, 3.6.2005
*/

class tloadcase
{
 public:
  tloadcase (void);
  ~tloadcase (void);
  void read (FILE *in);
  void assemble (long lcid,double *rhs,long n, double t);
  void assemble (double *rhs, double *lhs);
  void compute_reactions (long lcid);
  void seisminit (double *seism);
  
  ///  type of time load
  timeload ttl;
  ///  directions of seismic loads
  dirdynload *direction;

  //  number of subload cases
  long nslc;
  
  ///  number of loaded nodes
  long nln;
  ///  number of loaded elements
  long nle;
  ///  number of prescribed displacements
  long npd;
  
  dloadn  *lon;
  dloadpd *pd;
  
  //  subload cases
  loadcase *slc;
  //  time functions
  gfunct *gf;
  //  array for seismic loads
  double *seism;
};

#endif
