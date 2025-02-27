#ifndef FUZZYNUM_H
#define FUZZYNUM_H

#include "iotools.h"


/**
   class fuzzynum
   
   
   JK, 18.8.2005
*/
class fuzzynum
{
 public:
  fuzzynum (void);
  ~fuzzynum (void);
  
  void read (XFILE *in);
  void initiate (long n,double *al);
  
  double give_min (double alpha);
  double give_min (long acutid);
  double give_max (double alpha);
  double give_max (long acutid);
  
  void save_min (double alpha,double min);
  void save_min (long acutid,double min);
  void save_alp_min (long acutid,double al,double min);
  void save_max (double alpha,double max);
  void save_max (long acutid,double max);
  void save_alp_max (long acutid,double al,double max);
  
  void print (FILE *out);
  void print_minmax (FILE *out);
  
  double give_val (long i);
  //void onearray ();
  
  void save_value (double val);
  void minmax_init ();

  
  ///  number of alpha-cuts
  long nalph;
  ///  array of alpha-cut values
  double *alph;
  ///  minimum values of alpha-cuts
  double *xmin;
  ///  maximum values of alpha-cuts
  double *xmax;
  ///  number of values on axis
  long nval;
  ///  values of alpha-cuts in increasing ordering
  double *x;

  ///  required error/deviation
  double err;
    
};

#endif
