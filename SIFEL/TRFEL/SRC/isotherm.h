#ifndef ISOTHERM_H
#define ISOTHERM_H

#include "aliast.h"
#include "genfile.h"
#include "sorpisohansen.h"
#include "sorpisorootgrun.h"

/**
   class contains tools for any isotherms description and manipulation
   
   JK, 2.10.2013
*/
class isotherm
{
 public:
  isotherm ();    //  constructor
  ~isotherm ();   //  destructor
  
  void read (XFILE *in);
  void print (FILE *out);
  
  double isotherm_value (double in, long ipp=0L, long eid=0L);
  double inverse_isotherm_value (double in);
  double derivative_isotherm_value (double in);
  double derivative_inverse_isotherm_value (double in);

  
  ///  type of isotherm
  isotypet isothermtype;
  
  ///  general function for isotherm given by set of data
  gfunct data;
  
  ///  Hansen sorption isotherm
  sorpisohansen hanseni;
  ///  Root sorption isotherm
  sorpisorootgrun rooti;

};

#endif
