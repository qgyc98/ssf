#ifndef TEMPLOAD_H
#define TEMPLOAD_H
#include <stdio.h>
#include "iotools.h"
class tempload;



/**
  This class groups temperature load data and it is used for the mechprep preprocessor.
*/
class tempload
{
  public :
   long nlc;    ///< load case number
   long nslc;   ///< subload case number
   double val;  ///< tempreture value
   tempload();
   ~tempload();
   long read(XFILE *in, long lc, long *slc);
   void copy(tempload &tl);
};

#endif
