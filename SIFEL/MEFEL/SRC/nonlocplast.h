#ifndef NONLOCPLAST_H
#define NONLOCPLAST_H

#include "iotools.h"
#include "alias.h"
#include "vector.h"



/**
  The class defines general nonlocal formulation of plastic material models.

  Created by Tomas Koudelka,
*/
class nonlocplast
{
  public :
   nonlocplast(void);
   ~nonlocplast(void);
   /// reads data from the file in
   long read(XFILE *in);
   /// returns number of averaged quantities
   long give_num_averq (long ipp);
   /// function returns vector of averaged quantities
   void give_aver_quantv(long ipp,long im,long ido, vector &qv);
   /// averages values in the given integration point
   void average (long ipp, long ido);

   /// taken radius for ipp environment
   double r;
   /// flag for the averaged variable
   wavrg waf;
};
#endif
