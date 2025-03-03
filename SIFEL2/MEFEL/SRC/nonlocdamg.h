#ifndef NONLOCDAMG_H
#define NONLOCDAMG_H

#include "iotools.h"
#include "alias.h"
#include "vector.h"



/**
  The class defines general nonlocal formulation of damage material models.

  Created by Tomas Koudelka,
*/
class nonlocdamg
{
  public :
   nonlocdamg(void);
   ~nonlocdamg(void);
   /// reads data from the file in
   long read(XFILE *in);
   /// prints data to the file out
   void print(FILE *out);
   /// number of averaged quantities
   long give_num_averq (long ipp, long im);
   void give_aver_quantv(long ipp,long im,long ido, vector &qv);
   /// averages quantities in the given integration point
   void average (long ipp, long im, long ido);
   /// returns the length of fracture process zone
   double give_proczonelength(long ipp, long im, long ido);

   /// taken radius for ipp environment
   double r;
   /// flag for the averaged variable
   avrgf af;

};
#endif
