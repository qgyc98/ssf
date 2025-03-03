#ifndef BOCONT_H
#define BOCONT_H

#include <stdio.h>
#include "gfunct.h"
#include "iotools.h"

/**
  The class is used for storage of boundary conditions in 
  transport problems in the preprocessor.

  Created by TKo, 09.2010
*/
class bocont
{
  public:
   bocont();
   ~bocont();

   /// reads data from the opened text file
   long read(XFILE *in, long ndof);

   /// prints data to the opened text file in the TRFEL format
   long print(FILE *out);

   /// compares boundary conditions of the instance with the one defined in object tbc
   long compare(bocont &tbc);

   long     lcid;  ///< number of prescribed dofs
   double   iv;    ///< initial value
   double   con;   ///< prescribed value
   gfunct  *cgf;   ///< pointer to general function for BC in nonstationary problems
};

#endif
