#ifndef INICD_H
#define INICD_H

#include "iotools.h"
#include "alias.h"

/**
  The class holds data about initial condition in the given node.
*/
class inicd
{
  public :
   inicd (void);
   ~inicd(void);
   /// function reads initial condition from the opened text file
   void read(XFILE *in);

   /// function prints initial condition to the opened text file
   void print(FILE *out);

   /// copies data from the parameter to the actual object
   void copy(inicd &ic);

   /// merges function from the parameter with the actual object
   long merge(inicd &ic);

   inictype type;   ///< type of initial condition
   double  *val;    ///< values of initial conditions of given node
   long     nval;   ///< number of values
};

#endif
