#ifndef BOCON_H
#define BOCON_H

#include <stdio.h>
#include "iotools.h"

enum bctype {nobc = 0, sbc=1, gbc=2, tdbc = 3};

class gfunct;

class bocon
{
  public:
   bocon();
   ~bocon();
   void init(long ndofn);

   /// reads data from the opened text file
   long read(XFILE *in, long ndof);

   /// return pointer to new allocated object with copy of the actual instance
   bocon *copy(void);

   /// merges boundary conditions of the instance with the one defined in object tbc
   bocon *merge(bocon *tbc);

   long    nn;    ///< node number
   long    ndir;  ///< number of directions
   bctype  *dir;  ///< array flags for each direction
   long   *nspd;  ///< array with the numbers of nonzero prescribed displacements for particular directions in time independent load cases
   long   *ndpd;  ///< array with the numbers of nonzero prescribed displacements for particular directions in time dependent load cases
   double *con;   ///< array with prescribed values
   gfunct **gf;   ///< array with prescribed general functions
};

#endif
