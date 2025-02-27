#ifndef ENTITYTDLOAD_H
#define ENTITYTDLOAD_H

#include "parser.h"
#include "siftop.h"
#include "iotools.h"
#include "dloadel.h"
#include "alias.h"

#include <stdio.h>

/**
  This class is used for entity time dependnet load in the mechprep preprocessor.
  Entity load must be defined by general function.
  The general function has to be defined either with one parameter for time function independent on spatial coordinates 
  or four parameters parser function for x, y, z spatial coordinates and time t in order to compute load values 
  for given nodal coordinates.
  
  Created by Tomas Koudelka, 25.4.2016
*/
class entitytdload
{
  public :
   entitytdload();
   ~entitytdload();

   /// funtion reads entity load from the opened text file
   long read(XFILE *in, long lc);

   /// function computes values of the entity load at the given node
   long getval(snode &tn, gfunct *nv);

   /// function checks the equation for proper names of variables (x,y,z,t) and returns number of different spatial coordinates in the equation in the argument nspatc
   long var2coord(Equation *eq, long &nspatc);

   /// function returns pointer to a new dloadel structure with time dependent edge load created from the given entity load, element and edge property id
   dloadel *tdedge2dloadel(selement &el, long prop, snode *nodes, long *edgid, long nedgid);

   /// function returns pointer to a new dloadel structure with time dependent surface load created from the given entity load, element and surface property id
   dloadel *tdsurface2dloadel(selement &el, long prop, snode *nodes, long *surfid, long nsurfid);

   /// function returns pointer to a new dloadel structure with time dependent volume load created from the given entity load, element and volume property id
   dloadel *tdvol2dloadel(selement &el, long prop, snode *nodes);

   long    ncomp; ///< number of load components (directions)
   long    nlc;   ///< load case number
   long    lgcs;  ///< coordinate system on edge (local=1, global=2)
   gfunct *func;  ///< array of pointers to time functions in case of load defined by time function on the given entity
};

#endif
