#ifndef ENTITYLOAD_H
#define ENTITYLOAD_H

#include "parser.h"
#include "siftop.h"
#include "iotools.h"
#include "loadel.h"
#include "alias.h"

#include <stdio.h>

/**
  This class is used for entity load in the mechprep preprocessor.
  Entity load can be defined either as constant load or by parser function. 
  The parser function has to be defined with three parameters for
  x, y, z coordinate in order to compute load values for given nodal coordinates.
*/
class entityload
{
  public :
   entityload();
   ~entityload();

   /// funtion reads entity load from the opened text file
   long read(XFILE *in, long lc, long *slc);

   /// function computes values of the entity load at the given node
   long getval(snode &tn, double *nv);

   /// function returns indeces of x, y and z coordinate parameters in the load function
   long var2coord(Equation *eq, long &var_x, long &var_y, long &var_z);

   /// function returns pointer to a new loadel structure with edge load created from the given entity load, element and edge property id
   loadel *edge2loadel(selement &el, long prop, snode *nodes, long *edgid, long nedgid);

   /// function returns pointer to a new loadel structure with surface load created from the given entity load, element and surface property id
   loadel *surface2loadel(selement &el, long prop, snode *nodes, long *surfid, long nsurfid);

   /// function returns pointer to a new loadel structure with volume load created from the given entity load, element and volume property id
   loadel *vol2loadel(selement &el, long prop, snode *nodes);

   long ncomp;       ///< number of load components (directions)
   long nlc;         ///< load case number
   long nslc;        ///< subload case number
   generalfunct ft;  ///< type of load function
   long lgcs;        ///< coordinate system on edge (local=1, global=2)
   double *val;      ///< load values in case of constant load
   Equation **func;  ///< array of pointers to parsed functions in case of load defined by function on the edge
   char     **tfunc; ///< array of function defintion strings in case of load defined by function on the edge
};

#endif
