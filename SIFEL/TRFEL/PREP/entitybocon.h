#ifndef ENTITYBOCON_H
#define ENTITYBOCON_H

#include "parser.h"
#include "siftop.h"
#include "iotools.h"
#include "tablefunct.h"
#include "loadelt.h"
#include "aliast.h"

#include <stdio.h>

/**
  This class is used for entity boundary condition (BC) in the
  transprep preprocessor.  Entity BC can be defined either as constant, table
  or by parser function. The parser function may be defined
  with one up to four parameters for x, y, z coordinate and time (t) in
  order to compute BC values for given nodal coordinates. Parser expressions are evaluated for x,y,z 
  of the given node and converted to the new gfunct object: 
  'stat'             -> 'stat'
  'tab'(t)           -> 'tab'(t)
  'pars'(x,y,z)      -> 'stat
  'pars'(x,y,z,t)    -> 'pars'(t) 
  'pars_set'(x,y,z)  -> 'tab'(t) (interpoltype has to be specified) 
  'pars_set'(x,y,z,t)-> 'par_set'(t)

  Created by TKo, 09.2010
*/
class entitybocon
{
  public :
   entitybocon();
   ~entitybocon();

   /// funtion reads entity load from the opened text file
   long read(XFILE *in);

   /// function computes values of the entity BC at the given node
   long getval(snode &tn, gfunct &tgf);

   /// function checks x, y, z and t coordinate parameters in the BC function
   long checkvar2coord();

   /// function detects indeces of x, y, z and t coordinate parameters of the BC function for the i-th parsed expression
   long var2coord(long i);

   /// function sets internal variables of func to the coordinates of the given node
     void setvars(snode &tn, long i);

   /// function replaces occurence of x, y and z parameter by the true values of coordinates of the given node and returns the new resulting string
   char *substbcstr(const char *expr, snode &tn);

   gfunct       gf; ///< function describing the values of BC on the given entity
   interpoltype it; ///< interpolation type in the case that BC is defined by function on the given entity
   long         dc; ///< flag for dependency of gf on coordinates (x,y,z)
   long         dt; ///< flag for dependency of gf on time (t)
   

   long var_x;  ///< index of x-coordinate variable in func
   long var_y;  ///< index of y-coordinate variable in func
   long var_z;  ///< index of z-coordinate variable in func
};

#endif
