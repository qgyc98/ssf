#ifndef ADVECTVEL_H
#define ADVECTVEL_H

#include "parser.h"
#include "siftop.h"
#include "iotools.h"
#include "galias.h"
#include "gfunct.h"

#include <stdio.h>

/**
  This class is used for the definition of advection velocities at nodes in the mechprep preprocessor.
  Velocity components can be defined either as constant load or by parser function. 
  The parser function can be defined with up three parameters for
  x, y, z coordinate in order to compute velocity values for given nodal coordinates.
*/
class advectvel
{
  public :
   advectvel();
   ~advectvel();

   /// funtion reads entity load from the opened text file
   long read(XFILE *in, long ndof);

   /// function returns indeces of x, y and z coordinate parameters in the load function
   long var2coord(Equation *eq, long &var_x, long &var_y, long &var_z);

   /// function computes values of the entity load at the given node
   long getval(snode &tn, long compid, double &vc);

   long ncomp;   ///< number of load components (directions)
   long lcid;    ///< load case id, (=dof id = medium id)
   gfunct *v;    ///< array of velocity components
};
#endif
