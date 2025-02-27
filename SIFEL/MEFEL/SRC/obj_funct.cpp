#include "obj_funct.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/**
 This function is here only for linking purposes. Please don't use
 this file to put here your objective function. See
 SIFEL/PRG/GEFEL/OPTIM/INSTALL file for instructions, how to create
 your own optimization problem.

  created $Date: 2003/12/05 20:43:50 $, $Author: leps $
*/
void obj_funct::user_allocate( void )
{
  fprintf ( stderr, "You are using blank objective function! \n" ) ;
  fprintf ( stderr, "Create your own! \n" ) ;
  fprintf ( stderr, "Obj_funct::allocate() (%s, line %d)\n",\
	    __FILE__,__LINE__) ;
}

/**
 This function is here only for linking purposes. Please don't use
 this file to put here your objective function. See
 SIFEL/PRG/GEFEL/OPTIM/INSTALL file for instructions, how to create
 your own optimization problem.

  created $Date: 2003/12/05 20:43:50 $, $Author: leps $
*/
double obj_funct::user_value ( double * )
{
  return ( 0.0 ) ;
}

/**
 This function is here only for linking purposes. Please don't use
 this file to put here your objective function. See
 SIFEL/PRG/GEFEL/OPTIM/INSTALL file for instructions, how to create
 your own optimization problem.

  created $Date: 2003/12/05 20:43:50 $, $Author: leps $
*/
void obj_funct::user_deallocate ( void )
{
}

/**
 This function is here only for linking purposes. Please don't use
 this file to put here your objective function. See
 SIFEL/PRG/GEFEL/OPTIM/INSTALL file for instructions, how to create
 your own optimization problem.

  created $Date: 2003/12/05 20:43:50 $, $Author: leps $
*/
void obj_funct::user_evaluate ( double */*BSF*/ )
{
}
