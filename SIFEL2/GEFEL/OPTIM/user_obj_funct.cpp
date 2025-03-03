#include <obj_funct.h>

/**
 The function "user_allocate( void )" is here for user. Here can be
 allocated memory or defined some variables for user-defined
 objective function, like opening of needed files, definition of known
 optimum etc. Please, delete all allocated memory in user_deallocate()
 function.

  created $Date: 2004/03/19 13:58:35 $, $Author: leps $
*/
void obj_funct::user_allocate( void )
{
}

/**
 The "user_value ( double *oCH )" function returns value of objective
 function. Here the user can specify his own objective function, which
 must be returned be "return()" statement.

  @param oCH is pointer to actual set of optimized variables ( x[i] )

  @b Requirements :
    oCH must be different from NULL.

  @retval double precision number = objective value f( x[i] )

  created $Date: 2004/03/19 13:58:35 $, $Author: leps $
*/
double obj_funct::user_value ( double * /*oCH*/ )
{
  return ( 0.0 ) ;
}

/**
 In this function "obj_funct::user_deallocate ()" the user must specify
 deleting statements for variables allocated by him in
 obj_funct::user_allocate ().

  created $Date: 2004/03/19 13:58:35 $, $Author: leps $
*/
void obj_funct::user_deallocate ( void )
{
}

/**
 The "user_evaluate ( double *BSF )" function is called at the end of
 optimization process. It is called from Optimization algorithm and in
 BSF user can find best-so-far solution found during optimization. The
 usual use of this function is to get e.g. the stress-strain diagram of
 the optimal structure etc.

  @param oBSF is pointer to the best solution found

  created $Date: 2004/03/19 13:58:35 $, $Author: leps $
*/
void obj_funct::user_evaluate ( double * /*BSF*/ )
{
}

