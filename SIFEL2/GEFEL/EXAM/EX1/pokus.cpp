#include <stdio.h>

#include <obj_funct.h>

/**
  created $Date: 2003/12/23 16:15:00 $, $Author: leps $
*/
void obj_funct::user_allocate( void )
{
}

/**
   Example of objective function, returns sqruare root.

  created $Date: 2003/12/23 16:15:00 $, $Author: leps $
*/
double obj_funct::user_value ( double *oCH )
{
  return ( -oCH[0]*oCH[0] ) ;
}

/**
  created $Date: 2003/12/23 16:15:00 $, $Author: leps $
*/
void obj_funct::user_deallocate ( void )
{
}

/**
  created $Date: 2003/12/23 16:15:00 $, $Author: leps $
*/
void obj_funct::user_evaluate ( double *BSF )
{
  fprintf( stderr, "* The best value found was: %f \n", user_value(BSF) ) ;
  fprintf( stderr, "* At point:                 %f \n\n", BSF[0] ) ;

}

