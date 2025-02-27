#include "obj_funct.h"
#include <math.h>
#include <stdlib.h>

/**
 The constructor obj_funct ( long oDim ) allocates memmory for
 objective function. Each optimization variable is limited by three
 values in Domain: [min,max,precision]. Therefore the whole domain is
 defined by Domain[Dim][3]. User can also specify and allocate
 optimum, if it is known. If Return_to_domain is set to zero, the
 algorithm will have chance to find optimum even outside given
 boundaries. At the end, this function calls obj_funct::user_allocate (),
 placed in user_obj_funct.cpp, where the user can specify his own
 allocation statements.

  @param oDim is number of optimized variables

  @b Requirements :
    oDim must be non-negative.

  created $Date: 2003/12/05 20:31:07 $, $Author: leps $
*/
obj_funct::obj_funct ( long oDim )
{
    long i;
    Dim=oDim;
    Domain = new p_double [Dim];
    for ( i=0; i<Dim; i++ ) Domain[i] = new double [3];
    optimum = NULL;
    precision = 0.0;
    Return_to_domain=0;
    outfile = new char [256];

    this->user_allocate () ;
}

/**
 The destructor ~obj_funct () deallocates all basic dynamic
 variables. It also calls obj_funct::user_deallocate (), placed in
 user_obj_funct.cpp, where user must specify deleting statements for
 variables allocated by him in obj_funct::user_allocate ().

  created $Date: 2003/12/05 20:31:07 $, $Author: leps $
*/
obj_funct::~obj_funct ( void )
{
    this->user_deallocate () ;

    long i;
    for ( i=0; i<Dim; i++ ) delete [] Domain[i];
    delete [] Domain;
    if ( optimum ) delete optimum;
}

/**
 The "value ( double *oCH )" function returns value of objective
 function. Calls sifel_value(), where some preparations for
 appropriate problem (MEFEL,TRFEL,...) are placed and secondly, calls
 user_value(), which is placed in user_obj_funct.cpp, and where user
 can specify his own objective function.

  @param oCH is pointer to actual set of optimized variables ( x[i] )

  @b Requirements :
    oCH must be different from NULL.

  @retval double precision number = objective value f( x[i] )

  created $Date: 2003/12/05 20:31:07 $, $Author: leps $
*/
double obj_funct::value ( double *oCH )
{
  this->sifel_value( oCH ) ;
  return ( this->user_value( oCH ) ) ;
}

