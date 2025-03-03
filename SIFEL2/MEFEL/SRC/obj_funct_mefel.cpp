#include "obj_funct.h"
#include "solverm.h"
#include "lssolver.h"
#include "edsolver.h"
#include "fdsolver.h"
#include "nssolver.h"
#include "mtsolver.h"
#include "llssolver.h"
#include "epsolver.h"
#include "slsolver.h"
#include "global.h"

#include <stdlib.h>




/**
 The "sifel_value ( double * )" function prepares data for sifel-based
 optimization. Therefore all statements appropriate to the problem
 (MEFEL,TRFEL,...) are placed here. If the user want to use
 optimization without FEM computation, he must comment all lines in
 this function.

  @b Requirements :
    Global stochdriver St form "global.h" must be allocated.

  created $Date: 2003/12/05 17:14:53 $, $Author: leps $
*/
void obj_funct::sifel_value ( double * )
{
  long i=0 ;
  if ( St != NULL ) 
    {
      St->changevalues (i);
      solve_mefel_deterministic_problem ();
      St->extractor ();
      St->save_results (i);
    }
  else
    {
      fprintf ( stderr, "The stochastic driver is not allocated! \n" ) ;
      fprintf ( stderr, "obj_funct::sifel_value (%s, line %d)\n",\
		__FILE__,__LINE__) ;
      exit( 1 ) ;
    }
}
