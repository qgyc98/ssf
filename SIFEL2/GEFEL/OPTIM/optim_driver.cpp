#include <stdlib.h>

#include "cmlfile.h"
#include "optim_driver.h"


/**
  The constructor serves for allocating memmory

  @param oat is type of optimization algorithm
  @param oF is pointer to objective function
 
  created $Date: 2003/12/23 16:07:57 $, $Author: leps $
*/
optim_driver::optim_driver ( algorithm_type &oat,
			     obj_funct *oF  )
{
  if ( !oF ) {
      fprintf (stderr,"\n\n You must allocate obj. function first! Exit.");
      fprintf (stderr,"\n optim_driver (file %s, line %d).\n",
	       __FILE__,__LINE__);
      exit ( 1 ) ;
  }    

  at=oat ;
  switch (at)
    {
    case alg_sade:{
      Sade = new sade ; 
      Sade->F = oF ;
      break ; 
    }
    case alg_grade:{ 
      Grade = new grade ; 
      Grade->F = oF ;
      break ; 
    }
    default:
      fprintf (stderr,"\n\n unknown type of optimization algorithm in");
      fprintf (stderr,"\n optim_driver (file %s, line %d).\n",
	       __FILE__,__LINE__);
      exit ( 1 ) ;
    }
      
  cmlfile f( "parameters.cfg" ) ;

  f.set_labels( 10 ) ;
  f.set_label_string(  0,"SadePoolRate" ) ;
  f.set_label_string(  1,"SadeRadioactivity" ) ;
  f.set_label_string(  2,"SadeLocalRadioactivity" ) ;
  f.set_label_string(  3,"SadeGradient" ) ;
  f.set_label_string(  4,"SadeCrossRate" ) ;
  f.set_label_string(  5,"SadeMutationRate" ) ;
  f.set_label_string(  6,"SadeMutagenRate" ) ;
  
  f.set_label_string(  7,"GradePoolRate" ) ;
  f.set_label_string(  8,"GradeRadioactivity" ) ;
  f.set_label_string(  9,"GradeCrossLimit" ) ;

  switch (at)
    {
    case alg_sade:{
      f.get_value( 0,Sade->pool_rate ) ;
      f.get_value( 1,Sade->radioactivity ) ;
      f.get_value( 2,Sade->local_radioactivity ) ;
      f.get_value( 3,Sade->gradient ) ;
      f.get_value( 4,Sade->cross_rate ) ;
      f.get_value( 5,Sade->mutation_rate ) ;
      f.get_value( 6,Sade->mutagen_rate ) ;
      break;
    }
    case alg_grade: {
      f.get_value( 7,Grade->pool_rate ) ;
      f.get_value( 8,Grade->radioactivity ) ;
      f.get_value( 9,Grade->cross_limit ) ;
      break ;
    }
    default:
      fprintf (stderr,"\n\n unknown type of optimization algorithm in");
      fprintf (stderr,"\n optim_driver (file %s, line %d).\n",
	       __FILE__,__LINE__);
    }
  f.close() ;
}

/**
  Destructor

  created $Date: 2003/12/23 16:07:57 $, $Author: leps $
*/
optim_driver::~optim_driver ( void )
{
  if ( Sade ) delete Sade ;
  if ( Grade ) delete Grade ;
}

/**
  Inner optimization cycle

  created $Date: 2003/12/23 16:07:57 $, $Author: leps $
*/
void optim_driver::run ( long calls_limit )
{
  switch (at)
    {
    case alg_sade:{
      Sade->fitness_calls_limit=calls_limit ; 
      Sade->run() ; break ;
    }
    case alg_grade:{ 
      Grade->fitness_calls_limit=calls_limit ; 
      Grade->run() ; break ; 
    }
    default:
      fprintf (stderr,"\n\n unknown type of optimization algorithm in");
      fprintf (stderr,"\n optim_driver (file %s, line %d).\n",
	       __FILE__,__LINE__);
      exit ( 1 ) ;
    }
}

