/**
    class optim_driver

    File:         optim_driver.h
  
    Description:  this class controls optimization computations
  
    Author:       Matej Leps, leps@cml.fsv.cvut.cz
  
    $Id: optim_driver.h,v 1.4 2004/02/18 10:13:51 leps Exp $
*/  

# ifndef __optim_driver_h__
# define __optim_driver_h__

# include "obj_funct.h"
# include "sade.h"
# include "grade.h"

enum algorithm_type { alg_sade=1,
		      alg_grade=2 } ;

class optim_driver
{
    public:
        optim_driver ( algorithm_type &oat,
		       obj_funct *oF );       ///< constructor
        ~optim_driver( void );                ///< destructor
 
	void run ( long calls_limit );        ///< optimization cycle

    private:
	/// Type of algorithm
	algorithm_type at ;
	/// Individual optimization algorithms
	sade *Sade ;
	grade *Grade ;
};

# endif // __optim_driver_h__
