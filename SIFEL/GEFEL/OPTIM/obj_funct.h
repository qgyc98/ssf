/**
    class obj_funct

    File:         obj_funct.h
  
    Description:  General definition for all objective functions
                  that can appear in SIFEL project.
  
    Author:       Anna Kucerova, anicka@cml.fsv.cvut.cz
                  Matej Leps, leps@cml.fsv.cvut.cz
  
    $Id: obj_funct.h,v 1.6 2003/12/05 20:31:07 leps Exp $
*/  

# ifndef __obj_funct_h__
# define __obj_funct_h__

#if !defined ( __general_h_ )
#include "general.h"
#endif


class obj_funct
{
    public:
        obj_funct ( long oDim );             ///< constructor
        ~obj_funct ( void );                 ///< destructor

        void user_allocate ( void );         ///< user's allocation
	double value ( double *oCH );        ///< general objective value
	double user_value ( double *oCH );   ///< user's objective value
	void sifel_value ( double *oCH );    ///< SIFEL preparation
        void user_deallocate ( void );       ///< user's deallocation
        void user_evaluate ( double *oBSF ); ///< final evaluation

        long Dim;                   ///< number of optimization variables
        p_double *Domain;           ///< definition of searched domain
        double *optimum;            ///< if known, here is optimum
        double precision;           ///< if optimum, this is distance to it
        char *outfile;              ///< name of output file
        char name[256];             ///< name of function
        long Return_to_domain;      ///< solutions can be outside domain
    private:

};

# endif // __obj_funct_h__
