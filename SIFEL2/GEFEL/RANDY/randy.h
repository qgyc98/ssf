/**
    class randy

    File:         randy.h
  
    Description:  C++ interface to the random number generator class
                  with some useful functions
  
    Author:       Matej Leps, leps@cml.fsv.cvut.cz
  
    $Id: randy.h,v 1.4 2004/03/19 14:29:04 leps Exp $
*/  

# ifndef __randy_h__
# define __randy_h__

# ifdef __make_fast__
#   define INLINE inline
# else
#   define INLINE /* */
# endif

# ifndef NearZero
# define NearZero 2.2204460492503131e-16
# endif

class randy
{ 
  public:
    randy ( void ) {} ;                    ///< constructor
   ~randy ( void ) {} ;                    ///< destructor
    INLINE void init ( long orseed=-1 ) ;  ///< initialization of rseed
  public:
    INLINE long give_long ( long omax ) ;
    INLINE long give_long ( long olow , long ohigh ) ;
    INLINE long give_long ( long olow, long ohigh, long oprecision ) ;
    INLINE double give_double ( long oa=0, long ob=0 ) ;
    INLINE double give_double ( double olow, double ohigh, long oa=0, long ob=0) ;
    INLINE double give_double ( double olow, double ohigh, double oprecision, long oa=0, long ob=0 ) ; 
  public:
    INLINE long is_less ( double oa ) ;
    INLINE long round ( long oa, long oprecision ) ;
    INLINE long move_to_bounds ( long oa, long olow, long ohigh ) ;
    INLINE double round ( double oa, double oprecision ) ;
    INLINE double move_to_bounds ( double oa, double olow, double ohigh ) ;
  public:
    INLINE long give_gauss_long ( double omi=0.0, double oss=1.0 ) ;
    INLINE double give_gauss_double ( double omi=0.0, double oss=1.0 ) ;
  public:
    long rseed ;        ///< random seed
};

# ifdef __make_fast__
#   include "randy.cpp_h"
# endif

# endif

//
//  History:
//
//  0001:         first version of random class
//  0002:         added rounding functions
//  0003:         new Makefile with install-fast condition
//                for all-  version
//  0004:         added (long) version of functions
//  0005:         gauss random number generator
//  0006:         some new notes in code
//  0007:         new features in get_real
//  0008:         name changed to randy, move "general.h" into self project
//                CVS managing
//  0009:         added comments for Doxygen, deleted all int's
//                changed long generation, get changed to give
//                and many other small changes (see CVS status)
//  0010:         added README, INSTALL; deleted links to <general.h>
//  0011:         Fixed error for implicit arguments in give_gauss_long()
//                => deleted give_gauss_long(long,long), now only (double,double)
//
