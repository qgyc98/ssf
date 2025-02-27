//
//  File:         general.h
//
//  Description:  File with general statement for EA's
//                for "wien" project
//
//  Author:       Matej Leps
//
//  Date:         23/April/2001
//
//  Version:      0001
//

#ifndef __general_h__
#define __general_h__

#include <limits.h>

#if !defined ( __GENERIC_CONSTANTS__ )
#define        __GENERIC_CONSTANTS__


# define MinDouble  -1.7976931348623150e+308
# define MaxDouble  +1.7976931348623150e+308
# define LowDouble  MinDouble
# define HighDouble MaxDouble
# define NoDouble -7.87813047009E+177
# ifndef NearZero
# define NearZero 2.2204460492503131e-16
# endif

# endif // __GENERIC_CONSTANTS__

/* Data from '/usr/include/limits.h' */
# ifdef LONG_MAX
#   define MaxLong LONG_MAX
# else
#   define MaxLong 2147483647L
# endif
# define MinLong (-MaxLong - 1L)
# define LowLong  MinLong
# define HighLong MaxLong


# define Minimum 0
# define Maximum 1
# define Yes 1
# ifndef No
#   define No 0
# endif

typedef double* p_double ;
typedef int* p_int ;
typedef char* p_char ;
typedef long* p_long ;

#ifndef sqr
#define sqr( x ) ( (x)*(x) )
#endif 

template <class ItemType>
inline void SWAP( ItemType a, ItemType b )
{
  ItemType temp = a ;
  a = b ;
  b = temp ;
}


#endif // __general_h__
