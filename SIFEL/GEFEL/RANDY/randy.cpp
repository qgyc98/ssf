/**
    class randy
  
    File:         randy.cpp
  
    Description:  C++ interface to the random number generator class
                  with some useful functions
  
    Author:       Matej Leps, leps@cml.fsv.cvut.cz
  
    $Id: randy.cpp,v 1.4 2004/03/19 14:32:14 leps Exp $
*/  

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef __make_fast__
#include "randy.h"
#endif

# ifndef M_PI
#   define M_PI 3.14159265358979323846  /* pi */
# endif

# ifndef RAND_MAX
#   define RAND_MAX 32767
# endif


/**
 The rnd2 ( &idum ) function is random number generator from
  Numerical Recipes in C++, version 3?

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double ran2 ( long &idum )
{
  const long IM1=2147483563 ;
  const long IM2=2147483399 ;
  const long IA1=40014 ;
  const long IA2=40692 ;
  const long IQ1=53668 ;
  const long IQ2=52774 ;
  const long IR1=12211 ;
  const long IR2=3791 ;
  const long NTAB=32 ;
  const long IMM1=IM1-1 ; 
  const long NDIV=1+IMM1/NTAB ;
  const double EPS=3.0e-16 ;
  const double RNMX=1.0-EPS ;
  const double AM=1.0/(double)IM1 ;
  static long idum2=123456789 ;
  static long iy=0 ;
  static long iv[NTAB] ;

  long j,k ;
  double temp;

  if (idum <= 0) {
    idum=(idum==0 ? 1 : -idum) ;
    idum2=idum ;
    for (j=NTAB+7;j>=0; j--) { 
      k=idum/IQ1 ;
      idum=IA1*(idum-k*IQ1)-k*IR1 ;
      if (idum < 0) idum += IM1 ;
      if (j < NTAB) iv[j] = idum ;
    } 
    iy=iv[0] ;
  } 
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

/**
 The rnd function returns random number from < 0, RAND_MAX >.
  Is here only for debug reasons.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long rnd ( void )
{
#ifdef __TEST_BOUNDS__

  const long pole[4]={ 0L, 1L, RAND_MAX-1L, RAND_MAX } ;
  static long i=-1 ;

  i++ ;
  //printf ( "%12ld \n", pole[i%4] ) ;
  return ( pole[i%4] ) ;

#else

  return ( rand() ) ;

#endif
}

/**
 The init function initialize random seed

  @param orseed is user defined seed.
         If there is nothing, process will be pseudo-random

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
void randy::init ( long orseed )
{
  //
  // this code is running only under Linux :
  //
  // struct timeval tv;
  // struct timezone tz;
  // gettimeofday( &tv,&tz ) ;
  // rseed=tv.tv_usec+( tv.tv_sec%4194304 )*1000 ;
  //

  if ( orseed==-1 )
    {
      time_t t;
      t = time(&t);
      rseed=long(t) ;
    }
  else 
    {
      rseed=orseed ;
    }
  
  srand( rseed ) ;
}

/**
 The give_long ( long omax ) returns random number from < 0; omax-1 >

  @param omax is upper bound of the generated number

  @b Requirements :
   Number omax has to be greater than 0.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::give_long ( long omax )
{
  if ( omax==1 )
    return ( 0 ) ;
  else if ( omax>0 )
    {
      // higher bits are used
      return ( (long)floor((double)omax * ((double)(rnd() - 1)/(double)( RAND_MAX ))) ) ;
    }
  else
    {
      fprintf ( stderr, "\n\n Higher value is lower than zero. \n" ) ;
      fprintf ( stderr, " randy (%s, line %d)\n",__FILE__,__LINE__) ;
      exit ( 1 ) ;
    }
  return( 0 ) ;
}

/**
 The give_long ( long olow, long ohigh ) returns random number from < olow; ohigh >.
  If the value from ( olow; ohigh ) is requested, the function have to be called as
  give_long ( olow+1 ; ohigh-1 ) ;

  @param olow is lower bound of the generated number
  @param ohigh is upper bound of the generated number

  @b Requirements :
   Number olow has to be smaller than ohigh.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::give_long ( long olow, long ohigh )
{
  long div = ohigh-olow ;
  if ( div>0 )
    return ( this->give_long( div+1 )+olow ) ;
  else if ( div==0 )
    return ( olow ) ;
  else
    {
      fprintf ( stderr, "\n\n Higher value is smaller than lower value. \n" ) ;
      fprintf ( stderr, " randy (%s, line %d)\n",__FILE__,__LINE__) ;
      exit ( 1 ) ;
    }
  return( 0 ) ;
}

/**
 The give_long ( long olow, long ohigh, long oprecision ) returns random number
  from < olow; ohigh > with the given precision oprecision.
  If the value from ( olow; ohigh ) is requested, the function have to be called as
  give_long ( olow+1 ; ohigh-1 ) ;

  @param olow is (inclusive) lower bound of the generated number
  @param ohigh is (inclusive) upper bound of the generated number
  @param oprecision is precision of the generated number

  @b Requirements :
   Number olow has to be smaller than ohigh.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::give_long ( long olow, long ohigh, long oprecision )
{
  return( this->round( this->give_long( olow, ohigh ), oprecision )) ;
}

/**
 The give_double ( oa, ob ) returns random number 
  from < 0 ; 1 > for a,b = {0,0} (default settings),
  from < 0 ; 1 ) for a,b = {0,1},
  from ( 0 ; 1 > for a,b = {1,0},
  from ( 0 ; 1 ) for a,b = {1,1}.

  @param oa is tag for exclusion of zero
  @param ob is tag for exclusion of one

  @b Requirements :
   Numbers oa, ob has to be zeros or ones.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double randy::give_double ( long oa, long ob )
{
  if (( oa==0 || oa==1 ) && ( ob==0 || ob==1 ))
    {
      // re-type long to double to overcome overloading
      double tmp = ((double)( rnd() ) + (double)( oa )) / ((double)( RAND_MAX ) + (double)( oa + ob )) ;
      return ( tmp  ) ;
    }
  else
    {
      fprintf ( stderr, "\n\n Values oa and ob must be zeros or ones.. \n" ) ;
      fprintf ( stderr, " randy (%s, line %d)\n",__FILE__,__LINE__) ;
      exit ( 1 ) ;
    }

  return( 0. ) ;
}

/**
 The give_double ( olow, ohigh, oa, ob ) returns random number 
  from < olow ; ohigh > for a,b = {0,0} (default settings),
  from < olow ; ohigh ) for a,b = {0,1},
  from ( olow ; ohigh > for a,b = {1,0},
  from ( olow ; ohigh ) for a,b = {1,1}.

  @param olow is lower bound of the generated number
  @param ohigh is upper bound of the generated number
  @param oa is tag for exclusion of zero
  @param ob is tag for exclusion of one

  @b Requirements :
   Number olow has to be smaller than ohigh.
   Numbers oa, ob has to be zeros or ones.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double randy::give_double ( double olow, double ohigh, long oa, long ob )
{
  double div = ohigh-olow ;
  if ( div>NearZero )
    return ( this->give_double(oa, ob)*div+olow ) ;
  else if ( div>=0. )
    return ( olow ) ;
  else
    {
      fprintf ( stderr, "\n\n Higher value is smaller than lower value. \n" ) ;
      fprintf ( stderr, " randy (%s, line %d)\n",__FILE__,__LINE__) ;
      exit ( 1 ) ;
    }
  return( 0. ) ;
}

/**
 The give_double ( olow, ohigh, oprecision, oa, ob ) returns random number 
  from < olow ; ohigh > for a,b = {0,0} (default settings),
  from < olow ; ohigh ) for a,b = {0,1},
  from ( olow ; ohigh > for a,b = {1,0},
  from ( olow ; ohigh ) for a,b = {1,1}
  with the given precision oprecision.

  @param olow is lower bound of the generated number
  @param ohigh is upper bound of the generated number
  @param oa is tag for exclusion of zero
  @param ob is tag for exclusion of one

  @b Requirements :
   Number olow has to be smaller than ohigh.
   Numbers oa, ob has to be zeros or ones.

  @retval random number

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double randy::give_double ( double olow, double ohigh, double oprecision, long oa, long ob )
{
  return( this->round( this->give_double( olow, ohigh, oa, ob ), oprecision )) ;
}

/**
 The is_less ( oa ) function returns 0 if random number form <0,1) is
  less than oa, returns 1 otherwise.

  @param oa is probability mark, e.g. from roullete wheel selection

  @b Requirements :
   Number oa is hoped to be from <0,1>.

  @retval 0L or 1L

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::is_less ( double oa )
{
  return( this->give_double(0L,1L)<oa ) ;
}

/**
 The round ( oa, oprecision ) function rounds number oa to oprecision precision.

  @param oa is rounded number
  @param oprecision is precision

  @b Requirements :
   Number oprecision has to be positive.

  @retval oa if oprecision is 0, rounded oa otherwise

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::round ( long oa, long oprecision )
{
  if ( oprecision==0 )
    {
      return ( oa ) ;
    }
  else if ( oprecision > 0 )
    {
      return( (long)floor( this->round( (double) oa, (double) oprecision )) ) ; 
    }
  else
    {
      fprintf ( stderr, "\n\n Precision has to be non-negative. \n" ) ;
      fprintf ( stderr, " randy (%s, line %d)\n",__FILE__,__LINE__) ;
      exit ( 1 ) ;
    }
    return oa;
}

/**
 The move_to_bounds ( oa, olow, ohigh ) move oa to interval <olow,ohigh>.

  @param oa is input nuber
  @param olow is lower bound
  @param ohigh is upper bound

  @retval number from <olow,ohigh>

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::move_to_bounds ( long oa, long olow, long ohigh )
{
  if( oa<olow )
    return( olow ) ;
  else if( oa>ohigh )
    return( ohigh ) ;
  else
    return( oa ) ;
}

/**
 The round ( oa, oprecision ) function rounds number oa to oprecision precision.

  @param oa is rounded number
  @param oprecision is precision

  @b Requirements :
   Number oprecision has to be positive.

  @retval oa if oprecision is 0, rounded oa otherwise

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double randy::round ( double oa, double oprecision )
{
  if ( fabs(oprecision) < NearZero )
    {
      return ( oa ) ;
    }
  else if ( oprecision > 0 )
    {
      return( ( floor(( oa/oprecision )+0.5 ))*oprecision ) ;
    }
  else
    {
      fprintf ( stderr, "\n\n Precision has to be non-negative. \n" ) ;
      fprintf ( stderr, " randy (%s, line %d)\n",__FILE__,__LINE__) ;
      exit ( 1 ) ;
    }
  return oa;
}

/**
 The move_to_bounds ( oa, olow, ohigh ) move oa to interval <olow,ohigh>.

  @param oa is input nuber
  @param olow is lower bound
  @param ohigh is upper bound

  @retval number from <olow,ohigh>

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double randy::move_to_bounds ( double oa, double olow, double ohigh )
{
  if( oa<olow )
    return( olow ) ;
  else if( oa>ohigh )
    return( ohigh ) ;
  else
    return( oa ) ;
}

/**
 The give_gauss_long ( omi, oss ) function returns random numbers with
  Gaussian normal distribution.

  @param omi is mean value
  @param oss is std. dev.

  @b Requirements :
   Number oss has to be non-negative.

  @retval long

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
long randy::give_gauss_long ( double omi, double oss )
{
  return( (long)round( give_gauss_double( omi,oss ), 1.0 )) ;
}

/**
 The give_gauss_double ( omi, oss ) function returns random numbers with
  Gaussian normal distribution.

  @param omi is mean value
  @param oss is std. dev.

  @b Requirements :
   Number oss has to be non-negative.

  @retval double

  created $Date: 2004/03/19 14:32:14 $, $Author: leps $
*/
double randy::give_gauss_double ( double omi, double oss )
{
  // generate 2 indep numbers in every 2nd call
  static bool value_stored=false ;  
  static double value=0.0 ;
  double fac, v1, v2, rsq ;

  if ( value_stored ) 
    {
      value_stored=false ;
      return value;
    }
  else 
    {
    // v1 v2 in unitcircle
    do 
      {
        v1=2.*give_double( 0L,1L )-1. ;
        v2=2.*give_double( 0L,1L )-1. ;
        rsq=v1*v1 + v2*v2 ;
      } 
    while ( rsq>=1. || rsq==0. ) ;      

    fac=sqrt(-2. * log(rsq)/rsq * oss) ;
    value=v1*fac + omi ;
    value_stored=true ;
    return ( v2*fac + omi ) ;
  }
}

