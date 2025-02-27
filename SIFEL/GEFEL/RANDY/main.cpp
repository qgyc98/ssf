/*
   Description:  test file for random generator
   
   Author:       Matej Leps, leps@cml.fsv.cvut.cz
    
   $Id: main.cpp,v 1.2 2004/02/18 10:00:16 leps Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include "randy.h"


randy RND ;

//
// For printing basic data
//
void print_basic_data ( void ){
  printf( "RAND_MAX = %d \n", RAND_MAX ) ;
}

//
// For testing boundaries, is better to add -D__TEST_BOUNDS__ into Makefile
//
void test_double ( void )
{
  const double eps = NearZero ;
  const long max_iter = 4 ;
  long i ;

  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double() < eps )
	{
	  printf ( " < 0 ; 1 > : Zero was reached.    \t\t OK. \n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < 0 ; 1 > : Zero was not reached. \t\t BAD!!!!!!\n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (1.0-RND.give_double()) < eps )
	{
	  printf ( " < 0 ; 1 > : One was reached.    \t\t OK. \n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < 0 ; 1 > : One was not reached. \t\t BAD!!!!!!\n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double(0L,1L) < eps )
	{
	  printf ( " < 0 ; 1 ) : Zero was reached.    \t\t OK. \n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < 0 ; 1 ) : Zero was not reached. \t\t BAD!!!!!!\n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (1.0-RND.give_double(0L,1L)) < eps )
	{
	  printf ( " < 0 ; 1 ) : One was reached.    \t\t BAD!!!!!!\n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < 0 ; 1 ) : One was not reached. \t\t OK. \n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double(1L,0L) < eps )
	{
	  printf ( " ( 0 ; 1 > : Zero was reached.    \t\t BAD!!!!!!\n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( 0 ; 1 > : Zero was not reached. \t\t OK. \n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (1.0-RND.give_double(1L,0L)) < eps )
	{
	  printf ( " ( 0 ; 1 > : One was reached.    \t\t OK.\n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( 0 ; 1 > : One was not reached. \t\t BAD!!!!!!\n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double(1L,1L) < eps )
	{
	  printf ( " ( 0 ; 1 ) : Zero was reached.    \t\t BAD!!!!!!\n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( 0 ; 1 ) : Zero was not reached. \t\t OK. \n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (1.0-RND.give_double(1L,1L)) < eps )
	{
	  printf ( " ( 0 ; 1 ) : One was reached.    \t\t BAD!!!!!!\n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( 0 ; 1 ) : One was not reached. \t\t OK. \n" ) ;
}

//
// For testing boundaries, is better to add -D__TEST_BOUNDS__ into Makefile
//
void test_double_from ( void )
{
  const double LOW=0. ;
  const double HIGH=1. ;

  const double eps = NearZero ;
  const long max_iter = 4 ;
  long i ;

  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double( LOW, HIGH, 0L, 0L ) - LOW < eps )
	{
	  printf ( " < %10.5f ; %10.5f > : %10.5f was reached.    \t\t OK. \n",
		   LOW, HIGH, LOW ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < %10.5f ; %10.5f > : %10.5f was not reached. \t\t BAD!!!!!!\n",
		   LOW, HIGH, LOW ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (HIGH-RND.give_double( LOW, HIGH, 0L, 0L )) < eps )
	{
	  printf ( " < %10.5f ; %10.5f > : %10.5f was reached.   \t\t OK. \n",
		   LOW, HIGH, HIGH ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < %10.5f ; %10.5f > : %10.5f was not reached.\t\t BAD!!!!!!\n",
		   LOW, HIGH, HIGH) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double( LOW, HIGH, 0L,1L)- LOW < eps )
	{
	  printf ( " < %10.5f ; %10.5f ) : %10.5f was reached.    \t\t OK. \n",
		   LOW, HIGH, LOW ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < %10.5f ; %10.5f ) : %10.5f was not reached. \t\t BAD!!!!!!\n",
		   LOW, HIGH, LOW ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (HIGH-RND.give_double( LOW, HIGH, 0L,1L)) < eps )
	{
	  printf ( " < %10.5f ; %10.5f ) : %10.5f was reached.   \t\t BAD!!!!!!\n",
		   LOW, HIGH, HIGH ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " < %10.5f ; %10.5f ) : %10.5f was not reached.\t\t OK. \n",
		   LOW, HIGH, HIGH ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double( LOW, HIGH, 1L,0L)- LOW < eps )
	{
	  printf ( " ( %10.5f ; %10.5f > : %10.5f was reached.    \t\t BAD!!!!!!\n",
		   LOW, HIGH, LOW ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( %10.5f ; %10.5f > : %10.5f was not reached. \t\t OK. \n",
		   LOW, HIGH, LOW ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (HIGH-RND.give_double( LOW, HIGH, 1L,0L)) < eps )
	{
	  printf ( " ( %10.5f ; %10.5f > : %10.5f was reached.   \t\t OK.\n",
		   LOW, HIGH, HIGH ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( %10.5f ; %10.5f > : %10.5f was not reached.\t\t BAD!!!!!!\n",
		   LOW, HIGH, HIGH ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_double( LOW, HIGH, 1L,1L) - LOW < eps )
	{
	  printf ( " ( %10.5f ; %10.5f ) : %10.5f was reached.    \t\t BAD!!!!!!\n",
		   LOW, HIGH, LOW ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( %10.5f ; %10.5f ) : %10.5f was not reached. \t\t OK. \n",
		   LOW, HIGH, LOW ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (HIGH-RND.give_double( LOW, HIGH, 1L,1L)) < eps )
	{
	  printf ( " ( %10.5f ; %10.5f ) : %10.5f was reached.   \t\t BAD!!!!!!\n",
		   LOW, HIGH, HIGH ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " ( %10.5f ; %10.5f ) : %10.5f was not reached.\t\t OK. \n",
		   LOW, HIGH, HIGH ) ;
}

void test_long ( void )
{
  const long NUMBER=5 ;

  const long max_iter = 4 ;
  long i ;


  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_long( NUMBER ) == 0 )
	{
	  printf ( "    0 was reached.    \t\t OK.\n" ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( "    0 was not reached. \t\t BAD!!!!\n" ) ;
  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_long( NUMBER ) == NUMBER )
	{
	  printf ( " %4ld was reached.    \t\t BAD!!!!\n", NUMBER ) ;
	  break ;
	}
    }
  if ( i==max_iter ) printf ( " %4ld was not reached. \t\t OK.\n", NUMBER ) ;
} 

void test_long_from ( void )
{
  const long LOW=2 ;
  const long HIGH=4 ;


  const long max_iter = 4 ;
  long i ;

  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( RND.give_long( LOW, HIGH ) - LOW == 0 )
	{
	  printf ( " < %ld ; %ld > : %ld was reached.    \t\t OK. \n",
		   LOW, HIGH, LOW) ;
	  break ;
	}
    }
  if ( i==max_iter ) 
    printf ( " < %ld ; %ld > : %ld was not reached. \t\t BAD!!!!!! \n",
	     LOW, HIGH, LOW) ;

  for ( i=0 ; i<max_iter ; i++ )
    {
      if ( (HIGH - RND.give_long( LOW, HIGH )) == 0 )
	{
	  printf ( " < %ld ; %ld > : %ld was reached.    \t\t OK. \n",
		   LOW, HIGH, HIGH) ;

	  break ;
	}
    }
  if ( i==max_iter ) 
    printf ( " < %ld ; %ld > : %ld was not reached. \t\t BAD!!!!!! \n",
	     LOW, HIGH, HIGH) ;

} 



void test_gauss_double ( void )
{
  FILE *f ;
  f=fopen ( "random.out", "at" ) ;
  
  for ( long i=0 ; i<100000 ; i++ )
  {
    fprintf ( f, "%f \n", RND.give_gauss_double() ) ;
  }
  
  fclose ( f ) ;
}




void test_gauss_long ( void )
{
  printf ( "1-1 %ld \n", RND.give_long( 1,1 ) ) ;
  printf ( "1-2 %ld \n", RND.give_long( 1,2 ) ) ;
  //printf ( "1-0 %ld \n", RND.give_long( 1,0 ) ) ;
} 

int main ( int , char *[] )
{

  RND.init() ;

  print_basic_data() ;
  test_double() ;
  test_double_from() ;
  test_long() ;
  test_long_from() ;
  //test_gauss_long() ;
  //test_gauss_double() ;


  printf( "** everything done ** \n" ) ;
  return 0 ;
}
