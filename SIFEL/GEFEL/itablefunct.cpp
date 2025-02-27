/*
    File:                    itablefunct.cpp
    Author:                  
    Purpose:                 Class computes function value at exact point, where the function is given by the table x,y
                             Implemented interpolations: dirac <), lagrange, piecewiselinear 

*/

#include <stdlib.h>
#include "itablefunct.h"
#include "intools.h"



/**
   The constructor defines arrays for interpolation.
   Parameters:
   @param asize is the size of vectors x, y
   @param x is the vector where time values are stored
   @param y is the vector where the function values are stored
 
   10.2.2006, JK, TKo
*/  
itablefunct::itablefunct()
{
  asize=0; x=NULL; y=NULL;
}



/**
   The destructor deallocates the memory where the vectors x, y are stored
   Parameters:

   10.2.2006, JK, TKo
*/
itablefunct::~itablefunct()
{
  delete[] x; delete[] y;
}



/**
   function reads data
   
   Example:
   
   function is equal to y1 for x from -infinity to x less than x1
   function is equal to y2 for x from x1 to x less than x2
   function is equal to y3 for x from x2 to x less than x3
   function is equal to y4 for x from x3 to infinity
   
   asize is equal to 4 for this example
   
   
   @param in is input file
   
   10.2.2006, JK, TKo
*/
void itablefunct::read(XFILE *in)
{
  long i;
  xfscanf (in,"%k%ld","nitab_items",&asize);
  x=new double[asize];
  y=new long[asize];
  for (i=0;i<asize;i++){
    xfscanf (in,"%lf",&x[i]);
    xfscanf (in,"%ld",&y[i]);
  }
}



/**
   The function writes input parameters (asize, x, y) for interpolation

   Parameters: 
   @param out - output file
   
   10.2. 2006, JK, TKo
*/
void itablefunct::print(FILE *out)
{
  long i;
  fprintf (out,"\n%ld",asize);
  for (i=0;i<asize;i++){
    fprintf (out,"\n%le ",x[i]);
    fprintf (out,"%ld",y[i]);
  }
}



/**
   Parameters:
   @param in is input file

   10.2.2006, TKo
*/
void itablefunct::read_prop(FILE *in)
{
  long i;
  getlong(in, asize);
  x=new double[asize];
  y=new long[asize];

  for (i=0;i<asize;i++)
  {
    getdouble (in, x[i]);
    getlong (in, y[i]);
  }
}



/**
   The function returns the interpolated value at exact point (temp) by linear interpolation,
   if the point lies out of the range, extrapolation from the limit intervals takes effect.

   Parameters:
   @param temp is x-coordinate, where the function value will be found
   
   10.20.2006, JK, TKo
   TKo+JK 4.7.2013 changed evaluation of piecewise const function 
                   from (x_i; x_{i+1}> => y[i+1]
                   to   <x_i; x_{i+1}) => y[i]
*/
long itablefunct::getval (double temp)
{
  long i,ret;
  
  if (temp<x[0]){
    return y[0];
  }
  if (x[asize-1]<=temp){
    return y[asize-1];
  }
  
  ret = 0;
  for (i=0;i<asize-1;i++){
    if (x[i]<=temp && temp<x[i+1]){
      ret = y[i];
      break;
    }
  }
  
  return ret;
}



/**
   The function copies all data from subroutine argument to the actual object.
   
   @param tf - object of tablefunct for comparison

   @return The function returns copy of the actual object in the parameter tf.

   JK, 15.11.2008
*/
void itablefunct::copy (itablefunct &tf)
{
  long i;
  
  //  number of function values
  asize = tf.asize;
  
  if (x!=NULL)
    delete [] x;
  x = new double [asize];
  
  if (y!=NULL)
    delete [] y;
  y = new long [asize];
  
  for (i=0;i<asize;i++){
    x[i]=tf.x[i];
    y[i]=tf.y[i];
  }
}



/**
   The function compares all data from subroutine argument with the actual object
   
   @param tf - object of tablefunct for comparison

   @retval 0 - objects are identical
   @retval 1 - objects differs  

   TKo, 09.2010
*/
long itablefunct::compare(itablefunct &tf)
{
  long i;
  
  //  check for number of function values
  if (asize != tf.asize)
    return 1;
  
  for (i=0;i<asize;i++){
    if (x[i] != tf.x[i])
      return 1;
    if (y[i] != tf.y[i])
      return 1;
  }
  return 0;
}
