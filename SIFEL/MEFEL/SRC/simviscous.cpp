#include <math.h>
#include "simviscous.h"
#include "global.h"
#include "matrix.h"
#include "vector.h"
#include "intpoints.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK,
*/
simviscous::simviscous (void)
{
  eta=0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK,
*/
simviscous::~simviscous (void)
{

}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by JK,
*/
void simviscous::read (XFILE *in)
{
  xfscanf (in,"%k %lf", "eta", &eta);
}



/**
  The function prints model parameters in to the opened text file.

  @param out - pointer to the opened text file for output

  @return The function does not return anything.

  Created by TKo, 4.2016
*/
void simviscous::print (FILE *out)
{
  fprintf (out, "%le", eta);
}



/**
  The function returns value of viscous function g with respect to the
  actual value of yield function f.

  @param f - actual value of yield function with respect to attained stress state

  @return The function returns computed value of viscous function.

  Created by JK,
*/
double simviscous::gfun (double f)
{
  double g;
  if (f>0.0)  g=f*eta;
  else        g=0.0;
  return g;
}



/**
  The function returns value of derivative of viscous function g with respect to the
  actual value of yield function f, i.e. dg/df.

  @param f - actual value of yield function with respect to attained stress state

  @return The function returns computed value of derivetive of viscous function dg/df.

  Created by JK,
*/
double simviscous::dergfun (double f)
{
  double dg;
  if (f<0.0)  dg=0.0;
  else        dg=eta;
  return dg;
}

