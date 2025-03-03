#include "dloadpd.h"
#include "global.h"
#include "intools.h"



/**
  The constructor initializes class attributes to the zero values.

  Created by Tomas Koudelka,
*/
dloadpd::dloadpd()
{
}



/**
  The destructor deallocates used memory.

  Created by Tomas Koudelka,
*/
dloadpd::~dloadpd()
{
}



/**
  The function reads data for the dynamical prescribed displacements from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @retval 0 - on success.
  @retval 1 - error parsing string expression

  Created by Tomas Koudelka,
  Rewritten by Tomas Koudelka, 7.7.2014
*/
long dloadpd::read(XFILE *in)
{
  return gf.read(in);
}


/**
  The function prints data for the dynamical prescribed displacements into the opened text file
  given by the parameter out.

  @param out - pointer to the opened text file

  @retval 0 - on success.
  @retval 1 - error parsing string expression

  TKr, 08/02/2013 according to read(XFILE *in)
  Rewritten by Tomas Koudelka, 7.7.2014
*/
long dloadpd::print(FILE *out)
{
  gf.print(out);
  return(0);
}



/**
  The function return value of the dynamical prescribed displacement
  for the time given by the parameter t.

  @param t - desired value of time.

  @return The function returns value of the dynamical prescribed displacement for the given time.

  Created by Tomas Koudelka,
  Rewritten by Tomas Koudelka, 7.7.2014
*/
double dloadpd::getval(double t)
{
  double ret;
  ret = gf.getval(t);
  return(ret);
}
