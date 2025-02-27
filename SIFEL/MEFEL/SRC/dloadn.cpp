#include "dloadn.h"
#include "iotools.h"
#include "gfunct.h"
#include "global.h"
#include "mechtop.h"
#include <string.h>
#include <stdlib.h>



/**
  The constructor initializes class attributes to the zero values

  Created by Tomas Koudelka
*/
dloadn::dloadn()
{
  idn = 0L;
  gf = NULL;
}



/**
  The destructor deallocates used memory.

  Created by Tomas Koudelka
*/
dloadn::~dloadn()
{
  delete [] gf;  
}



/**
  The function reads data for the dynamical nodal load from the opened text file
  given by the parameter in.

  @param in - pointer to the opened text file

  @retval 0 - on success.
  @retval 1 - error parsing string expression

  Created by JK
  Modified by TKo
  Rewritten by Tomas Koudelka, 7.7.2014
*/
long dloadn::read(XFILE *in)
{
  long i,ndofn;

  xfscanf(in, "%ld",&idn);
  if ((idn<1) || (idn>Mt->nn))
  {
    print_err("number of loaded node is out of range <1,%ld>", __FILE__, __LINE__, __func__, Mt->nn);
    abort();
  }
  idn--;
  ndofn = Mt->give_ndofn (idn);
  gf = new gfunct[ndofn];

  for (i=0; i<ndofn; i++)
    gf[i].read(in);

  return(0);
}


/**
  The function prints data for the dynamical nodal load into the opened text file
  given by the parameter out.

  @param out - pointer to the opened text file

  @retval 0 - on success.
  @retval 1 - error parsing string expression

  TKr, 08/02/2013 according to read(XFILE *in)
  Rewritten by Tomas Koudelka, 7.7.2014
*/
long dloadn::print(FILE *out)
{
  long i,ndofn;

  fprintf(out, "\n\n  %ld",idn+1);

  ndofn = Mt->give_ndofn (idn);

  for(i=0; i<ndofn; i++)
    gf[i].print(out);

  return(0);
}



/**
  The function return value of the dynamical nodal load for the time given by the parameter t.

  @param t - required time.
  @param id - required load component id

  @return The function returns value of the dynamical nodal load for the given time.

  Created by JK, Tomas Koudelka,
  Rewritten by Tomas Koudelka, 7.7.2014
*/
double dloadn::getval(double t, long id)
{
  long ndofn = Mt->give_ndofn(idn);
  double ret;

  if ((id < 0) || (id >= ndofn))
  {
    print_err("required load component id=%ld is out of range <0;%ld>", 
              __FILE__, __LINE__, __func__, id, ndofn);
    abort();
  }

  ret = gf[id].getval(t);
  
  return(ret);
}
