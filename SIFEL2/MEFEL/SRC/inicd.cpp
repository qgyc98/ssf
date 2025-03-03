#include <string.h>
#include "inicd.h"



/**
  The constructor inializes attributes to zero values.

  Created by Tomas Koudelka
*/
inicd::inicd (void)
{
  val = NULL;
  nval = 0L;
  type = none;
}



/**
  The destructor deallocates used memory.

  Created by Tomas Koudelka
*/
inicd::~inicd(void)
{
  delete [] val;
}



/**
  The function reads initial conditions form the opened text file.

  @param in - pointer to the opened text file

  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void inicd :: read(XFILE *in)
{
  // in the following format, the %d must be used 
  // in order to allow for combinations of values of the inictype
  xfscanf(in, "%k%+m", "ini_cd_type" , &inictype_kwdset, &type);
  xfscanf(in, "%k%ld", "nval", &nval);

  val = new double[nval];
  memset(val, 0, sizeof(*val)*nval);

  for (long i = 0; i < nval; i++)
    xfscanf(in, "%le", val+i);
}



/**
  The function prints initial conditions to the opened text file

  @param out - pointer to the opened text file
 
  @return The function does not return anything.

  Created by Tomas Koudelka
*/
void inicd::print(FILE *out)
{
  fprintf(out, " %d", type);
  fprintf(out, " %ld", nval);
  for (long i = 0; i < nval; i++)
    fprintf(out, " %e", val[i]);
  fprintf(out, "\n");
}



/**
  The function copies initial condition form the parameter to the given object.

  @param ic - copied initial condition (output)
   
  @return The function returns copy of the object in the parameter ic.

  Created by Tomas Koudelka, 6.2009
*/
void inicd::copy(inicd &ic)
{
  long i;

  if (val != NULL)
    delete [] val;
  type = ic.type;
  nval = ic.nval;
  val = new double[nval];
  for (i=0; i<nval; i++)
    val[i] = ic.val[i];  
}



/**
  The function merges initial condition form the parameter to the given object.

  @param ic - merged initial condition
   
  @retval 0 - on success
  @retval 1 - in case of incompatible conditions

  Created by Tomas Koudelka, 6.2009
*/
long inicd::merge(inicd &ic)
{
  long i;
  if (ic.type == none)
    return 0; 
  if (type == none)
  {
    type = ic.type;
    nval = ic.nval;
    val = new double[nval];
    for (i=0; i<nval; i++)
      val[i] = ic.val[i];
    return 0;
  }
  if ((ic.type == type) && (ic.nval == nval))
  {
    // in the case of the initial condition of the same type
    for (i=0; i<nval; i++)
    {
      // check that values of initial conditions are identical
      if (val[i] == ic.val[i])
        continue;
      else
        return 1;
    }
    return 0;
  }
  // initial condition types are not identical => error
  return 1;
}
