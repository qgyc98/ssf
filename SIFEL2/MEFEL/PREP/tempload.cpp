#include "tempload.h"
#include "intools.h"



/**
  This constructor initializes class attributes to the zero values.
*/
tempload::tempload()
{
  nlc = 0L;
  nslc = 0L;
  val = 0.0;
}



/**
  Destructor
*/
tempload::~tempload()
{
}



/**
  The function reads temperature load data from the MEFEL preprocessor file given by the
  parameter in.
  Parameters :
  @param in - pointer to the opened text file
  @param lc - total number of load cases
  @param slc - pointer to array with numbers of subload cases, default value is NULL

  Returns :
  @retval 0 - on success
  @retval 1 - invalid load case number
  @retval 2 - error reading temperature value
*/
long tempload::read(XFILE *in, long lc, long *slc=NULL)
{
  xfscanf(in, "%k%ld", "lc_id", &nlc);
  if ((nlc < 1) && (nlc > lc))
    return 1L;
  if (slc)
  {
    xfscanf(in, "%k%ld", "slc_id", &nslc);
    if ((nslc < 1) || (nslc > slc[nlc-1]))
      return(1);
  }
  xfscanf(in, "%k%le", "temperature", &val);
  return 0L;
}



/**
  The function copies temperature load from the parameter tl to the given object.

  Parameters :
  @param tl - copied tempearture load
*/
void tempload::copy(tempload &tl)
{
  nlc = tl.nlc;
  nslc = tl.nslc;
  val = tl.val;
}
