#include <string.h>
#include <math.h>
#include <stdio.h>
#include "outdriverc.h"
#include "couptop.h"
#include "probdescc.h"
#include "intools.h"
#include "globalc.h"

/**
  Constructor initializes data members to zero values.
*/
outdriverc::outdriverc()
{
  outfn[0] = 0x0;
}



/**
  Destructor deallocates used memory.
*/
outdriverc::~outdriverc()
{
}



/**
  Function reads description of required output from the opened text file in.

  @param in - pointer to opened text file

  @retval 0 - on success
  @retval 1 - error reading output text filename
  @retval 2 - error reading nodal output description
  @retval 3 - error reading element output description
  @retval 4 - error reading user defined point output description
  @retval 5 - error reading output graphics description
  @retval 6 - error reading output diagram description
*/
long outdriverc::read(XFILE */*in*/)
{
  // xfscanf(in, " %1000a", outfn);
  return 0;
}


/**
  Function prints data with output description to the text file given by out.
  
  @param out - pointer to opened text file for output    
*/
void outdriverc::print(FILE */*out*/)
{
  //  fprintf(out, "\n%s\n", outfn);
}


/**
  Function prints header to the output text file.
  
  @param out - pointer to the opened text output file
*/
void outdriverc::print_header(FILE *out)
{
  fprintf(out, "%15s ****  *  ****  ****  *\n", " ");
  fprintf(out, "%15s *     *  *     *     *\n", " ");
  fprintf(out, "%15s  *    *  ***   ***   *\n", " ");
  fprintf(out, "%15s   *   *  *     *     *\n", " ");
  fprintf(out, "%15s****   *  *     ****  ****  METR OUTPUT\n", " ");

  fprintf(out, "\n%s\n", Cp->name);
  fprintf(out, "\n\n\n\n\n");
}



/**
  Function prints step number to the output text file.
  
  @param out  - pointer to the opened text output file
  @param step - integer step id
  @param time - time or load step
*/
void outdriverc::print_newstep(FILE *out, long step, double time)
{
  long i;
  for (i=0; i<53; i++)
    fprintf(out, "*");
  fprintf(out, "\n%10sStep number=% ld, time/load step=% g\n", " ", step, time);
  for (i=0; i<53; i++)
    fprintf(out, "*");
  fprintf(out, "\n\n\n\n");
}


/**
  Function prints required output values to the output text file.
  
  @param out - pointer to the opened text file 
  @param lcid - load case id
  @param istep - step id
  
*/
void outdriverc::print_out(FILE */*out*/, long /*lcid*/, long /*istep*/, double /*time*/)
{
}



/**
  Function prints diagrams.
  
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
*/
void outdriverc::print_diags(long /*lcid*/, double /*lambda*/, long /*istep*/, double */*fi*/)
{
}



/**
  Function prints required value to the graphics files.
  
  @param out - pointer to the opened text file 
  @param lcid - required load case id
  @param lambda - actual load parameter/actual time
  @param istep - actual step id
  @param fi - array with values of load vector
  
*/
void outdriverc::print_graphics(FILE */*out*/, long /*lcid*/, double /*lambda*/, long /*istep*/, double */*fi*/)
{
}
