#include "seebeckheat.h"

seebeckheat::seebeckheat (void)
{

}

seebeckheat::~seebeckheat (void)
{

}

/**
   function reads data needed for description of heat source
   
   @param in - input file
   
   TKr+JM, 15.6.2016
*/
void seebeckheat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf",&as,&sigma);

}

/**
   function prints data needed for description of heat source
   
   @param out - output file\
   
   TKr+JM, 15.6.2016
*/
void seebeckheat::print (FILE *out)
{
  fprintf (out,"%lf %lf\n",as,sigma);
}

/**
   function computes value for actual time
   
   @param t - time
   
   TKr+JM, 15.6.2016
*/
double seebeckheat::give_value ()
{
  double z;
  
  z = as*as*sigma;
  
  return z;
}

