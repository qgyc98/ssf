#include "hydrationheat.h"

hydrationheat::hydrationheat (void)
{

}

hydrationheat::~hydrationheat (void)
{

}

/**
   function reads data needed for description of heat source
   
   @param in - input file
   
   JK, 29.12.2009
*/
void hydrationheat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf %lf %lf %lf %lf",&a,&b,&c,&d,&e,&f,&r);
  
  //  auxiliary variable
  aa=a*f*r;
  bb=b/pow(c*3600.0,d);
  cc=bb*d*e;
}

/**
   function prints data needed for description of heat source
   
   @param out - output file
   
   JK, 29.12.2009
*/
void hydrationheat::print (FILE *out)
{
  fprintf (out,"%lf %lf %lf %lf %lf %lf %lf\n",a,b,c,d,e,f,r);
}

/**
   function computes value for actual time
   
   @param t - time
   
   JK, 29.12.2009
*/
double hydrationheat::give_value (double t)
{
  double z,denom,nom1,nom2;
  
  denom=e+bb*pow(t,d);
  nom1=bb*pow(t,d);
  nom2=cc*pow(t,d-1);
  
  z=aa*pow(nom1/denom,f-1)*nom2/denom/denom;
  
  return z;
}



/**
  The function compares content of the actual objects with 
  the content of hh.

  @param hh - compared object of hydrationheat

  Returns:
   @retval 0 - objects are identical
   @retval 1 - objects differ

  Created by TKo 09.2010
*/
long hydrationheat::compare(hydrationheat &hh)
{
  if (a != hh.a)
    return 1;
  if (b != hh.b)
    return 1;
  if (c != hh.c)
    return 1;
  if (d != hh.d)
    return 1;
  if (e != hh.e)
    return 1;
  if (f != hh.f)
    return 1;
  if (r != hh.r)
    return 1;
  if (aa != hh.aa)
    return 1;
  if (bb != hh.bb)
    return 1;
  if (cc != hh.cc)
    return 1;

  return 0;
}
