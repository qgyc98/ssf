#include "sorpisohansen.h"

sorpisohansen::sorpisohansen (void)
{
  //  maximum hygroscopically bound water by adsorption
  uh=0.0;
  //  empirical fixed exponent
  n=0.0;
  //  coefficient A = u_n/u_h (u_n - non-evaporable water content)
  a=0.0;;
}
sorpisohansen::~sorpisohansen (void)
{
}

/**
   function reads material parameters
   
   @param in - input file
   
   19. 11. 2012, JK
*/
void sorpisohansen::read (XFILE *in)
{
  xfscanf(in,"%le %le %le",&uh,&a,&n);
}

/**
   function prints material parameters

   @param out - output file

   19. 11. 2012, JM
*/
void sorpisohansen::print (FILE *out)
{
  fprintf(out, "\n %le %le %le",uh,a,n);

}

/**
   function computes derivative of the Hansen sorption isotherm
   with respect to the relative humidity
   
   @param rh - the relative humidity
   
   JK, 19. 11. 2012
*/
double sorpisohansen::derivative_relhum (double rh)
{
  if (rh > 1.0)
    rh = 1.0;
  if (rh < 0.0)
    rh = 0.0;
  
  return (1000 * (uh/(a*rh*n))*pow((1-(log(rh))/a),(-1-(1/n))));
  //return ((uh/(a*rh*n))*pow((1-(log(rh))/a),(-1-(1/n))));
}

/**
   function computes value of the Hansen sorption isotherm
   
   @param rh - the relative humidity
   
   JK, 19. 11. 2012
*/
double sorpisohansen::hansen_sorption_isotherm (double rh)
{
  return (uh*pow((1-(log(rh))/a),(-1/n)));
}

/**
   Function calculates relative humidity from water content with the help
   of inverse Hansen sorption isotherm
   
   @param w - water content
   
   19. 11. 2012, JK
*/
double sorpisohansen::hansen_inverse_sorption_isotherm (double w)
{
  //if (w > w_cap)//moisture content control (transient region II)
  //w = w_h_sorp;

  // relative humidity
  return (exp(a*(1.0-pow((uh/w),n))));
}
