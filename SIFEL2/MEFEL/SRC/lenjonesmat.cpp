#include "lenjonesmat.h"
#include "matrix.h"
#include "vector.h"
#include "stochdriver.h"
#include "global.h"
#include "intpoints.h"

lenjonesmat::lenjonesmat (void)
{
  eps = 0.0;  sig = 0.0;  mindist=0.0;  eqdist=0.0;
}

lenjonesmat::~lenjonesmat (void)
{

}

void lenjonesmat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf",&eps,&sig);
  
  //  computation of equilibrium distance between two atoms and for given constants of potential
  //r=sig*exp(log(sig)+1.0/6.0*log(2.0));
  eqdist = sig*1.122462048;
}

/**
   function computes magnitude of interatomic force
   
   @param r - distance between atoms
   
   JK, 19.6.2005
*/
double lenjonesmat::compute_force (double r)
{
  double f,c1,c2;
  
  if (r<mindist){
    fprintf (stderr,"\n\n distance between atoms (%le) is smaller than minimum defined distance (%le)",r,mindist);
    fprintf (stderr,"\n see file %s, line %d",__FILE__,__LINE__);
  }

  c1=24.0*eps*sig*sig*sig*sig*sig*sig/r/r/r/r/r/r/r;
  c2=2.0*sig*sig*sig*sig*sig*sig/r/r/r/r/r/r;
  
  f=c1*(1.0-c2);
  
  return f;
}

/**
   function returns equilibrium distance between two atoms for given constants of potential
   
   JK, 19.6.2005
*/
double lenjonesmat::equilib_distance ()
{
  return eqdist;
}

/**
   
*/
double lenjonesmat::first_derivative (double r)
{
  double f,c1,c2,tmp;
  
  if (r<mindist){
    fprintf (stderr,"\n\n distance between atoms (%le) is smaller than minimum defined distance (%le)",r,mindist);
    fprintf (stderr,"\n see file %s, line %d",__FILE__,__LINE__);
  }

  tmp = sig/r;
  c1=24.0*eps*tmp*tmp*tmp*tmp*tmp*tmp/r;
  c2=2.0*tmp*tmp*tmp*tmp*tmp*tmp;
  
  f=c1*(1.0-c2);
  
  return f;
}
/**
   
*/
double lenjonesmat::second_derivative (double r)
{
  double f,c1,c2,tmp;
  
  if (r<mindist){
    fprintf (stderr,"\n\n distance between atoms (%le) is smaller than minimum defined distance (%le)",r,mindist);
    fprintf (stderr,"\n see file %s, line %d",__FILE__,__LINE__);
  }

  tmp = sig/r;
  c1=tmp*tmp*tmp*tmp*tmp*tmp/r/r;
  c2=tmp*tmp*tmp*tmp*tmp*tmp;
  
  f=4.0*eps*c1*(12.0*13.0*c2-42.0);
  
  return f;
}

