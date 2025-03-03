#include "relaxeuroc.h"
#include "stochdriver.h"
#include "global.h"
#include "probdesc.h"


/**
  Constructor initializes data members to zero or default values.

  Created by JK, 11. 6. 2013
*/
relaxeuroc::relaxeuroc (void)
{
  //  initial prestress \sigma_{pm0}
  siginit=0.0;
  //  characteristic strength f_{pk}
  fpk=0.0;
  //  model coefficient A (determined by regression)
  a=0.0;
  //  model coefficient B (determined by regression)
  b=0.0;
  //  Young modulus of elasticity
  e=0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK, 11. 6. 2013
*/
relaxeuroc::~relaxeuroc (void)
{

}



/**
   Function reads material parameters from the opened text file.
   
   @param in - pointer to the opened XFILE
   
   @return The function does not return anything.
   
   Created by JK, 11. 6. 2013
*/
void relaxeuroc::read (XFILE *in)
{
  //  initial prestress \sigma_{pm0}
  //  characteristic strength f_{pk}
  //  model coefficient A (determined by regression)
  //  model coefficient B (determined by regression)
  //  Young modulus of elasticity
  xfscanf (in,"%k%lf %k%lf %k%lf %k%lf %k%lf", 
           "sig_ini", &siginit, "f_pk", &fpk, "a", &a, "b", &b, "e", &e);
}

/**
   Function prints material parameters into the opened text file.
   
   @param out - pointer to the opened FILE
   
   @return The function does not return anything.
   
   Created by JK,  11. 6. 2013
*/
void relaxeuroc::print (FILE *out)
{
  //  initial prestress \sigma_{pm0}
  //  characteristic strength f_{pk}
  //  model coefficient A (determined by regression)
  //  model coefficient B (determined by regression)
  //  Young modulus of elasticity
  fprintf (out,"%le %le %le %le %le",siginit,fpk,a,b,e);
}

/**
   function computes stress decrement
   time has to be measured in hours
   
   @param time - actual time
   
   JK, 11. 6. 2013
*/
double relaxeuroc::stress_decrement (void)
{
  double c,d,f,time;
  
  // actual time
  time = Mp->time/3600;

  //  auxiliary coefficients
  c=b*siginit/fpk;
  d=0.75*(1.0-siginit/fpk);
  
  f = siginit*a*exp(c)*pow(time/1000.0,d)/100000.0;

  return f;
}

/**
   function computes stress decrement
   time has to be measured in hours
   
   @param time - actual time
   
   JK, 11. 6. 2013
*/
void relaxeuroc::stress (vector &sig,vector &eps,strastrestate ssst)
{
  double c,d,f,time;
  
  // actual time
  time = Mp->time/3600.0;
  
  //  auxiliary coefficients
  c=b*siginit/fpk;
  d=0.75*(1.0-siginit/fpk);
  
  f = -(siginit - siginit*a*exp(c)*pow(time/1000.0,d)/100000.0); // this is contribution to the right hand side
  
  switch (ssst){
  case bar:{
    sig[0]=f;
    eps[0]=f/e; // eigenstrain
    break;
  }
  case spacestress:{
    sig[0]=f;
    sig[1]=f;
    sig[2]=f;
    sig[3]=0.0;
    sig[4]=0.0;
    sig[5]=0.0;

    // eigenstrain components
    eps[0]=f/e;
    eps[1]=f/e;
    eps[2]=f/e;
    eps[3]=0.0;
    eps[4]=0.0;
    eps[5]=0.0;
    break;
  }
  default:{
    print_err("unknown number of components of stress tensor is required", __FILE__, __LINE__, __func__);
  }
  }

}
