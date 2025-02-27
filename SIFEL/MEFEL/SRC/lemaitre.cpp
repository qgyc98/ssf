#include <math.h>
#include "lemaitre.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "matrix.h"
#include "vector.h"
#include "intpoints.h"



/**
  Constructor initializes data members to zero or default values.

  Created by JK, MM
*/
lemaitre::lemaitre (void)
{
  eta=0.0;  m = 0.0;  n = 0.0;
}



/**
  Destructor is defined only for the formal purposes.

  Created by JK, MM
*/
lemaitre::~lemaitre (void)
{

}



/**
  Function reads material parameters from the opened text file.
   
  @param in - pointer to the opened XFILE

  @return The function does not return anything.

  Created by JK, MM
*/
void lemaitre::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&eta,&m,&n);
}



/**
  The function prints model parameters in to the opened text file.

  @param out - pointer to the opened text file for output

  @return The function does not return anything.

  Created by TKo, 4.2016
*/
void lemaitre::print (FILE *out)
{
  fprintf (out, "%lf %lf %lf", eta, m, n);
}



/**
  The function returns value of viscous function g with respect to the
  actual value of yield function f and cumulativ strain.

  @param f - actual value of yield function with respect to attained stress state
  @param cs - cumulative strain

  @return The function returns computed value of viscous function.

  Created by JK, MM
*/
double lemaitre::gfun (double f,double cs)
  //  ramp function
  //  f - yield function
  //  cs - cumulative irreversible strain
  //  28.10.2001
{
  double g;
  
  if (f>0.0){
    if (cs<Mp->zero)   g=eta*pow(f,n);
    else               g=eta*pow(cs,m)*pow(f,n);
  }
  else{
    g=0.0;
  }
  
  return g;
}

/*
void lemaitre::nlstresses (long ipp,long ido)
  //  29.10.2001
{
  long i,ncomp=Mm->ip[ipp].ncompstr,ncomph=1;
  double dt;
  vector sig(ncomp),q(ncomph),epsn(ncomp),epso(ncomp),epsp(ncomp);
  matrix d(ncomp,ncomp);

  if (Mp->phase==1){
    // *********************************
    //  right hand side computation  //
// ********************************* 
    //  stress component
    for (i=0;i<ncomp;i++){
      sig[i]=Mm->ip[ipp].other[ido+i];
    }
    //  internal variable
    q[0]=Mm->ip[ipp].other[ido+3*ncomp+1];
    
    //dt=Mp->incrvs;
    //dt=Mp->timefun.getval(Mp->time);
    dt=Mp->timecon.actualforwtimeincr ();
    Mm->stiff_deps_vispl (ipp,0,0,sig,q,dt);
  }
  
  
  if (Mp->phase==2){
    for (i=0;i<ncomp;i++){
      //  new total strain
      epsn[i] = Mm->ip[ipp].strain[i];
      //  previous total strain
      epso[i] = Mm->ip[ipp].other[ido+i+ncomp*2];
      //  irreversible increment of strain
      epsp[i] = Mm->ip[ipp].other[ido+i+ncomp];
    }
    
    //  internal variable
    //q[0] = Mm->ip[ipp].other[ido+3];
    
    //  elastic strain
    subv (epsn,epso,epso);
    subv (epso,epsp,epsp);
    //  stiffness matrix of material
    Mm->matstiff(d,ipp);
    //  stress components
    mxv (d,epsp,sig);
    
    //  new data storage
    for (i=0;i<ncomp;i++){
      Mm->ip[ipp].other[ido+i] += sig[i];
      Mm->ip[ipp].other[ido+i+ncomp*2] = epsn[i];
    }
  }

}
*/

/**
   function returns tangent stiffness matrix of material
   
   @param ipp - integration point pointer
   @param d - elastic stiffness matrix

   12.8.2001
*/
/*
void lemaitre::matstiff (matrix &d,long ipp,long ido)
{
  Mm->elmatstiff (d,ipp);
}

*/
