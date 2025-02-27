#include "visplast.h"
#include "probdesc.h"
#include "mechmat.h"
#include "matrix.h"
#include "vector.h"
#include "elastisomat.h"
#include "global.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>

/**
  This constructor initializes attributes to zero values.
*/
visplast::visplast (void)
{

}

/**
  This destructor is only for the formal purposes.
*/
visplast::~visplast (void)
{

}

/**
  This function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number

*/
void visplast::matstiff (matrix &d,long ipp,long /*im*/,long /*ido*/)
{
  switch (Mp->nlman->stmat){
  case initial_stiff:{
    Mm->elmatstiff (d,ipp);
    break;
  }
  case tangent_stiff:{
    Mm->elmatstiff (d,ipp);
    break;
  }
  default:{
    print_err("unknown type of stifness matrix is required",__FILE__,__LINE__,__func__);
  }
  }
  
}



/**
   29.4.2008
*/
void visplast::givestressincr (long ipp,long ido,long fi,vector &sig)
{
  long i,j,ncomp;
  ncomp=Mm->ip[ipp].ncompstr;
  
  if (ncomp!=sig.n){
    print_err("number of components of sig is not equal to ncomp",__FILE__,__LINE__,__func__);
    abort ();
  }

  for (i=fi,j=0; j<sig.n; i++,j++){
    sig[j]=Mm->ip[ipp].eqother[ido+i];
  }
}

/**
   29.4.2008
*/
void visplast::storestressincr (long ipp,long ido,vector &sig)
{
  long i,ncomp;
  ncomp=Mm->ip[ipp].ncompstr;
  
  if (ncomp!=sig.n){
    print_err("number of components of sig is not equal to ncomp",__FILE__,__LINE__,__func__);
    abort ();
  }
  
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].eqother[ido+i]=sig[i];
  }
}

/**
   function computes increments of stress

   @param ipp - number of integration point
   @param im - id of material type
   @param ido - id of instance
   
   JK, 29.4.2008
*/
void visplast::nlstressesincr (long ipp,long im,long ido)
{
  long i,ncompo,ncomp=Mm->ip[ipp].ncompstr,ncompq=1;
  double dt;
  vector sig(ncomp),q(ncompq),epsn(ncomp),epso(ncomp),epsp(ncomp);
  matrix d(ncomp,ncomp);
  
  //  the number of components of viscous material model
  ncompo=Mm->givencompother(ipp,im+1);
  
  //  stress component
  for (i=0;i<ncomp;i++){
    sig[i]=Mm->ip[ipp].stress[i];
  }
  //  internal variable
  Mm->give_interparam (ipp,im+2,ido+ncompo,q);
  
  //  time increment
  dt=Mp->timecon.actualforwtimeincr ();
  //  computation of stress increment
  Mm->stiff_deps_vispl (ipp,im,ido,sig,q,dt);
  
  for (i=0;i<ncomp;i++){
    Mm->ip[ipp].eqother[ido+i]=sig[i];
  }
  
}

/**
   function computes stresses

   @param ipp - number of integration point
   @param ido - 

   JK, 27.10.2001
*/
void visplast::nlstresses (long ipp,long ido)
{
  long i,ncomp=Mm->ip[ipp].ncompstr,ncomph=1;
  vector sig(ncomp),sigp(ncomp),q(ncomph),epsn(ncomp),epso(ncomp),deps(ncomp);
  matrix d(ncomp,ncomp);
  
  for (i=0;i<ncomp;i++){
    //  new total strain
    epsn[i] = Mm->ip[ipp].strain[i];
    //  previous total strain
    epso[i] = Mm->ip[ipp].eqother[ido+ncomp+i];
    //  total previous stress
    sigp[i] = Mm->ip[ipp].eqother[ido+2*ncomp+i];
  }
  
  //  total strain increment
  subv (epsn,epso,deps);
  //  stiffness matrix of material
  Mm->matstiff(d,ipp);
  //  D deps (part of stress components)
  mxv (d,deps,sig);
  
  //  new data storage
  for (i=0;i<ncomp;i++){
    //  total stress
    Mm->ip[ipp].stress[i]=sig[i]+sigp[i]-Mm->ip[ipp].eqother[ido+i];
  }
}

/**
   function copies array other to the array eqother
   
   @param ipp - number of integration point
   @param im - index of material type
   @param ido - index in array other
   
*/
void visplast::updateval (long ipp, long /*im*/, long ido)
{
  long i;
  long ncomp;
  
  //  the number of components
  ncomp=Mm->ip[ipp].ncompstr;
  
  for (i=0;i<ncomp;i++){
    //  the actual total strains are stored for the next time step
    Mm->ip[ipp].eqother[ido+ncomp+i] = Mm->ip[ipp].strain[i];
    //  total stress
    Mm->ip[ipp].eqother[ido+2*ncomp+i]=Mm->ip[ipp].stress[i];
  }
}
