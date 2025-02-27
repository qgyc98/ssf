#include "boermat.h"
#include "iotools.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "mechcrsec.h"
#include "mechbclc.h"
#include "intpoints.h"
#include "matrix.h"
#include "vecttens.h"
#include "vector.h"
#include <math.h>



#define nijac 20
#define limit 1.0e-8

/**
  This constructor inializes attributes to zero values; exponent n
  is initialized with Groen constant.
*/
boermat::boermat (void)
{
  phi=0.0;  c=0.0;  psi=0.0;
  n=1.0/0.229; // Groen constant

}

/**
  This destructor is only for the formal purposes.
*/
boermat::~boermat (void)
{

}


/**
  This function reads material parameters from the opened text file given
  by the parameter in. Then it computes material constants alpha, alpha1, alpha2,
  beta and beta1

  @param in - pointer to the opned text file

*/
void boermat::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&phi,&c,&psi);
  sra.read (in);
  
  a=pow((3.0+sin(phi))/(3.0-sin(phi)),n);
  delta=(a-1.0)/(a+1.0);

  alpha=6.0*sin(phi)/(3.0-sin(phi));
  alpha1=alpha*pow(1.0-delta,1.0/n);
  alpha2=6.0*sin(psi)/(3.0-sin(psi));

  beta=6.0*cos(phi)/(3.0-sin(phi));
  beta1=beta*pow(1.0-delta,1.0/n);

}



/**
   This function computes the value of yield functions.

   @param sig - stress tensor

   @retval The function returns value of yield function for the given stress tensor

   25.3.2002
*/
double boermat::yieldfunction (vector &sig)
{
  double i1s,i3s,j2s,a1,a2,f,sin3m;

  i1s = first_invar (sig)/3.0;

  i3s = third_stress_invar (sig);
  j2s = j2_stress_invar (sig);

  if (fabs(j2s)>Mp->zero)
    sin3m = - 3.0*sqrt(3.0)/2.0 * i3s/sqrt(j2s*j2s*j2s);
  else
    sin3m = 0.0;

  a1 = sqrt(3.0*j2s) * pow((1.0-delta*sin3m), 1.0/n);
  a2 = alpha1 * i1s;

  f = a1 + a2 - beta1 * c;

  return f;
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector sigma.

   @param sig  - stress vector in Voigt notation
   @param dfds - %matrix where the resulting derivatives are stored


   4.1.2002
*/
void boermat::deryieldfsigma (vector &sig, vector &dfds)
{
  double j2s,j3s,sin3m, a1, a2;
  vector dev(ASTCKVEC(6));
  vector da1ds(ASTCKVEC(6)), da2ds(ASTCKVEC(6)), da3ds(ASTCKVEC(6)), dj3sds(ASTCKVEC(6));

  deviator (sig,dev);

  j2s = second_stress_invar (dev);
  j3s = third_stress_invar (dev);
  if (fabs(j2s)>Mp->zero)
    sin3m = - 3.0*sqrt(3.0)/2.0 * j3s/sqrt(j2s*j2s*j2s);
  else
    sin3m = 0.0;
  a1 = sqrt(3.0*j2s);
  a2 = pow(1.0-delta*sin3m, 1.0/n);
  // derivative d_a1/d_sigma
  copyv(dev,da1ds);
  cmulv(0.5*sqrt(3.0/j2s), da1ds);
  // derivative d_a2/d_sigma
  //  d I3/d sigma_ij
  dj3sds(0) = (sig(1)*sig(2)-sig(3)*sig(3))*(2.0/3.0); // d/d sigma_11
  dj3sds(1) = (sig(0)*sig(2)-sig(4)*sig(4))*(2.0/3.0); // d/d sigma_22
  dj3sds(2) = (sig(0)*sig(1)-sig(5)*sig(5))*(2.0/3.0); // d/d sigma_33
  dj3sds(3) = sig(5)*sig(4)-sig(0)*sig(3); // d/d sigma_23
  dj3sds(4) = sig(5)*sig(3)-sig(1)*sig(4); // d/d sigma_13
  dj3sds(5) = sig(3)*sig(4)-sig(5)*sig(2); // d/d sigma_12

  cmulv(sqrt(j2s*j2s*j2s), dj3sds);
  copyv(dev, da2ds);
  cmulv(3.0/2.0*sqrt(j2s)*j3s, da2ds);
  subv(dj3sds, da2ds, da2ds);
  cmulv(3.0*sqrt(3.0)/2.0/(j2s*j2s*j2s), da2ds);
  cmulv(delta/n*pow(1-delta*sin3m, 1-1.0/n), da2ds);
  // derivative d_a3/d_sigma
  fillv(0.0, da3ds);
  da3ds(0) = da3ds(1) = da3ds(2) = 1.0/3.0*alpha1;

  addv(da1ds, da2ds, dfds);
  addv(dfds, da3ds, dfds);

  c = -1.0*sqrt(3.0*j2s)*(1.0/n)*pow(1.0-delta*sin3m,1.0/n-1)*delta;
  cmulv (c,dfds);

  c = sqrt(3.0)/2.0/sqrt(j2s)*pow(1.0-delta*sin3m,n);
  cmulv (c,dev);

  addv (dev,dfds,dfds);

  dfds(0) += alpha1/3.0;
  dfds(1) += alpha1/3.0;
  dfds(2) += alpha1/3.0;
  // double offdiagonal components due to Voigt notation
  dfds(3) *= 2.0;
  dfds(4) *= 2.0;
  dfds(5) *= 2.0;
}



/**
   This function computes derivatives of plastic potential function
   with respect of vector sigma.

   @param sig - stress tensor in Voigt notation
   @param dgds - %vector where the resulting derivatives are stored in Voigt notation
*/
void boermat::derpotsigma (vector &sig, vector &dgds)
{
  double c;
  vector dev(ASTCKVEC(6));

  deviator(sig, dev);
  normed_stress_tensor(dev, dgds);
  c=sqrt(6.0)/2.0;
  cmulv(c, dgds);

  dgds(0)+=alpha1/3.0;
  dgds(1)+=alpha1/3.0;
  dgds(2)+=alpha1/3.0;
  // double offdiagonal components due to Voigt notation
  dgds(3) *= 2.0;
  dgds(4) *= 2.0;
  dgds(5) *= 2.0;
}



/**
  This function computes material stiffnes %matrix.

  @param d - allocated %matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

*/
void boermat::matstiff (matrix &d, long ipp, long ido)
{
  if (Mp->nlman->stmat==0)
  {
    //  initial elastic matrix
    Mm->elmatstiff (d, ipp, ido);
  }
  if (Mp->nlman->stmat==1)
  {
    //  tangent stiffness matrix
    //  not completed
    Mm->elmatstiff (d, ipp, ido);
  }
}



/**
  This function computes stresses at given integration point ipp,
  depending on the reached strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point is actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
*/
void boermat::nlstresses (long ipp, long im, long ido)
  //
{
  long i,ni,nc=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(ASTCKVEC(nc)), epsp(ASTCKVEC(nc)), q(0);

  //  initial values
  for (i=0; i<nc; i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  gamma=Mm->ip[ipp].eqother[ido+nc];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    fprintf (stderr,"\n\n wrong type of stress return algorithm is required in nlstresses (file %s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }
  
  //  new data storage
  for (i=0;i<nc;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+nc]=gamma;

}

/**
  This function computes stresses at given integration point ipp,
  depending on the reached averaged nonlocal strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point is actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
*/
void boermat::nonloc_nlstresses (long ipp,long im, long ido)
  //
{
  long i,ni,nc=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(ASTCKVEC(nc)), epsp(ASTCKVEC(nc)), q(0);

  //  initial values
  for (i=0; i<nc; i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].nonloc[i];
  }
  gamma=Mm->ip[ipp].eqother[ido+nc];

  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    fprintf (stderr,"\n\n wrong type of stress return algorithm is required in nlstresses (file %s, line %d).\n",__FILE__,__LINE__);
    abort ();
  }
  
  //  new data storage
  for (i=0;i<nc;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+nc]=gamma;
}

/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

*/
void boermat::updateval (long ipp,long im,long ido)
{
  long i,n = Mm->givencompeqother(ipp, im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  Returns vector of irreversible strains via parameter epsp
*/
void boermat::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



/**
  This function extracts consistency parametr gamma for the reached equilibrium state
  from the integration point other array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

*/
double boermat::give_consparam (long ipp,long ido)
{
  long ncompstr;
  double gamma;

  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].other[ido+ncompstr];

  return gamma;
}

void boermat::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      phi=val[i];
      break;
    }
    case 1:{
      c=val[i];
      break;
    }
    case 2:{
      psi=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }
  
  a=pow((3.0+sin(phi))/(3.0-sin(phi)),n);
  delta=(a-1.0)/(a+1.0);

  alpha=6.0*sin(phi)/(3.0-sin(phi));
  alpha1=alpha*pow(1.0-delta,1.0/n);
  alpha2=6.0*sin(psi)/(3.0-sin(psi));

  beta=6.0*cos(phi)/(3.0-sin(phi));
  beta1=beta*pow(1.0-delta,1.0/n);
  
}
