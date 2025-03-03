#include "mohrcparab.h"
#include "global.h"
#include "mechmat.h"
#include "probdesc.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "intpoints.h"
#include "vecttens.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define nijac 20
#define limit 1.0e-8



/**
  The constructor inializes attributes to zero values.
*/
mohrcoulombpar::mohrcoulombpar (void)
{
  phi=0.0;  c=0.0;  psi=0.0;
  alphaphi = 0.0; betaphi = 0.0; sigcphi = 0.0;
  alphapsi = 0.0; betapsi = 0.0; sigcpsi = 0.0;
}



/**
  The destructor is defined only for the formal purposes.
*/
mohrcoulombpar::~mohrcoulombpar (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in. Internal parameters alpha, beta and sigc
  for phi and psi are evaluated.

  @param in - pointer to the opned text file

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mohrcoulombpar::read (XFILE *in)
{
  xfscanf (in,"%lf %lf %lf",&phi,&c,&psi);
  sra.read (in);

  alphaphi = sin(phi)/sqrt(9+3*sin(phi)*sin(phi));
  betaphi  = alphaphi / sqrt(3*(3*c*c - alphaphi*alphaphi));
  sigcphi  = sqrt(3*(c*c - alphaphi*alphaphi/3.0));

  alphapsi = sin(psi)/sqrt(9+3*sin(psi)*sin(psi));
  betapsi  = alphapsi / sqrt(3*(3*c*c - alphapsi*alphapsi));
  sigcpsi  = sqrt(3*(c*c - alphapsi*alphapsi/3.0));
}



/**
  The function computes the value of yield functions.

  @param sig - stress tensor

  @return The function returns value of yield function for the given stress tensor.

  Created by Tomas Koudelka, 25.3.2002
*/
double mohrcoulombpar::yieldfunction (vector &sig)
{
  double f, tmp;
  double j2, i1;

  i1 = first_invar(sig);
  j2 = j2_stress_invar (sig);
  tmp = 3 * j2 + sqrt(3.0) * betaphi * sigcphi * i1;
  f = sqrt(tmp) - sigcphi;

  return f;
}



/**
  The function computes derivatives of as-th yield function
  with respect of %vector sigma.

  @param sig - stress components in Voigt notation
  @param dfds - %vector of derivatives of the yield function (output)

  @return The function returns resulting derivatives in the parameter dfds.

  Created by Tomas Koudelka, 4.1.2002
*/
void mohrcoulombpar::deryieldfsigma (vector &sig, vector &dfds)
{
  double j2, i1, k1, k2;
  vector dev(6);

  deviator(sig, dev);
  i1 = first_invar(sig);
  j2 = j2_stress_invar(sig);
  k1 = 0.5 / (sqrt(j2 + sqrt(3.0)*betaphi*sigcphi*i1));
  k2 = sqrt(3.0) * betaphi * sigcphi;
  dfds(0) = k1 * (3.0 * dev(0) + k2);
  dfds(1) = k1 * (3.0 * dev(1) + k2);
  dfds(2) = k1 * (3.0 * dev(2) + k2);
  // offdiagonal components must be scaled by 2 due to Voigt notation
  dfds(3) = 2.0 * k1 * (3.0 * dev(3));
  dfds(4) = 2.0 * k1 * (3.0 * dev(4));
  dfds(5) = 2.0 * k1 * (3.0 * dev(5));
  return;
}



/**
  The function computes derivatives of plastic potential function
  with respect of stress tensor sigma.

  @param sig - stress tensor
  @param dgds - %matrix of derivatives of the plastic potential function (output)
  
  @return The function returns resulting derivatives in the parameter dgds.

  Created by Tomas Koudelka, 4.1.2002
*/
void mohrcoulombpar::derplaspotsigma (vector &sig, vector &dfds)
{
  double j2, i1, k1, k2;
  vector dev(6);

  deviator(sig, dev);
  i1 = first_invar(sig);
  j2 = j2_stress_invar(sig);
  k1 = 0.5 / (sqrt(j2 + sqrt(3.0)*betapsi*sigcpsi*i1));
  k2 = sqrt(3.0) * betapsi * sigcpsi;
  dfds(0) = k1 * (3.0 * dev(0) + k2);
  dfds(1) = k1 * (3.0 * dev(1) + k2);
  dfds(2) = k1 * (3.0 * dev(2) + k2);
  // offdiagonal components must be scaled by 2 due to Voigt notation
  dfds(3) = 2.0 * k1 * (3.0 * dev(3));
  dfds(4) = 2.0 * k1 * (3.0 * dev(4));
  dfds(5) = 2.0 * k1 * (3.0 * dev(5));
  return;
}



/**
  The function computes material stiffnes matrix.

  @param d - allocated matrix structure for material stiffness %matrix (output)
  @param ipp - integration point number

  @return The function returns resulting stiffnes %matrix in the parameter d.

  Created by Tomas Koudelka, 4.1.2002
*/
void mohrcoulombpar::matstiff (matrix &d, long ipp, long /*ido*/)
{
  if (Mp->nlman->stmat==0)
  {
    //  initial elastic matrix
    Mm->elmatstiff (d,ipp);
  }
  if (Mp->nlman->stmat==1)
  {
    //  tangent stiffness matrix
    //  nedodelano
    Mm->elmatstiff (d,ipp);
  }
}


/**
  The function computes stresses at given integration point ipp,
  depending on the reached strains. The cutting plane algorithm is used. 
  The stresses and other array of given integration point are actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mohrcoulombpar::nlstresses (long ipp, long im, long ido)
  //
{
  long i,ni, nc=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(nc),epsp(nc),q(0);

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
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }

  //  new data storage
  for (i=0;i<nc;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+nc]=gamma;
}



/**
  The function computes stresses at given integration point ipp,
  depending on the reached averaged nonlocal strains.
  The cutting plane algorithm is used. The stresses and other array of
  given integration point are actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mohrcoulombpar::nonloc_nlstresses (long ipp, long im, long ido)
  //
{
  long i,ni, nc=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(nc),epsp(nc),q(0);

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
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
  
  //  new data storage
  for (i=0;i<nc;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+nc]=gamma;
}



/**
  The function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mohrcoulombpar::updateval (long ipp, long im, long ido)
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
  @param epsp  - %vector of irreversible strains (output)
 
  @return The function returns %vector of irreversible strains in the parameter epsp.

  Created by Tomas Koudelka,
*/
void mohrcoulombpar::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



/**
  The function extracts consistency parametr gamma for the reached equilibrium state
  from the other array of the given integration point.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns value of consistency parameter.

  Created by Tomas Koudelka,
*/
double mohrcoulombpar::give_consparam (long ipp, long ido)
{
  long ncompstr;
  double gamma;

  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].eqother[ido+ncompstr];

  return gamma;
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mohrcoulombpar::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++){
    switch (atm.atrib[i]){
    case 0:{
      phi=val[0];
      break;
    }
    case 1:{
      c=val[1];
      break;
    }
    case 2:{
      psi=val[2];
      break;
    }
    default:
      print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
  }
  
  alphaphi = sin(phi)/sqrt(9+3*sin(phi)*sin(phi));
  betaphi  = alphaphi / sqrt(3*(3*c*c - alphaphi*alphaphi));
  sigcphi  = sqrt(3*(c*c - alphaphi*alphaphi/3.0));

  alphapsi = sin(psi)/sqrt(9+3*sin(psi)*sin(psi));
  betapsi  = alphapsi / sqrt(3*(3*c*c - alphapsi*alphapsi));
  sigcpsi  = sqrt(3*(c*c - alphapsi*alphapsi/3.0));
}
