#include "mohrc.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "matrix.h"
#include "vector.h"
#include "intpoints.h"
#include "tensor.h"
#include "vecttens.h"
#include "elastisomat.h"
#include "galias.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define nijac 20
#define limit 1.0e-8




/**
  The constructor inializes attributes to zero values.
  
  Created by Tomas Koudelka,
*/
mohrcoulomb::mohrcoulomb (void)
{
  phi=0.0;  c=0.0;  psi=0.0;
}



/**
  The destructor is defined only for the formal purposes.
  
  Created by Tomas Koudelka,
*/
mohrcoulomb::~mohrcoulomb (void)
{

}



/**
  The function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void mohrcoulomb::read (XFILE *in)
{
  xfscanf (in,"%k%lf%k%lf%k%lf","phi",&phi,"coh",&c,"psi",&psi);
  sra.read (in);
}



/**
  This function prints material parameters to the opened text file given
  by the parameter in.

  @param out - pointer to the opned text file

  8.6.2015  by Tomas Koudelka
*/
void mohrcoulomb::print (FILE *out)
{
  fprintf (out,"%le %le %le ", phi, c, psi);
  sra.print (out);
}



/**
  Function computes the value of yield functions

  @param sig - stress components

  @return The function returns actual value of yield function.

  Created by Tomas Koudelka, 10.11.2001
*/
double mohrcoulomb::yieldfunction (vector &psig)
{
  double f,k1,k2, s1, s2;

  s1 = psig[0];
  s2 = psig[2];
  //double sigma   =  (s1 + s2) / 2.0;
  //double tau     = -(s1 - s2) / (2.0 * cos(phi));
  //double maxsigt =  c / tan(phi);

  k1 = -1.0 + sin(phi);  k2 = 1.0 + sin(phi);
  f = k1*s1 + k2*s2 - 2.0*c*cos(phi);
  return f;
}



/**
  The function computes derivatives of yield function
  with respect of the principal stresses %vector
   
  @param dfds - %vector for resulting derivatives (output)

  @return The function returns %vector of required derivatives.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::deryieldfsigma (vector &dfds)
{
  dfds[0] = -1.0 + sin(phi);
  dfds[1] =  0.0;
  dfds[2] =  1.0 + sin(phi);
  return;
}



/**
  The function computes derivatives of plastic potential function
  with respect of the principal stresses %vector.
   
  @param dgds - %vector for resulting derivatives (output)

  @return The function returns %vector of required derivatives.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::derplaspotsigma (vector &dgds)
{
  dgds[0] = -1.0 + sin(psi);
  dgds[1] =  0.0;
  dgds[2] =  1.0 + sin(psi);
  return;
}



/**
  The function returns material stiffness %matrix for given integration point.
  Result is stored in the parameter d.

  @param d   - material stiffness %matrix (output)
  @param ipp - integration point pointer
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns stiffness %matrix of the material in the parameter d.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::matstiff (matrix &d, long ipp,long ido)
{
  matrix ad(d.m,d.n);
  switch (Mp->nlman->stmat)
  {
    case initial_stiff:
      //  initial elastic matrix
      Mm->elmatstiff (d,ipp);
      break;
    case tangent_stiff:
      Mm->elmatstiff (ad,ipp);
      tangentstiff (ad,d,ipp,ido);
      break;
    default:
      print_err("unknown type of the stiffness matrix is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function returns tangent material stiffness %matrix for given integration point.
  Result is stored in the parameter d.

  @param d   - elastic material stiffness %matrix
  @param td  - tangent material stiffness %matrix
  @param ipp - integration point pointer
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns tangent stiffness %matrix of the material in the parameter d.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::tangentstiff (matrix &d,matrix &td,long ipp,long ido)
{
  long ncomp=Mm->ip[ipp].ncompstr;
  double denom,gamma;
  vector str,av(d.m),q(1),psig(3), b(3);
  matrix sig(3,3),am(d.m,d.n),pvect(3,3);
  
  gamma=Mm->ip[ipp].eqother[ido+ncomp];
  if (gamma<1.0e-10){
    copym (d,td);
  }
  else{
    
    reallocv (ncomp,str);
    
    Mm->givestress (0,ipp,str);
    vector_tensor (str,sig,Mm->ip[ipp].ssst,stress);
    princ_val (sig, psig, pvect, nijac, limit, Mp->zero, 3,1);
    deryieldfsigma (psig);
    mtxv(pvect, psig, b);
    if ((Mm->ip[ipp].ssst==planestress) || (Mm->ip[ipp].ssst==planestrain)){
      vector auxstr(3);
      auxstr[0]=str[0];auxstr[1]=str[1];auxstr[2]=str[2];
      reallocv (d.m,str);
      str[0]=auxstr[0];str[1]=auxstr[1];str[2]=auxstr[2];
    }
    
    mxv (d,b,av);
    scprd (av,b,denom);
    
    
    //    q[0] = Mm->ip[ipp].eqother[ido+ncomp+1];
    //    denom+= plasmodscalar(q);
    
    if (fabs(denom)<1.0e-10){
      copym (d,td);
    }
    else{
      vxv (b,b,am);
      mxm (d,am,td);
      mxm (td,d,am);
      
      cmulm (1.0/denom,am);
      
      subm (d,am,td);
    }
  }
  
}



/**
  The function returns material stiffness %matrix for given integration point for the 
  principal directions. Result is stored in the parameter d.

  @param ipp - integration point pointer
  @param d   - material stiffness %matrix (output)
  @param e   - Young modulus
  @param nu  - Poissons ratio

  @return The function returns stiffness %matrix in the parameter d.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::pelmatstiff (long ipp, matrix &d, double e, double nu)
{
  double b;

  switch (Mm->ip[ipp].ssst)
  {
    case planestress:
      b = e/(1.0-nu*nu);
      d[0][0] = b;     d[0][1] = b*nu;  d[0][2] = b*nu;
      d[1][0] = b*nu;  d[1][1] = b;     d[1][2] = b*nu;
      d[2][0] = b*nu;  d[2][1] = b*nu;  d[2][2] = b;
      break;
    case planestrain:
    case spacestress:
      b = e /((1.0 + nu)*(1.0 - 2.0*nu));
      d[0][0] = d[1][1] = d[2][2] = 1.0 - nu;
      d[0][1] = d[0][2] = d[1][0] = d[1][2] = d[2][0] = d[2][1] = nu;
      cmulm(b, d);
      break;
    default:
      print_err("unknown type of stress/strain state is required", __FILE__, __LINE__, __func__);
      break;
  }
}



/**
  The function computes correct stresses in the integration point and stores
  them into ip stress array.

  @param ipp - integration point pointer
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 10.11.2001
*/
void mohrcoulomb::nlstresses (long ipp, long ido)
  //
{
  long i,ni,n=Mm->ip[ipp].ncompstr,mu;
  double gamma,err;
  vector epsn(ASTCKVEC(n)), epsp(ASTCKVEC(n)), q(0);

  //  initial values
  for (i=0; i<n; i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  gamma=Mm->ip[ipp].eqother[ido+n];

  //  try single vector stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    mu = cutting_plane (ipp,gamma,epsn,epsp,q,ni,err);
  }
  else{
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
  
  // check singularity
  if (mu > 0)
  {
    // singularity was found, use multisurface stress return algorithm
    gamma=Mm->ip[ipp].eqother[ido+n];
    
    ni=sra.give_ni ();
    err=sra.give_err ();
    mc_msurf_cp(ipp, gamma, epsn, epsp, q, mu,ni,err);
  }

  //  new data storage
  for (i=0;i<n;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+n]=gamma;
}



/**
  This function computes stresses at given integration point ipp,
  depending on the reached averaged nonlocal strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point is actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function does not return anything.

  Created by Tomas Koudelka, 10.11.2001
*/
void mohrcoulomb::nonloc_nlstresses (long ipp,long ido)
{
  long i,ni,n=Mm->ip[ipp].ncompstr,mu;
  double gamma,err;
  vector epsn(n),epsp(n),q(0);

  //  initial values
  for (i=0; i<n; i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].nonloc[i];
  }
  gamma=Mm->ip[ipp].eqother[ido+n];

  //  try single vector stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    mu = cutting_plane (ipp,gamma,epsn,epsp,q,ni,err);
  }
  else{
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
  
  // check singularity
  if (mu > 0){
    // singularity was found, use multisurface stress return algorithm
    ni=sra.give_ni ();
    err=sra.give_err ();
    mc_msurf_cp(ipp, gamma, epsn, epsp, q, mu,ni,err);
  }
  
  //  new data storage
  for (i=0;i<n;i++){
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+n]=gamma;

}




/**
  Function returns stresses on sufrace of plasticity
  cutting plane method in principal stresses is used

  Parameters gamma, epsp and q will be replaced by new values

  @param ipp   - integration point pointer
  @param gamma - consistency parameter (input/output)
  @param epsn  - total strain components
  @param epsp  - plastic strain components (input/output)
  @param q     - hardening parameters (input/output)
  @param ni    - maximum number of iterations
  @param err   - required error
   
  @return Function returns indicator of active yield surfaces in case detection of the singularity
          or zero in case of successfull stress return. In addition, the function
          updates parameters gamma, epsp and q.

  Created by Tomas Koudelka, 4.8.2001
*/
long mohrcoulomb::cutting_plane(long ipp, double &gamma, vector &epsn, vector &epsp, vector &q, long ni, double err)
{
  long i,j,ncomp=epsn.n,idem,mu,zps;
  double f,denom,dgamma,e,nu;
  vector epsa(ncomp), sig(ncomp), epsptr(ncomp);
  vector psig(3), peps(3), pdepsp(3), pepsa(3), dfds(3), dgds(3), depsp(ncomp);
  matrix d(ncomp,ncomp), sigt(3,3), epsat(3,3), epspt(3,3);
  matrix pd(3,3), pvect(3,3);

  // Initialization
  // checking for the correct type of the elastic material
  idem = Mm->ip[ipp].gemid();
  if ((Mm->ip[ipp].tm[idem] != elisomat) && (Mm->ip[ipp].tm[idem] != elisopdmat))
  {
    print_err("Mohr-Coulomb material is combined with\n"
              " the unsupported elastic material.\n"
              " only elisomat type is supported", __FILE__, __LINE__, __func__);
    abort ();
  }

  e = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  // elastic stiffness matrix computation
  Mm->elmatstiff(d, ipp);
  pelmatstiff (ipp, pd, e, nu);
  // computation complementary values for plain stress/strain
  if (Mm->ip[ipp].ssst == planestrain)
  {
     d[0][3] = d[0][1]; d[1][3] = d[1][0];
     d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
  if (Mm->ip[ipp].ssst == planestress)
    epsn[3] = -nu / (1.0 - nu) * (epsn[0]+epsn[1]);
  copyv(epsp, epsptr);

  // Computing principal directions - they do not change during stress return
  // elastic strain
  subv (epsn, epsp, epsa);
  // trial stress computation
  mxv (d, epsa, sig);
  // principal stresses and their directions
  vector_tensor (sig, sigt, Mm->ip[ipp].ssst, stress);
  princ_val (sigt, psig, pvect, nijac, limit, Mp->zero, 3,1);
  // transformation elastic strain into the computed principal directions
  vector_tensor (epsa, epsat, Mm->ip[ipp].ssst, strain);
  glmatrixtransf(epsat, pvect);
  peps[0] = epsat[0][0];
  peps[1] = epsat[1][1];
  peps[2] = epsat[2][2];
  if (Mm->ip[ipp].ssst == planestress)
  {
    zps = checkzeropsig(psig);
    fillrow(0.0, zps, pd);
    fillcol(0.0, zps, pd);
  }
  //  Main iteration loop
  for (i=0;i<ni;i++){
    //  elastic strain
    subv (peps, pdepsp, pepsa);
    //  trial stress computation
    mxv (pd, pepsa, psig);
    // checking ordering of the reached principal stresses and vertex singularity
    mu = checkpsig(psig);
    if (mu > 0)
      return mu;
    // checking yield function
    f = yieldfunction (psig);
    if (f<err)  break;
    if (i==ni-1 && f>err)
    {
      print_err("yield surface was not reached for ipp = %ld", __FILE__, __LINE__, __func__, ipp);
      exit(0);
    }
    // compute necessary derivations of yield function and plastic potential function.
    deryieldfsigma (dfds);
    derplaspotsigma (dgds);

    mxv (pd, dgds, psig);
    scprd (dfds, psig, denom);
    denom += plasmodscalar(ipp, q);
    //  new increment of consistency parameter
    dgamma = f/denom;
    //  new internal variables
    gamma+=dgamma;
    // compute new principal plastic strain increments and transform them into global c.s.
    for (j=0; j<3; j++)
    {
      epspt[j][j]+=dgamma*dgds[j];
      pdepsp[j]+=dgamma*dgds[j];
    }
    lgmatrixtransf(epspt, pvect);
    tensor_vector(depsp, epspt, Mm->ip[ipp].ssst, strain);
    addv(epsp, depsp, epsptr);
    // update hardening parameter
    updateq(ipp, epsptr, q);
  }
  // updating plastic strains
  copyv(epsptr, epsp);
  //  elastic strain
  subv (epsn, epsp, epsa);
  //  trial stress computation
  mxv (d, epsa, sig);
  //  storing resulting stresses
  Mm->storestress (0, ipp, sig);
  return 0;
}



/**
  Function returns plastic modulus in case that hardening/softening rule is adopted

  @param ipp   - integration point pointer
  @param q     - hardening parameters

  @return Function returns value of plastic modulus.

  Created by Tomas Koudelka, 4.11.2003
*/
double mohrcoulomb::plasmodscalar(long /*ipp*/, vector &/*q*/)
{  
  return 0.0;
}



/**
  Function updates hardening/softening parameter in case that hardening/softening rule is adopted

  @param ipp   - integration point pointer
  @param epsp  - plastic strain components
  @param q     - hardening parameters (output)

  @return The function returns udated hardening parameters in the parameter q.

  Created by Tomas Koudelka, 4.11.2003
*/
void mohrcoulomb::updateq(long /*ipp*/, vector &/*epsp*/, vector &/*q*/)
{
  return;
}



/**
  Function checks ordering of principal stresses and vertex singularity.

  @param psig  - principal stress vector sorted in this way psig[0] < psig[1] < psig[2]
  
  @retval 12 - in case that the psig[0] < psig[1] condition is not statisfied
  @retval 23 - in case that the psig[1] < psig[2] condition is not statisfied
  @retval 1  - in case that the vertex singularity occurs i.e. tensile strength is exceeded
  @retval 0  - in other cases

  Created by Tomas Koudelka, 4.11.2003
*/
long mohrcoulomb::checkpsig(vector &psig)
{
  double s1, s2, sigma, tau, maxsigt;

  // vertex return
  s1 = psig[0];
  s2 = psig[2];
  sigma   =  (s1 + s2) / 2.0;
  tau     = -(s1 - s2) / (2.0 * cos(phi));
  maxsigt =  c / tan(phi);
  if ((phi != 0.0) && (tau < ((sigma - maxsigt)/tan(phi))))
    return 1;

  // double vector return
  if (psig[0] > psig[1])
    return 12;
  if (psig[1] > psig[2])
    return 23;

  // normal return
  return 0;
}



/**
  Function checks ordering of principal stresses and vertex singularity.

  @param psig  - principal stress vector sorted in this way psig[0] < psig[1] < psig[2]

  @return Function returns index of zero psig.

  Created by Tomas Koudelka,4.11.2003
*/
long mohrcoulomb::checkzeropsig(vector &psig)
{
  long i;

  for (i=0; i<3; i++)
  {
    if (psig[i] == 0.0)
      return i;
  }
  return -1;  
}



/**
  Function returns stresses on sufrace of plasticity
  multisurface cutting plane method at principal stresses is used.
  Parameters gamma,epsp and q will be replaced by new values.

  @param ipp - integration point pointer
  @param gamma - consistency parameter (input/output)
  @param epsn - total strain components
  @param epsp - plastic strain components (input/output)
  @param q - hardening parameters (input/output)
  @param mu - possibly active yield surfaces indicator
  @param ni - maximum number of iterations
  @param err - required error

  @return The function returns actual values of consistency parameter,
          plastic strains and hardening parameters in the parameters gamma, epsp and q.
   
  Created by Tomas Koudelka, 4.8.2001
*/
void mohrcoulomb::mc_msurf_cp(long ipp, double &gamma, vector &epsn, vector &epsp, vector &q, long mu,long ni,double err)
{
  long i,j,ncomp=epsn.n,idem,zps;
  double e,nu;
  long stat[2], nas;
  vector epsa(ncomp), sig(ncomp), epsptr(ncomp), depsp(ncomp);
  vector psig(3),peps(3),pdepsp(3),pepsa(3),dpepsp(3), f(2);
  matrix d(ncomp,ncomp),sigt(3,3),epsat(3,3),epspt(3,3);
  matrix pd(3,3),pvect(3,3);
  matrix dfds, dgds, am, cpm, hcpm;
  vector dgamma, ff;   

  // Initialization
  // checking for the correct type of the elastic material
  idem = Mm->ip[ipp].gemid();
  if ((Mm->ip[ipp].tm[idem] != elisomat) && (Mm->ip[ipp].tm[idem] != elisopdmat))
  {
    print_err("Mohr-Coulomb material is combined with\n"
              " the unsupported elastic material.\n"
              " only elisomat type is supported", __FILE__, __LINE__, __func__);
    abort ();
  }

  e = Mm->give_actual_ym(ipp);
  nu = Mm->give_actual_nu(ipp);
  // elastic stiffness matrix computation
  Mm->elmatstiff(d, ipp);
  pelmatstiff (ipp, pd, e, nu);
  // computation complementary values for plain stress/strain
  if (Mm->ip[ipp].ssst == planestrain)
  {
     d[0][3] = d[0][1]; d[1][3] = d[1][0];
     d[3][0] = d[1][0]; d[3][1] = d[1][0]; d[3][3] = d[1][1];
  }
/*  if (Mm->ip[ipp].ssst == planestress)
    epsn[3] = -nu / (1.0 - nu) * (epsn[0]+epsn[1]);*/
  copyv(epsp, epsptr);

  // Computing principal directions - they do not change during stress return
  // elastic strain
  subv (epsn, epsp, epsa);
  // trial stress computation
  mxv (d, epsa, sig);
  // principal stresses and their directions
  vector_tensor (sig, sigt, Mm->ip[ipp].ssst, stress);
  princ_val (sigt, psig, pvect, nijac, limit, Mp->zero, 3,1);
  // transformation elastic strain into the computed principal directions
  vector_tensor (epsa, epsat, Mm->ip[ipp].ssst, strain);
  glmatrixtransf(epsat, pvect);
  peps[0] = epsat[0][0];
  peps[1] = epsat[1][1];
  peps[2] = epsat[2][2];
  if (Mm->ip[ipp].ssst == planestress)
  {
    zps = checkzeropsig(psig);
    fillrow(0.0, zps, pd);
    fillcol(0.0, zps, pd);
  }

 // Main iteration loop
 // number of yield surfaces is always 2, for every mu value
 for (i=0; i<ni; i++)
 {
    //  elastic strain
    subv (peps, pdepsp, pepsa);
    //  trial stress computation
    mxv (pd, pepsa, psig);
    //  yield function control
    yieldfunction(psig, mu, f);
    // detection of active yield surfaces
    nas=0;
    for (j=0; j<2; j++)
    {
      if (f[j] >= err)
      {
        stat[j] = 1;
        nas++;
      }
      else
        stat[j] = 0;
    }
    if (nas==0)  // no active yield surface
      break;
    
    // allocation of necessarry matrices and vectors        
    reallocm (nas, 3, dfds); reallocm (nas, 3, dgds);
    reallocm (nas, 3, am);   reallocm (nas, nas, cpm);
    reallocm (nas, nas, hcpm);
    reallocv (nas, dgamma);
    reallocv (nas, ff);
    yieldfunction(psig, mu, stat, ff);
    // computing matrices of derivations
    dfdsigma(stat, mu, dfds);
    dgdsigma(stat, mu, dgds);
    plasmodscalar(ipp, q, mu, hcpm);
    // assembling resulting matrix for the dgamma solution
    mxm(dfds, pd, am);
    mxmt(am, dgds, cpm);
    addm(cpm, hcpm, cpm);
    // solution of the equation system
    gemp(cpm.a, dgamma.a, ff.a, nas, 1, Mp->zero, 1);
    // new plastic strain vector increments
    mtxv(dgds, dgamma, dpepsp);
    // compute new principal plastic strain increments and their transformed values into g.c.s.
    for (j=0; j<3; j++)
    {
      epspt[j][j]+=dpepsp[j];
      pdepsp[j]+=dpepsp[j];
    }
    for (j=0; j<nas; j++)
      gamma += dgamma[j];
    lgmatrixtransf(epspt, pvect);
    tensor_vector(depsp, epspt, Mm->ip[ipp].ssst, strain);
    addv(epsp, depsp, epsptr);
    // update hardening parameter
    updateq(ipp,epsptr, q);
  }
  // updating plastic strains
  copyv(epsptr, epsp);
  //  elastic strain
  subv (epsn, epsp, epsa);
  //  trial stress computation
  mxv (d, epsa, sig);
  //  storing resulting stresses
  Mm->storestress (0, ipp, sig);
  return;
}



/**
  The function assembles %matrix of plastic moduli used in the multisurface 
  cutting-plane algorithm.

  @param ipp  - integration point pointer
  @param q    - hardening parameters
  @param mu   - possibly active yield surfaces indicator
  @param hcpm - %matrix of plastic moduli (output)
  
  @return The function returns resulting %matrix in the parameter hcpm.

  Created by Tomas Koudelka,
*/
void mohrcoulomb::plasmodscalar(long /*ipp*/, vector &/*q*/, long /*mu*/, matrix &hcpm)
{
  fillm(0.0, hcpm);
}  



/**
  The function computes actual values of TRUE ACTIVE yield functions for the multisurface
  cutting-plane algorithm.

  @param psig - %vector of principal stresses
  @param mu   - possibly active yield surfaces indicator
  @param stat - 2 component array with indicator of active yield surface
                Two surfaces should be active(stat[i]=1), but
                one can become inactive for some limit cases during stress return.
                Array is set in the function mc_msurf_cp depending on the results of 
                yieldfunction(vector&, long, vector&) call.
  @param f - %vector of actual values of active yield functions (output)

  @return The function returns actual values of active yield functions in the parameter f.

  Created by Tomas Koudelka, 4.8.2001
*/
void mohrcoulomb::yieldfunction(vector &psig, long mu, long *stat, vector &f)
{
  long i = 0;
  double k1, k2, s1, s2;

  s1 = psig[0];
  s2 = psig[2];
  k1 = -1.0 + sin(phi);  k2 = 1.0 + sin(phi);

  if (stat[0])
  {
    s1 = psig[0];
    s2 = psig[2];
    f[i] = k1*s1 + k2*s2 - 2.0*c*cos(phi);
    i++;
  }
  if (stat[1])
  {
    switch (mu)
    {
      case 12: // swap s1 and s2 at the condition s1 < s2 < s3
        s1 = psig[1];
        s2 = psig[2];
        f[i] = k1*s1 + k2*s2 - 2.0*c*cos(phi);
        break;
      case 23: // swap s2 and s3 at the condition s1 < s2 < s3
        s1 = psig[0];
        s2 = psig[1];
        f[i] = k1*s1 + k2*s2 - 2.0*c*cos(phi);
        break;
      case 1: // tensile strength cutoff
        f[i] = s1 + s2 - 2.0*c/tan(phi);
        break;
      default:
        break;
    }
  }
  return;
} 



/**
  The function computes actual values of POSSIBLY ACTIVE yield functions for the multisurface
  cutting-plane algorithm.

  @param psig - %vector of principal stresses
  @param mu   - possibly active yield surfaces indicator
  @param stat - 2 component array with indicator of active yield surface
                Two surfaces should be active(stat[i]=1), but
                one can become inactive for some limit cases during stress return.
                Array is set in the function mc_msurf_cp depending on the results of 
                yieldfunction(vector&, long, vector&) call.
  @param f - %vector of actual values of active yield functions (output)

  @return The function returns actual values of active yield function in the parameter f.

  Created by Tomas Koudelka, 4.8.2001
*/
void mohrcoulomb::yieldfunction(vector &psig, long mu, vector &f)
{
  double k1, k2, s1, s2;

  k1 = -1.0 + sin(phi);  k2 = 1.0 + sin(phi);
  s1 = psig[0];
  s2 = psig[2];

  f[0] = k1*s1 + k2*s2 - 2.0*c*cos(phi);
  switch (mu)
  {
    case 12: // swap s1 and s2 at the condition s1 < s2 < s3
      s1 = psig[1];
      s2 = psig[2];
      f[1] = k1*s1 + k2*s2 - 2.0*c*cos(phi);
      break;
    case 23: // swap s2 and s3 at the condition s1 < s2 < s3
      s1 = psig[0];
      s2 = psig[1];
      f[1] = k1*s1 + k2*s2 - 2.0*c*cos(phi);
      break;
    case 1: // tensile strength cutoff
      f[1] = s1 + s2 - 2.0*c/tan(phi);
      break;
    default:
      break;
  }
  return;
}
 


/**
  The function computes derivatives of TRUE ACTIVE yield functions
  with respect of the principal stresses %vector
   
  @param mu   - possibly active yield surfaces indicator
  @param stat - 2 component array with indicator of active yield surface
                Two surfaces should be active(stat[i]=1), but
                one can become inactive for some limit cases during stress return.
                Array is set in the function mc_msurf_cp depending on the results of 
                yieldfunction(vector&, long, vector&) call.
  @param dfds - %vector for resulting derivatives (output)

  @return The function returns %vector of required derivatives.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::dfdsigma(long *stat, long mu, matrix &dfds)
{
  long i = 0;
  if (stat[0])
  {
    dfds[i][0] = -1.0 + sin(phi);
    dfds[i][1] =  0.0;
    dfds[i][2] =  1.0 + sin(phi);
    i++;
  }
  if (stat[1])
  {
    switch (mu)
    {
      case 12: // swap s1 and s2 at the condition s1 < s2 < s3
        dfds[i][0] =  0.0;
        dfds[i][1] = -1.0 + sin(phi);
        dfds[i][2] =  1.0 + sin(phi);
        break;
      case 23: // swap s2 and s3 at the condition s1 < s2 < s3
        dfds[i][0] = -1.0 + sin(phi);
        dfds[i][1] =  1.0 + sin(phi);
        dfds[i][2] =  0.0;
        break;
      case 1: // tensile strength cutoff
        dfds[i][0] = 1.0;
        dfds[i][1] = 0.0;
        dfds[i][2] = 1.0;
        break;
      default:
        break;
    }
  }
  return;
}



/**
  The function computes derivatives of TRUE ACTIVE plastic potential functions
  with respect of the principal stresses %vector
   
  @param mu   - possibly active yield surfaces indicator
  @param stat - 2 component array with indicator of active yield surface
                Two surfaces should be active(stat[i]=1), but
                one can become inactive for some limit cases during stress return.
                Array is set in the function mc_msurf_cp depending on the results of 
                yieldfunction(vector&, long, vector&) call.
  @param dgds - %vector for resulting derivatives (output)

  @return The function returns %vector of required derivatives.

  Created by Tomas Koudelka, 10.11.2003
*/
void mohrcoulomb::dgdsigma(long *stat, long mu, matrix &dgds)
{
  long i = 0;
  if (stat[0])
  {
    dgds[i][0] = -1.0 + sin(psi);
    dgds[i][1] =  0.0;
    dgds[i][2] =  1.0 + sin(psi);
    i++;
  }
  if (stat[1])
  {
    switch (mu)
    {
      case 12: // swap s1 and s2 at the condition s1 < s2 < s3
        dgds[i][0] =  0.0;
        dgds[i][1] = -1.0 + sin(psi);
        dgds[i][2] =  1.0 + sin(psi);
        break;
      case 23: // swap s2 and s3 at the condition s1 < s2 < s3
        dgds[i][0] = -1.0 + sin(psi);
        dgds[i][1] =  1.0 + sin(psi);
        dgds[i][2] =  0.0;
        break;
      case 1: // tensile strength cutoff
        dgds[i][0] = 1.0;
        dgds[i][1] = 0.0;
        dgds[i][2] = 1.0;
        break;
      default:
        break;
    }
  }
  return;
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
void mohrcoulomb::updateval (long ipp,long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp, im);

  for (i=0; i<n; i++)
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];

  n = Mm->ip[ipp].ncompstr;
  vector epsn(ASTCKVEC(n)), epsp(ASTCKVEC(n)), epse(ASTCKVEC(n)), fepsp(ASTCKVEC(6));
  vector sig(ASTCKVEC(n)), psig(ASTCKVEC(3));
  matrix d(ASTCKMAT(n,n));
  matrix sigt(ASTCKMAT(3,3)), pvect(ASTCKMAT(3,3));
  //  initial values
  for (i=0; i<n; i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  
  give_full_vector(fepsp, epsp, Mm->ip[ipp].ssst);
  Mm->ip[ipp].other[ido+n+1]=sqrt(2.0/3.0)*tensor_strain_norm(fepsp); // equivalent plastic strain

  Mm->elmatstiff(d, ipp, ido);
  subv(epsn, epsp, epse);
  mxv(d, epse, sig);
  vector_tensor(sig, sigt, Mm->ip[ipp].ssst, stress);
  princ_val(sigt, psig, pvect, 30, 1.0e-6, 1.0e-15, 3, 1);  
  Mm->ip[ipp].other[ido+n+2]=0.5*(psig(0)-psig(2)); // tau_max
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  @return The function returns %vector of irreversible strains in the parameter epsp.

  Created by Tomas Koudelka,
*/
void mohrcoulomb::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



/**
  The function returns actual value of consistency parameter gamma.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array

  @return The function returns consistency parameter.

  Created by Tomas Koudelka,
*/
double mohrcoulomb::give_consparam (long ipp,long ido)
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
void mohrcoulomb::changeparam (atsel &atm,vector &val)
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
    default:
      print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
  }
}
