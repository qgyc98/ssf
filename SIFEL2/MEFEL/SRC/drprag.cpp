#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "drprag.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"



#define nijac 20
#define limit 1.0e-8




/**
  This constructor inializes attributes to zero values.
*/
drprag::drprag (void)
{
  //  friction angle
  phi=0.0;
  //  cohesion
  c=0.0;
  //  dilatation
  psi=0.0;
  // angle of linear hardening/softening
  theta = 0.0;
  // limit cohesion
  clim = 0.0;
}



/**
  This destructor is only for the formal purposes.
*/
drprag::~drprag (void)
{

}



/**
  This function reads material parameters from the opened text file given
  by the parameter in. Then it computes material constants alpha, alpha1 and
  beta

  @param in - pointer to the opned text file

  25.3.2002  by Tomas Koudelka
*/
void drprag::read (XFILE *in)
{
  //  friction angle
  //  cohesion
  //  dilatation
  //  angle of linear hardening/softening
  //  limit cohesion
  xfscanf (in,"%k%lf%k%lf%k%lf%k%lf%k%lf","phi",&phi,"coh",&c,"psi",&psi,"theta",&theta,"clim",&clim);
  //  stress return algorithm
  sra.read (in);
  
  // outter coincidence with Mohr-Coulomb - compression cone
  // alpha=6.0*sin(phi)/sqrt(3.0)/(3.0-sin(phi));
  // inner coincidence with Mohr-Coulomb - extension cone
  // alpha=6.0*sin(phi)/sqrt(3.0)/(3.0+sin(phi));
  // identical collapse as Mohr-Coulomb in plain strain conditions
  alpha=3.0*tan(phi)/sqrt(9.0 + 12.0*tan(phi)*tan(phi));

  // alpha1=6.0*sin(psi)/sqrt(3.0)/(3.0-sin(psi));
  // alpha1=6.0*sin(psi)/sqrt(3.0)/(3.0+sin(psi));
  alpha1=3.0*tan(psi)/sqrt(9.0 + 12.0*tan(psi)*tan(psi));

  // beta=6.0*cos(phi)/sqrt(3.0)/(3.0-sin(phi));
  // beta=6.0*cos(phi)/sqrt(3.0)/(3.0+sin(phi));
  beta=3.0/sqrt(9.0 + 12.0*tan(phi)*tan(phi));
}



/**
  This function prints material parameters to the opened text file given
  by the parameter in.

  @param out - pointer to the opned text file

  25.3.2015  by Tomas Koudelka
*/
void drprag::print (FILE *out)
{
  fprintf (out,"%le %le %le %le %le ", phi, c, psi, theta, clim);
  sra.print (out);
}



/**
   function computes the value of yield functions
   
   @param sig - engineering components of stress, it contains 6 components
          sig[0] = sigma_x
          sig[1] = sigma_y
          sig[2] = sigma_z
          sig[3] = tau_yz
          sig[4] = tau_zx
          sig[5] = tau_xy
   @param q - internal variables (hardening)

   25.3.2002  by Tomas Koudelka
*/
double drprag::yieldfunction (vector &sig, vector &q)
{
  double p,j2,sqj2,f,h,sigm_max,alphav;
  
  //  mean stress (one third of the first invariant of the stress tensor)
  p = first_invar (sig)/3.0;
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  //  components of the deviator are not needed
  j2 = j2_stress_invar (sig);
  sqj2 = sqrt(j2);
  
  
  sigm_max = beta * cohesion(q)/alpha;
  
  if (p > sigm_max){
    h = p/alpha - sigm_max/alpha;
  }
  else{
    h = 0.0;
  }
  
  if (sqj2 >= h){
    // cone surface return
    f = sqj2 + alpha*p - beta*cohesion(q);
  }
  else{
    // appex return => use artifical slope of the cone alphav which is given 
    // by normal vector to the cone surface defined by appex and current stress point

    alphav = (p - sigm_max)/sqj2;
    f = sqj2 + alphav*p - alphav*sigm_max;
  }
  
  return f;
}



/**
   This function computes derivatives of as-th yield function
   with respect to stress components

   @param sig - stress tensor
   @param dfds - %vector where the resulting derivatives are stored
   @param q - %vector of hardening parameters

   @return The function returns resulting %vector of 
           derivatives in the parameter dfds.

   4.4.2002 by Tomas Koudelka
*/
void drprag::dfdsigma (vector &sig, vector &dfds, vector &q)
{
  double p,j2,h,sigm_max,alphav;
  vector dev(ASTCKVEC(6));
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  j2 = j2_stress_invar (sig);
  j2 = sqrt(j2);
  
  //  deviator of stress tensor
  deviator (sig,dev);

  //  actual mean stress
  p = first_invar(sig)/3.0;

  //  distance of the appex from the origin of coordinate system
  //  it is the maximum possible value of the mean stress
  sigm_max = beta*cohesion(q)/alpha;

  //  limit value of the square root of J2 evaluated from the mean stress
  //  if J2 > h^2, the stress is returned to the surface
  //  if J2 < h^2, the stress is returned to the appex
  h = p/alpha - sigm_max/alpha;
  
  if (j2 >= h){
    // cone surface return
    
    dfds[0] = dev[0]/j2/2.0 + alpha/3.0;
    dfds[1] = dev[1]/j2/2.0 + alpha/3.0;
    dfds[2] = dev[2]/j2/2.0 + alpha/3.0;
    dfds[3] = dev[3]/j2;
    dfds[4] = dev[4]/j2;
    dfds[5] = dev[5]/j2;
  }
  
  else{
    // appex return
    
    alphav = (p - sigm_max)/j2;

    dfds[0] = dev[0]/j2/2.0 + alphav/3.0;
    dfds[1] = dev[1]/j2/2.0 + alphav/3.0;
    dfds[2] = dev[2]/j2/2.0 + alphav/3.0;
    dfds[3] = dev[3]/j2;
    dfds[4] = dev[4]/j2;
    dfds[5] = dev[5]/j2;
  }
  
}

/**
   function assembles second derivatives of the yield function
   with respect to the stress tensor
   
   @param sig - %vector of stress components
   @param dfdsds - %matrix of the derivatives
   
   20. 4. 2015, JK
*/
void drprag::dfdsigmadsigma (vector &sig,matrix &dfdsds)
{
  double j2,j23;
  vector dev(ASTCKVEC(6));
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  j2 = j2_stress_invar (sig);
  j23 = j2*j2*j2;
  j2 = sqrt(j2);
  j23 = sqrt(j23);
  
  //  deviator of stress tensor
  deviator (sig,dev);
  
  //  only one half of components is evaluated
  //  due to symmetry

  dfdsds[0][0] = 0.0 - dev[0]*dev[0]/4.0/j23 + 1.0/3.0/j2;
  dfdsds[0][1] = 0.0 - dev[0]*dev[1]/4.0/j23 - 1.0/6.0/j2;
  dfdsds[0][2] = 0.0 - dev[0]*dev[2]/4.0/j23 - 1.0/6.0/j2;
  dfdsds[0][3] = 0.0 - dev[0]*sig[3]/2.0/j23;
  dfdsds[0][4] = 0.0 - dev[0]*sig[4]/2.0/j23;
  dfdsds[0][5] = 0.0 - dev[0]*sig[5]/2.0/j23;

  dfdsds[1][1] = 0.0 - dev[1]*dev[1]/4.0/j23 + 1.0/3.0/j2;
  dfdsds[1][2] = 0.0 - dev[1]*dev[2]/4.0/j23 - 1.0/6.0/j2;
  dfdsds[1][3] = 0.0 - dev[1]*sig[3]/2.0/j23;
  dfdsds[1][4] = 0.0 - dev[1]*sig[4]/2.0/j23;
  dfdsds[1][5] = 0.0 - dev[1]*sig[5]/2.0/j23;

  dfdsds[2][2] = 0.0 - dev[2]*dev[2]/4.0/j23 + 1.0/3.0/j2;
  dfdsds[2][3] = 0.0 - dev[2]*sig[3]/2.0/j23;
  dfdsds[2][4] = 0.0 - dev[2]*sig[4]/2.0/j23;
  dfdsds[2][5] = 0.0 - dev[2]*sig[5]/2.0/j23;

  dfdsds[3][3] = 0.0 - sig[3]*sig[3]/j23 + 1.0/j2;
  dfdsds[3][4] = 0.0 - sig[3]*sig[4]/j23;
  dfdsds[3][5] = 0.0 - sig[3]*sig[5]/j23;

  dfdsds[4][4] = 0.0 - sig[4]*sig[4]/j23 + 1.0/j2;
  dfdsds[4][5] = 0.0 - sig[4]*sig[5]/j23;

  dfdsds[5][5] = 0.0 - sig[5]*sig[5]/j23 + 1.0/j2;
  

  //  second half of components are copied from
  //  the first half

  dfdsds[1][0] = dfdsds[0][1];

  dfdsds[2][0] = dfdsds[0][2];
  dfdsds[2][1] = dfdsds[1][2];

  dfdsds[3][0] = dfdsds[0][3];
  dfdsds[3][1] = dfdsds[1][3];
  dfdsds[3][2] = dfdsds[2][3];

  dfdsds[4][0] = dfdsds[0][4];
  dfdsds[4][1] = dfdsds[1][4];
  dfdsds[4][2] = dfdsds[2][4];
  dfdsds[4][3] = dfdsds[3][4];

  dfdsds[5][0] = dfdsds[0][5];
  dfdsds[5][1] = dfdsds[1][5];
  dfdsds[5][2] = dfdsds[2][5];
  dfdsds[5][3] = dfdsds[3][5];
  dfdsds[5][4] = dfdsds[4][5];

}


/**
   This function computes derivatives of plastic potential function
   with respect to stress components

   @param sig - %vector of stress components
   @param dgds - %vector of derivatives
   @param q - %vector of hardening parameters

   @return The function returns resulting %vector of derivatives in the 
           parameter dgds.

   4.4.2002 by Tomas Koudelka
*/
void drprag::dgdsigma (vector &sig, vector &dfds, vector &/*q*/)
{
  double j2;
  vector dev(ASTCKVEC(6));
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  j2 = j2_stress_invar (sig);
  j2 = sqrt(j2);
  
  //  deviator of stress tensor
  deviator (sig,dev);
  
  dfds[0] = dev[0]/j2/2.0 + alpha/3.0;
  dfds[1] = dev[1]/j2/2.0 + alpha/3.0;
  dfds[2] = dev[2]/j2/2.0 + alpha/3.0;
  dfds[3] = dev[3]/j2;
  dfds[4] = dev[4]/j2;
  dfds[5] = dev[5]/j2;

}



/**
   This function computes the second derivatives of plastic potential function
   with respect to %vector sigma.

   @param sig - stress tensor
   @param dgds - %matrix where the resulting derivatives are stored
   @param q - %vector of hardening parameters

   @return The function returns resulting %matrix of derivatives in the 
           parameter dgds.

   4.4.2002 by Tomas Koudelka
*/
void drprag::dgdsigmadsigma (vector &sig,matrix &dgdsds, vector &/*q*/)
{
  double j2,j23;
  vector dev(ASTCKVEC(6));
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  j2 = j2_stress_invar (sig);
  j23 = j2*j2*j2;
  j2 = sqrt(j2);
  j23 = sqrt(j23);
  
  //  deviator of stress tensor
  deviator (sig,dev);
  
  //  only one half of components is evaluated
  //  due to symmetry

  dgdsds[0][0] = 0.0 - dev[0]*dev[0]/4.0/j23 + 1.0/3.0/j2;
  dgdsds[0][1] = 0.0 - dev[0]*dev[1]/4.0/j23 - 1.0/6.0/j2;
  dgdsds[0][2] = 0.0 - dev[0]*dev[2]/4.0/j23 - 1.0/6.0/j2;
  dgdsds[0][3] = 0.0 - dev[0]*sig[3]/2.0/j23;
  dgdsds[0][4] = 0.0 - dev[0]*sig[4]/2.0/j23;
  dgdsds[0][5] = 0.0 - dev[0]*sig[5]/2.0/j23;

  dgdsds[1][1] = 0.0 - dev[1]*dev[1]/4.0/j23 + 1.0/3.0/j2;
  dgdsds[1][2] = 0.0 - dev[1]*dev[2]/4.0/j23 - 1.0/6.0/j2;
  dgdsds[1][3] = 0.0 - dev[1]*sig[3]/2.0/j23;
  dgdsds[1][4] = 0.0 - dev[1]*sig[4]/2.0/j23;
  dgdsds[1][5] = 0.0 - dev[1]*sig[5]/2.0/j23;

  dgdsds[2][2] = 0.0 - dev[2]*dev[2]/4.0/j23 + 1.0/3.0/j2;
  dgdsds[2][3] = 0.0 - dev[2]*sig[3]/2.0/j23;
  dgdsds[2][4] = 0.0 - dev[2]*sig[4]/2.0/j23;
  dgdsds[2][5] = 0.0 - dev[2]*sig[5]/2.0/j23;

  dgdsds[3][3] = 0.0 - sig[3]*sig[3]/j23 + 1.0/j2;
  dgdsds[3][4] = 0.0 - sig[3]*sig[4]/j23;
  dgdsds[3][5] = 0.0 - sig[3]*sig[5]/j23;

  dgdsds[4][4] = 0.0 - sig[4]*sig[4]/j23 + 1.0/j2;
  dgdsds[4][5] = 0.0 - sig[4]*sig[5]/j23;

  dgdsds[5][5] = 0.0 - sig[5]*sig[5]/j23 + 1.0/j2;
  

  //  second half of components are copied from
  //  the first half

  dgdsds[1][0] = dgdsds[0][1];

  dgdsds[2][0] = dgdsds[0][2];
  dgdsds[2][1] = dgdsds[1][2];

  dgdsds[3][0] = dgdsds[0][3];
  dgdsds[3][1] = dgdsds[1][3];
  dgdsds[3][2] = dgdsds[2][3];

  dgdsds[4][0] = dgdsds[0][4];
  dgdsds[4][1] = dgdsds[1][4];
  dgdsds[4][2] = dgdsds[2][4];
  dgdsds[4][3] = dgdsds[3][4];

  dgdsds[5][0] = dgdsds[0][5];
  dgdsds[5][1] = dgdsds[1][5];
  dgdsds[5][2] = dgdsds[2][5];
  dgdsds[5][3] = dgdsds[3][5];
  dgdsds[5][4] = dgdsds[4][5];

}



/**
  This function computes value of cohesion which depends on the
  hardening parameter.

  @param qtr - %vector of hardening parameters

  @return This function returns value of the cohesion coresponding to the 
          given hardening parameter.

  4.4.2002 by Tomas Koudelka
*/
double drprag::cohesion(vector &qtr)
{
  double tgt = tan(theta);
  double cq = c + tgt*qtr[0];
  if (((cq < clim) && (c < clim)) ||
      ((cq > clim) && (c > clim)))
    return cq;
  return clim;
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector of hradening parameters.

   @param qtr - %vector of the hardening parameters
   @param dfds - %matrix where the resulting derivatives are stored

   @return The function returns resulting %vector of derivatives in 
           the parameter dfq.

   4.4.2002 by Tomas Koudelka
*/
void drprag::deryieldfq(vector &qtr, vector &dfq)
{
  double tgt = tan(theta);
  double cq = c + tgt*qtr[0];
  if (((cq < clim) && (c < clim)) ||
      ((cq > clim) && (c > clim)))
    dfq[0] = -beta*tgt;
  else
    dfq[0] = 0.0;
  return;
}



/**
   This function computes derivatives of hardening parameters
   with respect of consistency parameter gamma.

   @param dqdg - %matrix where the resulting derivatives are stored


   @return The function returns resulting %vector of derivatives in the parameter dqdg

   4.4.2002 by Tomas Koudelka
*/
void drprag::der_q_gamma(vector &dqdg)
{
  dqdg[0] = sqrt(1.0/3.0 + 2/9*alpha1*alpha1);
}



/**
  This function computes plastic modulus.

  @param qtr -%vector of hardening parameters

  @return The function returns value of the plastic modulus.

  4.4.2002 by Tomas Koudelka
*/
double drprag::plasmodscalar(vector &qtr)
{
  double ret;
  vector dfq(qtr.n);
  vector dqg(qtr.n);

  deryieldfq(qtr, dfq);
  der_q_gamma(dqg);
  scprd(dfq, dqg, ret);
//  ret = 0.0;

  return -ret;
}



/**
  This function computes new value of teh hardening parameter q.

  @param ipp  - integration point pointer
  @param epsp - %vector of the reached plastic strains
  @param q    - %vector of the hardening parameters

  @return The function updates components of the vector of 
          hardening parameters q.

  4.4.2002 by Tomas Koudelka
*/
void drprag::updateq(long ipp, vector &epsp, vector &q)
{
  matrix epst(3,3);
  strastrestate ssst=Mm->ip[ipp].ssst;

  vector_tensor (epsp, epst, ssst, strain);
  q(0) = sqrt(2.0/3.0)*norm(epst);
}



/**
  This function computes material stiffnes %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns calculated stiffness %matrix in the parameter d.

  4.4.2002 by Tomas Koudelka
*/
void drprag::matstiff (matrix &d, long ipp, long ido)
{
  if (Mp->nlman->stmat==initial_stiff)
  {
    //  initial elastic matrix
    Mm->elmatstiff(d, ipp, ido);
  }
  if (Mp->nlman->stmat==tangent_stiff)
  {
    //  tangent stiffness matrix
    matrix ad(d.m,d.n);
    tangentstiff(d, ipp, ido);
  }
}



/**
  This function computes elastic-plastic material stiffnes %matrix.

  @param td - allocated matrix structure for material tangent stiffness %matrix 
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns calculated tangent stiffness %matrix in the parameter td.

  4.4.2002 by Tomas Koudelka
*/
void drprag::tangentstiff (matrix &td,long ipp,long ido)
{
  long ncomp=Mm->ip[ipp].ncompstr;
  matrix d(ASTCKMAT(ncomp, ncomp));
  double gamma;
  
  Mm->elmatstiff(d, ipp);
  gamma = Mm->ip[ipp].other[ido+ncomp];

  if (gamma < 1.0e-10){
    // store block of the resulting matrix in the output argument
    // (it is due to plane stress/strain problems where element matrix has lower dimension)
    extractm(td, d, 0, td.m);
  }
  else{
    double denom;
    vector sig(ASTCKVEC(ncomp)), dfds(ASTCKVEC(ncomp)), av(ASTCKVEC(d.m)),q(ASTCKVEC(1));
    matrix am(ASTCKMAT(d.m, d.n));
    matrix am2(ASTCKMAT(d.m, d.n));
    
    Mm->givestress(0, ipp, sig);
    dfdsigma(sig, dfds, q);
    
    mxv(d, dfds, av);
    scprd(av, dfds, denom);        
    q(0) = Mm->ip[ipp].other[ido+ncomp+1];
    denom += plasmodscalar(q);
    
    if (fabs(denom) < 1.0e-10){
      // store block of the resulting matrix in the output argument
      // (it is due to plane stress/strain problems where element matrix has lower dimension)
      extractm(td, d, 0, td.m);
    }
    else{
      vxv(dfds, dfds, am);
      mxm(d, am, am2);
      mxm(am2, d, am);
      
      cmulm(1.0/denom, am);
      
      subm(d, am, d);
      // store block of the resulting matrix in the output argument
      // (it is due to plane stress/strain problems where element matrix has lower dimension)
      extractm(td, d, 0, td.m);
    }
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

  @return The function updates the stress array of the given integration point.

  4.4.2002 by Tomas Koudelka
*/
void drprag::nlstresses (long ipp, long im, long ido)
{
  long i,ni,n=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(n),epsp(n),q(1);
  
  //  initial values
  for (i=0; i<n; i++){
    //  total strains
    epsn[i]=Mm->ip[ipp].strain[i];
    //  plastic strains
    epsp[i]=Mm->ip[ipp].eqother[ido+i];
  }
  //  consistency parameter
  gamma=Mm->ip[ipp].eqother[ido+n];
  //  hardening variable
  q[0] = Mm->ip[ipp].eqother[ido+n+1];
  
  //  stress return algorithm
  switch (sra.tsra){
  case cp:{
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    //Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    Mm->cutting_plane3 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    break;
  }
  case gsra:{
    ni=sra.give_ni ();
    err=sra.give_err ();
    Mm->newton_stress_return (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    break;
  }
  default:{
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
  }
  
  //  new data storage
  for (i=0;i<n;i++){
    //  plastic strains
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  //  consistency parameter
  Mm->ip[ipp].other[ido+n]=gamma;
  //  hardening variable
  Mm->ip[ipp].other[ido+n+1]=q[0];

}



/**
  This function computes stresses at given integration point ipp,
  depending on the reached averaged nonlocal strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point are actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function updates the stress array of the given integration point.

  4.4.2002 by Tomas Koudelka
*/
void drprag::nonloc_nlstresses (long ipp, long im, long ido)
{
  long i,ni,n=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(n),epsp(n),q(1);

  //  initial values
  for (i=0; i<n; i++){
    epsn[i]=Mm->ip[ipp].strain[i];
    epsp[i]=Mm->ip[ipp].nonloc[i];
  }
  gamma=Mm->ip[ipp].eqother[ido+n];
  q[0] = Mm->ip[ipp].eqother[ido+n+1];
  
  //  stress return algorithm
  if (sra.give_tsra () == cp){
    ni=sra.give_ni ();
    err=sra.give_err ();
    
    //Mm->cutting_plane (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
    Mm->cutting_plane3 (ipp,im,ido,gamma,epsn,epsp,q,ni,err);
  }
  else{
    print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
    abort ();
  }
}



/**
   function assembles the %vector of hardening h which governs
   the evolution of internal parameters
   
   17. 6. 2015
*/
void drprag::hardvect (vector &hv)
{
  hv[0]=1.0;
  //hv[0]=(0.0-k)/sqrt(3.0);
  //hv[0]=0.0-k;
}


/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  @return The function updates the eqother array of internal variables in the given integration point.

  4.4.2002 by Tomas Koudelka
*/
void drprag::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->givencompeqother(ipp,im);

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  @return The function returns the %vector of irreversible strains in the parameter epsp

  4.4.2002 by Tomas Koudelka
*/
void drprag::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



/**
  This function extracts consistency parametr gamma for the attained equilibrium state
  from the integration point eqother array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of consistency parameter.

  4.4.2002 by Tomas Koudelka
*/
double drprag::give_consparam (long ipp, long ido)
{
  long ncompstr;
  double gamma;

  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].eqother[ido+ncompstr];

  return gamma;
}



/**
  This function changes required material parameters and updates
  material constants alpha, alpha1 and beta. It is used for stochastic 
  caluclations.

  @param atm - structure with description of selected parameters
  @param val - %vector of new values for the required material parameters.

  @retval The function returns value of consistency parameter.

  4.4.2002 by Tomas Koudelka
*/
void drprag::changeparam (atsel &atm,vector &val)
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
    case 3:{
      theta=val[i];
      break;
    }
    case 4:{
      clim=val[i];
      break;
    }
    default:{
      fprintf (stderr,"\n\n wrong number of atribute in function changeparam (file %s, line %d).\n",__FILE__,__LINE__);
    }
    }
  }

  alpha=6.0*sin(phi)/(3.0-sin(phi));
  alpha1=6.0*sin(psi)/(3.0-sin(psi));

  beta=6.0*cos(phi)/(3.0-sin(phi));
}
