#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "camclay.h"
#include "global.h"
#include "probdesc.h"
#include "mechtop.h"
#include "mechmat.h"
#include "globmat.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "elastisomat.h"
#include "mathem.h"



/**
  This constructor inializes attributes to zero values.
*/
camclay::camclay (void)
{
  m = 0.0;
  lambda = 0.0;
  kappa = 0.0;
}



/**
  This destructor is only for the formal purposes.
*/
camclay::~camclay (void)
{

}



/**_
  This function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opned text file

*/
void camclay::read (XFILE *in)
{
  xfscanf (in, "%k%lf%k%lf%k%lf", "m", &m, "lambda", &lambda, "kappa", &kappa);
  sra.read (in);
}



/**
  The function prints material parameters into the opened text file given
  by the parameter out.

  @param out - pointer to the opened input text file
  
  @return The function does not return anything.

  Created by T. Koudelka 2.7.2015
*/
void camclay::print (FILE *out)
{
  fprintf (out, "%le %le %le ", m, lambda, kappa);
  sra.print (out);
}



/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array

   12/06/2012 TKo
*/
void camclay::initval(long ipp, long ido)
{
  double v_ini, v_pc0, v_lambda1;
  double i1s, j2s;
  double v_kappa1, p1, pc_0;
  long i, ncompstr = Mm->ip[ipp].ncompstr;
  vector sig(ASTCKVEC(ncompstr)), sigt(ASTCKVEC(6)), q(ASTCKVEC(1));
  double err=sra.give_err ();
  
  if (Mm->ip[ipp].eqother[ido+ncompstr+1] == 0.0)  // actual value of p_c was not set
  {
    // initial value of specific volume under small pressure p_1
    v_kappa1 = Mm->ic[ipp][0];
    // initial reference presure
    p1       = Mm->ic[ipp][1];
    // initial preconsolidation pressure
    pc_0     = Mm->ic[ipp][2];

    //  initial stresses
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
      for (i=0; i<ncompstr; i++)
        sig(i) = Mm->ip[ipp].stress[i] = Mm->eigstresses[ipp][i];
    }
    else{
      print_err("initial stresses (eigentsresses) are not defined on element %ld, ip=%ld,\n"
                " cannot determine specific volume for initial stress state", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
      abort();
    }

    give_full_vector(sigt, sig, Mm->ip[ipp].ssst);
    // initial value of total mean stress
    i1s = first_invar(sigt)/3.0;
    q(0) = pc_0;
    // test if the initial state is elastic, if not calculate new value of inital preconsolidation pressure pc_0
    if (yieldfunction(sigt, q) > err){
      j2s = j2_stress_invar (sigt);
      pc_0 = j2s/(i1s*m*m) + i1s;
    }

    // specific volume for initial preconsolidation pressure
    v_pc0 = v_kappa1 - kappa*log(pc_0/p1);
    // actual value of p_c = initial preconsolidtaion pressure
    Mm->ip[ipp].other[ido+ncompstr+1] = Mm->ip[ipp].eqother[ido+ncompstr+1] = pc_0;
    // original specific volume before any loading
    Mm->ip[ipp].other[ido+ncompstr+2] = Mm->ip[ipp].eqother[ido+ncompstr+2] = v_lambda1 = v_pc0 + lambda*log(pc_0/p1);
    // specific volume for initial stress state
    Mm->ip[ipp].other[ido+ncompstr+3] = Mm->ip[ipp].eqother[ido+ncompstr+3] = v_ini     = v_pc0 + kappa*log(pc_0/i1s);
    // inital value of mean stress
    Mm->ip[ipp].other[ido+ncompstr+4] = Mm->ip[ipp].eqother[ido+ncompstr+4] = i1s;
    // porosity n
    Mm->ip[ipp].other[ido+ncompstr+10] = Mm->ip[ipp].eqother[ido+ncompstr+10] = (v_ini - 1.0)/v_ini; // v_ini = 1.0 + e_ini
  }  

  check_math_errel(Mm->elip[ipp]);

  return;
}


/**
   This function computes the value of yield functions.

   @param sig - stress tensor
   @param q   - %vector of hardening parameter

   @retval The function returns value of yield function for the given stress tensor

   25.3.2002
*/
double camclay::yieldfunction (vector &sig, vector &q)
{
  double i1s,j2s,f;
  
  //  second invariant of deviator
  //  it is expressed with the help of the stress components
  //  components of the deviator are not needed
  j2s = j2_stress_invar (sig);

  i1s = first_invar (sig)/3.0;

  f = j2s/(m*m) + i1s * (i1s - q[0]);
  if (f > 0.0) 
    return f;
  
  return f;
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector sigma.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameter
   @param dfds - %matrix where the resulting derivatives are stored


   4.1.2002
*/
void camclay::deryieldfsigma (vector &sig, vector &q, vector &dfds)
{
  double volume,i1s;
  
  i1s = first_invar (sig)/3.0;
  
  deviator (sig,dfds);
  
  cmulv(1.0/(m*m), dfds);
  volume=1.0/3.0*(2.0*i1s - q[0]);
  dfds(0)+=volume;
  dfds(1)+=volume;
  dfds(2)+=volume;  
  dfds(3)*= 2.0;
  dfds(4)*= 2.0;
  dfds(5)*= 2.0;
}



/**
   The function computes the second derivatives of yield function
   with respect of stress tensor sigma.

   @param ddfds - tensor of the 4-th order where the are derivatives stored

   19.12.2002
*/
void camclay::dderyieldfsigma (matrix &ddfds)
{
  fillm(0.0, ddfds);
  
  ddfds[0][0] = ddfds[1][1] = ddfds[2][2] = 2.0/(3.0*m*m) + 2.0/9.0;
  ddfds[0][1] = ddfds[0][2] = ddfds[1][0] = ddfds[1][2] = -1.0/(3.0*m*m) + 2.0/9.0;
  ddfds[2][0] = ddfds[2][1] = ddfds[0][1];
  ddfds[3][3] = 1.0/(m*m);
  ddfds[4][4] = 1.0/(m*m);
  ddfds[5][5] = 1.0/(m*m);
/*
  switch (ssst)
  {
    case planestress:
    case planestrain:
      ddfds[0][0] = ddfds[1][1] = 2.0/(3.0*m*m) + 2.0/9.0;
      ddfds[0][1] = ddfds[1][0] = -1.0/(3.0*m*m) + 2.0/9.0;
      ddfds[2][2] = 1.0/(m*m);
      break;
    case axisymm:
      ddfds[0][0] = ddfds[1][1] = ddfds[2][2] = 2.0/(3.0*m*m) + 2.0/9.0;
      ddfds[0][1] = ddfds[0][2] = ddfds[1][0] = ddfds[1][2] = -1.0/(3.0*m*m) + 2.0/9.0;
      ddfds[2][0] = ddfds[2][1] = -1.0/(3.0*m*m) + 2.0/9.0;
      ddfds[3][3] = 1.0/(m*m);
      break;
    case spacestress:
      ddfds[0][0] = ddfds[1][1] = ddfds[2][2] = 2.0/(3.0*m*m) + 2.0/9.0;
      ddfds[0][1] = ddfds[0][2] = ddfds[1][0] = ddfds[1][2] = -1.0/(3.0*m*m) + 2.0/9.0;
      ddfds[2][0] = ddfds[2][1] = ddfds[0][1];
      ddfds[3][3] = ddfds[4][4] = ddfds[5][5] = 1.0/(m*m);
      break;
    default:
      print_err("unknown type of stress/strain state is required.", __FILE__, __LINE__, __func__);
      break;
  }*/
}



/**
   The function computes derivatives of plastic potential function
   with respect of vector sigma.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameter
   @param dgds - %matrix where the resulting derivatives are stored
*/
void camclay::derpotsigma (vector &sig, vector &q, vector &dgds)
{
  deryieldfsigma (sig, q, dgds);
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector of hradening parameters.

   @param sig - stress tensor
   @param dfds - %vector where the resulting derivatives are stored


   4.1.2002
*/
void camclay::deryieldfq(vector &sig, vector &dfq)
{
  dfq[0] = -first_invar (sig)/3.0;

//  pro fci plasticity prenasobenou m*m  
//    dfq[0] = first_invar (sig)/3.0;
//  dfq[0] *= -m*m;
  return;
}



/**
   This function computes the second derivatives of yield function
   with respect to hradening parameter.

   @param dfds - %matrix, where the resulting derivatives are stored


   12.4.2005
*/
void camclay::deryieldfdqdq(matrix &dfdqdq)
{  
  dfdqdq[0][0] = 0.0;
  return;
}



/**
   This function computes the second derivatives of yield function
   with respect to stresses.

   @param dfds - tensor, where the resulting derivatives are stored.
                 size of dfds = (6,number_of_hardening_param)


   12.4.2005
*/
void camclay::deryieldfdsdq(matrix &dfdsdqt)
{  
  dfdsdqt[0][0] = dfdsdqt[1][0] = dfdsdqt[2][0] = -1.0/3.0;
  return;
}



/**
   This function computes the derivatives of hardening parameters
   with respect to stresses.

  @param sigt - stress tensor
  @param qt   - %vector of hardening variables
  @param dqds - tensor, where the resulting derivatives are stored.
                 size of dqds = (3,3)


   12.4.2005
*/
void camclay::dqdsigma(vector &sigt, vector &qt, vector &dqds)
{  
  double i1s;

  i1s = first_invar (sigt)/3.0;

  fillv (0.0, dqds);
  dqds[0] = dqds[1] = dqds[2] = qt[0]*kappa/i1s*1/3.0;
  return;
}



/**
   This function computes the derivatives of hardening function
   with respect to hardening parameters.

   
   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array
   @param dgamma - increment of consistency parameter gamma
   @param qt   - %vector of hardening variables
   @param dqdq - tensor, where the resulting derivatives are stored.
                 size of dqdq = (number_of_hardening_param, number_of_hardening_param)


   12.4.2005
*/
void camclay::dhardfdq(long ipp, long ido, double dgamma, vector &qt, vector &dqdq)
{  
  double v_ini;
  long ncompstr = Mm->ip[ipp].ncompstr;

  fillv(0.0, dqdq);
  v_ini = Mm->ip[ipp].eqother[ido+ncompstr+3];
  dqdq[0] = -dgamma * v_ini * qt[0]/(kappa-lambda);
  return;
}



/**
   This function computes derivatives of hardening paramters
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array
   @param sig - stress tensor
   @param qtr  - %vector of hardening variables
   @param dqdg - %matrix where the resulting derivatives are stored


   4.1.2002
*/
void camclay::der_q_gamma(long ipp, long ido, vector &sig, vector &qtr, vector &dqdg)
{
  double v_ini, i1s;
  long ncompstr = Mm->ip[ipp].ncompstr;

  // original specific volume before any loading
  i1s = first_invar (sig)/3.0;
  v_ini = Mm->ip[ipp].eqother[ido+ncompstr+3];

  dqdg[0] = qtr[0] * (v_ini)/(kappa-lambda) * (2.0*i1s-qtr[0]);

  //  pro fci plasticity prenasobenou m*m  
  //  dqdg[0] = qtr[0] * (-v_ini)/(lambda-kappa)*m*m*(2.0*i1s-qtr[0]);
  return;
}



/**
  This function computes plastic modulus.

  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  @param sig - stress tensor
  @param qtr - %vector of hardening parameters

  @retval The function returns value of the plastic modulus.
*/
double camclay::plasmodscalar(long ipp, long ido, vector &sig, vector &qtr)
{
  double ret;
  vector dfq(1);
  vector dqg(1);

  deryieldfq(sig, dfq);
  der_q_gamma(ipp, ido, sig, qtr, dqg);
  scprd(dfq, dqg, ret);

  return -ret;
}



/**
  This function computes new value of the hardening parameter q.

  @param ipp  - integration point pointer
  @param ido  - index of internal variables for given material in the ipp other array
  @param eps  - %vector of the attained strains
  @param epsp - %vector of the attained plastic strains
  @param sig  - attained stress tensor
  @param ssst - stress/strain state parameter (for the used vector_tensor function)
  @param q    - %vector of the hardening parameters

*/
void camclay::updateq(long ipp, long ido, vector &/*eps*/, vector &epsp, vector &/*sig*/, vector &q)
{
  vector epst(ASTCKVEC(6));
  vector depsp(ASTCKVEC(epsp.n));
  double v_kappa1, v_lambda1, p1, depsvp, v_ini, pc_new;
  strastrestate ssst = Mm->ip[ipp].ssst;
  long ncompstr = Mm->ip[ipp].ncompstr;

  // initial value of specific volume under small pressure p_1
  v_kappa1  = Mm->ic[ipp][0];
  // initial reference presure
  p1        = Mm->ic[ipp][1];
  // original specific volume before any loading
  v_lambda1 = Mm->ip[ipp].eqother[ido+ncompstr+2];
  // specific volume for initial stress state
  v_ini     = Mm->ip[ipp].eqother[ido+ncompstr+3];

//  Original version of pc calculation - problems with tension
//  vk = v_ini*(1+epsv)+kappa*log(i1s/p1);
//  pc_new = p1 * exp((vk-v_lambda1)/(kappa-lambda));
//
//  Alternative version of pc calculation - 'less' accurate in benchmarks
//  pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*epsv)/(kappa-lambda));
//  vk = v_ini*(1+epsv)+kappa*log(i1s/p1);

  for (long i=0; i<ncompstr; i++)  // add plastic strains and eigenstrains (inital value of elastic strain)
    depsp[i] = epsp[i];// + Mm->ic[ipp][3+i];

  give_full_vector (epst, depsp, ssst);
  depsvp = first_invar(epst);

  pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*depsvp)/(kappa-lambda));
  
  //if (pc_new < q(0))
  //  q(0) = pc_new;    
  q(0) = pc_new;    

  return;
}



/**
  The function returns actual Young modulus value.

  @param ipp - integration point number
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

*/
double camclay::give_actual_ym(long ipp, long im, long ido)
{
  long idem    = Mm->ip[ipp].gemid();
  long ncompo  = Mm->givencompeqother(ipp, im);
  double e;

  // initial Young modulus
  double e_ini = Mm->give_initial_ym(ipp, idem, ido+ncompo);
  // actual Poisson's ratio
  double nu = give_actual_nu(ipp, im, ido);

  // actual values of strains
  long ncomp=Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  //vector epsp(ASTCKVEC(ncomp));
  vector epst(ASTCKVEC(6));
  vector eps(ASTCKVEC(ncomp));
  vector sig(ASTCKVEC(ncomp));
  vector q(ASTCKVEC(1));

  // actual total volumetric strain
  Mm->givestrain(0, ipp, eps);
  give_full_vector(epst, eps, ssst);
  double epsv = first_invar(epst);
  // actual plastic volumteric strain
  //Mm->giveother(ipp, ido, ncomp, epsp);
  // total volumetric strain from the last converged step
  double epsvo = Mm->ip[ipp].eqother[ido+ncomp+6];

  // actual value of p_c
  //updateq(ipp, ido, eps, epsp, sig, q);
  //double pc = q(0);
  double pc = Mm->ip[ipp].other[ido+ncomp+1];

  /*
  //  initial value of elastic strain
  for (long i=0; i<ncomp; i++)
    eps(i) = Mm->ic[ipp][3+i];
  give_full_vector(auxt, eps, ssst);
  double epsv_ini = first_invar(auxt);
  if (epsv >= epsv_ini)
    return e_ini;*/
  // actual value of of specific volume
  double v_ini = Mm->ip[ipp].eqother[ido+ncomp+3];
  //double v     = v_ini*(1.0+epsv);

  //double p = pc/exp((v-v_pc)/kappa);
  double po = Mm->ip[ipp].eqother[ido+ncomp+4];
  if ((po == 0.0) || (epsv == 0.0))
    return e_ini;
  if (po > pc){
    // for the unloading: \Delta eps^{el}_V = \Delta eps_V
    // actual p according actual swelling line
    double p = po*exp(-(epsv-epsvo)*v_ini/kappa);
    if (epsv == epsvo)
      e = 3.0*(1.0-2.0*nu)*v_ini*(-p)/kappa;
    else 
      e = 3.0*(1.0-2.0*nu)*(p-po)/(epsv-epsvo);
  }
  else{    
    e = e_ini;
  }

  //e = 3.0*(1.0-2.0*nu)*v*(-p)/kappa;

  return e;
}



/**
  The function returns actual Poisson's ratio.

  @param ipp - integration point number
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

*/
double camclay::give_actual_nu(long ipp, long im, long ido)
{
 long idem = Mm->ip[ipp].gemid();
 long ncompo = Mm->givencompeqother(ipp, im);
 double nu = Mm->give_initial_nu(ipp, idem, ido+ncompo);

 return nu;
}



/**
  This function computes material stiffness %matrix.

  @param d - allocated matrix structure for material stiffness %matrix
  @param ipp - integration point number
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

*/
void camclay::matstiff (matrix &d, long ipp, long /*im*/, long /*ido*/)
{
  if (Mp->nlman->stmat==1)
  {
    //  initial elastic matrix
    Mm->elmatstiff (d,ipp);
  }
/*
  if (Mp->nlman->stmat > 1)
  {
    //  consistent tangent stiffness matrix
    long   n = Mm->ip[ipp].ncompstr, i, j;    
    double gamma = Mm->ip[ipp].other[ido+n];
    strastrestate ssst = Mm->ip[ipp].ssst;
    double zero = 1.0e-20;
 
    matrix de(ASTCKMAT(6,6)), ddfdst(ASTCKMAT(6,6)), tmp(ASTCKMAT(6,6)), auxm(ASTCKMAT(6,6));
    vector sig(ASTCKVEC(n)), dfds(ASTCKVEC(n)), dgds(ASTCKVEC(n));
    vector q(ASTCKVEC(1));

    matrix rm(ASTCKMAT(n,n));
    vector auxv(ASTCKVEC(n)), tmpv(ASTCKVEC(n));
    double aux1;

    // calculation of auxiliary matrix R
    // \mbf{R} = (\mbf{I} + \gamma \mbf{D}_e dd_f/dd_sigma)^(-1) \mbf{D}_e

    Mm->elmatstiff (de, ipp, spacestress);         // de = \mbf{D}_e

    if (Mp->stressstate == 0)
      nlstresses(ipp, im, ido);

    Mm->givestress(0, ipp, sig);                   // actual stresses
    q(0) = Mm->ip[ipp].eqother[ido+n+1];           // actual hardening parameter
    double gammap = Mm->ip[ipp].eqother[ido+n+5];  // previous consistency parameter
    double dgamma = gamma - gammap;                // consistency parameter increment  
 
    aux1 = yieldfunction(sig, q);    
    if ((aux1 < -1.0) || (dgamma == 0.0))
    {
      copym(de, d);
      return;
    }
    else
    {
      dderyieldfsigma(ddfdst); // ddfdst = dd_f/dd_sigma
      cmulm(gamma, de, tmp);    // tmp = \gamma \mbf{D}_e
      mxm(tmp, ddfdst, auxm);    // auxm = tmp dd_f/dd_sigma
      for (i = 0; i < n; i++)
        auxm(i,i) += 1.0;      // auxm = (\mbf{I} + \gamma \mbf{D}_e dd_f/dd_sigma)
      invm(auxm, tmp, zero);    // tmp = auxm^(-1)
      mxm(tmp, de, auxm);       // auxm = \mbf{R} = tmp \mbf{D}_e
      //  matrix representation of the fourth order tensor ddfdst
      tensor4_matrix(rm, auxm, ssst);
    }

    // transform derivatives from the tensor to the vector form
    deryieldfsigma(sig, q, auxt);
    tensor_vector(dfds, auxt, ssst, stress);
    copyv(dfds, dgds);

    //                              \mbf{R} (d_g/d_sigma) (d_f/d_sigma)^T \mbf{R}
    //D_{epc} = \mbf{D}_e - -------------------------------------------------------------------
    //                       (d_f/d_sigma)^T \mbf{R} (d_g/d_sigma) - (d_f/d_q)^T (d_q/d_gamma)
    
    // nominator calculation
    mxv(auxm, dgds, auxv); // auxv = \mbf{R} (d_g/d_sigma)
    vxm(dfds, auxm, tmpv); // tmpv = (d_f/d_sigma) \mbf{R}
    vxv(auxv, tmpv, tmp);  // tmp  = auxv tmpv

    // denominator calculation
    mxv(auxm, dgds, auxv);   // auxv =  \mbf{R} (d_g/d_sigma)
    scprd(dfds, auxv, aux1); // aux1 =  (d_f/d_sigma)^T auxv
    aux1 += plasmodscalar(ipp, ido, sig, q); // aux1 += -(d_f/d_q)^T (d_q/d_gamma)

    cmulm(1.0/aux1, tmp, tmp); // tmp = tmp / aux1

    // d = D_{epc} = \mbf{D}_e - tmp
    // this matrix subtraction must be calculated with respect to dimensions
    // resulting matrix because in the case of plane-strain stress state, the
    // d matrix has dimensions (3x3) but de and tmp have dimensions (4x4).
    tensor4_matrix(rm, de, ssst);
    for (i=0; i<d.m; i++)
    {
      for(j=0; j<d.n; j++)
        d(i,j) = rm(i,j) - tmp(i,j);
    }
  }
  */
}



/**
  This function computes stresses at given integration point ipp,
  depending on the reached strains.
  The cutting plane algorithm is used. The stress and the other attribute of
  given integration point is actualized.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  Created by Tomas Koudelka,  2005
*/
void camclay::nlstresses (long ipp, long im, long ido)
{
  int errr;
  long i,ni,ncomp=Mm->ip[ipp].ncompstr;
  double gamma,err;
  vector epsn(ASTCKVEC(ncomp)),epsp(ASTCKVEC(ncomp)),q(ASTCKVEC(1));
  vector sig(ASTCKVEC(ncomp));
  vector epst(ASTCKVEC(6)), epspt(ASTCKVEC(6)),sigt(ASTCKVEC(6));
  double i1s, j2s, epsv, epsvp;
  double dtime = Mp->timecon.actualforwtimeincr();

  //  initial values
  for (i=0; i<ncomp; i++)
  {
    epsn(i) = Mm->ip[ipp].strain[i];
    epsp(i) = Mm->ip[ipp].eqother[ido+i];
    // It is necessary to substract initial strains in order to invoke initial nonzero stress 
    // state for zero displacements (i.e. zero strains)
    // sig = E*(eps - eps_p), sig_0 = E*eps_0. For eps = 0: sig = sig_0  => eps_p = -eps_0
    // epsp[i] -= Mm->ic[ipp][3+i];
  }

  gamma = Mm->ip[ipp].eqother[ido+ncomp];
  q(0)  = Mm->ip[ipp].eqother[ido+ncomp+1];

  //  stress return algorithm
  switch(sra.give_tsra ())
  {
    case cp:
      ni=sra.give_ni ();
      err=sra.give_err ();
      Mm->cutting_plane (ipp, im, ido, gamma, epsn, epsp, q, ni, err);
      break;
    case gsra:
      ni=sra.give_ni ();
      err=sra.give_err ();
      Mm->newton_stress_return (ipp, im, ido, gamma, epsn, epsp, q, ni, err);
      break;
    default:
      print_err("wrong type of stress return algorithm is required", __FILE__, __LINE__, __func__);
      abort ();
  }
  // Restore plastic strains related to the initial strain state
  for (i=0; i<ncomp; i++)
  {
    //epsp[i] += Mm->ic[ipp][3+i];    
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  
  //  new data storage
  Mm->ip[ipp].other[ido+ncomp]=gamma;
  Mm->ip[ipp].other[ido+ncomp+1]=q[0];
  // other[ido+ncomp+2] = v_lambda1 is the initial value and therefore it does not require actualization
  // other[ido+ncomp+3] = v_ini is the initial value and therefore it does not require actualization
  Mm->givestress(0, ipp, sig);
  give_full_vector (sigt,sig,Mm->ip[ipp].ssst);
  give_full_vector (epst,epsn,Mm->ip[ipp].ssst);
  give_full_vector (epspt,epsp,Mm->ip[ipp].ssst);

  i1s = first_invar (sigt)/3.0;
  //  the second invariant of deviator
  //  it is expressed with the help of the stress components
  //  components of the deviator are not needed
  j2s = sqrt(j2_stress_invar (sig));

  epsv = first_invar (epst);
  epsvp = first_invar (epspt);
  Mm->ip[ipp].other[ido+ncomp+4] = i1s;  
  Mm->ip[ipp].other[ido+ncomp+5] = j2s;
  //Mm->ip[ipp].other[ido+ncomp+5] = give_actual_ym(ipp, im, ido);  
  Mm->ip[ipp].other[ido+ncomp+6] = epsv;  
  Mm->ip[ipp].other[ido+ncomp+7] = epsvp;
  Mm->ip[ipp].other[ido+ncomp+8] = 1.0;
  // volumetric strain rate
  Mm->ip[ipp].other[ido+ncomp+9] = (epsv - Mm->ip[ipp].eqother[ido+ncomp+6])/dtime;
  // the volumetric strain rate is computed via generalized trapesoidal rule
  //Mm->ip[ipp].other[ido+ncomp+9] = 0.5*((epsv - Mm->ip[ipp].eqother[ido+ncomp+6])/dtime + Mm->ip[ipp].eqother[ido+ncomp+9]);
  // reference presure = p1
  double p1 = Mm->ic[ipp][1];
  // specific volume at the reference pressure p_1 = v_lambda1 
  double v_lambda1 = Mm->ip[ipp].eqother[ido+ncomp+2];
  // specific volume for the actual preconsolidation pressure = v_pc
  double v_pc = v_lambda1 - lambda * log(q[0]/p1);
  // specific volmue for the actual pressure p
  double va;
  if (i1s < 0.0)
    va = v_pc + kappa*log(q[0]/i1s);
  else
    va = v_pc + kappa*log(-q[0]);
  // actual porosity n
  Mm->ip[ipp].other[ido+ncomp+10] = (va - 1.0)/va;  // va = 1.0 + ea

  errr = check_math_errel(Mm->elip[ipp]);
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

*/
void camclay::updateval (long ipp,long ido)
{
  long i,n = Mm->ip[ipp].ncompstr;

  Mm->ip[ipp].eqother[ido+n+5] = Mm->ip[ipp].eqother[ido+n];     // store previous gamma;  

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
  Mm->ip[ipp].eqother[ido+n]    = Mm->ip[ipp].other[ido+n];       // gamma
  Mm->ip[ipp].eqother[ido+n+1]  = Mm->ip[ipp].other[ido+n+1];     // pc (hardening) 
  //  Mm->ip[ipp].eqother[ido+n+2] = Mm->ip[ipp].other[ido+n+2];     // v_lambda1 - initial value
  //  Mm->ip[ipp].eqother[ido+n+3] = Mm->ip[ipp].other[ido+n+3];     // v_ini - initial value
  Mm->ip[ipp].eqother[ido+n+4]  = Mm->ip[ipp].other[ido+n+4];     // i1s;  
  Mm->ip[ipp].eqother[ido+n+5]  = Mm->ip[ipp].other[ido+n+5];     // j2s;  
  Mm->ip[ipp].eqother[ido+n+6]  = Mm->ip[ipp].other[ido+n+6];     // epsv;  
  Mm->ip[ipp].eqother[ido+n+7]  = Mm->ip[ipp].other[ido+n+7];     // epsvp;    
  Mm->ip[ipp].eqother[ido+n+8]  = Mm->ip[ipp].other[ido+n+8];     // degree of saturation Sr = 1.0, it should be constant in this model
  Mm->ip[ipp].eqother[ido+n+9]  = Mm->ip[ipp].other[ido+n+9];     // depsv/dt  
  Mm->ip[ipp].eqother[ido+n+10] = Mm->ip[ipp].other[ido+n+10];     // porosity n = Vp/V
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  Returns vector of irreversible strains via parameter epsp
*/
void camclay::giveirrstrains (long ipp, long ido, vector &epsp)
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

*/
double camclay::give_consparam (long ipp,long ido)
{ 
  long ncompstr;
  double gamma;

  ncompstr = Mm->ip[ipp].ncompstr;
  gamma    = Mm->ip[ipp].eqother[ido+ncompstr];

  return gamma;
}



/**
  The function extracts saturation degree s from the integration point eqother array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of saturation degree.

  Created by Tomas Krejci, 14/12/2014
*/
double camclay::give_saturation_degree(long ipp, long ido)
{
  long ncompstr;
  double s;

  ncompstr = Mm->ip[ipp].ncompstr;
  s        = Mm->ip[ipp].other[ido+ncompstr+8];
  
  return s;
}

/**
  The function extracts preconsolidation pressure p_c for the attained equilibrium state
  from the integration point eqother array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of preconsolidation pressure.

  Created by Tomas Koudelka, 8.10.2013  
*/
double camclay::give_preconspress(long ipp, long ido)
{
  long ncompstr;
  double pc;

  ncompstr = Mm->ip[ipp].ncompstr;
  pc       = Mm->ip[ipp].other[ido+ncompstr+1];
  
  return pc;
}



/**
  The function extracts virgin void ratio e_lambda1 for the last equilibrium state 
  from the integration point eqother array. This value remains constant from initialization of the model.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of virgin porosity e_lambda1.

  Created by Tomas Koudelka, 8.10.2013  
*/
double camclay::give_virgporosity(long ipp, long ido)
{
  long ncompstr;
  double e_lambda1;

  ncompstr  = Mm->ip[ipp].ncompstr;
  e_lambda1 = Mm->ip[ipp].eqother[ido+ncompstr+2]-1.0;
  
  return e_lambda1;
}



/**
  The function extracts initial porosity n_ini for the last equilibrium  state
  from the integration point other array. This value remains constant from initialization of the model.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of initial porosity e_ini.

  Created by Tomas Koudelka, 8.10.2013  
*/
double camclay::give_iniporosity(long ipp, long ido)
{
  long ncompstr;
  double v_ini, n_ini;

  ncompstr = Mm->ip[ipp].ncompstr;
  v_ini    = Mm->ip[ipp].eqother[ido+ncompstr+3];
  n_ini    = (v_ini - 1.0)/v_ini;
  
  return n_ini;
}



/**
  The function extracts porosity for the actual state from the integration point other array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of actual porosity n.

  Created by Tomas Koudelka, 06.2020  
*/
double camclay::give_porosity(long ipp, long ido)
{
  long ncompstr;
  double n;

  ncompstr = Mm->ip[ipp].ncompstr;
  n        = Mm->ip[ipp].other[ido+ncompstr+10];
  
  return n;
}



/**
  The function returns rate of the volumetric strain for the actual state 
  at the given integration point.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns rate of the volumetric strain.
  
  Created by Tomas Koudelka 06.2020
*/
double camclay::give_strain_vol_rate(long ipp, long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;

  return Mm->ip[ipp].other[ido+ncomp+9];
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by Tomas Koudelka,
*/
void camclay::changeparam (atsel &atm,vector &val)
{
  long i;
  
  for (i=0;i<atm.num;i++)
  {
    switch (atm.atrib[i])
    {
      case 0:
        m=val[i];
        break;
      case 1:
        lambda=val[i];
        break;
      case 2:
        kappa=val[i];
        break;
      default:
        print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
  }
}
