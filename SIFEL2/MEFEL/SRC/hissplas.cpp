#include "hissplas.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "iotools.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "vecttens.h"
#include "intpoints.h"
#include <stdlib.h>
#include <math.h>



/**
  Default constructor - sets material parameters to the 
  initial (zero) values
  
  Created by TKo 01.2010
*/
hissplas::hissplas ()
{
  beta = gamma_max = n = r = pa = alpha1 = eta1 = ksi_lim = k = 0.0;  
}



/**
  Default destructor
  
  Created by TKo 01.2010
*/
hissplas::~hissplas ()
{
}



/**
  The function reads material parameters of HISS plasticity from the 
  opened text file given by parameter in.

  @param in - pointer to the opened XFILE

  @return The function does not return anything.
  
  Created by TKo 02.2010
*/
void hissplas::read (XFILE *in)
{
  xfscanf(in, "%k%lf %k%lf %k%lf% k%lf% k%lf% k%lf% k%lf %k%lf", 
              "beta", &beta, "n", &n, "r", &r, "pa", &pa, 
          "alpha1", &alpha1, "eta1", &eta1, "ksi_lim", &ksi_lim, "gamma_max", &gamma_max);
  sra.read (in);
}



/**
  The function returns value of the yield function
  for the given stress state defined by stress invariants.
  Yield function is defined as f = \frac{J_2}{{p_a}^2} - F_b F_c = 0,
  where F_b = \left[ -\alpha\left(\frac{I_1+R}{p_a}\right)^n + \gamma \left(\frac{I_1+R}{p_a}\right)^2 \right]
        F_c = \left(1-\beta \cos 3\theta \right)^{-1/2}
        \cos 3\theta = \frac{3\sqrt{3}}{2} \frac{J_3}{{J_2}^{3/2}}

  @param i1   - first invariant of stress tensor I1
  @param j2   - second invariant of stress deviator J2
  @param cos3t - Lode's angle (cos3t = 3*sqrt(3)/2 * J3/sqrt(J2^2)
  @param q    - %vector of hardening parameters (i.e. alpha)

  @retval The function returns the value of the yield function

  Created by TKo 01.2010
*/
double hissplas::yieldfunction (double i1, double j2, double cos3t, vector &q)
{
  double f, fb, fc;
  double aux, alpha, gamma;

  // actual values of hardening parameters
  alpha = q[0];
  gamma = q[1];

  aux = (i1+r)/pa;
  fb = -alpha*pow(aux, n) + gamma*aux*aux;

  aux = 1.0-beta*cos3t;
  fc = 1.0/sqrt(aux);
  fc = 1.0;

  f = j2/(pa*pa)-fb*fc;  

  return f;
}



/**
  The function returns value of the yield function
  for the given stress state defined by stress tensor.
  Yield function is defined as f = \frac{J_2}{{p_a}^2} - F_b F_c = 0,
  where F_b = \left[ -\alpha\left(\frac{I_1+R}{p_a}\right)^n + \gamma \left(\frac{I_1+R}{p_a}\right)^2 \right]
        F_c = \left(1-\beta \cos 3\theta \right)^{-1/2}
        \cos 3\theta = \frac{3\sqrt{3}}{2} \frac{J_3}{{J_2}^{3/2}}

  @param sig - stress components in Voigt notation
  @param q    - %vector of hardening parameters (i.e. alpha)

  @retval The function returns the value of the yield function

  Created by TKo 01.2010
*/
double hissplas::yieldfunction (vector &sig, vector &q)
{

  double i1, j2, j3, cos3t, f;
  vector dev(ASTCKVEC(6));

  //
  // Compuation of stress invariants
  //
  i1    = first_invar(sig);
  deviator(sig, dev);
  j2    = j2_stress_invar(sig);
  j3    = third_stress_invar(dev);
  cos3t = 3.0*sqrt(3.0)/2.0 * j3 /pow(j2, 3.0/2.0);

  f = yieldfunction(i1, j2, cos3t, q);

  return f;
}



/**
  The function returns value of expression g1 = df/d[sqrt(J2)] * depsvp - df/dI1 * depsdp.
  The above expression represents condition of modified flow rule.

  @param depsvp - increment of volumetric plastic strain
  @param depsdp - increment of equivalent deviatoric plastic strain
  @param i1 - first invariant of stress tensor
  @param j2 - second invariant of trial stress tensor
  @param alpha - hardening parameter
  @param gamma - softening parameter
*/
double hissplas::g1function(double depsvp, double /*depsdp*/, double i1, double j2, double alpha, double gamma)
{
  double ret, aux;

  // computation of g1 function
  ret  = 2.0/(pa*pa) * sqrt(j2) * depsvp;
  aux  = (i1+r)/pa;
  ret -= -(1.0/pa)*(-n*alpha*pow(aux, n-1.0) + 2.0*gamma*aux);

  return ret;
}



/**
  The function computes a derivative of the yield function with respect to the plastic 
  volumetric strain for the given stress state.

  @param ipp  - integration point pointer
  @param i1   - first invariant of stress tensor
  @param alpha - hardening parameter
  @param gamma - softening parameter
  @param ksi   - norm of plastic strain tensor
  @param depsp - tensor of plastic strain increments for the actual step in Voigt notation

  @retval The function returns a derivative of yield function with respect to plastic volumetric strain.

  Created by TKo 01.2010
*/
double hissplas::dfdepsvp (long ipp, double i1, double alpha, double gamma, double ksi, vector &depsp)
{
  //  df/ddeps_vp = (df/dFa*dFa/ddeps_vp) + (df/dFb*dFb/ddeps_vp)
  //                     X1                       X2
  //
  //  X1 = df/dFa*dFa/d[sqrt(J2)]*d[sqrt(J2)]/ddeps_vp = 0
  //  X2 = df/dFb*(dFb/dI1*dI1/ddeps_vp + dFb/dalpha*dalpha/ddeps_vp)
  //                     A                         B  
  //  
  //  F_b = (-alpha*((I1+R)/p_a)^n + gamma*((I1+R)/p_a)^2)
  //  alpha = alpha1*(ksi_lim/(1+ksi^eta1)-ksi/(1+ksi^eta1))  according to Erkens and Liu

  double aux, auxA, auxB, ret;
  // Young modulus
  double e = Mm->give_actual_ym(ipp);
  // Poissons ratio
  double nu = Mm->give_actual_nu(ipp);
  // bulk modulus
  double k = e/(3.0*(1-2*nu));

  auxA = (1.0/pa)*(-n*alpha*pow((i1+r)/pa, n-1.0)+2.0*gamma*(i1+r)/pa);  // dFb/dI1
  auxA *= -9.0*k;   // dI1/ddeps_vp

  auxB = -pow((i1+r)/pa, n);  // dFb/dalpha
  // dalpha/ddeps_vp = dalpha/dksi * dksi/ddeps_p * ddeps_p/ddeps_vp
  aux = (1+pow(ksi,eta1)); 
  auxB *= -alpha1/(aux*aux)*(pow(ksi,eta1-1.0)*(eta1*ksi_lim + ksi*(1-eta1))+1.0); // dalpha/dksi
  auxB *= first_invar(depsp)/tensor_strain_norm(depsp);            // dksi/ddeps_p * ddeps_p/ddeps_vp

  aux = -1.0;  // df/dFb = -Fc
  ret = aux*(auxA + auxB); // df/dFb * (A + B)

  return ret;
}



/**
  The function returns derivatives of the yield function with respect to the plastic 
  deviatoric strain for the given stress state.

  @param ipp  - integration point pointer
  @param i1 - first invariant of trial stress tensor
  @param j2 - second invariant of trial stress tensor
  @param dev - deviator of trial stress second order tensor in Voigt notation
  @param ksi   - norm of plastic strain tensor
  @param depsp - tensor of plastic strain increments for the actual step in Voigt notation

  @retval The function returns a derivative of yield function with respect to plastic deviatoric strain.

  Created by TKo 01.2010
*/
double hissplas::dfdepspd (long ipp, double i1, double j2, vector &dev, double ksi, vector &depsp)
{
  // df/ddeps_dp = (df/dFa*dFa/ddeps_dp) + (df/dFb*dFb/ddeps_dp)
  //                     X1                       X2
  //
  //  X1 = df/dFa*dFa/d[sqrt(J2)]*d[sqrt(J2)]/ddeps_dp
  //  X2 = df/dFb*(dFb/dI1*dI1/ddeps_dp + dFb/dalpha*dalpha/ddeps_dp) =
  //     = (df/dFb)*(dFb/dalpha)*(dalpha/dksi)*(dksi/ddeps_p)*(ddeps_p/ddeps_dp)


  double aux, aux1, aux2, ret;
  // Young modulus
  double e = Mm->give_actual_ym(ipp);
  // Poissons ratio
  double nu = Mm->give_actual_nu(ipp);
  // shear modulus
  double g = e/(2.0*(1.0+nu));

  // df/dFa
  aux1 = 1.0;
  // dFa/d[sqrt(J2)]
  aux1 *= 2.0*sqrt(j2)/(pa*pa);
  // d[sqrt(J2)]/ddeps_dp
  aux1 *= -g/2.0*tensor_dbldot_prod(dev, dev, 2.0)/j2;

  // df/dFb
  aux2 = -1.0; 
  // dFb/dalpha
  aux2 *= -pow((i1+r)/pa, n);
  // dalpha/ddeps_dp = dalpha/dksi * dksi/ddeps_p * ddeps_p/ddeps_dp
  aux = (1+pow(ksi,eta1)); 
  aux2 *= -alpha1/(aux*aux)*(pow(ksi,eta1-1.0)*(eta1*ksi_lim + ksi*(1-eta1))+1.0); // dalpha/dksi
  // dksi/ddepsp = depsp/sqrt(depsp:depsp)
  // ddepsp/ddeps_dp = 1/(2*sqrt(j2))*dev
  aux2 *= tensor_dbldot_prod(depsp, dev, 1.0);
  aux2 /= 2.0*sqrt(j2)*tensor_strain_norm(depsp);
  
  ret = aux1+aux2;

  return ret;
}



/**
  The function returns derivatives of g1 = (df/d[sqrt(j2)])*deps_vp - df/dI1*deps_dp with 
  respect to increment of plastic volumetric strain.

  @param ipp  - integration point pointer
  @param i1 - first invariant of trial stress tensor
  @param j2 - second invariant of trial stress tensor
  @param alpha - hardening parameter
  @param gamma - softening parameter
  @param dev - deviator of trial stress second order tensor in Voigt notation
  @param depsdp - increment of equivalent deviatoric plastic strain in Voigt notation
  @param ksi   - norm of plastic strain tensor
  @param depsp - tensor of plastic strain increments for the actual step in Voigt notation

  @retval The function returns a derivative of yield function with respect to plastic deviatoric strain.

  Created by TKo 01.2010
*/
double hissplas::dg1depsvp (long ipp, double i1, double j2, double alpha, double gamma, vector &/*dev*/, double depsdp, double ksi, vector &depsp)
{
  // dg1/ddeps_vp = d/ddeps_vp{(df/d[sqrt(J2)])*deps_vp} - d/ddeps_vp{(df/dI1*deps_dp)}
  //                         X1                                        X2
  //
  //  X1 = df/d[sqrt(J2)] + deps_vp * d/ddeps_vp{df/d[sqrt(J2)]}
  //     = 2.0/pa^2 * sqrt(J2) +      0
  // 
  //  X2 = d/ddeps_vp{df/dI1 * deps_dp}
  //     = -deps_dp/pa^n * d/deps_vp{-n*alpha*(I1+r)^(n-1) + 2.0*gamma*pa^(n-2)*(I1+r)}
  //                                            A                        B 


  double ddfdsqrtJ2ddeps_vp;
  double aux, auxA, auxB, ret;
  // Young modulus
  double e = Mm->give_actual_ym(ipp);
  // Poissons ratio
  double nu = Mm->give_actual_nu(ipp);
  // bulk modulus
  double k = e/(3.0*(1-2*nu));

  // d/ddeps_vp{df/d[sqrt(J2)] * depsp_vp}
  ddfdsqrtJ2ddeps_vp = 2.0/(pa*pa)*sqrt(j2);

  // A = d/deps_vp{-n*alpha*(I1+r)^(n-1)} =
  //   = -n*{alpha*(I1+r)^(n-1)*dalpha/ddeps_vp + alpha * d/ddeps_vp{(I1+r)^(n-1)}
  //                     A1                                 A2

  // computation of A1
  auxA = pow(i1+r, n-1.0);
  // dalpha/ddeps_vp = dalpha/dksi * dksi/ddeps_p * ddeps_p/ddeps_vp
  aux = (1+pow(ksi,eta1)); 
  auxA *= -alpha1/(aux*aux)*(pow(ksi,eta1-1.0)*(eta1*ksi_lim + ksi*(1-eta1))+1.0); // dalpha/dksi
  // dksi/ddepsp = depsp/sqrt(depsp:depsp)
  // ddepsp/ddeps_vp = delta_ij
  auxA *= first_invar(depsp)/tensor_strain_norm(depsp);  // dksi/ddeps_p * ddeps_p/ddeps_vp

  // computation of A2
  auxA += -alpha*9.0*k*(n-1.0)*pow(i1+r, n-2.0);
  
  auxA *= -n;

  // B = d/deps_vp{2.0*gamma*pa^(n-2)*(I1+r)} =
  //   = 2.0*gamma * d/ddeps_vp{(I1+r)} = -18*gamma*k*pa^(n-2)
  auxB = -18.0*gamma*k*pow(pa, n-2.0);

  ret = ddfdsqrtJ2ddeps_vp-depsdp/pow(pa, n)*(auxA+auxB);

  return ret;
}



/**
  The function returns derivatives of g1 = (df/d[sqrt(j2)])*deps_vp - df/dI1*deps_dp with 
  respect to increment of plastic deviatoric strain.

  @param ipp  - integration point pointer
  @param i1 - first invariant of trial stress tensor
  @param j2 - second invariant of trial stress tensor
  @param alpha - hardening parameter
  @param gamma - softening parameter
  @param dev - deviator of trial stress second order tensor in Voigt notation
  @param depsvp - increment of volumetric plastic strains
  @param depsdp - increment of deviatoric 
  @param ksi   - norm of plastic strain tensor
  @param depsp - tensor of plastic strain increments for the actual step in Voigt notation

  @retval The function returns a derivative of yield function with respect to plastic deviatoric strain.

  Created by TKo 01.2010
*/
double hissplas::dg1depsdp (long ipp, double i1, double j2, double alpha, double gamma, vector &dev, double depsvp, double depsdp, double ksi, vector &depsp)
{
  // dg1/ddeps_dp = d/ddeps_dp{(df/d[sqrt(J2)])*deps_vp - df/dI1*deps_dp} =
  //              = d/ddeps_dp{df/d[sqrt(J2)]*deps_vp} - d/ddeps_dp{df/dI1*deps_dp}
  //                         A                                   B
  // A = d/ddeps_dp{df/d[sqrt(J2)]*deps_vp} =
  //   = deps_vp * d/ddeps_dp{df/d[sqrt(J2)]} + df/d[sqrt(J2)]*0 =
  //   = deps_vp * 2.0/(pa*pa)*d[sqrt(J2)]/ddeps_dp =
  //   = deps_vp * 2.0/(pa*pa)*(-G/2*(dev:dev)/j2
  //
  // B = d/ddeps_dp{df/dI1*deps_dp} =
  //   = df/dI1 + deps_dp * d/ddeps_dp{df/dI1} = 
  //   = df/dI1 + deps_dp * (-1/pa^n * (-n * d/ddeps_dp{alpha*(I1+r)^(n-1)} + 2*gamma*d/ddeps_dp{I1+r})) =
  //   =   X1   + deps_dp * (n/pa^n  * d/ddeps_dp{alpha*(I1+r)^(n-1)}       +        0) =
  //   =   X1   + deps_dp * (n/pa^n  * (I1+r)^(n-1)* dalpha/dksi * dksi/ddepsp * ddepsp/ddeps_dp)

  double aux, auxA, auxB, ret;

  // Young modulus
  double e = Mm->give_actual_ym(ipp);
  // Poissons ratio
  double nu = Mm->give_actual_nu(ipp);
  // shear modulus
  double g = e/(2.0*(1.0+nu));

  auxA  = depsvp*2.0/(pa*pa); // df/d[sqrt(J2)]  
  auxA *= -g/2.0*tensor_dbldot_prod(dev, dev, 2.0)/j2;

  auxB  = n/pow(pa, n);
  auxB *= pow(i1+r, n-1.0);
  aux = (1+pow(ksi,eta1)); 
  auxB *= -alpha1/(aux*aux)*(pow(ksi,eta1-1.0)*(eta1*ksi_lim + ksi*(1-eta1))+1.0); // dalpha/dksi
  auxB *= tensor_dbldot_prod(depsp, dev, 1.0)/2.0/sqrt(j2);
  auxB /= tensor_strain_norm(depsp);
  auxB *= depsdp;
  auxB += -1.0/pa*(-n*alpha*pow((i1+r)/pa, n-1) + 2.0*gamma*(i1+r)/pa); //X1

  ret = auxA-auxB;

  return ret;
}



/**
  The function returns material siffness %matrix.

  @param d[out]   - material stiffness %matrix 
  @param ipp[in] - integration point id
  @param ido[in] - index of internal variables for given material in the ipp other array

  @retval The function returns resulting stiffness %matrix in the output parameter d.

  Created by TKo 01.2010
*/
void hissplas::matstiff (matrix &d, long ipp, long ido)
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
//    Mm->elmatstiff (ad,ipp);
//    tangentstiff (ad,d,ipp,ido);
    Mm->elmatstiff(d, ipp, ido);
  }
}



/**
  The function computes stresses in the given integration point.
  
  @param ipp - integration point pointer
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns stresses computed from the actual strains 
          and nonlocal averaged values.

  Created by TKo 01.2010
*/
void hissplas::nonloc_nlstresses (long /*ipp*/, long /*ido*/)
{
}



/**
  @param ipp - integration point pointer
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function computes stresses for the actual strains and stores them in the stress 
          array at the given integration point.

  Created by TKo 01.2010
*/
void hissplas::nlstresses (long ipp,long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;  // number of stress/strain components
  strastrestate ssst = Mm->ip[ipp].ssst; // stress/strain state
  long eqstate;          // equilibrium state flag
  long i;
  vector epsn(ASTCKVEC(ncomp));   // total strains
  vector epsp(ASTCKVEC(ncomp));   // plastic strains
  vector depsp(ASTCKVEC(ncomp));  // plastic strain increments
  vector ddepsp(ASTCKVEC(ncomp)); // correction of plastic strain increments
  vector epsa(ASTCKVEC(ncomp));   // auxiliary/elastic strain vector
  vector sig(ASTCKVEC(ncomp));    // stress vector 
  vector dsig(ASTCKVEC(ncomp));   // correction of stresses
  vector qo(ASTCKVEC(4));         // vector of previous values of hardening parameters
  vector q(ASTCKVEC(4));          // vector of hardening parameters
  vector epspt(ASTCKVEC(6));      // tensor of plastic strains
  vector depspt(ASTCKVEC(6));     // tensor of plastic strain increments
  vector ddepspt(ASTCKVEC(6));    // tensor of corrections of plastic strain increments
  vector sigt(ASTCKVEC(6));       // stress tensor
  vector dsigt(ASTCKVEC(6));      // tensor of correction of stresses 
  vector dev(ASTCKVEC(6));        // trial stress deviator
  matrix d(ASTCKMAT(ncomp,ncomp)); // material stiffness matrix
  matrix h(ASTCKMAT(2,2));         // Hessian matrix
  vector e(ASTCKVEC(2));           // vector of unknowns
  vector r(ASTCKVEC(2));           // vector of residua
  double f;              // actual value of yield function
  double g1;             // actual value of g1 function
  double lambda;         // actual value of consistency parameter
  double norsig;         // norm of stress tensor 
  double nordsig;        // norm of tensor of stress corrections
  double norepsp;        // norm of plastic strains
  double norddepsp;      // norm of corrections of plastic strain increments
  double pnordsig;       // proportional norm of stress tensor corrections  (nordsig/norsig)
  double pnorddepsp;     // proportional norm of plastic strain increments corrections (norddepsp/norepsp)
  double i1;             // first invariant of stresses 
  double j2, j3, cos3t;  // deviatoric stress invariants 
  double depsvp;         // increment of volumetric plastic strain
  double depsdp;         // increment of equivalent deviatoric plastic strain
  double &alpha = q[0];  // reference to actual value of alpha
  double &gamma = q[1];  // reference to actual value of gamma
  double &ksi   = q[2];  // reference to actual value of ksi
  double &wpf   = q[3];  // reference to actual value of Wpf
  

  //
  //  Initial values
  //
  for (i=0; i<ncomp; i++)
  {
    epsn[i]=Mm->ip[ipp].strain[i];             // total strains
    epsp[i]=Mm->ip[ipp].eqother[ido+i];        // previous plastic strains
  }
  lambda = Mm->ip[ipp].eqother[ido+ncomp];     // initial value of consistency parameter
  qo[0]    = Mm->ip[ipp].eqother[ido+ncomp+1]; // previous value of hardening parameter alpha
  qo[1]    = Mm->ip[ipp].eqother[ido+ncomp+2]; // previous value of hardening parameter gamma
  qo[2]    = Mm->ip[ipp].eqother[ido+ncomp+3]; // previous value of cumulative plastic strain ksi
  qo[3]    = Mm->ip[ipp].eqother[ido+ncomp+4]; // previous value of cumulative post-failure plastic work
  copyv(qo, q); 

  //
  // Computation of trial stress state
  //
  Mm->elmatstiff(d, ipp);
  subv(epsn, epsp, epsa);
  mxv(d,epsa,sig);
  give_full_vector(sigt, sig, ssst);

  //
  // Compuation of stress invariants
  //
  i1    = first_invar(sigt);
  deviator(sigt, dev);
  j2    = j2_stress_invar(sigt);
  j3    = third_stress_invar(dev);
  cos3t = 3.0*sqrt(3.0)/2.0 * j3 /pow(j2, 3.0/2.0);

  //
  // Test of yield function
  //
  eqstate = 1;
  f = yieldfunction(i1, j2, cos3t, q);
  if (f > sra.err)
  {
    //
    // Initialization of plastic strain dependent quantities
    //
    give_full_vector(epspt, epsp, ssst);
    depsvp = 0.0;
    depsdp = 0.0;

    eqstate = 0;    

    if (alpha > 0.0)    
    {
      //
      // Hardening phase
      //
      g1 = g1function(depsvp, depsdp, i1, j2, alpha, gamma);
      for (i=0; i<sra.ni; i++)
      {
        // assembling Hessian matrix
        h[0][0] = dg1depsvp(ipp, i1, j2, alpha, gamma, dev, depsdp, ksi, depspt);
        h[0][1] = dg1depsdp(ipp, i1, j2, alpha, gamma, dev, depsvp, depsdp, ksi, depspt);
        h[1][0] = dfdepsvp (ipp, i1, alpha, gamma, ksi, depspt);
        h[1][1] = dfdepspd (ipp, i1, j2, dev, ksi, depspt);
        
        // assembling vector of unknowns
        e[0] = depsvp; 
        e[1] = depsdp;
        // assembling vector of residua
        r[0] = 0.0 - g1;
        r[1] = 0.0 - f; 
      
        // solve system of equations
        gause(h,h,e,r,1.0e-10);
      
        // norm of stresses and plastic strains
        norsig  = tensor_stress_norm(sigt);
        norepsp = tensor_strain_norm(epsp);

        // new corrections of plastic strain increments
        copyv(dev, ddepspt);
        cmulv(r[1]/(2.0*sqrt(j2)), ddepspt);
        ddepspt(0) = r[0];
        ddepspt(1) = r[0];
        ddepspt(2) = r[0];
        give_red_vector(ddepspt, ddepsp, ssst);

        // new plastic strain increments
        depsvp += r[0];
        depsdp += r[1];
        depspt(0) = depspt(1) = depspt(2) = depsvp;
        cmulv(depsdp/(2.0*sqrt(j2)), dev);
        addv(depspt, dev, depspt);
        give_red_vector(depspt, depsp, Mm->ip[ipp].ssst);

        // new plastic strains
        addv(epsp, depsp, epsp); 
        // new correction of stresses
        mxv(d, ddepsp, dsig);
        // new stresses
        subv(epsn, epsp, epsa);
        mxv(d, epsa, sig);
        give_full_vector(sigt, sig, ssst);
        // new stress invariants
        i1 = first_invar(sigt);
        deviator(sigt, dev);
        j2 = j2_stress_invar(sigt);
        j3 = third_stress_invar(dev);
        cos3t = 3.0*sqrt(3.0)/2.0 * j3 /pow(j2, 3.0/2.0);
        // update of hardening parameters
        copyv(qo, q);
        updateq(ipp, depspt, sigt, depspt, q);

        // new value of yieldfunction
        f = yieldfunction(i1, j2, cos3t, q);
        // new value of g1 function
        g1 = g1function(depsvp, depsdp, i1, j2, alpha, gamma);
        // norm of correction of plastic strain increment      
        norddepsp = tensor_strain_norm(ddepspt);
        // norm of stress increment      
        give_full_vector(dsigt, dsig, ssst);  
        nordsig = tensor_stress_norm(dsigt);

        // test of convergence
        pnordsig  = nordsig/norsig;
        pnorddepsp = norddepsp/norepsp;
        if ((pnordsig < sra.err) && (pnorddepsp < sra.err) && (f < sra.err) && (fabs(g1) < sra.err))
        {
          // new consistency parameter computed from the modified flow rule 
          // for deviatoric equivalent plastic strain
          // depsdp = dlambda * df/d[sqrt(J2)] = dlambda * 2/pa^2 * sqrt(J2)
          // dlambda = depsdp * pa^2 / (2*sqrt(J2))
          lambda = depsdp/(2.0*sqrt(j2)) * pa * pa;
          eqstate = 1;
          break;
        }
      }
    }
    else
    {
      //
      // alpha == 0 => Softening phase (Degradation)     
      //
    }
  }
  if (!eqstate) // convergency problem
  {
    print_err("Stress return algorithm does not converge on element %ld (ipp %ld)", 
              __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
  }

  //
  // Storage of attained internal variables
  //
  // storage of resulting stresses
  Mm->storestress(0,ipp,sig);
  // storage of plastic strains
  for (i=0; i<ncomp; i++)
    Mm->ip[ipp].other[ido+i] = epsp[i];

  // storage of consistency parameter
  Mm->ip[ipp].other[ido+ncomp+0]   = lambda;
  // storage of hardening parameters
  Mm->ip[ipp].other[ido+ncomp+1] = alpha;
  Mm->ip[ipp].other[ido+ncomp+2] = gamma;
  Mm->ip[ipp].other[ido+ncomp+3] = ksi;
  Mm->ip[ipp].other[ido+ncomp+4] = wpf;
}



/**
  The function updates hardening parameters according to actual values of plastic strains
  and other parameters used for the hardening.

  @param ipp    - integration point number in the mechmat ip array.
  @param depspt - tensor of actual plastic strain increments in Voigt notation
  @param sigt   - stress tensor in Voigt notation
  @param depst  - tensor of actual strain increments in Voigt notation
  @param q      - %vector of updated values of hardening parameters (input/output)

  @retval The function returns updated %vector of hardening parameters 
          in the output parameter q.

  Created by TKo 01.2010
*/
void hissplas::updateq (long /*ipp*/, vector &depspt, vector &sigt, vector &depst, vector &q)
{
  double alpha, gamma, ksi, wpf;

  // hardening 
  ksi  = q[2]; // previous value of ksi
  // ksi = \int_0^t{sqrt(depsp:depsp)}dt
  ksi += tensor_strain_norm(depspt);
  if (ksi_lim >= ksi)
    alpha = alpha1*(ksi_lim - ksi)/(1+pow(ksi, eta1));
  else
    alpha = 0.0;

  // softening (degradation)
  if (q[0] == 0.0) // previuos value of alpha == 0.0
  {
    wpf = q[3]; // previous value of Wpf
    // Wpf = \int_0^t{sigma : depsp}dt 
    wpf += tensor_dbldot_prod(sigt, depst, 1.0); // post-failure work of plastic strain increments
    gamma = sqrt(gamma_max*gamma_max-k*pow(wpf,3.0/2.0));
  }
  else{
    gamma = gamma_max;
    wpf = q[3];
  }

  // update of hardening parameters
  q[0] = alpha;
  q[1] = gamma;
  q[2] = ksi;
  q[3] = wpf;
}



/**
  The function updates internal variables according to actual values of plastic strains
  and other parameters used for the hardening.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function does not return anything.

  Created by TKo 01.2010
*/
void hissplas::updateval (long ipp, long ido)
{
  long i,n = Mm->ip[ipp].ncompstr;

  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
  Mm->ip[ipp].eqother[ido+n]   = Mm->ip[ipp].other[ido+n+0];     // lambda
  Mm->ip[ipp].eqother[ido+n+1] = Mm->ip[ipp].other[ido+n+1];     // alpha
  Mm->ip[ipp].eqother[ido+n+2] = Mm->ip[ipp].other[ido+n+2];     // gamma
  Mm->ip[ipp].eqother[ido+n+3] = Mm->ip[ipp].other[ido+n+3];     // ksi
  Mm->ip[ipp].eqother[ido+n+4] = Mm->ip[ipp].other[ido+n+4];     // wpf
}



/**
  The function returns actual value of the irrevesible(plastic)
  strains.

  @param ipp  - integration point number in the mechmat ip array.
  @param ido  - index of internal variables for given material in the ipp other array
  @param epsp - %vector actual values of plastic strains

  @retval The function returns %vector of actual values of the irreversible strains in the output 
          parameter epsp.
   
  Created by TKo 01.2010
*/
void hissplas::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  long ncomp = Mm->ip[ipp].ncompstr;

  if (epsp.n != ncomp)
  {
    print_err("wrong number of components of plastic strain vector", __FILE__, __LINE__, __func__);
    abort();
  }
  for(i=0; i<epsp.n; i++)
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
}



/**
  The function returns actual value of the consistency parameter.

  @param ipp  - integration point number in the mechmat ip array.
  @param ido  - index of internal variables for given material in the ipp other array

  @retval The function returns actual value of the consistency parameter.
   
  Created by TKo 01.2010
*/
double hissplas::give_consparam (long ipp,long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;
  return Mm->ip[ipp].eqother[ido+ncomp];
}



/**
  The function changes values of selected material parameters to the new ones. It is used 
  for the stochastic calculations.

  @param atm  - selection of material parameters
  @param val  - %vector of new values

  @retval The function changes value of the material parameters. It does not return anything.
   
  Created by TKo 01.2010
*/
void changeparam (atsel &/*atm*/,vector &/*val*/)
{
}


