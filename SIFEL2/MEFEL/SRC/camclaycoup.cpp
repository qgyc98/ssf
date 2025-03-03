#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "camclaycoup.h"
#include "global.h"
#include "probdesc.h"
#include "mechmat.h"
#include "vecttens.h"
#include "intpoints.h"
#include "matrix.h"
#include "vector.h"
#include "tensor.h"
#include "elastisomat.h"
#include "mathem.h"
#include "nonlinman.h"


/**
  This constructor inializes attributes to zero values.
*/
camclaycoup::camclaycoup (void)
{
  m = 0.0;
  lambda = 0.0;
  kappa = 0.0;
  c1_tilde = 0.0;
  c2 = 0.0;
  ks = 0.0;
}



/**
  This destructor is only for the formal purposes.
*/
camclaycoup::~camclaycoup (void)
{

}



/**
  This function reads material parameters from the opened text file given
  by the parameter in.

  @param in - pointer to the opned text file

*/
void camclaycoup::read (XFILE *in)
{
  xfscanf (in, "%k%lf%k%lf%k%lf", "m", &m, "lambda", &lambda, "kappa", &kappa);
  xfscanf (in, "%k%lf%k%lf%k%lf", "c1_tilde", &c1_tilde, "c2", &c2, "ks", &ks);
  sra.read (in);
}

/**
   This function initializes material model data
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param im  - index of material type for given ip
   @param ido - index of internal variables for given material in the ipp other array

   12/06/2012 TKo
*/
void camclaycoup::initval(long ipp, long im, long ido)
{
  double v_ini, v_pc0, v_lambda1;
  double i1s, j2s, s;
  double v_kappa1, p1, bar_pc_0, pc_0;
  long i, ncompstr = Mm->ip[ipp].ncompstr;
  vector sig(ASTCKVEC(ncompstr)),sigt(ASTCKVEC(6)), q(ASTCKVEC(2));
  double err=sra.give_err ();
  // variables used for the calculation of initial value of preconsolidation pressure for the fully saturated state
  double a, b, fs, ksi, sr, qq, p_atm = 101325.0;

  if (Mm->ip[ipp].eqother[ido+ncompstr+1] == 0.0)  // actual value of p_c was not set
  {
    // initial value of specific volume under small pressure p_1 on the normal consolidation line (v_lambda1 = 1+N)
    v_lambda1 = Mm->ic[ipp][0];
    // initial reference presure
    p1        = Mm->ic[ipp][1];
    // initial effective preconsolidation pressure \bar{p}_{c0}}
    bar_pc_0  = Mm->ic[ipp][2];

    //  initial effective stresses
    if ((Mp->eigstrains == 4) || (Mp->eigstrains == 5)){
      for (i=0; i<ncompstr; i++)
        sig(i) = Mm->ip[ipp].stress[i] = Mm->eigstresses[ipp][i];
    }
    else{
      print_err("initial effective stresses (eigentsresses) are not defined on element %ld, ip=%ld,\n"
                " cannot determine specific volume for initial stress state", __FILE__, __LINE__, __func__, Mm->elip[ipp]+1, ipp+1);
      abort();
    }

    //vector_tensor (sig,sigt,Mm->ip[ipp].ssst,stress);
    give_full_vector (sigt, sig, Mm->ip[ipp].ssst);
    // initial value of mean stress
    i1s = first_invar(sigt)/3.0;
    q(0) = bar_pc_0;
    s  = Mm->givenonmechq(suction, ipp);  // suction pressure s
    sr = Mm->givenonmechq(saturation_degree, ipp);  // degree of saturation S_r
    //limitation of positive suction
    if(s >=0.0)
    s=0.0;
    s = fabs(s);
    q(1) = ks*s;
    // test if the initial state is elastic, if not calculate new value of inital effective preconsolidation pressure pc_0 so the OCR=1
    if (yieldfunction(sigt, q) > err){
      j2s = j2_stress_invar (sigt);
      bar_pc_0 = j2s/(m*m*(i1s-q(1))) + i1s;
    }

    //
    // Calculation of inital value of preconsolidation pressure for the fully saturated state pc_0
    //
    // suction function:
    //
    //                  s / p_{atm}
    // f(s) = 1 + -----------------------
    //            10.7 + 2.4(s / p_{atm})
    //
    fs = 1.0 + (s/p_atm) / (10.7 + 2.4*(s/p_atm)); 
    if (fs < 1.0)
      fs = 1.0;
    
    // bonding variable:
    //
    // \ksi = f(s) (1 - S_r)
    //
    ksi = (1-sr)*fs;

    // according to Gallipoli
    qq = 1.0 - c1_tilde*(1.0 - exp(c2*ksi)); // e/e_s
    a = (qq - 1.0)/(lambda*qq - kappa)*(v_lambda1 - 1.0);
    b = (lambda - kappa)/(lambda*qq - kappa);
    pc_0 = -pow(-bar_pc_0 /exp(a), 1.0/b);

    
    // compute initial value of specific volume
    v_pc0 = v_lambda1 - lambda*log(pc_0/p1);               

    // actual value of p_c = initial effective preconsolidtaion pressure
    Mm->ip[ipp].other[ido+ncompstr+1] = Mm->ip[ipp].eqother[ido+ncompstr+1] = bar_pc_0;
    Mm->ip[ipp].other[ido+ncompstr+2] = Mm->ip[ipp].eqother[ido+ncompstr+2] = q(1);
    // initial value of specific volume under small pressure p_1 on the swelling line (v_{\kappa l}) at the fully saturated state
    Mm->ip[ipp].other[ido+ncompstr+3] = Mm->ip[ipp].eqother[ido+ncompstr+3] = v_kappa1 = v_pc0 + kappa*log(pc_0/p1);
    // specific volume for intial stress state
    Mm->ip[ipp].other[ido+ncompstr+4] = Mm->ip[ipp].eqother[ido+ncompstr+4] = v_ini     = v_pc0 + kappa*log(pc_0/i1s);
    // intial mean stress
    Mm->ip[ipp].other[ido+ncompstr+5] = Mm->ip[ipp].eqother[ido+ncompstr+5] = i1s;
    // saturation degree
    Mm->ip[ipp].other[ido+ncompstr+9] = Mm->ip[ipp].eqother[ido+ncompstr+9] = Mm->givenonmechq(saturation_degree, ipp);  // degree of saturation S_r
    // initial depsv_dt
    Mm->ip[ipp].other[ido+ncompstr+10] = Mm->ip[ipp].eqother[ido+ncompstr+10] = 0.0;
    // porosity n for intial stress state
    Mm->ip[ipp].other[ido+ncompstr+11] = Mm->ip[ipp].eqother[ido+ncompstr+11] = (v_ini-1.0)/v_ini; // v_ini = 1.0 + e_ini
    // suction
    Mm->ip[ipp].other[ido+ncompstr+12] = Mm->ip[ipp].eqother[ido+ncompstr+12] = Mm->givenonmechq(suction, ipp);  // suction
  }
  // set inital value of Young's modulus
  long idem = Mm->ip[ipp].gemid();
  long ncompo  = Mm->givencompeqother(ipp, im);
  double e_ini = Mm->give_initial_ym(ipp, idem, ido+ncompo);
  Mm->ip[ipp].other[ido+ncompstr+13] = Mm->ip[ipp].eqother[ido+ncompstr+13] = e_ini;  
  // time step scaling factor
  Mm->ip[ipp].other[ido+ncompstr+14] = Mm->ip[ipp].eqother[ido+ncompstr+14] = 1.0;
  // set inital stiffness matrix  
  matrix d;
  makerefm(d,  Mm->ip[ipp].other+ido+ncompstr+15, ncompstr, ncompstr);
  Mm->elmatstiff(d, ipp, ido);
  copym(d, Mm->ip[ipp].eqother+ido+ncompstr+15);
  check_math_errel(Mm->elip[ipp]);

  return;
}


/**
   This function computes the value of yield functions.

   @param sig - full stress tensor in Voigt notation (6 components)
   @param q   - %vector of hardening parameter

   @retval The function returns value of yield function for the given stress tensor

   25.3.2002
*/
double camclaycoup::yieldfunction (vector &sig, vector &q)
{
  double i1s,j2s,f;
  
  i1s = first_invar (sig)/3.0;
  j2s = j2_stress_invar (sig);

  f = j2s/(m*m) + (i1s - q[1]) * (i1s - q[0]);
  
  return f;
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector sigma.

   @param sig - stress tensor in Voigt notation (6 components)
   @param q   - %vector of the hardening parameter
   @param dfds - %vector where the resulting derivatives are stored (6 components)


   4.1.2002
*/
void camclaycoup::deryieldfsigma (vector &sig, vector &q, vector &dfds)
{
  double volume,i1s;
  
  i1s = first_invar (sig)/3.0;
  
  deviator (sig,dfds);

  // q[0] = \bar{p}_c
  // q[1] = p_s = k * s
  cmulv(1.0/(m*m), dfds);
  volume=1.0/3.0*(2.0*i1s - q[0] - q[1]);
  dfds(0) += volume;
  dfds(1) += volume;
  dfds(2) += volume;  
  dfds(3) *= 2.0;
  dfds(4) *= 2.0;
  dfds(5) *= 2.0;
}



/**
   The function computes the second derivatives of yield function
   with respect of stress tensor sigma.

   @param ddfds - tensor of the 4-th order where the are derivatives stored

   19.12.2002
*/
void camclaycoup::dderyieldfsigma (matrix &ddfds)
{
  fillm(0.0, ddfds);
  ddfds[0][0] = ddfds[1][1] = ddfds[2][2] = 2.0/(3.0*m*m) + 2.0/9.0;
  ddfds[0][1] = ddfds[0][2] = ddfds[1][0] = ddfds[1][2] = -1.0/(3.0*m*m) + 2.0/9.0;
  ddfds[2][0] = ddfds[2][1] = ddfds[0][1];
  // doubled shear components due to tensor-> vector transformation
  ddfds[3][3] = 2.0/(m*m); 
  ddfds[4][4] = 2.0/(m*m);
  ddfds[5][5] = 2.0/(m*m);
}



/**
   The function computes derivatives of plastic potential function
   with respect of vector sigma.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameters
   @param ddgds - %matrix where the resulting derivatives are stored
*/
void camclaycoup::derpotsigma (vector &sig, vector &q, vector &ddgds)
{
  deryieldfsigma (sig, q, ddgds);
}



/**
   This function computes derivatives of as-th yield function
   with respect of vector of hradening parameters.

   @param sig - stress tensor
   @param q   - %vector of the hardening parameters
   @param dfds - %vector where the resulting derivatives are stored


   4.1.2002
*/
void camclaycoup::deryieldfq(vector &sig, vector &q, vector &dfq)
{
  //     df
  // ---------- = dfq[0]
  // d\bar{p}_c
  dfq[0] = -first_invar (sig)/3.0 + q[1];
  // df
  // --- = dfq[1]
  // p_s
  dfq[1] = -first_invar (sig)/3.0 + q[0];

  return;
}



/**
  The function computes the second derivatives of yield function
  with respect to hradening parameters.

  @param ddfdq - %matrix, where the resulting derivatives are stored


  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::deryieldfdqdq(matrix &ddfdq)
{  
  ddfdq[0][0] = 0.0;
  ddfdq[0][1] = 1.0;
  ddfdq[1][0] = 1.0;
  ddfdq[1][1] = 0.0;

  return;
}



/**
  The function computes the second derivatives of yield function
  with respect to stresses.

  @param dfds - tensor, where the resulting derivatives are stored.
                size of dfds = (6,number_of_hardening_param)


  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::deryieldfdsdq(matrix &dfdsdqt)
{  
  dfdsdqt[0][0] = dfdsdqt[1][0] = dfdsdqt[2][0] = -1.0/3.0;
  dfdsdqt[0][1] = dfdsdqt[1][1] = dfdsdqt[2][1] = -1.0/3.0;
  return;
}



/**
   The function computes derivatives of hardening paramters
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array
   @param sig - full stress tensor
   @param qtr  - %vector of hardening variables
   @param epsp - %vector of attained plastic strains
   @param dqdg - %vector where the resulting derivatives are stored in Voigt notation


  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::der_q_gamma(long ipp, long ido, vector &sig, vector &qtr, vector &epsp, vector &dqdg)
{
  double v_ini, v_lambda1, v_kappa1, p1, i1s, epsvp;
  double p_atm = 101325.0;
  double s, sr, fs, b, c, c1, ksi;
  long ncompstr = Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector epst(ASTCKVEC(6));
  
  // mean stress
  i1s = first_invar (sig)/3.0;
  // initial value of specific volume under small pressure p_1 for the normal consolidation line (v_{\lambda 1} = N)
  v_lambda1  = Mm->ic[ipp][0];
  // initial reference presure
  p1         = Mm->ic[ipp][1];
  // initial value of specific volume under small pressure p_1 for the swelling line (v_{\kappa 1})
  v_kappa1 = Mm->ip[ipp].eqother[ido+ncompstr+3];
  // specific volume for intial stress state
  v_ini = Mm->ip[ipp].eqother[ido+ncompstr+4];

  s  = Mm->givenonmechq(suction, ipp);  // suction pressure s
  //limitation of positive suction
  if(s >=0.0)
    s=0.0;
  s = fabs(s);
  sr = Mm->givenonmechq(saturation_degree, ipp);  // degree of saturation S_r

  // actual value of preconsolidation pressure for the fully saturated state
  give_full_vector (epst, epsp, ssst);
  epsvp = first_invar(epst);
  //double pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*epsvp)/(kappa-lambda));
  
  // specific volume for the actual preconsolidation pressure:
  //
  // v_{pc} = v_{\lambda1} - \lambda \ln(p_c/p_1)
  // double v_pc = v_lambda1 - lambda * log(pc_new/p1);

  // suction function:
  //
  //                  s / p_{atm}
  // f(s) = 1 + -----------------------
  //            10.7 + 2.4(s / p_{atm})
  //
  fs = 1.0 + (s/p_atm) / (10.7 + 2.4*(s/p_atm)); 
  if (fs < 1.0)
    fs = 1.0;

  // bonding variable:
  //
  // \ksi = f(s) (1 - S_r)
  //
  ksi = (1-sr)*fs;    

  //            \tilde{c}_1
  // c_1 = --------------------
  //           1
  //        -------------- + 1
  //        v_{pc}-1

  //  c1 = c1_tilde/(1.0/(v_pc-1.0)+1.0); // according to Borja
  c1 = c1_tilde;  // according to Gallipoli


  //
  // c(\ksi) = 1 - c_1 [1-\exp(c_2 \ksi)]
  //
  c = 1.0 - c1*(1.0-exp(c2*ksi));

  //
  //                \lambda - \kappa
  // b(\ksi) = --------------------------
  //            \lambda c(\ksi) - \kappa
  //
  b = (lambda - kappa)/(lambda*c - kappa);


  // q[0] = \bar{p}_c = -\exp[a(\ksi)] (-p_c)^{b(\ksi)}
  //               
  // p_c  = p_1 exp[(v_{\kappa_1} - v_{\lambda_1} + eps_{vp} v_{ini})/(\kappa - \lambda)]
  //
  // p_1           = const
  // p_{c0}        = const
  // \sigma_{m0}   = const
  // v_{\kappa_1}  = const
  // v_{p_{c0}}    = v_{\kappa_1} - \kappa  \ln(p_{c0}/p_1) = const
  // v_{ini}       = v_{p_{c0}} + \kappa \ln(p_{c0}/\sigma_{m0}) = const
  // v_{\lambda_1} = v_{p_{c0}}   - \lambda \ln(p_{c0}/p_1) = const
  // 
  // d \bar{p}_c                                              d p_c                             (-p_c)^{b(\ksi)}    d p_c  
  // ----------- = -\exp[a(\ksi)] b(\ksi)(-p_c)^{b(\ksi)-1} ---------- = -\exp[a(\ksi)] b(\ksi) ----------------  ---------- , 
  //  d \gamma                                               d \gamma                                p_c          d \gamma
  //                             
  //   d p_c               v_ini
  // ---------- = p_c ----------------- (2p - \bar{p}_c - p_s),
  //  d \gamma        \kappa - \lambda
  //
  //  d \bar{p}_c                b(\ksi) v_ini                              
  //  ----------- = -\bar{p}_c  ---------------- (2p - \bar{p}_c - p_s)
  //   d \gamma                 \kappa - \lambda                       
  //
  //
  dqdg[0] = -qtr[0] * (b*v_ini)/(kappa-lambda) * (2.0*i1s-qtr[0]-qtr[1]);


  // q[1] = p_s = k*s
  //
  //    d p_s
  // ----------- = 0.0 because p_s is constant for the given strain level
  //  d \gamma
  //
  dqdg[1] = 0.0;

  return;
}



/**
   The function computes derivatives of hardening paramters
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array
   @param sig - full stress tensor
   @param qtr    - %vector of hardening variables
   @param epsp   - %vector of attained plastic strains
   @param dqdgds - %matrix where the resulting derivatives are stored in Voigt notation, 
                   its dimensions must be qtr.n x 6


  Created by Tomas Koudelka,  05.2022
*/
void camclaycoup::dqpardgammadsigma(long ipp, long ido, vector &sig, vector &qtr, vector &epsp, matrix &dqdgds)
{
  double v_ini, v_lambda1, v_kappa1, p1, i1s, epsvp;
  double p_atm = 101325.0;
  double s, sr, fs, b, c, c1, ksi;
  long ncompstr = Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector epst(ASTCKVEC(6));
  
  // mean stress
  i1s = first_invar (sig)/3.0;
  // initial value of specific volume under small pressure p_1 for the normal consolidation line (v_{\lambda 1} = N)
  v_lambda1  = Mm->ic[ipp][0];
  // initial reference presure
  p1         = Mm->ic[ipp][1];
  // initial value of specific volume under small pressure p_1 for the swelling line (v_{\kappa 1})
  v_kappa1 = Mm->ip[ipp].eqother[ido+ncompstr+3];
  // specific volume for intial stress state
  v_ini = Mm->ip[ipp].eqother[ido+ncompstr+4];

  s  = Mm->givenonmechq(suction, ipp);  // suction pressure s
  //limitation of positive suction
  if(s >=0.0)
    s=0.0;
  s = fabs(s);
  sr = Mm->givenonmechq(saturation_degree, ipp);  // degree of saturation S_r

  // actual value of preconsolidation pressure for the fully saturated state
  give_full_vector (epst, epsp, ssst);
  epsvp = first_invar(epst);
  //double pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*epsvp)/(kappa-lambda));
  
  // specific volume for the actual preconsolidation pressure:
  //
  // v_{pc} = v_{\lambda1} - \lambda \ln(p_c/p_1)
  //double v_pc = v_lambda1 - lambda * log(pc_new/p1);

  // suction function:
  //
  //                  s / p_{atm}
  // f(s) = 1 + -----------------------
  //            10.7 + 2.4(s / p_{atm})
  //
  fs = 1.0 + (s/p_atm) / (10.7 + 2.4*(s/p_atm)); 
  if (fs < 1.0)
    fs = 1.0;

  // bonding variable:
  //
  // \ksi = f(s) (1 - S_r)
  //
  ksi = (1-sr)*fs;    

  //            \tilde{c}_1
  // c_1 = --------------------
  //           1
  //        -------------- + 1
  //        v_{pc}-1

  //  c1 = c1_tilde/(1.0/(v_pc-1.0)+1.0); // according to Borja
  c1 = c1_tilde;  // according to Gallipoli


  //
  // c(\ksi) = 1 - c_1 [1-\exp(c_2 \ksi)]
  //
  c = 1.0 - c1*(1.0-exp(c2*ksi));

  //
  //                \lambda - \kappa
  // b(\ksi) = --------------------------
  //            \lambda c(\ksi) - \kappa
  //
  b = (lambda - kappa)/(lambda*c - kappa);


  // q[0] = \bar{p}_c = -\exp[a(\ksi)] (-p_c)^{b(\ksi)}
  //               
  // p_c  = p_1 exp[(v_{\kappa_1} - v_{\lambda_1} + eps_{vp} v_{ini})/(\kappa - \lambda)]
  //
  // p_1           = const
  // p_{c0}        = const
  // \sigma_{m0}   = const
  // v_{\kappa_1}  = const
  // v_{p_{c0}}    = v_{\kappa_1} - \kappa  \ln(p_{c0}/p_1) = const
  // v_{ini}       = v_{p_{c0}} + \kappa \ln(p_{c0}/\sigma_{m0}) = const
  // v_{\lambda_1} = v_{p_{c0}}   - \lambda \ln(p_{c0}/p_1) = const
  // 
  // d \bar{p}_c                                              d p_c                             (-p_c)^{b(\ksi)}    d p_c  
  // ----------- = -\exp[a(\ksi)] b(\ksi)(-p_c)^{b(\ksi)-1} ---------- = -\exp[a(\ksi)] b(\ksi) ----------------  ---------- , 
  //  d \gamma                                               d \gamma                                p_c          d \gamma
  //                             
  //   d p_c               v_ini
  // ---------- = p_c ----------------- (2p - \bar{p}_c - p_s),
  //  d \gamma        \kappa - \lambda
  //
  //  d \bar{p}_c                b(\ksi) v_ini                              
  //  ----------- = -\bar{p}_c  ---------------- (2p - \bar{p}_c - p_s)
  //   d \gamma                 \kappa - \lambda                       
  //
  //
  // dqdg(0) = -qtr[0] * (b*v_ini)/(kappa-lambda) * (2.0*i1s-qtr[0]-qtr[1]);
  dqdgds(0,0) = dqdgds(0,1) = dqdgds(0,2) = -qtr[0] * (b*v_ini)/(kappa-lambda) * 2.0/3.0;


  // q[1] = p_s = k*s
  //
  //    d p_s
  // ----------- = 0.0 because p_s is constant for the given strain level
  //  d \gamma
  //
  dqdgds(1,0) = dqdgds(1,1) = dqdgds(1,2) = dqdgds(1,3) = dqdgds(1,4) = dqdgds(1,5) = 0.0;

  return;
}



/**
   The function computes derivatives of hardening paramters
   with respect of consistency parameter gamma.

   @param ipp - integration point number
   @param ido - index of internal variables for given material in the ipp other array
   @param sig - full stress tensor
   @param qtr    - %vector of hardening variables
   @param epsp   - %vector of attained plastic strains
   @param dqdgds - %matrix where the resulting derivatives are stored,  its dimensions must be qtr.n x qtr.n


  Created by Tomas Koudelka,  05.2022
*/
void camclaycoup::dqpardgammadqpar(long ipp, long ido, vector &sig, vector &qtr, vector &epsp, matrix &dqdgdq)
{
  double v_ini, v_lambda1, v_kappa1, p1, i1s, epsvp;
  double p_atm = 101325.0;
  double s, sr, fs, b, c, c1, ksi;
  long ncompstr = Mm->ip[ipp].ncompstr;
  strastrestate ssst = Mm->ip[ipp].ssst;
  vector epst(ASTCKVEC(6));
  
  // mean stress
  i1s = first_invar (sig)/3.0;
  // initial value of specific volume under small pressure p_1 for the normal consolidation line (v_{\lambda 1} = N)
  v_lambda1  = Mm->ic[ipp][0];
  // initial reference presure
  p1         = Mm->ic[ipp][1];
  // initial value of specific volume under small pressure p_1 for the swelling line (v_{\kappa 1})
  v_kappa1 = Mm->ip[ipp].eqother[ido+ncompstr+3];
  // specific volume for intial stress state
  v_ini = Mm->ip[ipp].eqother[ido+ncompstr+4];

  s  = Mm->givenonmechq(suction, ipp);  // suction pressure s
  //limitation of positive suction
  if(s >=0.0)
    s=0.0;
  s = fabs(s);
  sr = Mm->givenonmechq(saturation_degree, ipp);  // degree of saturation S_r

  // actual value of preconsolidation pressure for the fully saturated state
  give_full_vector (epst, epsp, ssst);
  epsvp = first_invar(epst);
  //double pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*epsvp)/(kappa-lambda));
  
  // specific volume for the actual preconsolidation pressure:
  //
  // v_{pc} = v_{\lambda1} - \lambda \ln(p_c/p_1)
  //double v_pc = v_lambda1 - lambda * log(pc_new/p1);

  // suction function:
  //
  //                  s / p_{atm}
  // f(s) = 1 + -----------------------
  //            10.7 + 2.4(s / p_{atm})
  //
  fs = 1.0 + (s/p_atm) / (10.7 + 2.4*(s/p_atm)); 
  if (fs < 1.0)
    fs = 1.0;

  // bonding variable:
  //
  // \ksi = f(s) (1 - S_r)
  //
  ksi = (1-sr)*fs;    

  //            \tilde{c}_1
  // c_1 = --------------------
  //           1
  //        -------------- + 1
  //        v_{pc}-1

  //  c1 = c1_tilde/(1.0/(v_pc-1.0)+1.0); // according to Borja
  c1 = c1_tilde;  // according to Gallipoli


  //
  // c(\ksi) = 1 - c_1 [1-\exp(c_2 \ksi)]
  //
  c = 1.0 - c1*(1.0-exp(c2*ksi));

  //
  //                \lambda - \kappa
  // b(\ksi) = --------------------------
  //            \lambda c(\ksi) - \kappa
  //
  b = (lambda - kappa)/(lambda*c - kappa);


  // q[0] = \bar{p}_c = -\exp[a(\ksi)] (-p_c)^{b(\ksi)}
  //               
  // p_c  = p_1 exp[(v_{\kappa_1} - v_{\lambda_1} + eps_{vp} v_{ini})/(\kappa - \lambda)]
  //
  // p_1           = const
  // p_{c0}        = const
  // \sigma_{m0}   = const
  // v_{\kappa_1}  = const
  // v_{p_{c0}}    = v_{\kappa_1} - \kappa  \ln(p_{c0}/p_1) = const
  // v_{ini}       = v_{p_{c0}} + \kappa \ln(p_{c0}/\sigma_{m0}) = const
  // v_{\lambda_1} = v_{p_{c0}}   - \lambda \ln(p_{c0}/p_1) = const
  // 
  // d \bar{p}_c                                              d p_c                             (-p_c)^{b(\ksi)}    d p_c  
  // ----------- = -\exp[a(\ksi)] b(\ksi)(-p_c)^{b(\ksi)-1} ---------- = -\exp[a(\ksi)] b(\ksi) ----------------  ---------- , 
  //  d \gamma                                               d \gamma                                p_c          d \gamma
  //                             
  //   d p_c               v_ini
  // ---------- = p_c ----------------- (2p - \bar{p}_c - p_s),
  //  d \gamma        \kappa - \lambda
  //
  //  d \bar{p}_c                b(\ksi) v_ini                              
  //  ----------- = -\bar{p}_c  ---------------- (2p - \bar{p}_c - p_s)
  //   d \gamma                 \kappa - \lambda                       
  //
  //
  // dqdg(0) = -qtr[0] * (b*v_ini)/(kappa-lambda) * (2.0*i1s-qtr[0]-qtr[1]);
  // wrong //dqdgdq(0,0) =  2.0*qtr[0] * (b*v_ini)/(kappa-lambda);
  dqdgdq(0,0) =  -(b*v_ini)/(kappa-lambda) * (2.0*i1s-2.0*qtr[0]-qtr[1]);
  dqdgdq(0,1) =  qtr[0] * (b*v_ini)/(kappa-lambda);


  // q[1] = p_s = k*s
  //
  //    d p_s
  // ----------- = 0.0 because p_s is constant for the given strain level
  //  d \gamma
  //
  dqdgdq(1,0) = dqdgdq(1,1) = 0.0;

  return;
}



/**
  This function computes plastic modulus.

  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array
  @param sig - stress tensor in Voigt notation
  @param qtr - %vector of hardening parameters
  @param epsp - %vector of attained plastic strains

  @retval The function returns value of the plastic modulus.

  Created by Tomas Koudelka,  09.2012
*/
double camclaycoup::plasmodscalar(long ipp, long ido, vector &sig, vector &qtr, vector &epsp)
{
  double ret;
  vector dfq(ASTCKVEC(2));
  vector dqg(ASTCKVEC(2));

  deryieldfq(sig, qtr, dfq);
  der_q_gamma(ipp, ido, sig, qtr, epsp, dqg);
  scprd(dfq, dqg, ret);

  return -ret;
}



/**
  This function computes new value of the hardening parameter q.

  @param ipp  - integration point pointer
  @param ido  - index of internal variables for given material in the ipp other array
  @param eps  - %vector of the attained strains in Voigt notation
  @param epsp - %vector of the attained plastic strains in Voigt notation
  @param sig  - attained stress tensor in Voigt notation
  @param ssst - stress/strain state parameter (for the used vector_tensor function)
  @param q    - %vector of the hardening parameters

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::updateq(long ipp, long ido, vector &/*eps*/, vector &epsp, vector &/*sig*/, vector &q)
{
  vector epst(ASTCKVEC(6));
  vector depsp(epsp.n);
  double v_kappa1, v_lambda1, p1, depsvp, v_ini, pc_new;
  double v_pc, p_atm = 101325.0;
  double s, sr, fs, a, b, c, c1, ksi;
  strastrestate ssst = Mm->ip[ipp].ssst;
  long ncompstr = Mm->ip[ipp].ncompstr;
  long i;
  
  // initial value of specific volume under small pressure p_1 for the normal consolidation line (v_lambda1 = 1+N)
  v_lambda1 = Mm->ic[ipp][0];
  // initial reference presure
  p1 = Mm->ic[ipp][1];
  // original specific volume before any loading
  v_kappa1 = Mm->ip[ipp].eqother[ido+ncompstr+3];
  // specific volume for intial stress state
  v_ini     = Mm->ip[ipp].eqother[ido+ncompstr+4];

//  Original version of pc calculation - problems with tension
//  vk = v_ini*(1+epsv)+kappa*log(i1s/p1);
//  pc_new = p1 * exp((vk-v_lambda1)/(kappa-lambda));
//
//  Alternative version of pc calculation - 'less' accurate in benchmarks
//  pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*epsv)/(kappa-lambda));
//  vk = v_ini*(1+epsv)+kappa*log(i1s/p1);

  for (i=0; i<ncompstr; i++)
    depsp[i] = epsp[i];// + Mm->ic[ipp][3+i];

  give_full_vector(epst, depsp, ssst);
  depsvp = first_invar(epst);

  pc_new = p1 * exp((v_kappa1 - v_lambda1 + v_ini*depsvp)/(kappa-lambda));

  //  For the partially saturated soil, 
  //  the effective preconsolidation pressure \bar{p}_c is related to the saturated 
  //  preconsolidation pressure p_c according to bonding variable \ksi
  //
  // \bar{p}_c = -exp[a(\ksi)] (-p_c)^{b(\ksi)}
  //
  //
 
  // specific volume for the actual preconsolidation pressure:
  //
  // v_{pc} = v_{\lambda1} - \lambda \ln(p_c/p_1)
  v_pc = v_lambda1 - lambda * log(pc_new/p1);

  s  = Mm->givenonmechq(suction, ipp);  // suction pressure s
  //limitation of positive suction
  if(s >=0.0)
    s=0.0;
  s = fabs(s);
  sr = Mm->givenonmechq(saturation_degree, ipp);  // degree of saturation S_r

  // suction function:
  //
  //                  s / p_{atm}
  // f(s) = 1 + -----------------------
  //            10.7 + 2.4(s / p_{atm})
  //
  fs = 1.0 + (s/p_atm) / (10.7 + 2.4*(s/p_atm)); 
  if (fs < 1.0)
    fs = 1.0;

  // bonding variable:
  //
  // \ksi = f(s) (1 - S_r)
  //
  ksi = (1-sr)*fs;    

  //            \tilde{c}_1
  // c_1 = --------------------
  //           1
  //        -------------- + 1
  //        v_{pc}-1
  // c1 = c1_tilde/(1.0/(v_pc-1.0)+1.0); // according to Borja
  c1 = c1_tilde; // according to Gallipoli


  //
  // c(\ksi) = 1 - c_1 [1-\exp(c_2 \ksi)]
  //
  c = 1.0 - c1*(1.0-exp(c2*ksi));

  //
  // a(\ksi) = v_{lambda_1} [c(\ksi)-1]/[\lambda c(\ksi) - \kappa)
  //
  a = (v_lambda1 - 1.0)*(c - 1.0)/(lambda*c - kappa);

  //
  //                \lambda - \kappa
  // b(\ksi) = --------------------------
  //            \lambda c(\ksi) - \kappa
  //
  b = (lambda - kappa)/(lambda*c - kappa);

  //
  // effective preconsolidation pressure:
  //
  // q[0] = \bar{p}_c = -exp[a(\ksi)] (-p_c)^{b(\ksi)}
  //
  q[0] = -exp(a)*pow(-pc_new, b);
  

  // apparent adhesion due to suction
  // (right intersection of the yield surface with p axis):
  //
  // q[1] = p_s = k s
  //
  q[1] = ks * s;
  return;
}



/**
  This function computes material stiffnes matrix.

  @param d - allocated matrix structure for material stiffness %matrix in the reduced format
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::matstiff (matrix &d, long ipp, long ido)
{
  double zero=1.0e-20;

  switch (Mp->nlman->stmat){
    case initial_stiff:
      //  initial elastic matrix
      Mm->elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff:
    case secant_stiff:
    case incr_tangent_stiff:
    case ijth_tangent_stiff:{
      //  consistent tangent stiffness matrix
      long n = Mm->ip[ipp].ncompstr;
      if (sra.give_tsra() == cp){
        // algorithmic consistent stiffness matrix for the cutting plane algorithm is computed in the course of stress return
        if (d.m != n){ // material stiffness matrix used for element in the case of plane stress/strain states
          for(long i=0; i<d.m; i++) // copy just the first d.n(=3) components at all matrix rows
            memcpy(&d(i,0), Mm->ip[ipp].other+ido+n+15+n*i, d.n*sizeof(*Mm->ip[ipp].other));
        }
        else{ // dimensions of material stiffness matrix matches the ones for element material stiffness matrix
          matrix auxm;
          makerefm(auxm, Mm->ip[ipp].other+ido+n+15, n, n);
          copym(auxm, d);
        }
      }
      else{
        //  compute algorithmic stiffness matrix for closest point projection stress return algorithm
        matrix de(ASTCKMAT(6, 6)), ddfdst(ASTCKMAT(6,6));
        matrix e(ASTCKMAT(6, 6)), auxd(ASTCKMAT(6,6)), tmp(ASTCKMAT(6, 6));
        double gamma = Mm->ip[ipp].other[ido+n];
        identm(e);
        Mm->elmatstiff(de, ipp, spacestress);
        dderyieldfsigma(ddfdst);
        cmulm(gamma, de, auxd);
        mxm(auxd, ddfdst, tmp);
        addm(e, tmp, auxd);
        invm(auxd, tmp, zero);
        mxm(tmp, de, auxd);
        //  matrix representation of the fourth order tensor
        tensor4_ematrix(d, auxd, Mm->ip[ipp].ssst);
      }
      break;
    }
  }
}



/**
  This function computes elasto-plastic material stiffnes matrix.

  @param d - allocated matrix structure for material stiffness %matrix in the reduced format
  @param ipp - integration point number
  @param ido - index of internal variables for given material in the ipp other array

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::depmatstiff (matrix &d, long ipp, long ido)
{
  switch (Mp->nlman->stmat){
    case initial_stiff:
      //  initial elastic matrix
      Mm->elmatstiff(d, ipp, ido);
      break;
    case tangent_stiff:
    case secant_stiff:
    case incr_tangent_stiff:
    case ijth_tangent_stiff:{
      //  consitent tangent stiffness matrix
      long n = Mm->ip[ipp].ncompstr;
      long ncompq = 2;
      matrix de(ASTCKMAT(n, n)), auxm;
      vector sig(ASTCKVEC(n)), q(ASTCKVEC(ncompq)), epsp(ASTCKVEC(n));
      vector dfds(ASTCKVEC(n)), dgds(ASTCKVEC(n));
      vector sigt(ASTCKVEC(6)), dfdst(ASTCKVEC(6)), dgdst(ASTCKVEC(6));
      vector dedgds(ASTCKVEC(n)), dedfds(ASTCKVEC(n));
      strastrestate ssst = Mm->ip[ipp].ssst;
      double denomh, denomd, denom, gamma, f, err=sra.give_err();

      // prepare actual values of needed quantities
      copyv(Mm->ip[ipp].stress, sig);     // stresses
      copyv(Mm->ip[ipp].other+ido, epsp); // plastic strains
      gamma = Mm->ip[ipp].other[ido+n+0]; // cumulated plastic multiplier
      q(0)  = Mm->ip[ipp].other[ido+n+1]; // hardening parameter \bar{p}_c
      q(1)  = Mm->ip[ipp].other[ido+n+2]; // hardening parameter k_s
      give_full_vector(sigt, sig, ssst);  // 6 component vector of stresses

      // check actual stress state: f <= 0.0
      f = yieldfunction(sigt, q);
      if ((f < -err) || (gamma == 0.0 && (f < err))){
        // stress state is in the elastic domain,
        // i.e. material has been subjected by elastic loading/unloading ->
        // -> return actual elastic stiffness matrix
        Mm->elmatstiff(d, ipp, ssst);
        return;
      }

      // stress state is on/close to  the yield surface
      // material has been subjected by plastic loading ->

      // compute elasto-plastic matrix

      //                 D_e \frac{dg}{d\sigma} (\frac{df}{d\sigma})^T D_e
      // D_{ep} = D_e -  ---------------------------------------------------
      //                  (\frac{df}{d\sigma})^T D_e \frac{dg}{d\sigma} + H
      //
      // H = -\frac{df}{dq} \frac{dq}{d\lambda}
      
      Mm->elmatstiff(de, ipp, ssst);
      deryieldfsigma(sigt, q, dfdst);
      derpotsigma(sigt, q, dgdst);
      //  conversion of full tensor notation to vector notation
      give_red_vector(dfdst, dfds, ssst);
      //  conversion of full tensor notation to vector notation
      give_red_vector(dgdst, dgds, ssst);
      mxv(de,dgds,dedgds);
      mxv(de,dfds,dedfds);
      denomh = plasmodscalar(ipp, ido, sig, q, epsp);
      scprd(dfds, dedgds, denomd);
      // preserve minimum stiffness for perfect plastic material, i.e. if the hardening modulus is 'close' to zero (H->0)
      //      if (fabs(denomh/denomd) < 1.0e-3)
      //        denomh = sgn(denomh)*1.0e-3*denomd;
      denom = denomd+denomh;

      if (n == d.m){
        // material stiffness matrix has identical dimensions with the ones of stiffness matrix used in
        // assembling of element stiffness matrix on finite element -> auxm will be reference on matrix d
        makerefm(auxm, d.a, d.m, d.n);
        nullm(auxm);
      }
      else{
        // material stiffness matrix has different dimensions than the stiffness matrix used in
        // assembling of element stiffness matrix on finite element -> auxm will be reference on matrix d
        reallocm(RSTCKMAT(n, n, auxm));
      }
      // calculate matrix in the nominator
      vxv(dedgds, dedfds, auxm);
      // scaling by the denominator
      cmulm(1.0/denom, auxm);
      // resulting material stiffness matrix
      subm(de, auxm, auxm);

      // resulting material stiffness matrix with dimensions adjusted for finite element computation
      if (n != d.m){
        // resulting material stiffness matrix dimensions must be adjusted for finite element computation
        extractblock(auxm, d, 0, 0);
      }
      else{
        // nothing is needed, matrix auxm is direct reference to the matrix in argument d
        // hence d contains resulting material stiffness matrix with appropriate dimensions for
        // the stiffness matrix on the finite element
      }
      /*
      double zero=1.0e-20;
      matrix de(ASTCKMAT(6, 6)), ddfdst(ASTCKMAT(6,6));
      matrix e(ASTCKMAT(6, 6)), auxd(ASTCKMAT(6,6)), tmp(ASTCKMAT(6, 6));
      long n = Mm->ip[ipp].ncompstr;
      double gamma = Mm->ip[ipp].other[ido+n];
      identm(e);
      Mm->elmatstiff(de, ipp, spacestress);
      dderyieldfsigma(ddfdst);
      cmulm(gamma, de, auxd);
      mxm(auxd, ddfdst, tmp);
      addm(e, tmp, auxd);
      invm(auxd, tmp, zero);
      mxm(tmp, de, auxd);
      //  matrix representation of the fourth order tensor
      tensor4_ematrix(d, auxd, Mm->ip[ipp].ssst);*/
      break;
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

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::nlstresses (long ipp, long im, long ido)
  //
{
  int errr;
  long i,ni,ncomp=Mm->ip[ipp].ncompstr, matcomp=0, ret;
  double gamma, err, pnewdt = 1.0;
  vector epsn(ASTCKVEC(ncomp)), epsp(ASTCKVEC(ncomp)), q(ASTCKVEC(2));
  vector epsa(ASTCKVEC(ncomp)), sig(ASTCKVEC(ncomp)), dfds(ASTCKVEC(ncomp)), dgds(ASTCKVEC(ncomp));
  matrix d;
  vector dfdst(ASTCKVEC(6)), dgdst(ASTCKVEC(6));
  vector epst(ASTCKVEC(6)), epsnt(ASTCKVEC(6)), epspt(ASTCKVEC(6)), sigt(ASTCKVEC(6));
  //vector dev(ASTCKVEC(6));
  double i1s, j2s, epsv, epsvp;
  double dtime = Mp->timecon.actualforwtimeincr();

  //  initial values
  for (i=0; i<ncomp; i++)
  {
    epsn[i] = Mm->ip[ipp].strain[i];
    epsp[i] = Mm->ip[ipp].eqother[ido+i];
    // It is necessary to substract initial strains in order to invoke initial nonzero stress 
    // state for zero displacements (i.e. zero strains)
    // sig = E*(eps - eps_p), sig_0 = E*eps_0. For eps = 0: sig = sig_0  => eps_p = -eps_0
    // epsp[i] -= Mm->ic[ipp][3+i];
  }
  gamma=Mm->ip[ipp].eqother[ido+ncomp];
  q[0] = Mm->ip[ipp].eqother[ido+ncomp+1];  // effective preconsolidation pressure \bar{p}_c
  q[1] = Mm->ip[ipp].eqother[ido+ncomp+2];  // apparent adhesion  p_s = k * s

  //  stress return algorithm
  switch(sra.give_tsra ())
  {
    case cp:
      ni=sra.give_ni ();
      err=sra.give_err ();
      matcomp = actualize_stiffmat(*Mp->nlman, Mp->istep, Mp->jstep);
      if (Mp->istep >= 265){
        matcomp = 0;
      }
      if (matcomp)   reallocm(RSTCKMAT(ncomp, ncomp, d));
      ret = Mm->cutting_plane (ipp, im, ido, gamma, epsn, epsp, q, ni, err, matcomp, d);
      if (ret)
        pnewdt = 0.5;
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
  //  new data storage
  for (i=0; i<ncomp; i++){
    //epsp[i] += Mm->ic[ipp][3+i];    
    Mm->ip[ipp].other[ido+i]=epsp[i];
  }
  Mm->ip[ipp].other[ido+ncomp]=gamma;
  Mm->ip[ipp].other[ido+ncomp+1]=q[0]; // the first hardening parameter \bar{p}_c
  Mm->ip[ipp].other[ido+ncomp+2]=q[1]; // the second hardening parameter p_s
  // other[ido+ncomp+3] = v_lambda1; it is the initial value and therefore it does not require actualization
  // other[ido+ncomp+4] = v_ini; it is the initial value and therefore it does not require actualization

  Mm->givestress(0, ipp, sig);

  give_full_vector(sigt,sig,Mm->ip[ipp].ssst);
  give_full_vector(epsnt,epsn,Mm->ip[ipp].ssst);
  give_full_vector(epspt,epsp,Mm->ip[ipp].ssst);

  i1s = first_invar (sigt)/3.0;
  //  the second invariant of deviator
  //  it is expressed with the help of the stress components
  //  components of the deviator are not needed
  j2s = sqrt(j2_stress_invar (sigt));

  epsv = first_invar (epsnt);
  epsvp = first_invar (epspt);
  Mm->ip[ipp].other[ido+ncomp+5] = i1s;  
  Mm->ip[ipp].other[ido+ncomp+6] = j2s;  
  Mm->ip[ipp].other[ido+ncomp+7] = epsv;  
  Mm->ip[ipp].other[ido+ncomp+8] = epsvp;
  Mm->ip[ipp].other[ido+ncomp+9] = Mm->givenonmechq(saturation_degree, ipp);
  // volumetric strain rate
  //Mm->ip[ipp].other[ido+ncomp+10] = (epsv - Mm->ip[ipp].eqother[ido+ncomp+7])/dtime;
  // the volumetric strain rate is computed via generalized trapesoidal rule
  Mm->ip[ipp].other[ido+ncomp+10] = 0.5*((epsv - Mm->ip[ipp].eqother[ido+ncomp+7])/dtime + Mm->ip[ipp].eqother[ido+ncomp+10]);

  // reference presure = p1
  double p1 = Mm->ic[ipp][1];
  // specific volume at the reference pressure p_1 = v_lambda1 
  double v_lambda1 =  Mm->ic[ipp][0];
  // specific volume for the actual preconsolidation pressure = v_pc
  double v_pc = v_lambda1 - lambda * log(q[0]/p1);
  // specific volume for the actual pressure p
  double v;
  if (i1s < 0.0)
    v = v_pc + kappa*log(q[0]/i1s);
  else
    v = v_pc + kappa*log(-q[0]);
  // actual porosity n
  Mm->ip[ipp].other[ido+ncomp+11] = (v-1.0)/v;
  // suction
  Mm->ip[ipp].other[ido+ncomp+12] = Mm->givenonmechq(suction, ipp);
  // Young modulus
  Mm->ip[ipp].other[ido+ncomp+13] = give_actual_ym(ipp, im, ido);
  // time step scaling factor
  Mm->ip[ipp].other[ido+ncomp+14] = pnewdt;

  if (matcomp)
    copym(d, Mm->ip[ipp].other+ido+ncomp+15);

  errr = check_math_errel(Mm->elip[ipp]);
}




/**
  Function returns flag whether actualize the stiffness matrix according to the setup of the nonlinear solver and
  attained number of steps.

  @param[in] smt - type of the stiffness matrix modification in NR method
  @param[in] istep - number of load increment step in NR method
  @param[in] jstep - number of inner iteration step in NR method

  @retval 0 - in the case that actualization of the stiffness matrix is NOT required
  @retval 1 - in the case that actualization of the stiffness matrix is required
*/
long camclaycoup::actualize_stiffmat(const nonlinman &nlman, long istep, long jstep)
{
  long ret = 0;
  switch (nlman.stmat){
    case initial_stiff:
      break;
    case tangent_stiff:
      ret = 1;
      break;
    case secant_stiff:
    case incr_tangent_stiff:{
      if (jstep <= 0)
        ret = 1;
      break;
    }
    case ijth_tangent_stiff:{
      if (((jstep <= 0) && (istep % nlman.ithstiffchange == 0)) ||
          ((jstep+1) % nlman.jthstiffchange == 0))
        ret = 1;
      break;
    }
    default:{
      print_err("unknown type of stiffness matrix is required",__FILE__,__LINE__,__func__);
    }
  }
  return ret;
}



/**
  This function updates values in the other array reached in the previous equlibrium state to
  values reached in the new actual equilibrium state.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::updateval (long ipp, long im, long ido)
{
  long i,n = Mm->ip[ipp].ncompstr;

  double e = comp_actual_ym(ipp, im, ido);
  
  for (i=0;i<n;i++){
    Mm->ip[ipp].eqother[ido+i]=Mm->ip[ipp].other[ido+i];
  }
  Mm->ip[ipp].eqother[ido+n]   = Mm->ip[ipp].other[ido+n];       // gamma
  Mm->ip[ipp].eqother[ido+n+1] = Mm->ip[ipp].other[ido+n+1];     // p_c (hardening parameter - effective preconsolidation pressure) 
  Mm->ip[ipp].eqother[ido+n+2] = Mm->ip[ipp].other[ido+n+2];     // p_s (hardening parameter - apparent adhesion) 
  //
  // Following three values of the specific volume are constant and they 
  // are stored directly in the eqother array during the intial step
  //
  //  Mm->ip[ipp].eqother[ido+n+3] = Mm->ip[ipp].other[ido+n+3];     // v_lambda1
  //  Mm->ip[ipp].eqother[ido+n+4] = Mm->ip[ipp].other[ido+n+4];     // v_ini
  //
  Mm->ip[ipp].eqother[ido+n+5] = Mm->ip[ipp].other[ido+n+5];     // i1s;  
  Mm->ip[ipp].eqother[ido+n+6] = Mm->ip[ipp].other[ido+n+6];     // j2s;  
  Mm->ip[ipp].eqother[ido+n+7] = Mm->ip[ipp].other[ido+n+7];     // epsv;  
  Mm->ip[ipp].eqother[ido+n+8] = Mm->ip[ipp].other[ido+n+8];     // epsvp;    
  // saturation degree
  Mm->ip[ipp].eqother[ido+n+9] = Mm->ip[ipp].other[ido+n+9];    // Sr
  // volumteric strain rate 
  Mm->ip[ipp].eqother[ido+n+10] = Mm->ip[ipp].other[ido+n+10];    // depsv/dt
  // actual porosity n
  Mm->ip[ipp].eqother[ido+n+11] = Mm->ip[ipp].other[ido+n+11]; 
  // suction
  Mm->ip[ipp].eqother[ido+n+12] = Mm->ip[ipp].other[ido+n+12];
  // Young modulus
  Mm->ip[ipp].eqother[ido+n+13] = Mm->ip[ipp].other[ido+n+13] = e;
  // time step scaling factor
  Mm->ip[ipp].eqother[ido+n+14] = Mm->ip[ipp].other[ido+n+14];
  // stiffness matrix
  for (i=0;i<n*n;i++){
    Mm->ip[ipp].eqother[ido+n+15+i] = Mm->ip[ipp].other[ido+n+15+i];
  }
}



/**
  The function returns actual Young modulus value.

  @param ipp - integration point number
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

*/
double camclaycoup::give_actual_ym(long ipp, long /*im*/, long ido)
{
  long ncomp=Mm->ip[ipp].ncompstr;  
  double e = Mm->ip[ipp].other[ido+ncomp+13];

  return e;
}



/**
  The function computes actual Young modulus value.

  @param ipp - integration point number
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array

*/
double camclaycoup::comp_actual_ym(long ipp, long im, long ido)
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
  vector epsp(ASTCKVEC(ncomp));
  vector sig(ASTCKVEC(ncomp));
  vector q(ASTCKVEC(1));

  // actual total volumetric strain
  Mm->givestrain(0, ipp, eps);
  give_full_vector(epst, eps, ssst);
  double epsv = first_invar(epst);
  // actual plastic volumteric strain
  //Mm->giveother(ipp, ido, ncomp, epsp);
  // total volumetric strain from the last converged step
  double epsvo = Mm->ip[ipp].eqother[ido+ncomp+7];
  // total volumetric plastic strain from the last converged step
  double epsvpo = Mm->ip[ipp].eqother[ido+ncomp+8];

  // actual plastic strain
  for(long i=0; i<ncomp; i++)
    epsp(i) = Mm->ip[ipp].other[ido];
  give_full_vector(epst, epsp, ssst);
  double epsvp = first_invar(epst);

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
  double v_ini = Mm->ip[ipp].eqother[ido+ncomp+4];
  //double v     = v_ini*(1.0+epsv);

  //double p = pc/exp((v-v_pc)/kappa);
  double po = Mm->ip[ipp].eqother[ido+ncomp+5];
  if ((po == 0.0) || (epsv == 0.0))
    return e_ini;
  if (po > pc){
    // for the unloading: \Delta eps^{el}_V = \Delta eps_V
    // actual p according actual swelling line
    double epsve = (epsv-epsvo) - (epsvp - epsvpo);
    double p = po*exp(-(epsve)*v_ini/kappa);
    if (epsve == 0.0)
      e = 3.0*(1.0-2.0*nu)*v_ini*(-p)/kappa;
    else 
      e = 3.0*(1.0-2.0*nu)*(p-po)/(epsve);
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
double camclaycoup::give_actual_nu(long ipp, long im, long ido)
{
 long idem = Mm->ip[ipp].gemid();
 long ncompo = Mm->givencompeqother(ipp, im);
 double nu = Mm->give_initial_nu(ipp, idem, ido+ncompo);

 return nu;
}


/**
  The function returns time step scaling factor required by the model. It is represented by ratio 
  of required time step size to the actual one or 1.0 for no change in time step size.

  @param ipp - integration point number in the mechmat ip array.
  @param im  - index of material type for given ip
  @param ido - index of internal variables for given material in the ipp other array
  
  @return The function returns required time step scaling factor provided by the camclaycoup model

  Created by Tomas Koudelka, 5.2022
*/
double camclaycoup::dstep_red(long ipp, long /*im*/, long ido)
{
  long ncomp  =  Mm->ip[ipp].ncompstr;
  return (Mm->ip[ipp].other[ido+ncomp+14]);
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the ipp other array
  @param epsp  - %vector of irreversible strains
 
  Returns vector of irreversible strains via parameter epsp

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::giveirrstrains (long ipp, long ido, vector &epsp)
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

  @retval The function returns value of consistency parameter.

  Created by Tomas Koudelka,  09.2012
*/
double camclaycoup::give_consparam (long ipp,long ido)
{ 
  long ncompstr;
  double gamma;

  ncompstr=Mm->ip[ipp].ncompstr;
  gamma = Mm->ip[ipp].eqother[ido+ncompstr];

  return gamma;
}



/**
  The function extracts saturation degree s from the integration point eqother array.
  This value remains constant from initialization of the model.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of saturation degree.

  Created by Tomas Krejci, 14/12/2014
*/
double camclaycoup::give_saturation_degree(long ipp, long ido)
{
  long ncompstr;
  double s;

  ncompstr=Mm->ip[ipp].ncompstr;
  s = Mm->ip[ipp].other[ido+ncompstr+9];
  
  return s;
}

/**
  The function extracts preconsolidation pressure p_c for the attained state
  from the integration point other array.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of preconsolidation pressure.

  Created by Tomas Koudelka, 8.10.2013  
*/
double camclaycoup::give_preconspress(long ipp, long ido)
{
  long ncompstr;
  double pc;

  ncompstr=Mm->ip[ipp].ncompstr;
  pc = Mm->ip[ipp].other[ido+ncompstr+1];
  
  return pc;
}



/**
  The function extracts virgin void ratio e_lambda1 for the attained state (it is constant initial value)
  from the integration point eqother array. This value remains constant from initialization of the model.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of virgin porosity e_lambda1.

  Created by Tomas Koudelka, 8.10.2013  
*/
double camclaycoup::give_virgporosity(long ipp, long ido)
{
  long ncompstr;
  double e_lambda1;

  ncompstr=Mm->ip[ipp].ncompstr;
  e_lambda1 = Mm->ip[ipp].eqother[ido+ncompstr+3]-1.0;
  
  return e_lambda1;
}



/**
  The function extracts initial porosity n_ini for the attained state (it is constant initial value)
  from the integration point eqother array. This value remains constant from initialization of the model.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of initial porosity e_ini.

  Created by Tomas Koudelka, 8.10.2013  
*/
double camclaycoup::give_iniporosity(long ipp, long ido)
{
  long ncompstr;
  double v_ini, n_ini;

  ncompstr=Mm->ip[ipp].ncompstr;
  v_ini = Mm->ip[ipp].eqother[ido+ncompstr+4]; // v_ini = 1.0 + e_ini
  n_ini = (v_ini-1.0)/v_ini;
  
  return n_ini;
}



/**
  The function extracts porosity for the actual state from the integration point other array. 

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @retval The function returns value of actual porosity e.

  Created by Tomas Koudelka, 06.2020  
*/
double camclaycoup::give_porosity(long ipp, long ido)
{
  long ncompstr;
  double n;

  ncompstr=Mm->ip[ipp].ncompstr;
  n = Mm->ip[ipp].other[ido+ncompstr+11];
  
  return (n);
}



/**
  The function returns rate of the volumetric strain at the given integration point.

  @param ipp - integration point number in the mechmat ip array.
  @param ido - index of internal variables for given material in the ipp other array

  @return The function returns rate of the volumetric strain.
  
  Created by Tomas Koudelka 06.2020
*/
double camclaycoup::give_strain_vol_rate(long ipp, long ido)
{
  long ncomp = Mm->ip[ipp].ncompstr;

  return Mm->ip[ipp].other[ido+ncomp+10];
}



/**
  The funtion marks required non-mechanical quantities in the array anmq.

  @param anmq - array with flags for used material types
                anmq[i] = 1 => qunatity type nonmechquant(i+1) is required
                anmq[i] = 0 => qunatity type nonmechquant(i+1) is not required

  @return The function does not return anything, but it may change content of anmq array.
*/
void camclaycoup::give_reqnmq(long *anmq)
{
  anmq[saturation_degree-1] = 1;
  anmq[suction-1] = 1;
}



/**
  The function changes material parameters for stochastic analysis.
   
  @param atm - selected material parameters (parameters which are changed)
  @param val - array containing new values of parameters
   
  @return The function does not return anything.

  Created by Tomas Koudelka,  09.2012
*/
void camclaycoup::changeparam (atsel &atm,vector &val)
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
      case 3:
        c1_tilde=val[i];
        break;
      case 4:
        c2=val[i];
        break;
      case 5:
        ks=val[i];
        break;
      default:
        print_err("wrong number of atribute", __FILE__, __LINE__, __func__);
    }
  }
}
